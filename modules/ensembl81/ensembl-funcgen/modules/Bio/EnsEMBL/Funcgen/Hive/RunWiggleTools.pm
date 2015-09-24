
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::RunWiggleTools

=head1 DESCRIPTION

A simple module to convert the bam files associated with a ResultSet into a BigWig

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools;

use warnings;
use strict;
use File::Temp                             qw( tempfile );
use Fatal                                  qw( close );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names 
                                               run_system_cmd
                                               run_backtick_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( write_chr_length_file);
use Bio::EnsEMBL::Utils::Exception         qw( throw );

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

# Todo
# 1 Move this to SeqTools so it can be re-used in other contexts?
#   Limited utility as it is already a very thin wrapper.
#   But convert_bam_to_bigwig would site well with other SeqTools
#   convert methods
#   That would spoil the generic functionality of RunWiggleTools
#   Other option would be to create a WiggleTools.pm wrapper
#   And have SeqTools::convert_bam_to_bigwig point to that?
# 2 Add file version support? Or handle this separately?
# 3 Change directory naming? There is not need for result_feature reference now?
#   Remember we can't change the structure on nfs live easily!
#   Is there any chance that get_alignment_files_by_ResultSet_formats


####################################################
## Preparations
####################################################


sub fetch_input {
  my $self = shift;
  $self->param('disconnect_if_idle', 1);  # Set before DB connection via SUPER::fetch_input
  $self->SUPER::fetch_input;
  $self->get_param_method('reduce_operator', 'silent', '');
  $self->get_param_method('map_operator', 'silent', '');
  $self->get_param_method('map_args', 'silent', '');
  my $mode = $self->get_param_method('mode', 'silent');

  # We could (cross-)validate the operators/mode here
  # but let's just let wiggletools fail for now.
  if($mode){
    my $mode_method = '_build_'.lc($mode).'_cmd';

    if(! $self->can($mode_method)){
      $self->throw_no_retry("Unsupported mode:\t".$mode);
    }
    
    $self->set_param_method('mode_method', $mode_method);
  }


  # Need an RPKM mode which does all the right things
  # samtools idxstats file.bam to get num reads



  if($self->get_param_method('set_type', 'silent')){  # Set input params
    my $rset = $self->fetch_Set_input('ResultSet');
    $self->helper->debug(1, "RunWiggleTools::fetch_input got ResultSet:\t".$rset->name);

    # check we haven't sepcified input_files
    if($self->get_param_method('input_files', 'silent')){
      $self->throw_no_retry("It is unsafe to specify the input_files are set parameters.\n".
        'Please remove input_files from input_id');
    }

    if($rset->has_status('IMPORTED')){
      # This is not handling assembly specific IMPROTED states!
      # But we have moved to separate DBs for idfferent assemblies
      # So we need to drop that?


      if($self->param_silent('force')){ #force/recover
        $rset->adaptor->revoke_status('IMPORTED', $rset);
      }
      else{
        $self->throw_no_retry("Cannot write bigWig as result set is already marked as IMPORTED:\t".
          $rset->name.'('.$rset->dbID.")\nTo over-ride, run reseed_jobs.pl -append '{force>1}'");
        # DebugJob -f will not do this, as this is a beekeeper option
      }
    }

    # This could potentially trigger a filter/sort step to produce the 
    # filtered bam. This would be unsafe if parallel jobs are using the same file
    # This normally is done directly after merge, so should be safe, unless the sorted/filtered
    # bam has been deleted by mistake
    $self->input_files([$self->get_alignment_files_by_ResultSet_formats($rset, ['bam'])->{bam}]);
    $self->set_param_method('output_format', 'BigWig');
    $self->get_param_method('output_prefix',
                            'silent',
                            $self->db_output_dir.'/result_feature/'.$rset->name);
  }
  else{  # File input params
    my $inp_files = $self->get_param_method('input_files', 'required');  # Let's be forgiving

    if(! ref($inp_files)){
      $self->input_files([$inp_files]);
    }
    else{
      assert_ref($inp_files, 'ARRAY', 'input_files');
    }

    $self->get_param_method('output_format', 'silent', 'wig');  # default to wig
    $self->get_param_method('output_prefix', 'required');

    $self->helper->debug(1, "RunWiggleTools::fetch_input got input_files:\n\t".
      join("\n\t", @{$self->input_files}));
  }


  # TODO
  # 

  # This may not actually be used if output_file is specified
  $self->get_output_work_dir_methods($self->db_output_dir.'/result_feature/');
  return;
}


####################################################
## Core function
####################################################

sub run {
  my $self        = shift;
  
  # DEFINE THE OUTPUT
  my $write_out    = '-';
  my $reformat_cmd = '';
  my ($output, $cmd, @tmpfiles);

  if(lc($self->output_format) eq 'wig'){
    $write_out = $self->output_prefix.'.wig';
    $output    = $write_out;
  }
  elsif(lc($self->output_format) eq 'bigwig'){
    push @tmpfiles, write_chr_length_file($self->slice_objects);

    # This is now in /tmp, so maybe add to DESTROY?
    $output       = $self->output_prefix.'.bw';
    $reformat_cmd = ' | wigToBigWig -fixedSummaries stdin '.$tmpfiles[0].' '.$output;
  }
  else{
    $self->throw_no_retry("Output format not supported:\t".$self->output_format);
    # Should really be done in fetch_input, but easier here.
  }

  # Build map command
  my $mode_method = $self->mode_method;

  if($mode_method){
    $cmd = $self->$mode_method();
  }
  else{
    $cmd = $self->reduce_operator.' '.$self->map_operator.' '.$self->map_args.' '.
     join(' ', @{$self->input_files});
  }

  # RUN THE COMMAND
  # Presence of index files is input format specific, so not explicitly tested/generated here
  # Always use write to avoid redirects
  $cmd = "wiggletools write $write_out ".$cmd.$reformat_cmd;
  $self->helper->debug(1, "Running:\n\t".$cmd);
  run_system_cmd($cmd);
  # Simple 'bam to wig' took about 20 mins with ~ 8GB mem for single bam

  # PUT OUT THE TRASH
  if(@tmpfiles){
    unlink @tmpfiles or
     warn "Could not unlinke tmp file:\t$!\n";
   }

  # UPDATE THE DB
  if($self->set_type){  # Update dbfile_registry.path and set ResultSet states
    my $rset = $self->ResultSet;
    $rset->adaptor->dbfile_data_root($self->db_output_dir);
    $rset->dbfile_path($output);
    $rset->adaptor->store_dbfile_path($rset);
    $rset->adaptor->set_imported_states_by_Set($rset);
  }

  return;
}


sub write_output {  return; }  # Nothing to do/flow here?
# Will still autoflow if wired on branch 1.


sub _build_rpkm_cmd{
  my $self = shift;
  my $cmd  = 'mean ';
  my ($total_mapped, $count_cmd);

  # This assumes bam input with index
  foreach my $file(@{$self->input_files}){
    
    $count_cmd = 'samtools idxstats '.$file.' | awk \'{total = total + $2} END{print total}\'';
    $self->helper->debug(2, "Running:\n\t".$count_cmd);
    $total_mapped = run_backtick_cmd($count_cmd);
    # pipe causes uncaught failures on absence of bai file, so test here?

    if(! $total_mapped){
      $self->throw_no_retry("Failed to get number of mapped reads from index of:\n\t".$file);
    }

    $cmd .= ' scale '.(10**9 / $total_mapped).' '.$file;
  }

  return $cmd;
}


1;
