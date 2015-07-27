=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::RunSWEmbl

=head1 DESCRIPTION

'RunSWEmbl' Runs SWEmbl peak caller and stores peaks as an annotated feature set.
Assumes Files are organized and named with a specific convention
($repository_folder)/experiment_name/cell-type_feature-type/
unless specific files are given

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::RunSWEmbl;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
#use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args generate_slices_from_names strip_param_flags);

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::SWEmbl');

#use Data::Dumper;

sub fetch_input {   # fetch and preprocess the input file plus do some more checking
  my $self = shift @_;

  $self->SUPER::fetch_input();

  # Consider using: Bio::EnsEMBL::Helper to debug and/or log to a file (given parameter)
  
  my $efgdba = $self->_efgdba();
  my $fsa = $efgdba->get_FeatureSetAdaptor();
  if(!$self->_feature_set($fsa->fetch_by_name($self->_feature_set_name()))){
    throw "Feature Set was not Created";
  }

  #if(!$self->_feature_set($fsa->fetch_by_name($self->_feature_set_name()))){
  #  throw "Feature set was not created. Maybe you're skipping the Setup Step of the Pipeline";
  #}

  my $analysis = $self->_analysis();
  #Use the definitions that are on the database
  $self->_command($analysis->program_file);
  $self->_command_options($analysis->parameters);
  
  my $cell_type = $self->_cell_type()->name;
  my $feature_type = $self->_feature_type()->name;
  my $experiment_name = $self->_experiment_name();
  
  my $file_type = $self->_file_type();

  my $output_dir = $self->_output_dir();  

  my $input_dir = $self->_input_dir();  
  
  my $file_name = $self->_input_file();
  my $input_file =  $input_dir."/".$file_name;
  if(-e $input_file){ 
    my $output_file = $output_dir."/".$file_name;
    #TODO Validate if existent file is ok. 
    #TODO Add flag to enable override existent file...
    if(!$self->param('reenter')){
      $self->_preprocess_file($input_file, $output_file, $file_type) || 
	throw "Error processing data file $input_file";
    } else {
      if(! -e $output_file){ warn "$output_file not found. May need to rerun from scratch"; }
    }
    $self->_input_file($output_file); 
    
    if(!$self->_results_file($self->param('results_file'))){ 
      $self->_results_file($output_file.".".$analysis->logic_name);
    }
  } else { throw "No valid data file was given: ".$input_file; }

  if(!$self->_skip_control()){
    my $control_file = $output_dir."/".$self->_control_file();
    if(! -e $control_file){ throw "No valid control file was given: ".$control_file; }
    $self->_control_file($control_file); 
  }

  return 1;
}




sub run {   # call SWEmbl
  my $self = shift @_;

  if($self->param('reenter')){ return 1; }

  #Don't leave this here... as it can go wrong!!
  #if( -e $self->_results_file()){
  #  return 1;
  #}

  my %fswitches = (
  		   bed   => '-B',
  		   sam   => '-S',
  		   bam   => '-F',
  		   #maq   => '-M',
  		   #eland => '-E',
  		  );
  
  #Ideally it will unify to sam or bam
  my $format_switch = $fswitches{$self->_file_type()};
  
  throw("Could not identify valid SWEmbl input format for ".$self->_file_type()) if(! $format_switch);
  
  my $command = $self->_bin_dir()."/".$self->_command() . " $format_switch -V -i " . $self->_input_file() . ' ' . 
    $self->_command_options() . ' -o ' . $self->_results_file();
  $command .= ' -z ' if ($self->_input_file() =~ m/\.gz$/);
  if($self->_control_file() && !$self->_skip_control()){  
    $command = $command." -r ".$self->_control_file(); 
  }
  
  #warn "Running analysis:\t$command";
  if(system($command)) { throw("FAILED to run $command"); }
  
  return 1;
}

sub write_output {  # Store SWEmbl results
  my $self = shift @_;
  
  #This is now handled in SetupPeaksPipeline...
  #$self->_add_experiment_to_db();

  #TODO Add an option to only process certain slices...
  #$self->_parse_result_file(@slices);

  $self->_parse_result_file();

  #idea of calling the load reads only after peak calling...
  #my $input_file = $self->_input_file();
  #my $nbr_reads = $self->_get_number_of_reads($input_file, $self->_file_type());
  #if($nbr_reads < 1) { throw "$input_file is empty"; }
  #Get relevant slices to avoid creating unnecessary jobs...
  #my @slices = $self->_get_reads_slices($input_file, $self->_file_type());

  #my @rep_out_ids;
  #Create the necessary LoadReads jobs
  #@slices = @{&generate_slices_from_names($slice_adaptor, @slices, undef, 1, 1, 1)};#toplevel, nonref, incdups
  #foreach my $slice (@slices){
  #  my $new_input_id = eval($self->input_id);
  #  $new_input_id->{"slice"} = $slice;
  #  $new_input_id->{"nbr_reads"} = $nbr_reads;
  #  push(@rep_out_ids,$new_input_id);
  #}

  #WrapUp... 
  #my ($funnel_job_id) = @{ $self->dataflow_output_id($new_input_id, 3, { -semaphore_count => scalar(@rep_out_ids) })  };
  #All the fanned jobs...
  #my $fan_job_ids = $self->dataflow_output_id(\@rep_out_ids, 2, { -semaphored_job_id => $funnel_job_id } );

  return 1;
}



sub _parse_result_file {
  
  my $self = shift @_;

  ### annotated features
  my $fset = $self->_feature_set();	
  
  my $efgdba = $self->_efgdba();
  my $sa = $efgdba->get_SliceAdaptor();
  
  #Cache slices and features...
  my %slice;
  my @af;
  

  #Parse the output file
  my $results_file = $self->_results_file();
  if(!open(FILE,$results_file)){ throw "Could not open results file : ".$results_file; }
  while(<FILE>){
    chomp;
    #Ignore headers
    if(/^#/ || /^Region/){ next; }
    #Parse SWEmbl output here... not much sense in making a parser, 
    # unless SWEmbl ouput is used elsewhere (or becomes more complex to parse);
    my ($seqid, $start, $end, undef, undef, undef, $score, undef, undef, $summit) = split(/\s+/);

    #This seqid may vary depending on the input given to SWEmbl... 
    # handle it manually at least for the moment... namely the sam seqid...
    #Make sure to test thoroughly to see if it works...
    #e.g. chromosome:GRCh37:15:1:102531392:1
    if($seqid =~ /^\S*:\S*:(\S+):\S+:\S+:\S/) { $seqid = $1; }
    #In case UCSC input is used... carefull names may not match with ensembl db!
    $seqid =~ s/^chr//i;

    if($self->param('slice') && ($seqid ne $self->param('slice'))){
      warn "Feature being ignored as it is not in specified slice ".$self->param('slice')." : Region:".
	$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
      next;
    }

    #Just in case some of the non-aligned were missed in a filtering step... though this seriously affects peak calling
    #if($seqid =~ /\*/){ next; }
   
    # skip mitochondria calls - remove this when we have pre-processing step to filter alignments using blacklist?
    #next if ($seqid =~ m/^M/);

    # filtering is done as a post-processing e.g. FilterBlacklist.pm
    $summit = int($summit);#Round up?	  
  
    unless (exists $slice{"$seqid"}) {
      $slice{"$seqid"} = $sa->fetch_by_region(undef, $seqid);
    }

    #Do not enter peaks that are invalid ENSEMBL slices...
    #See if this slows down the process significantly...
    #May do that as post-processing?? Cannot cache since each feature should be independent...
    
    #This step seems irrelevant as negative coordinates are still passing and errors are likely in further steps...
    #my $feature_slice = $sa->fetch_by_region(undef, $seqid, $start, $end);
    #if(!$slice{"$seqid"} || !$feature_slice){

    #Sometimes there are naming issues with the slices... e.g. special contigs... which are not "valid" slices in ENSEMBL
    #if(!$slice{"$seqid"} || !(($start =~ /^-?\d+$/) && ($end =~ /^\d+$/))){
  
    #  warn "Feature being ignored due to incorrect coordinates: Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
    #  next;
    #}

    if( ($start < 1) || ($end > $slice{"$seqid"}->end)){
      warn "Feature being ignored due to coordinates out of slice: Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
      next;
      }
  
    #Gracefully handle errors...
    my $af;
    eval{
      $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
    	(
	 -slice         => $slice{"$seqid"},
	 -start         => $start,
	 -end           => $end,
	 -strand        => 0,
	 -score         => $score,
	 -summit        => $summit,
	 -feature_set   => $fset,
	);
    }; 
    if($@) { warn($@); next; }
    if(!$af) { warn("Could not create feature - Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit); next; }
    
    #if($self->param('slice') && ($af->seq_region_name ne $self->param('slice'))){
    #  warn "Feature being ignored as it is not in specified slice ".$self->param('slice')." : Region:".
    #	$af->seq_region_name." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
    #  next;
    #}

    push(@af, $af);
  } 
  close FILE;
 
  # Batch store features...
  $efgdba->get_AnnotatedFeatureAdaptor->store(@af);
  
  #Do this on a wrapup runnable...so it will only be visible after reads are loaded...
  $fset->adaptor->set_imported_states_by_Set($fset);

}


#Private getter / setter to the feature set
sub _feature_set {
  return $_[0]->_getter_setter('feature_set',$_[1]);
}

#Private getter / setter to the command to be executed
sub _command {
  return $_[0]->_getter_setter('command',$_[1]);
}

#Private getter / setter to the command options 
sub _command_options {
  return $_[0]->_getter_setter('command_options',$_[1]);
}

#Private getter / setter to the results file 
sub _results_file {
  return $_[0]->_getter_setter('results_file',$_[1]);
}

1;

