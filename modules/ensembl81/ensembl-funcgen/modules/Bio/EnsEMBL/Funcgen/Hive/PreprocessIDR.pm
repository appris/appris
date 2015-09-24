
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR

=head1 DESCRIPTION

This is a very simple module which calculates the IDR threshold and number of peaks to filter before
submitting the individual IDR jobs. This enables the correct semaphore config (fan/funnel), which is currently 
impossible from the previous analyses.

# This could be omitted in multi-sempahores are supported e.g.
#    '2->A'    => ['DefineReplicateDataSet'], #fan
#    'A->3->B' => ['RunIDR'],                 #fan and immediate funnel
#    'B->4'    => ['PostProcessIDRReplicates],

=cut


#TODO
# 1 Enable theshold override
# 2 Enable dropping ResultSets (alter job id from env). Need to list rset dbids in output?
# 3 Fan filter/sort jobs, this will require a dummy analysis and maybe a job factory
# 4 Pick up peak QC failures from status table(no FeatureSet)? Or accu?
#   Or tracking table? Accu is probably enough here
#   But we do need to track the QC for individual reps, so we can mark them as dodgy   
# 5 Currently hardcoded for SWEmbl output with txt file suffix.

package Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( pre_process_IDR );                                          
                                               
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->check_analysis_can_run;
  $self->SUPER::fetch_input;
  
  if($self->param_required('set_type') ne 'ResultSet'){
    throw("Param set_type should be 'ResultSet', not:\t".$self->param('set_type'));  
  }
  
  my $rset_ids = $self->get_param_method('dbIDs',  'required');
  assert_ref($rset_ids, 'ARRAY', 'ResultSet dbIDs');
  
  
  #warn "pre fudge ids are ".join(' ', @$rset_ids);
  
  $rset_ids = [keys %{{ map { $_ => 1 } @$rset_ids }}];
   
  #warn "post fudge ids are ".join(' ', @$rset_ids);
  
  #Temporary bug fix
  #This accounts for duplicate result_set_ids which 
  #which are generated in DefineResultSet, due to redundant 
  #replicate numbers
  #We need to delete all downstream jobs for the affected PreprocessIDR jobs
  #also need to put a sanity check in DefineResultSets, or Identify, which
  #catches any redundant signal rep numbers!
  $self->dbIDs($rset_ids);
   
  $self->set_param_method('peak_analysis', 
                          &scalars_to_objects($self->out_db, 'Analysis',
                                              'fetch_by_logic_name',
                                              $self->param_required('permissive_peaks'))->[0]);
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self  = shift;  
  my $rsets = &scalars_to_objects($self->out_db, 'ResultSet',
                                                 'fetch_by_dbID',
                                                 $self->dbIDs);
  my $peak_analysis = $self->peak_analysis;                                                 
  my $batch_params  = $self->batch_params;
  my $exp_name      = $rsets->[0]->experiment->name;
  #This is also done in RunPeaks, so we really need a single method to do this?
  my $lname         =  $peak_analysis->logic_name;
  $self->get_output_work_dir_methods($self->db_output_dir.'/peaks/'.$exp_name.'/'.$lname);
  my $out_dir = $self->output_dir;
  my $max_peaks     = 300000;
 
  # Validate analysis  
  if($lname !~ /swembl/io){
    #$self->input_job->transient_error( 0 );
    $self->throw_no_retry('Pre-IDR peak filtering and reformating is currently only optimised for the SWEmbl analysis');
    #This maybe fine, but not for MACS2, see IDR docs here 
    #https://sites.google.com/site/anshulkundaje/projects/idr#TOC-CALL-PEAKS-ON-INDIVIDUAL-REPLICATES
  }
  elsif($lname =~ /macs/io){#future proofing as will currently never be tested
    $max_peaks = 100000;
    warn "Reseting max filtered peaks value to 100000 for MACS analysis:\t$lname\n";   
  } 
      
      
  my (@pre_idr_files, @rep_nums, $ctrl_ids, @bams);
  
  #sub _pre_process_IDR_ResultSets?
  
  foreach my $rset(@$rsets){
 
    if($exp_name ne $rset->experiment->name){
      $self->throw_no_retry("IDR Replicate ResultSets are not from the same experiment:\t".
        $rsets->[0]->name.' '.$rset->name);  
    }
    
    my $seen_rep = 0;
    
    if($rset->table_name ne 'input_subset'){
      $self->throw_no_retry('Expected input_subset support but found '.
        $rset->table_name." support for ResultSet:\t".$rset->name);  
    }
    
    my @issets = @{$rset->get_support};
    my @ctrl_ids;
    
    foreach my $isset(@issets){
      
      if($isset->is_control){
        push @ctrl_ids, $isset->dbID;;        
      }
      else{
        if($seen_rep){
          $self->throw_no_retry("Found more than 1 replicate (non-control) InputSet supporting an an IDR ResultSet:\n\t".
            join("\n\t", map {$_->name} @issets)."\n");  
        }  
        
        push @rep_nums, $isset->replicate;
        $seen_rep = 1;
      }
    }
    
    if(! defined $ctrl_ids){
      $ctrl_ids = join(' ', (sort {$a <=> $b} @ctrl_ids));  
    }
    
    my $current_ctrl_ids = join(' ', (sort {$a <=> $b} @ctrl_ids));  
    
    if($current_ctrl_ids ne $ctrl_ids){
      $self->throw_no_retry("Found mismatched controls between IDR ResultSets:\n\t".
        join("\t", (map {$_->name} @$rsets)));  
    }
    
    #todo validate ResultSet analysis here too?
    #bed file is currently defined by PeakCaller::init_files
    #based on the out_dir, outfile_prefix and the output_format
    #We should probably pass this whole path, such that we can centralise how the path is generated?
    #If we let the default PeakCaller formats be used, we can't know what the suffix will be at this point
    
    #There is currently an issues with set nameing, as this will integrate the alignment analysis
    #too. Although this maybe desirable to avoid clashes between features sets with different alignments
    #The API does not handle this yet.
    
    push @pre_idr_files, $out_dir.'/'.$rset->name.'.'.$lname.'.txt';
    #do read counts in RunIDR to parallelise
    push @bams,         $self->get_alignment_files_by_ResultSet_formats($rset, ['bam'])->{bam};
  }
  
  
  #Put batch_name code in BaseSequencing or Base? 
  my $batch_name                 = $exp_name.'_'.$lname.'_'.join('_', sort @rep_nums);
  my ($np_files, $threshold, $x_thresh_adjust);
  
  if(! eval { ($np_files, $threshold, $x_thresh_adjust) = 
                pre_process_IDR($out_dir, \@pre_idr_files, $batch_name, $max_peaks); 1}){
    $self->throw_no_retry("Failed to pre_process_IDR $batch_name\n$@");                
  } 
  
  # end _pre_process_IDR_ResultSets

  #Build 2 way rep combinations for IDR jobs
  my @idr_job_ids;
  my @out_prefixes;
  my $last_i = $#{$np_files} - 1;
  
  foreach my $i(0..$last_i){
    my $next_i = $i + 1;
  
    #2 way jobs 
    foreach my $j($next_i..$#{$np_files}){
      my @names;
    
      foreach my $file($np_files->[$i], $np_files->[$j]){
        (my $name = $file) =~ s/.*\///;  
        $name =~ s/(?:\.np_idr)\.txt$//;
        push @names, $name;
      }
      
      my $output_prefix = $names[0].'_VS_'.$names[1];
      push @out_prefixes, $output_prefix;
      
      my $job = {%$batch_params,
                 output_dir    => $out_dir,
                 idr_threshold => $threshold,
                 accu_idx      => $i,
                 output_prefix => $output_prefix,
                 batch_name    => $batch_name,
                 bed_files     => [$np_files->[$i], $np_files->[$j]]};
      
      if($x_thresh_adjust){
        $job->{bam_files} = [$bams[$i], $bams[$j]];  
      }
                                
      push @idr_job_ids, $job;
    }  
  }
  
  #Now we need to pool and produce pseudo reps? This should be done way before here?!!
  #In between MergeControlAlignments_and_QC and Submit_IDR
  $self->branch_job_group(2, \@idr_job_ids,
                          3, [{%$batch_params,
                               output_dir      => $out_dir,
                               batch_name    => $batch_name,
                               output_prefixes => \@out_prefixes,
                               dbIDs           => $self->dbIDs,
                               set_type        => 'ResultSet'}]);
      
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;