
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( post_process_IDR );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects
                                                    get_set_prefix_from_Set
                                                    run_system_cmd );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

#TODO 
# 1 ADD IDR based QC
# 2 Build output id

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  
  #From SubmitIDR
  #{%{$batch_params},
  #                            dbIDs    => \@rset_ids,
  #                            set_type => 'ResultSet'}

  #Can either get the accu idr peaks or read from file
  #based on the 
    
  #How are we going to differentiate pseudo reps?
  #SubmitIDR will know about these,as they will be passed explicitly from GeneratePseudoReplicates
  #as file names (not sets!?).
  $self->get_output_work_dir_methods;
  
  if($self->param_required('set_type') ne 'ResultSet'){
    throw("Param set_type should be 'ResultSet', not:\t".$self->param('set_type'));  
  }
  
  my $rset_ids = $self->get_param_method('dbIDs',  'required');
  assert_ref($rset_ids, 'ARRAY', 'ResultSet dbIDs');
  $self->get_param_method('permissive_peaks', 'required');
  $self->get_param_method('idr_peak_counts',  'required');
  $self->get_param_method('output_prefixes',  'required');
  $self->get_param_method('batch_name',       'required');
  assert_ref($self->output_prefixes, 'ARRAY', 'output_prefixes');
  
  #This is accumulated data from the RunIDR fan jobs, submitted & semaphored from PreprocessIDR
  my $idr_peaks = $self->get_param_method('idr_peak_counts', 'required');   
  assert_ref($idr_peaks, 'ARRAY', 'IDR peaks');
 
  return;
}



#To do
#1 Move most of this to SeqTools
#2 Add an IDR peak analysis, and record the max_peaks in the feature_set description? or tracking?

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self         = shift;  
  my $out_dir      = $self->output_dir;
  my $batch_params = $self->batch_params;
  my $batch_name   = $self->batch_name;
  my $rsets        = &scalars_to_objects($self->out_db,  'ResultSet',
                                         'fetch_by_dbID', $self->dbIDs);
                                          
  #Need to flow peak_analysis here, such that we over-ride the no IDR mode defaults
  #Otherwise DefineMergedDataSet will use the pre-IDR analysis.????
  #Is this correct? No, we always want to run permissive, so this entirely depends on no_idr
  my $max_peaks;
  
  if(! eval { $max_peaks = post_process_IDR($out_dir, $batch_name, $self->idr_peak_counts,
                                   {'-idr_files'     => $self->output_prefixes,
                                    '-debug'         => $self->debug}); 1 }){
    $self->throw_no_retry("Failed to post_process_IDR $batch_name\n$@");                                    
  }
  
  #Pooled file should already exist from Pseudo-rep creation
  #but handle absence, just incase we want to run without this step
 
  
  #Do the alignment merge here as we already have access to the rep ResultSets here?
  #or do it in 
  #controls should be same across all fsets
  my @controls   = grep {$_->is_control}   @{$rsets->[0]->get_support};
  my @signals    = grep {! $_->is_control} (map {@{$_->get_support}} @$rsets);  
  my $set_prefix = get_set_prefix_from_Set($rsets->[0]); 
  my %controls   = ();
  
  if(@controls){
    #Only flow controls if they are defined
    #mainly for tidyness, but just incase downstream analysis
    #doesn't handle emtpy controls
    $controls{controls} =  [ map {$_->dbID } @controls];  
  }
  
  #This is flowing to DefineMergedReplicateResultSet
  $self->branch_job_group(2, [{%{$batch_params},
                              max_peaks                  => $max_peaks, #This is batch flown
                              #peak_analysis              => $peak_analysis->logic_name,
                              merge_idr_replicates       => 1,
                              alignment_analysis         => $rsets->[0]->analysis->logic_name,
                             
                              input_subset_ids => 
                               {
                                %controls,
                                $set_prefix  => [ map {$_->dbID } @signals], }}]);
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;