
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

Bio::EnsEMBL::Funcgen::Hive::RunIDR

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunIDR;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref ); 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects 
                                                    run_system_cmd 
                                                    run_backtick_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( run_IDR );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

### CURRENT ISSUES 
#
# 1 This is currently tied to the SWEmbl format
#
# 2 SWEmbl has not been tested for IDing large numbers of peaks, so peak filtering may have issues like MACS2
#
# 3 IDR analysis is currently not stored in the DB. Where are we going to store the threshold and num peaks used for the IDR?
#   Should probably go in the FeatureSet.description rather than different analyses, or just omit for now? 
#IDR analysis should probably be specified as a default/pipeline_wide and batch flowable logic_name
#This will need to be defined in BaseSequenceAnalysis as it is needed by several confs
#This should be merged with the SWEMBL analysis e.g. 
  

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
 
  my $rset_ids = $self->get_param_method('bed_files',  'required');
  assert_ref($rset_ids, 'ARRAY', 'ResultSet dbIDs');
  

  $self->get_param_method('idr_threshold', 'required');
  
  #One of these is required, why do we have 2?
  #It's easier to generate this in the caller
  #Or do we just use the names to generate the output prefix?
  #and drop idr_name completely?
  #How do we get what matches between two string, then do a vs on the rest?
  
  $self->get_param_method('output_prefix', 'required');
  $self->get_param_method('batch_name',    'required');
  $self->param_required('accu_idx');
  $self->get_output_work_dir_methods;

  #temp solution to adjust maxpeaks for poor reps
  $self->get_param_method('bam_files', 'silent');  
  assert_ref($self->bam_files, 'ARRAY', 'bam_files') if $self->bam_files;
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self          = shift;
  
  #eval this and throw no retry
  my $num_peaks;
  
  if(! eval { $num_peaks = run_IDR(-out_dir       => $self->output_dir, 
                                   -output_prefix => $self->output_prefix, 
                                   -threshold     => $self->idr_threshold, 
                                   -bed_files     => $self->bed_files,
                                   -batch_name    => $self->batch_name, 
                                   -bam_files     => $self->bam_files); 1 }){
    $self->throw_no_retry("Failed to run_IDR ".$self->output_prefix."\n$@");                                   
  }
  
  $self->set_param_method('num_peaks', $num_peaks);
  
  #This probably needs writing to the accu table too! 
  #So the Postprocess jobs can pick up the min $num_peaks without reading the file

  #This is old 2 rep only method which sets a new threshold for a subsampled SWEmbl run 
  # as opposed to counting the number of peaks as per the IDR docs   
  #Compute cutoff value (min score) for IDR = 1%
  #my $factor = run_backtick_cmd("awk 'NF > 10 && $11 < 0.01 && (min == 0 || min > $5) {min = $5} END {print min}' ${out_dir}/${exp_name}-overlapped-peaks.txt");
  # This is when I realize that starting off in an actual scripting language would have been better... oh well
  #LOW_CUTOFF='0.0001'
  #IDR_CUTOFF='0.01'
  #new_cutoff="$(echo "$LOW_CUTOFF*$factor" | bc)"
  #my $low_cutoff = 0.0001;
  #my $new_cutoff = $factor * 
  # echo Subsample
  #samtools view -s 0.5 -b $bam1 > $sub_bam1
  #samtools view -s 0.5 -b $bam2 > $sub_bam2
  #samtools merge -f $merged_bam $sub_bam1 $sub_bam2
  #rm -f $sub_bam1 $sub_bam2
  #echo Final run with new -R param set to $new_cutoff
  #SWEMBL -F -V -i $merged_bam -f 150 -R $new_cutoff -o $output -r $control 
   


  #TO DO Capture IDR QC output and threshold
  #accu num_peaks depedant on pass/fail
  #WHere are the QC metrics????
  #em.sav and iri.sav are both internally gzip compressed but have no .(t)gz suffix
  #Its unclear what exactly can read these formats as zless does not
  

    
    
  return;
}



sub write_output {  # Create the relevant jobs
  my $self = shift;
  #This is insert into an accu array based on the value of accu_idx
  #This is picked up by PostprocessIDR analysis on another branch.
  $self->dataflow_output_id( {'idr_peak_counts'   => $self->num_peaks}, 2);
  return;
}

1;