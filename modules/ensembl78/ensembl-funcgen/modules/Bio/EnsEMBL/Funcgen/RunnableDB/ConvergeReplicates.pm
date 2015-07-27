=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::ConvergeReplicates

=head1 DESCRIPTION

'ConvergeReplicates' Just merges the replicates BAM and convert it to a SAM
into the alignments folder for the peaks pipeline...

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::ConvergeReplicates;

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
#use Data::Dumper;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Alignment');

#TODO... Maybe use and update the tracking database...
sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  
  $self->SUPER::fetch_input();

  my $species = $self->_species();
  my $gender = $self->_cell_type()->gender;
  $gender = $gender ? $gender : "male";
  my $assembly = $self->_assembly();

  my $sam_header = $self->_work_dir()."/sam_header/".$species."/".$species."_".$gender."_".$assembly."_unmasked.header.sam";
  $self->_sam_header($sam_header);

  my $repository = $self->_repository();
  if(! -d $repository){ 
    system("mkdir -p $repository") && throw("Couldn't create directory $repository");
  }

  #my $nbr_replicates = $self->param('nbr_replicates') || throw "Number of replicates not given";
  my $nbr_replicates = $self->param('nbr_replicates') || 2;
  $self->_nbr_replicates($nbr_replicates);

  return 1;
}

sub run {   
  my $self = shift @_;

  my $sam_header = $self->_sam_header();

  my $output_file_prefix = $self->_repository()."/".$self->_set_name();
  my $file_prefix = $self->_output_dir()."/".$self->_set_name();
  my $merge_cmd;
  if($self->_nbr_replicates()>1){
    $merge_cmd="samtools merge -h $sam_header ${file_prefix}.bam ${file_prefix}.[0-9]*.sorted.bam "; 
  }  else {
    $merge_cmd = "cp ${file_prefix}.[0-9]*.sorted.bam ${file_prefix}.bam";    
  }
  if(system($merge_cmd) != 0){ throw "Error merging replicate bam files: $merge_cmd"; }
  

  my $convert_cmd = "samtools view -h ${file_prefix}.bam > ${output_file_prefix}.samse.sam";
  if(system($convert_cmd) != 0){ throw "Error converting merged bam to sam: $convert_cmd"; }

  my $zip_cmd = "gzip ${output_file_prefix}.samse.sam";
  if(system($zip_cmd) != 0){ warn "Error zipping sam"; }

  #Merge the alignment logs too...
  my $alignment_log = $output_file_prefix.".alignment.log";
  my $log_cmd="echo \"Alignment QC - total reads as input: \" >> ${alignment_log}";
  $log_cmd="${log_cmd};samtools flagstat ${file_prefix}.bam | head -n 1 >> ${alignment_log}";
  $log_cmd="${log_cmd}; echo \"Alignment QC - mapped reads: \" >> ${alignment_log} ";
  $log_cmd="${log_cmd};samtools view -u -F 4 ${file_prefix}.bam | samtools flagstat - | head -n 1 >> ${alignment_log}";
  $log_cmd="${log_cmd}; echo \"Alignment QC - reliably aligned reads (mapping quality >= 1): \" >> ${alignment_log}";
  $log_cmd="${log_cmd};samtools view -u -F 4 -q 1 ${file_prefix}.bam | samtools flagstat - | head -n 1 >> ${alignment_log}";
    
  if(system($log_cmd) != 0){ warn "Error making the alignment statistics"; }

  my $rm_cmd="rm -f ${file_prefix}.bam";
  if(system($rm_cmd) != 0){ warn "Error removing temp files. Remove them manually: $rm_cmd"; }

  return 1;
}


sub write_output {  # Nothing to do here
  my $self = shift @_;

  return 1;

}

#Private getter / setter to the sam header
sub _sam_header {
  return $_[0]->_getter_setter('sam_header',$_[1]);
}

#Private getter / setter to the sam header
sub _nbr_replicates {
  return $_[0]->_getter_setter('nbr_replicates',$_[1]);
}

1;
