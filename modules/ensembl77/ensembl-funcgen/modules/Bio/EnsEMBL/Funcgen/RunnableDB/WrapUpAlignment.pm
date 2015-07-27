=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::WrapUpAlignment;

=head1 DESCRIPTION

'WrapUpAlignment' Merges all results from alignment jobs into a single file

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::WrapUpAlignment;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Alignment');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;
use File::Copy;

sub fetch_input {
  my $self = shift @_;

  $self->SUPER::fetch_input();

  #my $replicate = $self->param('replicate') || throw "No replicate given";
  #$self->_replicate($replicate);

  #my $nbr_subfiles = $self->param('nbr_subfiles') || throw "Number of subfiles not given";
  #if($nbr_subfiles<1){ throw "We need at least one subfile. Empty or invalid set given: $nbr_subfiles"; }
  #$self->_nbr_subfiles($nbr_subfiles);

  my $species=$self->_species();
  my $gender=$self->_cell_type->gender();
  $gender = $gender ? $gender : "male";
  my $assembly=$self->_assembly();
  my $sam_header = $self->_work_dir()."/sam_header/".$species."/".$species."_".$gender."_".$assembly."_unmasked.header.sam";
  $self->_sam_header($sam_header);


  my $repository = $self->_repository();
  if(! -d $repository){
    system("mkdir -p $repository") && throw("Couldn't create directory $repository");
  }

  return 1;
}

sub run {
  my $self = shift @_;

 # my $replicate = $self->_replicate();
  my $sam_header = $self->_sam_header();

  my $set_name = $self->_set_name();
  my $file_prefix = $self->_output_dir()."/".$set_name;
  my $input_file_pattern = "${file_prefix}_[0-9]*_[0-9]*_[0-9]*.sorted.bam";
  my $merge_cmd="samtools merge -f -h $sam_header ${file_prefix}.sorted.bam $input_file_pattern";

  #Maybe remove duplicates as default behavior ?.
  #my $merge_cmd="samtools merge -h $sam_header - ${file_prefix}.[0-9]*.sorted.bam | ";
  #$merge_cmd.=" samtools rmdup -s - ${file_prefix}.sorted.bam";

  if(system($merge_cmd) != 0){ throw "Error merging file: $merge_cmd"; }

  #keep end result as bam... it's more compact: only when passing to the peak caller we pass to sam or other format...
  #merge_cmd.=" samtools view -h - | gzip -c > ${file_prefix}${align_type}.sam.gz"

  my $alignment_log = $file_prefix.".alignment.log";

  my $log_cmd="echo \"Alignment QC - total reads as input: \" >> ${alignment_log}";
  $log_cmd.=";samtools flagstat ${file_prefix}.sorted.bam | head -n 1 >> ${alignment_log}";
  $log_cmd.=";echo \"Alignment QC - mapped reads: \" >> ${alignment_log} ";
  $log_cmd.=";samtools view -u -F 4 ${file_prefix}.sorted.bam | samtools flagstat - | head -n 1 >> ${alignment_log}";
  $log_cmd.="; echo \"Alignment QC - reliably aligned reads (mapping quality >= 1): \" >> ${alignment_log}";
  $log_cmd.=";samtools view -u -F 4 -q 1 ${file_prefix}.sorted.bam | samtools flagstat - | head -n 1 >> ${alignment_log}";
  #Maybe do some percentages?

  if(system($log_cmd) != 0){ warn "Error making the alignment statistics"; }
  print STDERR "logged\n";
  my $repository = $self->_repository();

  move("${file_prefix}.sorted.bam","${repository}/${set_name}.samse.bam");

  my $convert_cmd =  "samtools view -h ${repository}/${set_name}.samse.bam | gzip -c - > ${repository}/${set_name}.samse.sam.gz";


  #print STDERR "using convert cmd: ".$convert_cmd;

  if(system($convert_cmd) != 0){ warn "Error converting to zipped sam"; }
  #print STDERR "converted\n";
  my $rm_cmd="rm -f $input_file_pattern";
  if(system($rm_cmd) != 0){ warn "Error removing temp files. Remove them manually"; }
  #$rm_cmd="rm -f ${file_prefix}.sorted.bam";
  #if(system($rm_cmd) != 0){ warn "Error removing final bam files. Remove it manually"; }

  return 1;
}


sub write_output {
  my $self = shift @_;

  return 1;

}

#Private getter / setter to the sam header
sub _sam_header {
  return $_[0]->_getter_setter('sam_header',$_[1]);
}

#Private getter / setter to the replicate number
sub _replicate {
  return $_[0]->_getter_setter('replicate',$_[1]);
}

#Private getter / setter to the replicate number
sub _nbr_subfiles {
  return $_[0]->_getter_setter('nbr_subfiles',$_[1]);
}

1;
