=pod

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::SetupAlignmentPipeline

=head1 DESCRIPTION

'SetupAlignmentPipeline' Does all the setup before the Alignment is run
Checks for existence of input files, etc...
This Runnable CAN be run multiple times in parallell!

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::SetupAlignmentPipeline;

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(is_gzipped);
use Data::Dumper;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Alignment');

#TODO... Maybe use and update the tracking database...
sub fetch_input {   # fetch parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();

  #Magic default number...
  my $fastq_chunk_size = 8000000;
  if($self->param("fastq_chunk_size")){ $fastq_chunk_size = $self->param("fastq_chunk_size")};
  $self->_fastq_chunk_size($fastq_chunk_size);

  #Sets up the output dir for this experiment_name
  my $output_dir = $self->_output_dir();
  if(! -d $output_dir){
    system("mkdir -p $output_dir") && throw("Couldn't create output directory $output_dir");
  }

  my $input_dir = $self->_input_dir();
  if(! -d $input_dir ){ throw " Couldn't find input directory $input_dir"; }

  opendir(DIR,$input_dir);
  my @dirs = grep(/^\d/,readdir(DIR));
  closedir(DIR);

  if(scalar(@dirs)==0){ throw "No replicates found in $input_dir"; }

  my @input_files;
  my @replicates;
  foreach my $dir (@dirs){
    #TODO: maybe use some other code for replicates? (e.g. Rep\d )
    if($dir =~ /^(\d)$/){
      my $replicate = $1;

      opendir(DIR,$input_dir."/".$replicate);
      my @files = grep(/.fastq/,readdir(DIR));
      closedir(DIR);

      if(scalar(@files)==0){ throw "No files for replicate $replicate"; }

    my $file_count = 0;
    for my $file (@files){
    push @input_files, {
      path => $input_dir."/".$replicate."/".$file,
      replicate => $replicate,
      file_index => $file_count++,
    };
      }

    push @replicates, $replicate;
    } else { warn "Invalid replicate $dir ignored";   }
  }

  $self->_input_files(\@input_files);
  $self->_replicates(\@replicates);

  return 1;
}

sub run {
  my $self = shift @_;

  my $fastq_chunk_size = $self->_fastq_chunk_size();

  my @output_ids;
  my $set_name = $self->_set_name();

  foreach my $file_info (@{$self->_input_files()}){
  my $file = $file_info->{'path'};
  my $replicate = $file_info->{'replicate'};
  my $file_index = $file_info->{'file_index'};

  my $cmd;

  if($file =~ /^(.*.fastq).gz$/){
      $cmd = "gunzip -c";
  }
  elsif($file =~ /^(.*.fastq).bz2$/){
      $cmd = "bunzip2 -c"
  }
  else {
    $cmd = "cat";
  }

  $cmd .= ' '.$file.' | split -d -a 4 -l '.$fastq_chunk_size.' - '. $self->_output_dir().'/'.$set_name."_".$replicate.'_'.$file_index.'_';

    if(system($cmd) != 0){ throw "Problems running $cmd";  }
  }


  return 1;
}


sub write_output {  # Create the relevant job
  my $self = shift @_;

  my $set_name = $self->_set_name;

  my (@align_output_ids, @merge_output_ids);

  opendir(DIR,$self->_output_dir());
  for my $split_file ( grep(/^${set_name}_\d+_\d+_\d+$/,readdir(DIR)) ){
  my $output = eval($self->input_id);
  $output->{input_file} = $split_file;
  push @align_output_ids, $output;
  }
  closedir(DIR);

  # merge data for each replicate

  for my $rep (@{$self->_replicates}){
  my $output = eval($self->input_id);
  $output->{replicate} = $rep;
  push @merge_output_ids, $output;
  }


  # files to align
  $self->dataflow_output_id(\@align_output_ids, 1);

  # merge data acros replicates
  $self->dataflow_output_id($self->input_id, 2);#input_id
  return 1;

}

#Private getter / setter to the fastq chunk size
sub _fastq_chunk_size {
  return $_[0]->_getter_setter('fastq_chunk_size',$_[1]);
}

#Private getter / setter to the output_ids list
sub _output_ids {
  return $_[0]->_getter_setter('output_ids',$_[1]);
}

#Private getter / setter to the output_ids list
sub _replicates {
  return $_[0]->_getter_setter('replicates',$_[1]);
}

sub _input_files {
  return $_[0]->_getter_setter('input_files',$_[1]);
}

1;
