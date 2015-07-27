=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::RunBWA

=head1 DESCRIPTION

'RunBWA' runs BWA with an input file

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::RunBWA;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Alignment');

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  
  $self->SUPER::fetch_input();
  $self->_input_dir($self->_output_dir());

  my $input_file = $self->param('input_file') || throw "No input file given";
  $input_file = $self->_input_dir()."/".$input_file;
  if(! -e $input_file){ throw " Couldn't find $input_file"; }
  $self->_input_file($input_file);

  my $species = $self->_species();
  my $gender = $self->_cell_type()->gender();
  $gender = $gender ? $gender : "male";
  my $assembly = $self->_assembly();
  #TODO: Maybe check if the index file is really there? Eventually the bwa output will tell you though
  my $bwa_index = $self->_work_dir()."/bwa_indexes/".$species."/".$species."_".$gender."_".$assembly."_unmasked.fasta";
  $self->_bwa_index($bwa_index);

  my $bwa_bin = $self->param('bwa_bin') || $self->_bin_dir().'/bwa';
  $self->_bwa_bin($bwa_bin);

  return 1;
}

sub run {   
  my $self = shift @_;

  my $input_file = $self->_input_file();
  my $bwa_index = $self->_bwa_index();

  my $bwa_bin = $self->_bwa_bin();

  #TODO Pass the location of the binary to be sure we'll be running the right version?
#  my $bwa_cmd = "$bwa_bin aln $bwa_index $input_file";
  #Allow this to work with paired reads?? Maybe not for the moment...
  #in that case pass bwa algorithm as parameter...
  #If using -q make sure we've got the correct Sanger quality scores...
#  $bwa_cmd .= " | $bwa_bin samse $bwa_index - $input_file";
#  $bwa_cmd .= " | samtools view -uS - ";
#  $bwa_cmd .= " | samtools sort - ${input_file}.sorted";
#  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }

  #Silent errors seem to be passing... running one command at a time!?
  my $bwa_cmd = "$bwa_bin aln $bwa_index $input_file > ${input_file}.aln";
  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }
  $bwa_cmd = "$bwa_bin samse $bwa_index ${input_file}.aln $input_file > ${input_file}.aln.sam";
  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }
  if(system("rm ${input_file}.aln") != 0){ warn "Couldn't remove tmp file. Remove it manually."; }
  $bwa_cmd = $self->_bin_dir()."/samtools view -uS ${input_file}.aln.sam > ${input_file}.aln.bam";
  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }
  if(system("rm ${input_file}.aln.sam") != 0){ warn "Couldn't remove tmp file. Remove it manually."; }
  $bwa_cmd = $self->_bin_dir()."/samtools sort ${input_file}.aln.bam ${input_file}.sorted";
  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }
  if(system("rm ${input_file}.aln.bam") != 0){ warn "Couldn't remove tmp file. Remove it manually."; }

  if(system("rm $input_file") != 0){ warn "Couldn't remove tmp file. Remove it manually."; }

  return 1;
}


sub write_output {  # Nothing to write
  my $self = shift @_;

  return 1;

}


#Private getter / setter to the bwa indexes
sub _bwa_index {
  return $_[0]->_getter_setter('bwa_index',$_[1]);
}

#Private getter / setter to the bwa bin
sub _bwa_bin {
  return $_[0]->_getter_setter('bwa_bin',$_[1]);
}

1;
