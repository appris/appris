=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Alignment;

=head1 DESCRIPTION

'Alignment' Is a base class for other classes dealing with Alignment
It contains virtually nothing so it may disappear and just pass to Funcgen

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Alignment;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Funcgen');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;

sub fetch_input {   
  my $self = shift @_;

  $self->param("file_type","sam");

  $self->SUPER::fetch_input();
  
  #Just override input folders... maybe consider overriding output folders too?
  my $input_dir = $self->_work_dir()."/fastq/".$self->_species()."/".$self->_experiment_name()."/".
    $self->_cell_type()->name."_".$self->_feature_type()->name;
  $self->_input_dir($input_dir);

  #Folder where to send the final results...
  my $repository = $self->_work_dir()."/alignments/".$self->_species."/".$self->_assembly()."/".$self->_experiment_name();
  $self->_repository($repository);
  
  return 1;
}

#Private getter / setter to the input folder
sub _input_dir {
  return $_[0]->_getter_setter('input_dir',$_[1]);
}

#Private getter / setter to the repository folder
sub _repository {
  return $_[0]->_getter_setter('repository',$_[1]);
}

#Private getter / setter to an input file
sub _input_file {
  return $_[0]->_getter_setter('input_file',$_[1]);
}


1;

