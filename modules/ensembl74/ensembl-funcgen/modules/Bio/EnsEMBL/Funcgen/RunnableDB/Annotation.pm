=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Annotation

=head1 DESCRIPTION

'Annotation' is a base class for runnables running the Annotatin Pipeline
It performs common tasks such as connecting to the EFG DB etc...

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Annotation;

use warnings;
use strict;

use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;
use DBI;

#use base ('Bio::EnsEMBL::Hive::ProcessWithParams');
use base ('Bio::EnsEMBL::Hive::Process');

#This defines a set of parameters based on given parameters from the pipeline:
sub fetch_input {   # nothing to fetch... just the DB parameters...
  my $self = shift @_;

  my $dnadb_params = $self->param('dnadb') || throw "No parameters for Core DB";
  eval{  $self->_dnadba(Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{ $dnadb_params })); };
  if($@) { throw "Error creating the Core DB Adaptor: $@";  }    
  if(!$self->_dnadba()){ throw "Could not connect to Core DB"; }
  $self->_dnadb_params($dnadb_params);

  #Get efg connection, otherwise fail..
  my $efgdb_params = $self->param('efgdb') || throw "No parameters for EFG DB";

  eval{
       $self->_efgdba(Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
								   %{ $efgdb_params },
								   #Why is this not working???
								   #-dnadb => $self->_dnadba,
								   -dnadb_user => $self->_dnadba->dbc->username,
								   -dnadb_port => $self->_dnadba->dbc->port,
								   -dnadb_host => $self->_dnadba->dbc->host,
								   -dnadb_dbname => $self->_dnadba->dbc->dbname,
								  ));
  };

  if($@) { throw "Error creating the EFG DB Adaptor: $@";  }    
  if(!$self->_efgdba()){ throw "Could not connect to EFG DB"; }
  $self->_efgdb_params($efgdb_params);

  #Get work db connection, otherwise fail..
  my $workdb_params = $self->param('workdb') || throw "No parameters for Work server";
  $self->_workdb_params($workdb_params);

  my $species = $self->param('species') || throw "No species defined";
  $self->_species($species);

  #Not is use, currently...
  #my $release = $self->param('release') || throw "Release number not specified";
  #$self->_release($release);

  my $work_dir = $self->param('work_dir') || throw "'work_dir' is a required parameter"; 
  $self->_work_dir($work_dir);

  #Work with conventions here too?? work_dir/output/dbname ??
  my $output_dir = $self->param('output_dir') || throw "'output_dir' is a required parameter";
  $self->_output_dir($output_dir);

  return 1;
}


sub run {   
  my $self = shift @_;

  return 1;
}


sub write_output {  
  my $self = shift @_;
  
  return 1;

}


#Private Generic getter and setter
sub _getter_setter {
  my ($self, $param_name, $param_value) = @_;
  if(!$param_name){ return undef; }
  if(!$param_value){ 
    $param_value = $self->param($param_name);   
  } else {
    $self->param($param_name, $param_value);
  }
  return $param_value;
}

# Private getter / setters : Maybe do some validation in some cases...

#Private getter / setter to the Work DB Connection params
sub _workdb_params {
  return $_[0]->_getter_setter('workdb_params',$_[1]);
}

#Private getter / setter to the EFG DB Adaptor
sub _efgdba {
  return $_[0]->_getter_setter('efgdb',$_[1]);
}

#Private getter / setter to the EFG DB Connection params
sub _efgdb_params {
  return $_[0]->_getter_setter('efgdb_params',$_[1]);
}

#Private getter / setter to the Core DB Adaptor
sub _dnadba {
  return $_[0]->_getter_setter('dnadb',$_[1]);
}

#Private getter / setter to the Core DB Connection params
sub _dnadb_params {
  return $_[0]->_getter_setter('dnadb_params',$_[1]);
}

#Private getter / setter to the Species name
sub _species {
  return $_[0]->_getter_setter('species',$_[1]);
}

#Private getter / setter to the work folder
sub _work_dir {
  return $_[0]->_getter_setter('work_dir',$_[1]);
}

#Private getter / setter to the output folder
sub _output_dir {
  return $_[0]->_getter_setter('output_dir',$_[1]);
}

#Private getter / setter to the release number
sub _release {
  return $_[0]->_getter_setter('release',$_[1]);
}


1;
