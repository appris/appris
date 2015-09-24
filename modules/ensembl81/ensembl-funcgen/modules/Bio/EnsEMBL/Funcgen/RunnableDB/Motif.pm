=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Motif

=head1 DESCRIPTION

'Funcgen::Motif' 

TODO: Try to collate that to Funcgen

=cut


package Bio::EnsEMBL::Funcgen::RunnableDB::Motif;

use warnings;
use strict;

use base ('Bio::EnsEMBL::Hive::Process');

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;

sub fetch_input {   # nothing to fetch... just the parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();
  
  if(!$self->param('bin_dir')){ throw "No binary folder given"; }
  $self->_bin_folder($self->param('bin_dir'));

  if(!$self->param('species')){ throw "No species given"; }
  $self->_species($self->param('species'));

  if(!$self->param('bin_size')){ throw "No bin size given"; }
  $self->_bin_size($self->param('bin_size'));

  if(!$self->param('window_size')){ throw "No window size given"; }
  $self->_window_size($self->param('window_size'));

  #Get the core db... possibly make these non-mandatory
  if(!$self->param('dnadb_host')){ throw "No core host given"; }
  if(!$self->param('dnadb_port')){ throw "No core port given"; }
  if(!$self->param('dnadb_user')){ throw "No core user given"; }
  #if(!$self->param('dnadb_name')){ throw "No core dbname given"; }
  my $dba;
  eval{
       $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
						  -host => $self->param('dnadb_host'),
						  -port => $self->param('dnadb_port'),
						  -user => $self->param('dnadb_user'),
						  -dbname => $self->param('dnadb_name'),
						  -species => $self->_species,
						 );
  };
  if($@) { throw "Error creating the Core DB Adaptor: $@";  }    
  if(!$dba){ throw "Could not connect to Core DB"; }

  #Get the efg db... always read only user
  if(!$self->param('dbhost')){ throw "No host given"; }
  if(!$self->param('dbport')){ throw "No port given"; }
  if(!$self->param('dbuser')){ throw "No user given"; }
  if(!$self->param('dbname')){ throw "No dbname given"; }
  eval{
       $self->_efgdba(Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
								   -host => $self->param('dbhost'),
								   -port => $self->param('dbport'),
								   -user => $self->param('dbuser'),
								   -dbname => $self->param('dbname'),
								   -species => $self->_species,
								   -dnadb => $dba,
								  ));
  };
  if($@) { throw "Error creating the EFG DB Adaptor: $@";  }    
  if(!$self->_efgdba()){ throw "Could not connect to EFG DB"; }

  my $dnadbc = $self->_efgdba()->dnadb->dbc;
  warn $dnadbc->host." ".$dnadbc->port." ".$dnadbc->username." ".$dnadbc->dbname;


  if(!$self->param('feature_set')){ throw "No feature set given"; }
  my $fseta = $self->_efgdba()->get_FeatureSetAdaptor();
  my $fset = $fseta->fetch_by_name($self->param('feature_set'));
  if(!$fset){
    throw $self->param('feature_set')." is not a valid Feature Set";
  }
  $self->_feature_set($fset);

  if(!$self->param('output_dir')){ throw "No output dir given"; }
  $self->_output_dir($self->param('output_dir')."/".$self->_feature_set->name);

  return 1;
}

sub run {   
  my $self = shift @_;
  
  return 1;
}


sub write_output {  # Nothing is written at this stage (for the moment)

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

#Private getter / setter to the bin folder
sub _bin_folder {
  return $_[0]->_getter_setter('bin_folder',$_[1]);
}

#Private getter / setter to the EFG DB Adaptor
sub _efgdba {
  return $_[0]->_getter_setter('efgdb',$_[1]);
}

#Private getter / setter to the bin size
sub _bin_size {
  return $_[0]->_getter_setter('bin_size',$_[1]);
}

#Private getter / setter to the window size
sub _window_size {
  return $_[0]->_getter_setter('window_size',$_[1]);
}

#Private getter / setter to the output dir
sub _output_dir {
  return $_[0]->_getter_setter('output_dir',$_[1]);
}

#Private getter / setter to the feature set
sub _feature_set {
  return $_[0]->_getter_setter('feature_set',$_[1]);
}

#Private getter / setter to the species
sub _species {
  return $_[0]->_getter_setter('species',$_[1]);
}

1;
