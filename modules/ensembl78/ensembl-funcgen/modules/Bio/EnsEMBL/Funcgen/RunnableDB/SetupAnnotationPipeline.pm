=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::SetupAnnotationPipeline

=head1 DESCRIPTION

'SetupAnnotationPipeline' Checks cell types and creates annotation processes for each
This Runnable CANNOT be run multiple times in parallell!

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::SetupAnnotationPipeline;

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Annotation');

sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  
  $self->SUPER::fetch_input();

  #Sets up the output dir 
  my $output_dir = $self->_output_dir();
  if(! -d $output_dir){ 
    system("mkdir -p $output_dir") && throw("Couldn't create output directory $output_dir");
  }

  return 1;
}

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;

  my $efgdba = $self->_efgdba();
  #Check how many different cell types exist with regulatory features for current species
  #Creates the appropriate jobs
  my @reg_sets = @{$efgdba->get_FeatureSetAdaptor->fetch_all_by_type('regulatory')};
  my @cell_types;
  foreach my $set (@reg_sets){
    if($set->cell_type->name ne 'MultiCell'){
      push @cell_types, $set->cell_type->name;
    }
  }
  $self->_cell_types_to_run(\@cell_types);

  return 1;
}


sub write_output {  # Create the relevant job
  my $self = shift @_;

  foreach my $cell_type (@{$self->_cell_types_to_run()}){
    my $new_input_id = eval($self->input_id);
    $new_input_id->{"cell_type"} = $cell_type;
    $self->dataflow_output_id($new_input_id, 2, { } );
  }

  return 1;

}

#Private getter / setter to the cell_types_to_run list
sub _cell_types_to_run {
  return $_[0]->_getter_setter('cell_types_to_run',$_[1]);
}

1;
