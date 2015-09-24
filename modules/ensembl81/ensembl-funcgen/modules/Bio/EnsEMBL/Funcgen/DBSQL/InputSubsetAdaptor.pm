#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor
my $iss_a = $db->get_InputSubsetAdaptor;


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::InputSubset;
use Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor; #DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor);


=head2 fetch_all_by_InputSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::InputSet
  Example    : my @subsets = @{$iss_adaptor->fetch_all_by_InputSet($input_set)};
  Description: Retrieves InputSubset objects from the database based that are part
               of the given InputSet
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::InputSubset objects
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut


#Was stable in 74
#This will go away when we drop InputSets

sub fetch_all_by_InputSet {
  my ($self, $iset) = @_;

  my $params = {constraints => {input_sets => [$iset]}};
	my $objs = $self->generic_fetch($self->compose_constraint_query($params));
  $self->reset_true_tables;

  return($objs);
}


=head2 fetch_all_by_Experiments

  Arg [1]    : Arrayref of Bio::EnsEMBL::Funcgen::Experiment objects
  Example    : my @subsets = @{$iss_adaptopr->fetch_all_by_Experiments($exp_objs)};
  Description: Retrieves InputSubset objects from the database that belong to
               the given Experiments
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::InputSubset objects
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_Experiments {
  my ($self, $exps) = @_;

  my $params = {constraints => {experiments => $exps}};
	return $self->generic_fetch($self->compose_constraint_query($params));
}

#experiment is only a 2nd order index (unique)key and analysis_id is not indexe at all


=head2 fetch_by_name

  Arg [1]    : String - InputSubset name
  Arg [2]    : Bio::EnsEMBL::Funcgen::Experiment (optional)
  Example    : my $iss = $iss_a->fetch_by_name('Iss_name');
  Description: Retrieves a InputSubset object which matches the passed name 
               and optional Experiment
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : Throws if no name or defined
               Throws if no Experiment and cannot identify unique InputSubset
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name {
  my ($self, $name, $exp) = @_;

  my $params = {constraints => {name        => $name}};
  $params->{contraints}{experiments} = [$exp] if $exp;
  
  my $sets = $self->generic_fetch($self->compose_constraint_query($params));
  
  if(scalar(@$sets) > 1){
    throw("Could not identify a unique InputSubset with name:\t$name".
      "\nPlease specify an Experiment argument");  
  }
  
  return $sets->[0]; 
}



=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for creating the object
               and queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  my ($self) = @_;

  return qw(
      iss.input_subset_id
      iss.cell_type_id
      iss.experiment_id
      iss.feature_type_id
      iss.is_control
      iss.name
      iss.replicate
      iss.analysis_id
      );
}


=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  return ([ 'input_subset',  'iss' ]);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
               handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::InputSubset objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (@result,
      $iss_id,
      $ct_id,
      $exp_id,
      $ft_id,
      $is_control,
      $name,
      $replicate,
      $analysis_id
     );

#  my $is_adaptor = $self->db->get_InputSetAdaptor;

  #Fetching the InputSet will need to be removed once we implement the lin ktable as there maybe >!
  #so we would need a hash of input_sets => {dbDB1 => $iset1, dbID2 => $iset2}
  #we're initially the values are undef, and lazy loaded on request.

  $sth->bind_columns(
      \$iss_id,
      \$ct_id,
      \$exp_id,
      \$ft_id,
      \$is_control,
      \$name,
      \$replicate,
      \$analysis_id
      );

  my $ct_adaptor  = $self->db->get_CellTypeAdaptor;
  my $exp_adaptor = $self->db->get_ExperimentAdaptor;
  my $ft_adaptor  = $self->db->get_FeatureTypeAdaptor;
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor; 
  my (%ctypes, %ftypes, %exps);  
 
  while($sth->fetch){
    
    if(! exists $ctypes{$ct_id}){
      $ctypes{$ct_id} = $ct_adaptor->fetch_by_dbID($ct_id);
    
      if(! defined $ctypes{$ct_id}){
        throw("Could not fetch linked CellType (dbID: $ct_id) for InputSubset (dbID: $iss_id) ");
      }
    }
    
   if(! exists $exps{$exp_id}){
      $exps{$exp_id} = $exp_adaptor->fetch_by_dbID($exp_id);
    
      if(! defined $exps{$exp_id}){
        throw("Could not fetch linked Experiment (dbID: $exp_id) for InputSubset (dbID: $iss_id) ");
      }
    }
    
    if(! exists $ftypes{$ft_id}){
      $ftypes{$ft_id} = $ft_adaptor->fetch_by_dbID($ft_id);
    
      if(! defined $ftypes{$ft_id}){
        throw("Could not fetch linked FeatureType (dbID: $ft_id) for InputSubset (dbID: $iss_id) ");
      }
    }
    
    #Don't cache analyses, as the adaptor does this
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id) ||
      throw("Could not fetch linked Analysis (dbID: $analysis_id) for InputSubset (dbID: $iss_id)");

    push @result, Bio::EnsEMBL::Funcgen::InputSubset->new (
       -dbID         => $iss_id,
       -CELL_TYPE    => $ctypes{$ct_id},
       -EXPERIMENT   => $exps{$exp_id},
       -FEATURE_TYPE => $ftypes{$ft_id},
       -IS_CONTROL   => $is_control,
       -NAME         => $name,
       -REPLICATE    => $replicate,
       -ANALYSIS     => $analysis,
       -ADAPTOR      => $self,
      );
  }
  
  return \@result;
}


=head2 store

  Args       : Bio::EnsEMBL::Funcgen::InputSet
  Example    : $iss_adaptor->store(\@e_subsets);
  Description: Stores or updates previously stored InputSubset objects in the database.
  Returntype : None
  Exceptions : Throws if a List of InputSet objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut

sub store{
  my ($self, @subsets) = @_;

  if(scalar(@subsets == 0)){
    throw("Must provide a list of InputSubset objects")
  }

  my $db = $self->db();

  my $sth = $self->prepare("
        INSERT INTO
          input_subset (
            cell_type_id,
            experiment_id,
            feature_type_id,
            is_control,
            name,
            replicate,
            analysis_id
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ");


  #Store and set all previously unstored table_ids
  foreach my $subset(@subsets){
    if( ! ref $subset || ! $subset->isa('Bio::EnsEMBL::Funcgen::InputSubset') ) {
      throw("InputSubset object required, not '".ref($subset) ."'");
    }

    if($subset->is_stored($db)){
      throw( "InputSubset '" . $subset->name . "' already stored in the DB.".
            "InputSubsetAdaptor does not support updating yet.");
    }
    
    #Test object attrs are stored
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $subset->feature_type);
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType',    $subset->cell_type);
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment',  $subset->experiment);
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis',             $subset->analysis);
    

    $sth->bind_param(1, $subset->cell_type->dbID,     SQL_INTEGER);
    $sth->bind_param(2, $subset->experiment->dbID,    SQL_INTEGER);
    $sth->bind_param(3, $subset->feature_type->dbID,  SQL_INTEGER);
    $sth->bind_param(4, $subset->is_control,          SQL_INTEGER);
    $sth->bind_param(5, $subset->name,                SQL_VARCHAR);
    $sth->bind_param(6, $subset->replicate,           SQL_INTEGER);
    $sth->bind_param(7, $subset->analysis->dbID,      SQL_INTEGER);
    $sth->execute();

    $subset->dbID($self->last_insert_id);
    $subset->adaptor($self);
    $self->store_states($subset);
  }

  return \@subsets;
}


sub _final_clause {
  return ( 'group by iss.input_subset_id');
}



### CONTRAIN METHODS ###
# Used by BaseAdaptor::compose_query_contraint

sub _constrain_experiments {
  my ($self, $exps) = @_;

  my $constraint = ' iss.experiment_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exps, 'dbID')}).')';

  #{} = not futher constraint conf
  return ($constraint, {});

}

sub _constrain_input_sets {
  my ($self, $isets) = @_;

  my $constraint_conf = {
    tables => [
      ['input_set_input_subset', 'isiss']
      ]};

  my $constraint = '
    iss.input_subset_id = isiss.input_subset_id AND
    isiss.input_set_id  IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $isets, 'dbID')}).')';

  return ($constraint, $constraint_conf);
}





###DEPRECATED METHODS ###

#should have been fetch_by_name_and_Experiment
#Experiment is not needed in the majority of cases
#was also returning Arrayref

sub fetch_by_name_and_experiment { #Deprecated in 75
  my ($self, $name, $exp) = @_;

  deprecate('Please use the new fetch_by_name method. This method will be removed in release 79');

  return $self->fetch_by_name($name, $exp);
  #my $params =
  # {constraints => {experiments => [$exp],
  #                  name        => $name}};
  #return $self->generic_fetch($self->compose_constraint_query($params));
}

sub _constrain_archive_ids {
  my ($self, $archive_ids) = @_;

  if ( (ref($archive_ids) ne 'ARRAY') || scalar(@$archive_ids) == 0 ) {
    throw('Must pass an arrayref of archive IDs');
  }
  my $constraint;

  #{} = not futher constraint conf
  return (' iss.archive_id IN '.join(', ', (map {uc($_)} @$archive_ids)) , {});
}

1;

