#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::InputSetAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::InputSetAdaptor - A database adaptor for fetching and
storing InputSet objects.

=head1 SYNOPSIS

my $rset_adaptor = $db->get_InputSetAdaptor();

my @rsets = @{$rset_adaptor->fetch_all_InputSets_by_Experiment()};
my @displayable_rsets = @{$rset_adaptor->fetch_all_displayable_InputSets()};

#Other methods?
#by FeatureType, CellType all with displayable flag?


=head1 DESCRIPTION

The InputSetAdaptor is a database adaptor for storing and retrieving
InputSet objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::InputSetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor; #DBI sql_types import

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor );


=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : $exp_set = $eseta->fetch_by_Experiment($exp);
  Description: Retrieves a InputSet based on the given Experiment
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : Throws if no valid stored Experiment provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Experiment {
  my ($self, $exp) = @_;
  my $params = {constraints => {experiments => [$exp]}};
  return $self->generic_fetch($self->compose_constraint_query($params));
}


=head2 fetch_all_by_InputSubsets

  Arg [1]    : List of Bio::EnsEMBL::Funcgen::InputSubsets
  Example    :
  Description:
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions :
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_InputSubsets {
  my ($self, $subsets) = @_;

  my $dbIDs =
    $self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSubset', $subsets, 'dbID');

  my $constraint = "iss.input_set_id IN ($dbIDs)";
  return $self->generic_fetch($constraint);
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
  #can't have 'is' as an alias as it is reserved
  return (
    [ 'input_set',              'inp' ],
    [ 'input_set_input_subset', 'isiss'],
    [ 'input_subset',           'iss' ]);
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
    return qw(
      inp.input_set_id
      inp.analysis_id
      inp.cell_type_id
      inp.experiment_id
      inp.feature_type_id
      inp.name
      inp.replicate
      inp.type
         );
}


=head2 _left_join

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
               queries.
  Returntype : List
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

#rather than default_where so we still get back InputSets without InputSubsets
#Doc why left join is needed here!

sub _left_join {
  return (
      ['input_set_input_subset',  'inp.input_set_id    = isiss.input_set_id'],
      ['input_subset',            'iss.input_subset_id = isiss.input_subset_id']
      );
}


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
      $dbid,
      $anal_id,
      $anal,
      $ctype_id,
      $ctype,
      $exp_id,
      $exp,
      $ftype_id,
      $ftype,
      $name,
      $replicate,
      $type,
      $input_set,
      @input_sets,
      );

  $sth->bind_columns(
      \$dbid,
      \$anal_id,
      \$ctype_id,
      \$exp_id,
      \$ftype_id,
      \$name,
      \$replicate,
      \$type,
      );

  my $anal_adaptor  = $self->db->get_AnalysisAdaptor();
  my $ct_adaptor    = $self->db->get_CellTypeAdaptor();
  my $exp_adaptor   = $self->db->get_ExperimentAdaptor();
  my $ft_adaptor    = $self->db->get_FeatureTypeAdaptor();

  while ( $sth->fetch() ) {

    if(! $input_set || ($input_set->dbID() != $dbid)){

      push @input_sets, $input_set if $input_set;

      $anal = (defined $anal_id) ? $anal_adaptor->fetch_by_dbID($anal_id) : undef;
      throw("Could not fetch Analysis with dbID '$anal_id' for InputSet '$name'") if ! $anal;

      $ctype = (defined $ctype_id) ? $ct_adaptor->fetch_by_dbID($ctype_id) : undef;
      throw("Could not fetch CellType with dbID '$ctype_id' for InputSet '$name'") if ! $ctype;

      $exp = (defined $exp_id) ? $exp_adaptor->fetch_by_dbID($exp_id) : undef;
      throw("Could not fetch Experiment with dbID '$exp_id' for InputSet '$name'") if ! $exp;

      $ftype = (defined $ftype_id) ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
      throw("Could not fetch FeatureType with dbID '$ftype_id' for InputSet '$name'") if ! $ftype;


      $input_set = Bio::EnsEMBL::Funcgen::InputSet->new(
          -DBID          => $dbid,
          -ANALYSIS      => $anal,
          -CELL_TYPE     => $ctype,
          -EXPERIMENT    => $exp,
          -FEATURE_TYPE  => $ftype,
          -NAME          => $name,
          -REPLICATE     => $replicate,
          -FEATURE_CLASS => $type,
          -ADAPTOR       => $self,
          );
          
      push @input_sets, $input_set;
    }
  }
  
  return \@input_sets;
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::InputSet objects
  Example    : $rsa->store(@esets);
  Description: Stores or updates previously stored InputSet objects in the database.
  Returntype : None
  Exceptions : Throws if a List of InputSet objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut

sub store{
  my ($self, @input_sets) = @_;

  throw("Must provide a list of InputSet objects") if(scalar(@input_sets == 0));

  my $sth = $self->prepare('
      INSERT INTO
        input_set (
          analysis_id,
          experiment_id,
          cell_type_id,
          feature_type_id,
          name,
          replicate,
          type
          )
       VALUES (?, ?, ?, ?, ?, ?, ?)');

  my $sth_link_table = $self->prepare('
      INSERT INTO
        input_set_input_subset (
          input_set_id,
          input_subset_id
          )
      VALUES
        (?,?)
      ');

  my $db = $self->db();
  foreach my $set (@input_sets) {

    if( ! ref $set || ! $set->isa('Bio::EnsEMBL::Funcgen::InputSet') ) {
      throw('Must be Bio::EnsEMBL::Funcgen::InputSet, not '. ref($set));
    }
    if(! exists $set->{'subsets'}){
      throw("No subset(s) provided for " . $set->name);
    }
    # Subsets are not array, but hash with the name as key
    for my $name (keys %{$set->{'subsets'}}){
      my $ss = $set->{'subsets'}->{$name};
      if( ! ref $ss || ! $ss->isa('Bio::EnsEMBL::Funcgen::InputSubset') ) {
        throw("InputSubset object `required, not ".ref($ss)." InputSet: ".$ss->name);
        if(!$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $ss)){
          # Store or demand stored subset?
        }
      }
    }

    if ( $set->is_stored($db) ) {
      throw('InputSet [' . $set->dbID() . '] is already stored in the database - InputSetAdaptor does not yet accomodate updating InputSets');
      #would need to retrive stored result set and update table_ids
    }

### add ->dbID and mmove tp bind_param
    my $an_id = $set->analysis()->dbID;
    my $ct_id = $set->cell_type()->dbID;
    my $ex_id = $set->get_Experiment->dbID;
    my $ft_id = $set->feature_type()->dbID;

    $sth->bind_param(1, $an_id,               SQL_INTEGER);
    $sth->bind_param(2, $ex_id,               SQL_INTEGER);
    $sth->bind_param(3, $ct_id,               SQL_INTEGER);
    $sth->bind_param(4, $ft_id,               SQL_INTEGER);
    $sth->bind_param(5, $set->name,           SQL_VARCHAR);
    $sth->bind_param(6, $set->replicate,      SQL_INTEGER);
    $sth->bind_param(7, $set->feature_class,  SQL_VARCHAR);
    $sth->execute();

    $set->dbID( $self->last_insert_id );
    $set->adaptor($self);
    $self->store_states($set);

    for my $name (sort keys %{$set->{'subsets'}}){
      my $ss = $set->{'subsets'}->{$name};
      $sth_link_table->bind_param(1, $set->dbID, SQL_INTEGER);
      $sth_link_table->bind_param(2, $ss->dbID,  SQL_INTEGER);
      $sth_link_table->execute;
    }
  }

  return \@input_sets;
}

### GENERIC CONSTRAIN METHODS ###

#All these _constrain methods must return a valid constraint string,
#and a hashref of any other constraint config

#Need to bind param any of these which come from URL parameters and are not tested

sub _constrain_analyses {
  my ($self, $anals) = @_;

  my $constraint = 'inp.analysis_id IN ('.
      join(', ', @{$self->db_are_stored_and_valid('Bio::EnsEMBL::Analysis',$anals, 'dbID')}).')';

  return($constraint);
}

sub _constrain_cell_types {
  my ($self, $cts) = @_;

  #Don't need to bind param this as we validate
  my $constraint = ' inp.cell_type_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $cts, 'dbID')}).')';

  #{} = no further config
  return ($constraint, {});
}

sub _constrain_experiments {
  my ($self, $exps) = @_;


  #This needs updating to link through via input_subsets (ignoring controls)
  #and having a final clause to avoid the product!


  #Don't need to bind param this as we validate
  my $constraint = ' inp.experiment_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exps, 'dbID')}).')';

  #{} = not futher constraint conf
  return ($constraint, {});
}

sub _constrain_feature_types {
  my ($self, $fts) = @_;

  my @tables = $self->_tables;
  my (undef, $syn) = @{$tables[0]};

  #Don't need to bind param this as we validate
  my $constraint = " ${syn}.feature_type_id IN (".
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fts, 'dbID')}).')';

  #{} = not futher constraint conf
  return ($constraint, {});
}


######### Deprecated ###########
# remove format?

sub _constrain_format {
  throw('v74. Use analysis instead');
  my ($self, $format) = @_;

  #Is not currently enum'd so have to hardcode current values for now
  #likely to change
  #SEQUENCING EQTL

  my %valid_formats = (SEQUENCING=>1);
  #SEGMENTATION?

  if (! exists $valid_formats{uc($format)}) {
    throw("$format is not a valid InputSet format, please specify one of:\t".
          join(', ', keys %valid_formats));
  }

  my $constraint = ' inp.format="'.uc($format).'"';

  #{} = not futher constraint conf
  return ($constraint, {});
}

=head2 store_InputSubsets

  Args       : Bio::EnsEMBL::Funcgen::InputSet
  Example    : $esa->store_InputSubsets(\@e_subsets);
  Description: Convenience methods extracted from store to allow updating of InputSubset entries
               during inline result processing which would otherwise be troublesome due to the need
               for an InputSet
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : Throws if a stored InputSet object is not provided
               Throws if no InputSubsets present
  Caller     : General
  Status     : At Risk

=cut

#DEPRECATE THIS!!
sub store_InputSubsets{
  throw('v74. Use InputSubsetAdaptor->store instead');
  my ($self, $ssets) = @_;

  my $sth = $self->prepare("
        INSERT INTO input_subset (
            input_set_id, name, archive_id, display_url, replicate, is_control
        ) VALUES (?, ?, ?, ? ,?, ?)
    ");

  throw('Must provide at least one InputSubset') if(! @$ssets);

  #Store and set all previously unstored table_ids
  foreach my $sset(@$ssets){

    #use is_stored here?
    if($sset->dbID()){
      warn "Skipping InputSubset ".$sset->name()." - already stored in the DB";
      next;
    }


    $sth->bind_param(1, $sset->input_set->dbID, SQL_INTEGER);
    $sth->bind_param(2, $sset->name,            SQL_VARCHAR);
    $sth->bind_param(3, $sset->archive_id,      SQL_VARCHAR);
    $sth->bind_param(4, $sset->display_url,     SQL_VARCHAR);
    $sth->bind_param(5, $sset->replicate,       SQL_INTEGER);
    $sth->bind_param(6, $sset->is_control,      SQL_INTEGER);
    $sth->execute();

    $sset->dbID($self->last_insert_id);
    $sset->adaptor($self);
  }

  #don't really need to return as we're working on a ref
  return $ssets;
}

1;

