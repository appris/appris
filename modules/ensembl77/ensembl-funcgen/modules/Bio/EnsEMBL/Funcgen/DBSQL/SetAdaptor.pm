#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::SetAdaptor
#

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::DBSQL::Funcgen::SetAdaptor - A base adaptor for Input/Result/FeatureSetAdaptors

=head1 SYNOPSIS




=head1 DESCRIPTION

This base SetAdaptor provides generic methods applicable to Input/Result/FeatureSetAdaptors

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import
use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor );

use vars qw(@EXPORT); #require Exporter done in parent
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});


=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : String  (optional) - status e.g. 'DISPLAYABLE'
  Example    : my @sets = @{$set_adaptor->fetch_all_by_FeatureType($ftype)};
  Description: Retrieves Set objects from the database based on FeatureType
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Set objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
  my ($self, $ftype, $status) = @_;
  my $params = {constraints => {feature_types => [$ftype]}};
  $params->{constraints}{states} = [$status] if defined $status;
  my $results = $self->generic_fetch($self->compose_constraint_query($params));
  $self->reset_true_tables;  #As we may have added status
  return $results;
}

=head2 fetch_all_by_CellType

  Arg [1]    : Bio::EnsEMBL::Funcgen::CellType
  Arg [2]    : String  (optional) - status e.g. 'DISPLAYABLE'
  Example    : my @sets = @{$set_adaptor->fetch_all_by_CellType($ctype)};
  Description: Retrieves Set objects from the database based on CellType
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Set objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_CellType {
  my ($self, $ctype, $status) = @_;
  my $params = {constraints => {cell_types => [$ctype]}};
  $params->{constraints}{states} = [$status] if defined $status;
  my $results = $self->generic_fetch($self->compose_constraint_query($params));
  $self->reset_true_tables; #As we may have added status
  return $results;
}


=head2 fetch_all_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen::Analysis
  Example    : my @sets = @{$set_adaptopr->fetch_all_by_Analysis($analysis)};
  Description: Retrieves Set objects from the database based on Analysis
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Set objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Analysis {
  my ($self, $ctype, $status) = @_;
  my $params = {constraints => {analyses => [$ctype]}};
  $params->{constraints}{states} = [$status] if defined $status;
  my $results = $self->generic_fetch($self->compose_constraint_query($params));
  $self->reset_true_tables; #As we may have added status
  return $results;
}


=head2 fetch_all_by_FeatureType_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : Bio::EnsEMBL::Analysis
  Arg [3]    : (optional) Bio::EnsEMBL::Funcgen::CellType
  Example    : my @sets = $set_adaptopr->fetch_all_by_FeatureType_Analysis($ftype, $anal, $ctype);
  Description: Retrieves Set objects from the database based on FeatureType, Analysis and
               CellType if defined.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Set objects
  Exceptions : Throws if args 1 and 2 are not valid or stored
  Caller     : General
  Status     : At Risk

=cut

#Historical method from the FeatureSetAdaptor

sub fetch_all_by_FeatureType_Analysis {
  my ($self, $ftype, $anal, $ctype) = @_;
  my $params = {constraints =>
                {
                 feature_types => [$ftype],
                 analyses     => [$anal],
                }
               };

  $params->{constraints}{cell_types} = [$ctype] if $ctype;
  return $self->generic_fetch($self->compose_constraint_query($params));
}


#No fetch_all_by_name as this is not useful

=head2 fetch_by_name

  Arg [1]    : String - Set name
  Example    : my $iss = $set_a->fetch_by_name('Iss_name');
  Description: Retrieves a Set object which matches the specified name
  Returntype : Bio::EnsEMBL::Funcgen::Set
  Exceptions : Throws if more than one returned
  Caller     : General
  Status     : Stable

=cut

sub fetch_by_name {
  my ($self, $name) = @_;

  my $params = {constraints => {name => $name}};
  my $results = $self->generic_fetch($self->compose_constraint_query($params));

  if(scalar(@$results) >1){
    throw('The name specified is not unique, please use another fetch_by_name method');  
  }
  
  return $results->[0];
}



# can't have fetch_by_name as this name is not unique for ResultSets

### GENERIC CONSTRAIN METHODS ###

#All these _constrain methods must return a valid constraint string, and a hashref of any other constraint config

#Need to bind param any of these which come from URL parameters and are not tested

#These type constraints are generic to all sets (apart from DataSet)
#and so could go in a SetAdaptor (doesn't exist yet)
#currently have

sub _constrain_cell_types {
  my ($self, $cts) = @_;

  my $constraint = $self->_table_syn.'.cell_type_id IN ('.
        join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $cts, 'dbID')}
        ).')';

  #{} = no futher contraint config
  return ($constraint, {});
}


sub _constrain_feature_types {
  my ($self, $fts) = @_;

  #Don't need to bind param this as we validate
  my $constraint = $self->_table_syn.'.feature_type_id IN ('.
		join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fts, 'dbID')}).')';

  #{} = no futher constraint conf
  return ($constraint, {});
}


sub _constrain_analyses {
  my ($self, $anals) = @_;

  #Don't need to bind param this as we validate
  my $constraint = $self->_table_syn.'.analysis_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Analysis', $anals, 'dbID')}).')';

  return ($constraint, {});   #{} = no futher constraint conf
}

sub _constrain_name {
  my ($self, $name) = @_;

  if(! defined $name) {
    throw("Need to specify a name argument");
  }
  my $constraint = $self->_table_syn . ".name = '$name'";
  return ($constraint, {});
}

1;

