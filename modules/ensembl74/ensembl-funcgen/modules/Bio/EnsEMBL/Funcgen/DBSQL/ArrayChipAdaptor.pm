#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor - A database adaptor for fetching and
storing Funcgen ArrayChip objects.

=head1 SYNOPSIS

my $ac_adaptor = $db->get_ArrayChipAdaptor();

my @achips = @{$ec_adaptor->fetch_all_by_Array($array)};


=head1 DESCRIPTION

The ArrayChipAdaptor is a database adaptor for storing and retrieving
Funcgen ArrayChip objects.


=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_all_by_array_id

  Arg [1]    : int - dbID of Array
  Example    : my @ccs = @{$ec_a->fetch_all_by_array_dbID($array->dbID());
  Description: Does what it says on the tin
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ArrayChip
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_array_id {
  my ($self, $array_id) = @_;
  throw("Must specify an array dbID") if(! $array_id);

  $self->bind_param_generic_fetch($array_id, SQL_INTEGER);
  return $self->generic_fetch('ac.array_id= ?');
}


=head2 fetch_all_by_Array

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array
  Example    : my @ccs = @{$ec_a->fetch_all_by_Array($array);
  Description: Returns all ArrayChips which belong to the given Array
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::ArrayChip objects
  Exceptions : Throws if ARg not valid and stored
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Array {
  my ($self, $array) = @_;

  if(! (ref($array) && $array->isa->('Bio::EnsEMBL::Funcgen::Array') && $array->dbID())){
	throw("Must pass a valid stored Bio::EnsEMBL::Funcgen::Array");
  }

  return $self->fetch_all_by_array_id($array->dbID);
}


=head2 fetch_all_by_ExperimentalChips

  Arg [1]    : arrayref of - Bio::EnsEMBL::Funcgen::ExperimentalChips
  Example    : my @achips = @{$ec_a->fetch_all_by_ExperimentalChips($echips);
  Description: Gets a non-redundant list of the corresponding ArrayChips
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ArrayChips
  Exceptions : Throws if no ExperimentalChips passed or if ExperimentalChips are not stored or valid
  Caller     : General
  Status     : at risk

=cut

sub fetch_all_by_ExperimentalChips {
  my ($self, $echips) = @_;

  my %ac_ids;

  foreach my $echip(@$echips){

	if(! ($echip->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip') && $echip->dbID())){
	  throw('Must provide an arrayref of valid stored Bio::EnsEMBL::Funcgen::ExperimentalChips');
	}

	$ac_ids{$echip->array_chip_id} = 1;

  }

  if(! keys(%ac_ids)){
	throw('Must provide an arrayref of valid stored Bio::EnsEMBL::Funcgen::ExperimentalChips');
  }

  return $self->generic_fetch('ac.array_chip_id IN ('.join(', ', keys(%ac_ids)).')');
}


=head2 fetch_by_array_design_ids

  Arg [1]    : int - dbID of Array
  Arg [2]    : string - Design ID of ArrayChip
  Example    : my $ac = $ac_adaptpr->fetch_by_array_design_ids($array->dbID, $design_id);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::ArrayChip
  Exceptions : Throws if args not met.
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_array_design_ids{
    my ($self, $array_id, $design_id) = @_;

	if( ! ($array_id && $design_id)){
	  throw('You must pass an Array dbID and a ArrayChip design ID');
	}

    my $constraint = "ac.array_id='$array_id' and ac.design_id='$design_id'";

    my ($ac) = @{$self->generic_fetch($constraint)};
    #unique key means this always has just one element

    return $ac;
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
  return (['array_chip', 'ac']);
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _columns {
  return qw( ac.array_chip_id ac.design_id ac.array_id ac.name );
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ArrayChip objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (@result, $ac_id, $design_id, $array_id, $name);
  $sth->bind_columns(\$ac_id, \$design_id, \$array_id, \$name);

  while ( $sth->fetch() ) {
    my $array = Bio::EnsEMBL::Funcgen::ArrayChip->new(
						      -dbID      => $ac_id,
						      -design_id => $design_id,
						      -array_id  => $array_id,
						      -name      => $name,
						      -adaptor   => $self,
						     );

    push @result, $array;

  }
  return \@result;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ArrayChip objects
  Example    : $arraychip_adaptor->store($ac1, $ac2, $ac3);
  Description: Stores given ArrayChip objects in the database. Should only be
               called once per ArrayChip because no checks are made for duplicates.
  Returntype : ARRAYREF
  Exceptions : warns if passed non-ArrayChip arg
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;

  #my ($stored_ac);

  #Should we implement a throw here is the caller is not Array?
  #make private _store?


  my $sth = $self->prepare("
			INSERT INTO array_chip
			(design_id, array_id, name)
			VALUES (?, ?, ?)"
			  );



  foreach my $ac (@args) {
    if ( ! $ac->isa('Bio::EnsEMBL::Funcgen::ArrayChip') ) {
      warning('Can only store ExperimentalChip objects, skipping $ec');
      next;
    }

    throw("ArrayChip must have an array_id to be stored") if ! $ac->array_id();

    #check for array_id here? this is done by not null in sql


    #check for previously stored array_chips is done in Array via add_ArrayChip



    if (!( $ac->dbID() && $ac->adaptor() == $self )){

      my $s_ac = $self->fetch_by_array_design_ids($ac->array_id(), $ac->design_id());

      if(! $s_ac){

	$sth->bind_param(1, $ac->design_id(), SQL_VARCHAR);
	$sth->bind_param(2, $ac->array_id(),  SQL_INTEGER);
	$sth->bind_param(3, $ac->name(),      SQL_VARCHAR);
	$sth->execute();

	$ac->dbID($self->last_insert_id);
	$ac->adaptor($self);

      }else{
	  $ac = $s_ac;
	  #	my @states = @{$self->db->fetch_all_states('experimental_chip', $ec->dbID())};
	  #	warn("Using previously stored ExperimentalChip (".$ec->unique_id().") with states\t@states\n");
	}
    }
  }
  return \@args;
}


1;

