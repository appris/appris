#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor - A database adaptor for fetching and
storing Funcgen ExperimentalChip objects.

=head1 SYNOPSIS

my $ec_a = $db->get_ExperimentalChipAdaptor();

my @ecs = @{$ec_a->fetch_all_by_Experiment($exp)};


=head1 DESCRIPTION

The ExperimentalChipAdaptor is a database adaptor for storing and retrieving
Funcgen ExperimentalChip objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_all_by_experiment_dbID

  Arg [1]    : int - dbID of Experiment
  Example    : my @ecs = @{$ec_a->fetch_all_by_experiment_dbID($ac_dbid);
  Description: Does what it says on the tin
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_experiment_dbID {
    my $self = shift;
    my $e_dbid = shift;

	my ($ec_id, @results);

	throw("Must specify an experiemntal dbID") if(! $e_dbid);


	my $sth = $self->prepare("
		SELECT ec.experimental_chip_id
		FROM experimental_chip ec, experiment e
		WHERE ec.experiment_id = $e_dbid
	");



	#can we do a generic fetch here?


	$sth->execute();


	while ($ec_id = $sth->fetchrow()){
	  #warn("got ec id $ec_id\n");
	  push @results, $self->fetch_by_dbID($ec_id);
	}

	return \@results;
}

=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my @ecs = @{$ec_a->fetch_all_by_Experiment($exp)};
  Description: Does what it says on the tin
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalChips
  Exceptions : throws if valid stored Experiment not passed
  Caller     : General
  Status     : at risk

=cut


sub fetch_all_by_Experiment(){
  my ($self, $exp) = @_;

  if(! ($exp && $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') && $exp->dbID())){
	throw('Must provide a valid stored Bio::EnsEMBL::Funcgen::Experiment');
  }

  $self->generic_fetch("ec.experiment_id=".$exp->dbID());
}

=head2 fetch_all_by_ArrayChip

  Arg [1]    : Bio::EnsEMBL::Funcgen::ArrayChip
  Example    : my @ecs = @{$ec_a->fetch_all_by_ArrayChip($achip)};
  Description: Retrieves all ExperimentChips which have the corresponding ArrayChip design
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalChips
  Exceptions : throws if valid stored ArrayChip not passed
  Caller     : General
  Status     : at risk

=cut


sub fetch_all_by_ArrayChip{
  my ($self, $achip) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ArrayChip', $achip);

  $self->generic_fetch("ec.array_chip_id=".$achip->dbID());
}

=head2 fetch_by_unique_and_experiment_id

  Arg [2]    : int - unique_id
  Arg [1]    : int - dbID of Experiment
  Example    : my $ec = ec_a->fetch_by_unique_and_experiment_id($c_uid, $exp_dbid);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_unique_and_experiment_id {
  my ($self, $c_uid, $e_dbid) = @_;




  throw("Must provide and unique_id and and experiment_id") if(! $c_uid || ! $e_dbid);

  my $sth = $self->prepare("
		SELECT ec.experimental_chip_id
		FROM experimental_chip ec
		WHERE ec.unique_id ='$c_uid'
        AND ec.experiment_id = $e_dbid
	");


  $sth->execute();
  my ($ec_id) = $sth->fetchrow();


  return $self->fetch_by_dbID($ec_id) if $ec_id;
}


=head2 fetch_by_unique_id_vendor

  Arg [2]    : string - unique_id
  Arg [1]    : string - name of array vendor e.g. NIMBLEGEN
  Example    : my $ec = $ec_a->fetch_by_unique_id_vendor($c_uid, $vendor);
  Description: Fetches a unique ExperimentalChip for a given unique id and vendor
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None
  Caller     : General
  Status     : at Risk

=cut

sub fetch_by_unique_id_vendor {
  my ($self, $c_uid, $vendor) = @_;

  throw("Must provide and unique_id and and vendor") if(! $c_uid || ! $vendor);

  my ($ec_id, $ac_id, @ecids, $avendor);


  my $sth = $self->prepare("
		SELECT ec.experimental_chip_id, ec.array_chip_id
		FROM experimental_chip ec
		WHERE ec.unique_id ='$c_uid'
	");


  $sth->execute();

  while (($ec_id, $ac_id) = $sth->fetchrow()){

    my $sql = "SELECT a.vendor from array a, array_chip ac where a.array_id=ac.array_id and ac.array_chip_id=$ac_id;";
    ($avendor) = @{$self->db->dbc->db_handle->selectrow_arrayref($sql)};
    push @ecids, $ec_id if (uc($avendor) eq uc($vendor));
  }

  #This check shouldn't be necessary if this control is illicited on import
  #no unique key possible so just for safety until import fully tested.
  if(scalar(@ecids) > 1){
    throw("Found more than one ExperimentalChip with the same unique_id($c_uid) for $vendor");

  }

  return $self->fetch_by_dbID($ecids[0]) if $ecids[0];
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
  return (['experimental_chip', 'ec']);
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
  return qw( ec.experimental_chip_id  ec.unique_id
	   	     ec.experiment_id         ec.array_chip_id
		     ec.feature_type_id       ec.cell_type_id
		     ec.biological_replicate  ec.technical_replicate);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (@result, $ec_id, $c_uid, $exp_id, $ac_id, $ftype_id, $ctype_id, $brep, $trep, $ftype, $ctype);

  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
  my $ct_adaptor = $self->db->get_CellTypeAdaptor();


  $sth->bind_columns(\$ec_id, \$c_uid, \$exp_id, \$ac_id, \$ftype_id, \$ctype_id, \$brep, \$trep);

  while ( $sth->fetch() ) {

    $ftype = (defined $ftype_id) ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
    $ctype = (defined $ctype_id) ? $ct_adaptor->fetch_by_dbID($ctype_id) : undef;

    my $array = Bio::EnsEMBL::Funcgen::ExperimentalChip->new(
															 -dbID           => $ec_id,
															 -unique_id      => $c_uid,
															 -experiment_id  => $exp_id,
															 -array_chip_id  => $ac_id,
															 -feature_type   => $ftype,
															 -cell_type      => $ctype,
															 -biological_replicate => $brep,
															 -technical_replicate  => $trep,
															 -adaptor        => $self,
															);

    push @result, $array;

  }
  return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Example    : $oaa->store($ec1, $ec2, $ec3);
  Description: Stores given ExperimentalChip objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : ARRAYREF
  Exceptions : Throws if passed non-ExperimentalChip arg or if ExperimentalChip already stored but arg has no dbID
  Caller     : General
  Status     : Medium Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;

  my ($sarray);

  my $sth = $self->prepare("
			INSERT INTO experimental_chip
			(unique_id, experiment_id, array_chip_id, feature_type_id,
             cell_type_id, biological_replicate, technical_replicate)
			VALUES (?, ?, ?, ?, ?, ?, ?)");

  foreach my $ec (@args) {
    throw('Can only store ExperimentalChip objects') if ( ! $ec->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip') );

    if (!( $ec->dbID() && $ec->adaptor() == $self )){

      #Some slight jiggery pokery here to get the vendor as the EC does not yet have an adaptor
      #cache these?
      my $vendor = $self->db->get_ArrayChipAdaptor->fetch_by_dbID($ec->array_chip_id)->get_Array->vendor();
      my $s_ec = $self->fetch_by_unique_id_vendor($ec->unique_id(), $vendor);

      throw("ExperimentalChip already exists in the database with dbID:".$s_ec->dbID().
	    "\nTo reuse/update this ExperimentalChip you must retrieve it using the ExperimentalChipAdaptor".
	    "\nMaybe you want to use the -recover option?") if $s_ec;

      my $ftype_id = (defined $ec->feature_type()) ? $ec->feature_type->dbID() : undef;
      my $ctype_id = (defined $ec->cell_type()) ? $ec->cell_type->dbID() : undef;

      $sth->bind_param(1, $ec->unique_id(),            SQL_VARCHAR);
      $sth->bind_param(2, $ec->experiment_id(),        SQL_VARCHAR);
      $sth->bind_param(3, $ec->array_chip_id(),        SQL_VARCHAR);
      $sth->bind_param(4, $ftype_id,                   SQL_INTEGER);
      $sth->bind_param(5, $ctype_id,                   SQL_INTEGER);
      $sth->bind_param(6, $ec->biological_replicate(), SQL_VARCHAR);
	  $sth->bind_param(7, $ec->technical_replicate(),  SQL_VARCHAR);

      $sth->execute();
      $ec->dbID($self->last_insert_id);
      $ec->adaptor($self);

      #}
      #else{
      #	  $ec = $s_ec;

      #my @states = @{$self->db->fetch_all_states('experimental_chip', $ec->dbID())};
      #	  my @states = @{$self->db->get_StatusAdaptor->fetch_all_states($ec)};
      #	  warn("Using previously stored ExperimentalChip (".$ec->unique_id().") with states\t@states\n");
      #  }
    }else{
      #assume we want to update the states
      warn('You may want to use $exp_chip->adaptor->store_states($exp_chip)');
      $self->store_states($ec);
    }
  }

  return \@args;
}




sub update_replicate_types{
  my ($self, $echip) = @_;

  if(! ($echip && $echip->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip') && $echip->dbID())){
	throw('Must provide a valid store Bio::EnsEMBL::Funcgen::ExperimentalChip');
  }

  my $sql = 'UPDATE experimental_chip set biological_replicate="'.$echip->biological_replicate().
	'" where experimental_chip_id='.$echip->dbID();

  $self->db->dbc->do($sql);

  $sql = 'UPDATE experimental_chip set technical_replicate="'.$echip->technical_replicate().
	'" where experimental_chip_id='.$echip->dbID();
  $self->db->dbc->do($sql);


  if(defined $echip->cell_type()){
	$sql = 'UPDATE experimental_chip set cell_type_id="'.$echip->cell_type()->dbID.
	  '" where experimental_chip_id='.$echip->dbID();
	$self->db->dbc->do($sql);
  }

  if(defined $echip->feature_type()){
	$sql = 'UPDATE experimental_chip set feature_type_id="'.$echip->feature_type()->dbID.
	  '" where experimental_chip_id='.$echip->dbID();
	$self->db->dbc->do($sql);
  }



  return;
}


1;

