#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor - A database adaptor for fetching and
storing Funcgen Channel objects.

=head1 SYNOPSIS

my $chan_a = $db->get_ChannelAdaptor();

my @channels = @{$chan_a->fetch_all_by_ExperimentalChip($exp)};


=head1 DESCRIPTION

The ChannelAdaptor is a database adaptor for storing and retrieving
Funcgen Channel objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#sql_types barewords import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_by_type_experimental_chip_id

  Arg [1]    : string - type
  Arg [2]    : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Example    : my $chan = $chan_a->fetch_by_type_experimental_chip_dbID('EXPERIMENTAL', $ec_dbid);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::Channel object
  Exceptions : Throws if args not met
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_type_experimental_chip_id {
    my ($self, $type, $ec_id) = @_;

    throw("Must specify a channel type e.g. EXPERIMENTAL, TOTAL and an ExperimentalChip id") if(! ($type && $ec_id));


	my $constraint = 'c.experimental_chip_id = ? AND c.type = ?';

    $self->bind_param_generic_fetch($ec_id,     SQL_INTEGER);
    $self->bind_param_generic_fetch($type,      SQL_VARCHAR);


    return $self->generic_fetch($constraint);
  }


=head2 fetch_by_dye_experimental_chip_dbID

  Arg [1]    : string - dye
  Arg [1]    : int - dbID of Experiment
  Example    : my $chan = $chan_a->fetch_by_type_experimental_chip_dbID('Cy5', $ec_dbid);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::Channel object
  Exceptions : Throws is experiment dbID or dye not passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_dye_experimental_chip_dbID {
    my $self = shift;
    my $dye = shift;
    my $ec_dbid = shift;

    throw("not yet impemented");

    my ($chan_id, @results);

    throw("Must specify an experiemntal dbID") if(! $ec_dbid);
    #Need to validate dye against VendorDefs?  Or leave and just return whatever is in DB e.g. nothing if dye name is wrong.


    my $sth = $self->prepare("
		SELECT c.channel_id
		FROM experimental_chip ec, channel c
		WHERE c.experimental_chip_id = ec.experimental_chip_id
        AND ec.experimental_chip_id = $ec_dbid
	");



	#can we do a generic fetch here?


	$sth->execute();


	while ($chan_id = $sth->fetchrow()){
		push @results, $self->fetch_by_dbID($chan_id);
	}

	return \@results;
}


=head2 fetch_all_by_experimental_chip_dbID

  Arg [1]    : int - dbID of Experiment
  Example    : my @chans = @{$ec_a->fetch_all_by_experimental_chip_dbID($ac_dbid);
  Description: Does what it says on the tin
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel object
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_experimental_chip_dbID {
    my $self = shift;
    my $ec_dbid = shift;

	my ($chan_id, @results);

    throw("Must specify an experiemntal dbID") if(! $ec_dbid);


    #my $sth = $self->prepare("
    #SELECT c.channel_id
    #		FROM experimental_chip ec, channel c
    #		WHERE c.experimental_chip_id = ec.experimental_chip_id
    #        AND ec.experimental_chip_id = $ec_dbid
    #	");


    my $constraint = "c.experimental_chip_id=$ec_dbid";

	#can we do a generic fetch here?


    #$sth->execute();


    #while ($chan_id = $sth->fetchrow()){
    #	  push @results, $self->fetch_by_dbID($chan_id);
    #	      }

    return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_ExperimentalChip

  Arg [1]    : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Example    : my @chans = @{$ec_a->fetch_all_by_ExperimentalChip($echip);
  Description: Returns all channels associated with a given ExperimentalChip
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : Throws if no ExperimentalChip defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ExperimentalChip{
  my ($self, $exp) = @_;
  throw("Must provide an ExperimentChip object") if(! $exp->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip'));
  return $self->fetch_all_by_experimental_chip_dbID($exp->dbID());
}

=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty Array objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Array getters
  Status     : Medium Risk

=cut

sub fetch_attributes {
    my $self = shift;
    my $array = shift;

    my $tmp_array = $self->fetch_by_dbID( $array->dbID() );
    %$array = %$tmp_array;
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
  return (['channel', 'c']);
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
  return qw( c.channel_id c.experimental_chip_id c.sample_id c.dye c.type );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $ec_id, $chan_id, $sample_id, $type, $dye);

	$sth->bind_columns(\$chan_id, \$ec_id, \$sample_id, \$dye, \$type);

	while ( $sth->fetch() ) {
	  my $chan = Bio::EnsEMBL::Funcgen::Channel->new(
							 -DBID                 => $chan_id,
							 -EXPERIMENTAL_CHIP_ID => $ec_id,
							 -SAMPLE_ID            => $sample_id,
							 -TYPE                 => $type,
							 -DYE                  => $dye,
							 -ADAPTOR              => $self,
							);

	  push @result, $chan;

	}
	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Channel objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores given Channel objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if object is not a Bio::EnsEMBL::Funcgen::Channel
               Throws if object is already present in the DB but has no dbID
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;

  my $sth = $self->prepare("
			INSERT INTO channel
			(experimental_chip_id, sample_id, dye, type)
			VALUES (?, ?, ?, ?)");


  foreach my $chan (@args) {
    throw('Can only store Channel objects') if ( ! $chan->isa('Bio::EnsEMBL::Funcgen::Channel'));

    if (!( $chan->dbID() && $chan->adaptor() == $self )){#use is_stored?


      my $s_chan = $self->fetch_by_type_experimental_chip_id($chan->type(), $chan->experimental_chip_id());
      throw("Channel already exists in the database with dbID:".$s_chan->dbID().
	    "\nTo reuse/update this Channel you must retrieve it using the ChannelAdaptor".
	    "\nMaybe you want to use the -recover option?") if $s_chan;

      #if(! $s_chan){
      $sth->bind_param(1, $chan->experimental_chip_id(),  SQL_INTEGER);
      $sth->bind_param(2, $chan->sample_id(),             SQL_VARCHAR);
      $sth->bind_param(3, $chan->dye() ,                  SQL_VARCHAR);
      $sth->bind_param(4, $chan->type(),                  SQL_VARCHAR);

      $sth->execute();
      $chan->dbID($self->last_insert_id);
      $chan->adaptor($self);
      #}
      #else{
      #  #do some status checks here, check IMPORTED
      #  #Need to account for recover in Importer?
      #  $chan = $s_chan;

      #  my @states = @{$self->db->fetch_all_states('channel', $chan->dbID())};

      #  #need better id than dbID?
      #  warn("Using previously stored Channel (".$chan->experimental_chip_id().":".$chan->type().") with states\t@states\n");
      #}
    }else{
      #assume we want to update the states
      warn('You may want to use $chan->adaptor->store_states($chan)');
      $self->store_states($chan);
    }
  }

  return \@args;
}


1;

