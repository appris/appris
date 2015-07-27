#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalGroupAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalGroupAdaptor - A database adaptor for fetching and
storing Funcgen ExperimentalGroup objects.

=head1 SYNOPSIS

my $eg_adaptor = $db->get_ExperimentalGroupAdaptor();

my $group = $eg_adaptor->fetch_by_name("EBI");

=head1 DESCRIPTION

The ExperimentalGroupAdaptor is a database adaptor for storing and retrieving
Funcgen ExperimentalGroup objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::ExperimentalGroup

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalGroupAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Funcgen::ExperimentalGroup;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_by_name

  Arg [1]    : string - name of ExperimentalGroup
  Example    : my $group = $eg_adaptor->fetch_by_name('EBI');
  Description: Fetches the group with the given name
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalGroup object
  Exceptions :
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name{
  my ($self, $name) = @_;

  throw("Must specify a Group name") if(! $name);

  my $constraint = " name = ?";

  $self->bind_param_generic_fetch($name,   SQL_VARCHAR);

  my @fts = @{$self->generic_fetch($constraint)};

  return $fts[0];
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
  return (['experimental_group', 'eg']);
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
  return qw( eg.experimental_group_id eg.name eg.location
             eg.contact eg.url eg.description eg.is_project );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalGroup objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $eg_id, $name, $location, $contact, $url, $desc, $is_project);

	$sth->bind_columns(\$eg_id, \$name, \$location, \$contact, \$url, \$desc, \$is_project);

	while ( $sth->fetch() ) {
	  my $group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
								    -dbID        => $eg_id,
								    -NAME        => $name,
								    -LOCATION    => $location,
								    -CONTACT     => $contact,
								    -URL         => $url,
								    -DESCRIPTION => $desc,
								    -IS_PROJECT  => $is_project,
								    -ADAPTOR     => $self,
								   );

	  push @result, $group;

	}
	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ExperimentalGroup objects
  Example    : $eg_adaptor->store($g1, $g2, $g3);
  Description: Stores given ExperimentalGroup objects in the database.
  Returntype : None
  Exceptions : Throws if ExperimentalGroup not valid
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;


  my $sth = $self->prepare("
			INSERT INTO experimental_group
			(name, location, contact, url, description, is_project)
			VALUES (?, ?, ?, ?, ?, ?)");



  foreach my $group (@args) {

    if ( ! (ref($group) && $group->isa('Bio::EnsEMBL::Funcgen::ExperimentalGroup') )) {
      throw('Can only store ExperimentalGroup objects, skipping $group');
    }

    if (!( $group->dbID() && $group->adaptor() == $self )){

      #Check for previously stored FeatureType
      my $s_eg = $self->fetch_by_name($group->name());

      if(! $s_eg){
	$sth->bind_param(1, $group->name(),           SQL_VARCHAR);
	$sth->bind_param(2, $group->location(),       SQL_VARCHAR);
	$sth->bind_param(3, $group->contact(),        SQL_VARCHAR);
	$sth->bind_param(4, $group->url(),            SQL_VARCHAR);
	$sth->bind_param(5, $group->description(),    SQL_VARCHAR);
	$sth->bind_param(6, $group->is_project(),     SQL_BOOLEAN);

	$sth->execute();
	$group->dbID($self->last_insert_id);
	$group->adaptor($self);
      }
      else{
	$group = $s_eg;
	warn("Using previously stored ExperimentalGroup:\t".$group->name()."\n");
      }
    }
  }

  return \@args;
}



1;

