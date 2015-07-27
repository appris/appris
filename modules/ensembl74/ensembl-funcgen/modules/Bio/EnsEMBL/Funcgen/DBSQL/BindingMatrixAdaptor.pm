#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor - A database adaptor for fetching and
storing Funcgen BindingMatrix objects.

=head1 SYNOPSIS

my $matrix_adaptor = $db->get_BindingMatrixAdaptor();
my @matrices = @{$matrix_adaptor->fetch_all_by_name("MA0122.1")};

=head1 DESCRIPTION

The BindingMatrixAdaptor is a database adaptor for storing and retrieving
Funcgen BindingMatrix objects.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::BindingMatrix

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#sql_types barewords import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_all_by_name

  Arg [1]    : string - name of Matrix
  Arg [2]    : Bio::EnsEMBL::Analysis (optional) Analysis indicating Matrix origin
  Example    : my @matrices = @{$matrix_adaptor->fetch_all_by_name('MA0122.1')};
  Description: Fetches matrix objects given a name and an optional Analysis object.
               If both are specified, only one unique BindingMatrix will be returned
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Exceptions : Throws if no name if defined
  Caller     : General
  Status     : At risk - Change this to fetch_all_by_name_FeatureType

=cut

sub fetch_all_by_name{
  my ($self, $name, $analysis) = @_;

  throw("Must specify a BindingMatrix name") if(! $name);

  my $constraint = " bm.name = ? ";
  $constraint .= " AND bm.analysis_id = ?" if $analysis;

  $self->bind_param_generic_fetch($name,           SQL_VARCHAR);
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER) if $analysis;

  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_name_FeatureType

  Arg [1]    : string - name of Matrix
  Arg [2]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [3]    : Bio::EnsEMBL::Analysis (optional) Analysis indicating Matrix origin
  Description: Fetches matrix objects given a name and a FeatureType.
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Exceptions : Throws if no name if defined or if FeatureType is not valid
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_name_FeatureType{
  my ($self, $name, $ftype, $analysis) = @_;

  throw("Must specify a BindingMatrix name") if(! $name);
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);

  my $constraint = " bm.name = ? and bm.feature_type_id = ?";
  $constraint .= " AND bm.analysis_id = ?" if $analysis;

  $self->bind_param_generic_fetch($name,           SQL_VARCHAR);
  $self->bind_param_generic_fetch($ftype->dbID,    SQL_INTEGER);
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER) if $analysis;

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : Bio::EnsEMBL::Analysis (optional) Analysis indicating Matrix origin
  Example    : my @matrices = @{$matrix_adaptor->fetch_all_by_FeatureType($ftype)};
  Description: Fetches BindingMatrix objects given it's FeatureType
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : Throws if FeatureType is not valid
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_FeatureType{
  my ($self, $ftype, $analysis) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);

  my $constraint = " bm.feature_type_id = ?";
  $constraint .= " AND bm.analysis_id = ?" if $analysis;

  $self->bind_param_generic_fetch($ftype->dbID,    SQL_INTEGER);
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER) if $analysis;

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
  return (['binding_matrix', 'bm']);
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
  return qw( bm.binding_matrix_id bm.name bm.analysis_id bm.frequencies
             bm.description bm.feature_type_id bm.threshold);
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates objects from an executed DBI statement handle.
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $matrix_id, $name, $analysis_id, $freq, $desc, $ftype_id, $thresh);
	$sth->bind_columns(\$matrix_id, \$name, \$analysis_id, \$freq, \$desc, \$ftype_id, \$thresh);

	my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
	my %ftype_cache;

	my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
	my %analysis_cache;

	while ( $sth->fetch() ) {

	  if(! exists $ftype_cache{$ftype_id}){
		$ftype_cache{$ftype_id} = $ftype_adaptor->fetch_by_dbID($ftype_id);
	  }

	  if(! exists $analysis_cache{$analysis_id}){
		$analysis_cache{$analysis_id} = $analysis_adaptor->fetch_by_dbID($analysis_id);
	  }

	  my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new
		(
		 -dbID         => $matrix_id,
		 -NAME         => $name,
		 -ANALYSIS     => $analysis_cache{$analysis_id},
		 -FREQUENCIES  => $freq,
		 -DESCRIPTION  => $desc,
		 -FEATURE_TYPE => $ftype_cache{$ftype_id},
		 -THRESHOLD    => $thresh,
		 -ADAPTOR      => $self,
		);

	  push @result, $matrix;

	}

	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Example    : $matrix_adaptor->store($m1, $m2, $m3);
  Description: Stores given Matrix objects in the database.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;

  my $sth = $self->prepare("
			INSERT INTO binding_matrix
			(name, analysis_id, frequencies, description, feature_type_id, threshold)
			VALUES (?, ?, ?, ?, ?, ?)");

  my $s_matrix;

  foreach my $matrix (@args) {

    if ( ! $matrix->isa('Bio::EnsEMBL::Funcgen::BindingMatrix') ) {
      throw('Can only store BindingMatrix objects, skipping $matrix');
    }

	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $matrix->feature_type);



    if (!( $matrix->dbID() && $matrix->adaptor() == $self )){

      #Check for previously stored BindingMatrix
      ($s_matrix) = @{$self->fetch_all_by_name_FeatureType($matrix->name(), $matrix->feature_type, $matrix->analysis())};

      if(! $s_matrix){

		$sth->bind_param(1, $matrix->name(),               SQL_VARCHAR);
		$sth->bind_param(2, $matrix->analysis()->dbID(),   SQL_INTEGER);
		$sth->bind_param(3, $matrix->frequencies(),        SQL_LONGVARCHAR);
		$sth->bind_param(4, $matrix->description(),        SQL_VARCHAR);
		$sth->bind_param(5, $matrix->feature_type->dbID(), SQL_INTEGER);
		$sth->bind_param(6, $matrix->threshold(),          SQL_DOUBLE);

		$sth->execute();
		$matrix->dbID($self->last_insert_id);
		$matrix->adaptor($self);

		$self->store_associated_feature_types($matrix);
	  }
      else{
		$matrix = $s_matrix;
		warn("Using previously stored Matrix:\t".$matrix->name()."\n");
		#Could update associated FeatureTypes here
      }
    }
  }

  return \@args;
}



1;
