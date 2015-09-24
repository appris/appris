#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ExternalFeatureAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::ExternalFeatureAdaptor - A database adaptor for fetching and
storing ExternalFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_ExternalFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The ExternalFeatureAdaptor is a database adaptor for storing and retrieving
ExternalFeature objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ExternalFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::ExternalFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor; #DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor);



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
  return ([ 'external_feature', 'ef' ],
          [ 'feature_set',      'fs' ]);
}
#feature_set is required for fetching on analysis, external_db (cell_type or feature_type).
#and should be removed in favour of query composition
#could put some of these type constraint methods in SetFeatureAdaptor
#but would have to have another feature_type method in here
#as these are in the external_feature table rather than feature_set table


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
			ef.external_feature_id   ef.seq_region_id
			ef.seq_region_start      ef.seq_region_end
			ef.seq_region_strand     ef.display_label
			ef.feature_type_id       ef.feature_set_id
			ef.interdb_stable_id
	   );
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates ExternalFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::ExternalFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;

	#For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	#So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.


	my ($sa, @features, $seq_region_id, %fset_hash, %slice_hash, %sr_name_hash, %sr_cs_hash, %ftype_hash);
	$sa = $dest_slice->adaptor->db->get_SliceAdaptor() if($dest_slice);#don't really need this if we're using DNADBSliceAdaptor?
	$sa ||= $self->db->get_SliceAdaptor();

	my $fset_adaptor = $self->db->get_FeatureSetAdaptor();
	my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor();
	my (
	    $external_feature_id,  $efg_seq_region_id,
	    $seq_region_start,      $seq_region_end,
	    $seq_region_strand,     $fset_id,
		$display_label,         $ftype_id,
		$interdb_stable_id
	   );

	$sth->bind_columns(
					   \$external_feature_id,   \$efg_seq_region_id,
					   \$seq_region_start,      \$seq_region_end,
					   \$seq_region_strand,     \$display_label,
					   \$ftype_id,              \$fset_id,
					   \$interdb_stable_id
					  );


	#abstract to BaseFeatureAdaptor?
	my ($asm_cs, $cmp_cs, $asm_cs_name, $asm_cs_vers, $cmp_cs_name, $cmp_cs_vers);

	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my ($dest_slice_start, $dest_slice_end, $dest_slice_strand, $dest_slice_length, $dest_slice_sr_name);

	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}


  FEATURE: while ( $sth->fetch() ) {
	  #Need to build a slice adaptor cache here?
	  #Would only ever want to do this if we enable mapping between assemblies??
	  #Or if we supported the mapping between cs systems for a given schema_build,
	  #which would have to be handled by the core api

	  #get core seq_region_id
	  $seq_region_id = $self->get_core_seq_region_id($efg_seq_region_id);

	  if(! $seq_region_id){
		warn "Cannot get slice for eFG seq_region_id $efg_seq_region_id\n".
		  "The region you are using is not present in the current dna DB";
		next;
	  }


	  #Get the FeatureSet/Type objects
	  $fset_hash{$fset_id} = $fset_adaptor->fetch_by_dbID($fset_id) if(! exists $fset_hash{$fset_id});
	  $ftype_hash{$ftype_id} = $ftype_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});

	  # Get the slice object
	  my $slice = $slice_hash{'ID:'.$seq_region_id};

	  if (! $slice) {
		$slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
		$slice_hash{'ID:'.$seq_region_id} = $slice;
		$sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
		$sr_cs_hash{$seq_region_id}       = $slice->coord_system();
	  }

	  my $sr_name = $sr_name_hash{$seq_region_id};
	  my $sr_cs   = $sr_cs_hash{$seq_region_id};

	  #abstract to BaseFeatureAdaptor?
	  # Remap the feature coordinates to another coord system if a mapper was provided
	  if ($mapper) {

		throw("Not yet implmented mapper, check equals are Funcgen calls too!");

	      ($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand)
			= $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);

	      # Skip features that map to gaps or coord system boundaries
	      next FEATURE if !defined $sr_name;

	      # Get a slice in the coord system we just mapped to
	      if ( $asm_cs == $sr_cs || ( $cmp_cs != $sr_cs && $asm_cs->equals($sr_cs) ) ) {
		$slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"}
		  ||= $sa->fetch_by_region($cmp_cs_name, $sr_name, undef, undef, undef, $cmp_cs_vers);
	      } else {
		$slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"}
		  ||= $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef, $asm_cs_vers);
	      }
	    }

	    # If a destination slice was provided convert the coords
	    # If the destination slice starts at 1 and is forward strand, nothing needs doing
	    if ($dest_slice) {
	      unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {
		if ($dest_slice_strand == 1) {
		  $seq_region_start = $seq_region_start - $dest_slice_start + 1;
		  $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
		} else {
		  my $tmp_seq_region_start = $seq_region_start;
		  $seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
		  $seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
		  $seq_region_strand      *= -1;
		}
	      }

	      # Throw away features off the end of the requested slice
	      next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
		|| ( $dest_slice_sr_name ne $sr_name );

	      $slice = $dest_slice;
	    }

	  push @features, Bio::EnsEMBL::Funcgen::ExternalFeature->new_fast
		({
		  'start'          => $seq_region_start,
		  'end'            => $seq_region_end,
		  'strand'         => $seq_region_strand,
		  'slice'          => $slice,
		  'analysis'       => $fset_hash{$fset_id}->analysis(),
		  'adaptor'        => $self,
		  'dbID'           => $external_feature_id,
		  'display_label'  => $display_label,
		  'set'    => $fset_hash{$fset_id},
		  'feature_type'   => $ftype_hash{$ftype_id},
		  'interdb_stable_id',    => $interdb_stable_id,
		 });
	  }

	return \@features;
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ExternalFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given ExternalFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Arrayref of stored ExternalFeatures
  Exceptions : Throws if a list of ExternalFeature objects is not provided or if
               the Analysis, CellType and FeatureType objects are not attached or stored
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @efs) = @_;

	if (scalar(@efs) == 0) {
		throw('Must call store with a list of ExternalFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO external_feature (
			seq_region_id,   seq_region_start,
			seq_region_end,  seq_region_strand,
            display_label,   feature_type_id,
            feature_set_id
		) VALUES (?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();

  FEATURE: foreach my $ef (@efs) {

	  if(! (ref($ef) && $ef->isa('Bio::EnsEMBL::Funcgen::ExternalFeature') )) {
		throw('Feature must be an ExternalFeature object');
	  }

	  if ( $ef->is_stored($db) ) {
		#does not accomodate adding Feature to >1 feature_set
		warning('ExternalFeature [' . $ef->dbID() . '] is already stored in the database');
		next FEATURE;
	  }

	  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $ef->feature_set);

	  my $seq_region_id;
	  ($ef, $seq_region_id) = $self->_pre_store($ef);

	  $sth->bind_param(1, $seq_region_id,             SQL_INTEGER);
	  $sth->bind_param(2, $ef->start(),               SQL_INTEGER);
	  $sth->bind_param(3, $ef->end(),                 SQL_INTEGER);
	  $sth->bind_param(4, $ef->strand(),              SQL_TINYINT);
	  $sth->bind_param(5, $ef->display_label(),       SQL_VARCHAR);
	  $sth->bind_param(6, $ef->feature_type->dbID(),  SQL_INTEGER);
	  $sth->bind_param(7, $ef->feature_set->dbID(),   SQL_INTEGER);

	  $sth->execute();
	  $ef->dbID( $self->last_insert_id );
	  $ef->adaptor($self);
	  $self->store_associated_feature_types($ef) if (defined $ef->{'associated_feature_types'});
	}

  return \@efs;
}

=head2 fetch_by_interdb_stable_id

  Arg [1]    : Integer $stable_id - The 'interdb stable id' of the ExternalFeature to retrieve
  Example    : my $rf = $rf_adaptor->fetch_by_interdb_stable_id(1);
  Description: Retrieves a ExternalFeature via its interdb_stable id. This is really an internal
               method to facilitate inter DB linking.
  Returntype : Bio::EnsEMBL::Funcgen::ExternalFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_interdb_stable_id {
  my ($self, $stable_id) = @_;

  $self->bind_param_generic_fetch($stable_id, SQL_INTEGER);

  return $self->generic_fetch('ef.interdb_stable_id=?')->[0];
}

1;
