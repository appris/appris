#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::AnnotatedFeatureAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::AnnotatedFeatureAdaptor - A database adaptor for fetching and
storing AnnotatedFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_AnnotatedFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The AnnotatedFeatureAdaptor is a database adaptor for storing and retrieving
AnnotatedFeature objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::AnnotatedFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor; #import sql_types barewords

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
  return (['annotated_feature', 'af'], ['feature_set', 'fs']);
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
			af.annotated_feature_id  af.seq_region_id
			af.seq_region_start      af.seq_region_end
			af.seq_region_strand     af.feature_set_id
			af.display_label         af.score
			af.summit
			);
}



=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates AnnotatedFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::AnnotatedFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;

	my $sa           = $self->db->dnadb->get_SliceAdaptor;
	my $fset_adaptor = $self->db->get_FeatureSetAdaptor;

	my ($seq_region_id, @features, %fset_hash,
      %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
	    $annotated_feature_id,  $efg_seq_region_id,
	    $seq_region_start,      $seq_region_end,
	    $seq_region_strand,     $fset_id,
      $display_label,         $score,
      $summit
     );

	$sth->bind_columns(
                     \$annotated_feature_id,  \$efg_seq_region_id,
                     \$seq_region_start,      \$seq_region_end,
                     \$seq_region_strand,     \$fset_id,
                     \$display_label,         \$score,
                     \$summit
                    );


	my ($asm_cs, $cmp_cs, $asm_cs_name, $asm_cs_vers, $cmp_cs_name, $cmp_cs_vers);

	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my ($dest_slice_start, $dest_slice_end, $dest_slice_strand,
      $dest_slice_length, $dest_slice_sr_name, $slice, $sr_name, $sr_cs);

	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}


 FEATURE: while ( $sth->fetch ) {
    #get core seq_region_id
	  #This was failing if we are using a 'comparable' CoordSystem as we don't have a cache
	  #for the new DB. Fixed with the tmp seq_region_cache
	  $seq_region_id = $self->get_core_seq_region_id($efg_seq_region_id);

	  if (! $seq_region_id) {
      warn "Cannot get slice for eFG seq_region_id $efg_seq_region_id\n".
        "The region you are using is not present in the current seq_region_cache.\n".
          "Maybe you need to redefine the dnadb or update_DB_for_release?";
      next;
	  }

	  #Get the FeatureSet object
	  $fset_hash{$fset_id} = $fset_adaptor->fetch_by_dbID($fset_id) if(! exists $fset_hash{$fset_id});


	  # Get the slice object
	  $slice = $slice_hash{'ID:'.$seq_region_id};

	  if (! $slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{'ID:'.$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
	  }

	  $sr_name = $sr_name_hash{$seq_region_id};
    $sr_cs   = $sr_cs_hash{$seq_region_id};

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
		  if (! $self->force_reslice) {
        #force_reslice set by RegulatoryFeature::regulatory_attributes
        #so we don't lose attrs which are not on the dest_slice

        next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
          || ( $dest_slice_sr_name ne $sr_name );
      }

      $slice = $dest_slice;
		}

    push @features, Bio::EnsEMBL::Funcgen::AnnotatedFeature->new_fast
		  ({
        start          => $seq_region_start,
        end            => $seq_region_end,
        strand         => $seq_region_strand,
        slice          => $slice,
        adaptor        => $self,
        dbID           => $annotated_feature_id,
        score          => $score,
        summit         => $summit,
        display_label  => $display_label,
        set            => $fset_hash{$fset_id},
		   });
	}

	return \@features;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::AnnotatedFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given AnnotatedFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Listref of stored AnnotatedFeatures
  Exceptions : Throws if a list of AnnotatedFeature objects is not provided or if
               the Analysis, CellType and FeatureType objects are not attached or stored
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @pfs) = @_;

	if (scalar(@pfs) == 0) {
		throw('Must call store with a list of AnnotatedFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO annotated_feature (
			seq_region_id,   seq_region_start,
			seq_region_end,  seq_region_strand,
            feature_set_id,  display_label, score, summit
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
	");


	my $db = $self->db();

  FEATURE: foreach my $pf (@pfs) {

		if( !ref $pf || !$pf->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature') ) {
			throw('Feature must be an AnnotatedFeature object');
		}

		if ( $pf->is_stored($db) ) {
			#does not accomodate adding Feature to >1 feature_set
			warning('AnnotatedFeature [' . $pf->dbID() . '] is already stored in the database');
			next FEATURE;
		}


		#Only need to check FeatureSet here, as FeatureSet store will check analysis

		if (! $pf->feature_set->is_stored($db)) {
			throw('A stored Bio::EnsEMBL::Funcgen::FeatureSet must be attached to the AnnotatedFeature objects to be stored.');
		}


		my $seq_region_id;
		($pf, $seq_region_id) = $self->_pre_store($pf);

		$sth->bind_param(1, $seq_region_id,             SQL_INTEGER);
		$sth->bind_param(2, $pf->start(),               SQL_INTEGER);
		$sth->bind_param(3, $pf->end(),                 SQL_INTEGER);
		$sth->bind_param(4, $pf->strand(),              SQL_TINYINT);
		$sth->bind_param(5, $pf->feature_set->dbID(),   SQL_INTEGER);
		$sth->bind_param(6, $pf->display_label(),       SQL_VARCHAR);
		$sth->bind_param(7, $pf->score(),                SQL_DOUBLE);
		$sth->bind_param(8, $pf->summit,                SQL_INTEGER);

		$sth->execute();
		$pf->dbID( $self->last_insert_id );
		$pf->adaptor($self);
	}

  return \@pfs;
}


=head2 fetch_all_by_associated_MotifFeature

  Arg [1]    : Bio::EnsEMBL::Funcgen::MotifFeature
  Example    : my $assoc_afs = $af_adaptor->fetch_all_by_associated_MotifFeature($motif_feature);
  Description: Fetches all associated AnnotatedFeatures for a given MotifFeature.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::AnnotatedFeature objects
  Exceptions : Throws if arg is not valid or stored
  Caller     : General
  Status     : At risk

=cut

# Can we standardise this across AnnotatedFeatures and Regulatory/SegmentationFeatures?
#wrt regulatory_attribute vs associated_motif_feature

sub fetch_all_by_associated_MotifFeature{
  my ($self, $mf) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::MotifFeature', $mf);
  $self->_tables([['associated_motif_feature', 'amf']]);
  my $constraint = 'amf.annotated_feature_id=af.annotated_feature_id AND amf.motif_feature_id='.$mf->dbID;
  my $afs        =  $self->generic_fetch($constraint);
  $self->reset_true_tables;
  return $afs;
}


#TODO GENERIC CONSTRAIN METHODS
#1 motif_features



1;
