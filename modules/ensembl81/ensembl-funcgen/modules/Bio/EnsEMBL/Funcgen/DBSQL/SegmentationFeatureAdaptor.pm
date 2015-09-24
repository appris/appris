#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::SegmentationFeatureAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::SegmentationFeatureAdaptor - A database adaptor for fetching and
storing SegmentationFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_SegmentationFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The SegmentationFeatureAdaptor is a database adaptor for storing and retrieving
SegmentationFeature objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::SegmentationFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::SegmentationFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor;#DBI sql_types import

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
  my $self = shift;

  return (
		  ['segmentation_feature', 'sf'],
		  ['feature_set', 'fs']
		 );
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
			sf.segmentation_feature_id  sf.seq_region_id
			sf.seq_region_start      sf.seq_region_end
			sf.seq_region_strand     sf.feature_type_id
			sf.feature_set_id        sf.score
			sf.display_label
	   );
}



=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates SegmentationFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::SegmentationFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;
	#For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	#So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.
	my $sa = $dest_slice ? $dest_slice->adaptor->db->get_SliceAdaptor : 
                         $self->db->dnadb->get_SliceAdaptor;

	#Some of this in now probably overkill as we'll always be using the DNADB as the slice DB
	#Hence it should always be on the same coord system
	my $fset_adaptor  = $self->db->get_FeatureSetAdaptor;
	my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
	my (@features, %fset_hash, %slice_hash, %sr_name_hash, %sr_cs_hash, %ftype_hash);

	my ($seq_region_id,  $segmentation_feature_id, $efg_seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand,       $ftype_id,
      $fset_id,        $score,                   $display_label	);

	$sth->bind_columns(
    \$segmentation_feature_id, \$efg_seq_region_id,
    \$seq_region_start,        \$seq_region_end,
    \$seq_region_strand,       \$ftype_id,
    \$fset_id,        	   	  \$score,
    \$display_label
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

	my ($dest_slice_start,  $dest_slice_end, $dest_slice_strand, 
      $dest_slice_length, $dest_slice_sr_name);

	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}

  FEATURE: while ( $sth->fetch ) {
	  #Need to build a slice adaptor cache here?
	  #Would only ever want to do this if we enable mapping between assemblies??
	  #Or if we supported the mapping between cs systems for a given schema_build, which would have to be handled by the core api

	  #get core seq_region_id
	  #This fails if we are using a 'comparable' CoordSystem as we don't have a cache
	  #for the new DB. Wasn't this fixed with the tmp seq_region_cache?
	  $seq_region_id = $self->get_core_seq_region_id($efg_seq_region_id);

	  if(! $seq_region_id){
      warn "Cannot get slice for eFG seq_region_id $efg_seq_region_id\n".
       "The region you are using is not present in the current seq_region_cache.\n".
			 "Maybe you need to redefine the dnadb or update_DB_for_release?";
      next;
	  }

	  #Get the FeatureSet FeatureType & Slice objects
	  $fset_hash{$fset_id} = $fset_adaptor->fetch_by_dbID($fset_id) if ! exists $fset_hash{$fset_id};
	  $ftype_hash{$ftype_id} = $ftype_adaptor->fetch_by_dbID($ftype_id) if ! exists $fset_hash{$ftype_id};
	  my $slice = $slice_hash{'ID:'.$seq_region_id};

	  if (! $slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{'ID:'.$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name;
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system;
	  }

	  my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

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
		  if(! $self->force_reslice){
  			#force_reslice set by RegulatoryFeature::regulatory_attributes
	  		#so we don't lose attrs which are not on the dest_slice

        next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
         || ( $dest_slice_sr_name ne $sr_name );

        $slice = $dest_slice;
		  }
		}

	  push @features, Bio::EnsEMBL::Funcgen::SegmentationFeature->new_fast
     ({start          => $seq_region_start,
		   end            => $seq_region_end,
		   strand         => $seq_region_strand,
		   slice          => $slice,
		   #'analysis'       => $fset_hash{$fset_id}->analysis(),
		   #Need to pass this to keep Feature.pm happy
		   #Let's grab this from the FeatureSet in SetFeature new and pass
		   score          => $score,
		   adaptor        => $self,
		   dbID           => $segmentation_feature_id,
		   display_label  => $display_label,
		   set            => $fset_hash{$fset_id},
		   feature_type   =>  $ftype_hash{$ftype_id},
		  });
    }

	return \@features;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::SegmentationFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given SegmentationFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Listref of stored SegmentationFeatures
  Exceptions : Throws if a list of SegmentationFeature objects is not provided or if
               the Analysis, CellType and FeatureType objects are not attached or stored
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @segs) = @_;

	if (scalar(@segs) == 0) {
		throw('Must call store with a list of SegmentationFeature objects');
	}

	my $sth = $self->prepare('
		INSERT INTO segmentation_feature (
			seq_region_id,   seq_region_start,
			seq_region_end,  seq_region_strand,
      feature_type_id, feature_set_id,   
      score,           display_label
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?)');
	my $db = $self->db;


  FEATURE: foreach my $seg (@segs) {

		if( !ref $seg || !$seg->isa('Bio::EnsEMBL::Funcgen::SegmentationFeature') ) {
			throw('Feature must be an SegmentationFeature object');
		}

		if ( $seg->is_stored($db) ) {
			#does not accomodate adding Feature to >1 feature_set
			warning('SegmentationFeature [' . $seg->dbID() . '] is already stored in the database');
			next FEATURE;
		}

		$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet',  $seg->feature_set);
		$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $seg->feature_type);

		my $seq_region_id;
		($seg, $seq_region_id) = $self->_pre_store($seg);

		$sth->bind_param(1, $seq_region_id,           SQL_INTEGER);
		$sth->bind_param(2, $seg->start,              SQL_INTEGER);
		$sth->bind_param(3, $seg->end,                SQL_INTEGER);
		$sth->bind_param(4, $seg->strand,             SQL_TINYINT);
		$sth->bind_param(5, $seg->feature_type->dbID, SQL_INTEGER);
		$sth->bind_param(6, $seg->feature_set->dbID,  SQL_INTEGER);
		$sth->bind_param(7, $seg->score,               SQL_DOUBLE);
		$sth->bind_param(8, $seg->display_label,      SQL_VARCHAR);

		$sth->execute();
		$seg->dbID( $self->last_insert_id );
		$seg->adaptor($self);
	}

  return \@segs;
}


1;
