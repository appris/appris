#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::RegulatoryFeatureAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::RegulatoryFeatureAdaptor

=head1 SYNOPSIS

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;

my $reg = Bio::EnsEMBL::Registry->load_adaptors_from_db(-host    => 'ensembldb.ensembl.org',
                                                        -user    => 'anonymous');

my $regfeat_adaptor = $reg->get_adaptor($species, 'funcgen', 'RegulatoryFeature');


#Fetch MultiCell RegulatoryFeatures
my @features = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};

#Fetch cell type specific RegulatoryFeatures
my @ctype_features = @{$regfeat_adaptor->fetch_all_by_Slice_FeatureSets($slice, [$ctype_fset1, $ctype_fset2])};

#Fetch all cell type RegulatoryFeatures for a given stable ID
my @ctype_features = @{$regfeat_adaptor->fetch_all_by_stable_ID('ENSR00001348194')};


=head1 DESCRIPTION

The RegulatoryFeatureAdaptor is a database adaptor for storing and retrieving
RegulatoryFeature objects. The FeatureSet class provides convenient wrapper
methods to the Slice functionality within this adaptor.

=cut


package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor;#DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor);


my %valid_attribute_features = (
								'Bio::EnsEMBL::Funcgen::MotifFeature' => 'motif',
								'Bio::EnsEMBL::Funcgen::AnnotatedFeature' => 'annotated',
							   );



sub fetch_MultiCell_by_stable_ID_Slice {
  my ($self, $stable_id, $slice) = @_;
  my  $fset = $self->_get_current_FeatureSet;

  throw('You must provide a stable_id argument') if ! $stable_id;
  $stable_id =~ s/[A-Z0]+//;

  $self->bind_param_generic_fetch($stable_id,  SQL_INTEGER);
  $self->bind_param_generic_fetch($fset->dbID, SQL_INTEGER);
  return $self->fetch_all_by_Slice_constraint($slice, 'rf.stable_id=? and rf.feature_set_id=?')->[0];
}



=head2 fetch_all

  Example    : my $rfs = $rf_adaptor->fetch_all;
  Description: Over-ride generic fetch_all method to return only MultiCell features.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::RegulatoryFeature objects
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub fetch_all{
  my $self = $_[0];

  return $self->SUPER::generic_fetch(' rf.feature_set_id='.
                                     $self->_get_current_FeatureSet->dbID);
}


=head2 _get_current_FeatureSet

  Example    : my $regf_featureset = $self->_get_current_FeatureSet;
  Description: Convenience method to get and test the current
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws is FeatureSet is not available
  Caller     : general
  Status     : at risk - change to _get_core_FeatureSet?

=cut

sub _get_current_FeatureSet{
  my $self = shift;


  if(! $self->{'multicell_set'}){
	$self->{'multicell_set'} = $self->db->get_FeatureSetAdaptor->fetch_by_name('RegulatoryFeatures:MultiCell');

	if(! $self->{'multicell_set'}){
	  warn('Could not retrieve current default RegulatoryFeatures:MultiCell FeatureSet');
	}
  }

  return  $self->{'multicell_set'};
}


=head2 fetch_by_stable_id

  Arg [1]    : String $stable_id - The stable id of the regulatory feature to retrieve
  Arg [2]    : optional - Bio::EnsEMBL::FeatureSet
  Example    : my $rf = $rf_adaptor->fetch_by_stable_id('ENSR00000309301');
  Description: Retrieves a regulatory feature via its stable id.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id, $fset) = @_;

  $fset ||= $self->_get_current_FeatureSet;


  #remove this ternary operatory when changes to fetch_all_by_stable_id_FeatureSets done?
  #No, this should not die when no data is present!
  return (defined $fset) ? $self->fetch_all_by_stable_id_FeatureSets($stable_id, $fset)->[0] : undef;
}

=head2 fetch_all_by_stable_id_FeatureSets

  Arg [1]    : String $stable_id - The stable id of the regulatory feature to retrieve
  Arg [2]    : optional list of FeatureSets
  Example    : my $rf = $rf_adaptor->fetch_by_stable_id('ENSR00000309301');
  Description: Retrieves a regulatory feature via its stable id.
  Returntype : Array ref of Bio::EnsEMBL::Funcgen::RegulatoryFeature objects
  Exceptions : throws if no stable ID provided or FeatureSets aren't valid
               warns if not FeatureSets defined
  Caller     : general
  Status     : at risk

=cut


#change this to fetch_all_by_stable_id and remove method of same name below or vice versa?
#check usage of this method first


sub fetch_all_by_stable_id_FeatureSets {
  my ($self, $stable_id, @fsets) = @_;

  #Change this to arrayref of fsets

  #Standard implementation exposes logic name as a parameter
  #But it will always be RegulatoryFeature/Build

  throw('Must provide a stable ID') if ! defined $stable_id;

  $stable_id =~ s/ENS[A-Z]*R0*//;


  #Need to test stable_id here as there is a chance that this argument has been omitted and we are dealing with
  #a feature set object
  $self->bind_param_generic_fetch($stable_id, SQL_INTEGER);
  my $constraint = 'rf.stable_id=?';


  #Change this to use _generate_feature_set_id_clause

  if(@fsets){

	#need to catch empty array and invalid FeatureSets
	if(scalar(@fsets == 0)){
	  warning("You have not specified any FeatureSets to fetch the RegulatoryFeature from, defaulting to all");
	}
	else{

	  #validate FeatureSets
	  #Need to check $fset->feature_class eq 'regulatory' too?
	  map { $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $_)} @fsets;

	  if(scalar(@fsets) == 1){
      $constraint .= ' and rf.feature_set_id=?';
      $self->bind_param_generic_fetch($fsets[0]->dbID, SQL_INTEGER);
	  }else{
      #How can we bind param this?

      my @bind_slots;

      foreach my $dbid(map $_->dbID, @fsets){
        push @bind_slots, '?';
        $self->bind_param_generic_fetch($dbid, SQL_INTEGER);
      }

      $constraint .= ' AND rf.feature_set_id IN ('.join(', ', @bind_slots).')';
	  }
	}
}

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
  return (
		  [ 'regulatory_feature',   'rf' ],
		  [ 'feature_set',          'fs' ],
		  [ 'regulatory_attribute', 'ra' ],
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
  my $self = shift;

  return qw(
			rf.regulatory_feature_id 
      rf.seq_region_id
			rf.seq_region_start      
      rf.seq_region_end
			rf.seq_region_strand     
      rf.bound_start_length
			rf.bound_end_length      
      rf.display_label
			rf.feature_type_id       
      rf.feature_set_id
			rf.stable_id             
      rf.binary_string
			rf.projected             
      rf.has_evidence
			rf.cell_type_count       
      ra.attribute_feature_id
			ra.attribute_feature_table
	   );
}
#rf.bound_seq_region_start			rf.bound_seq_region_end

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

sub _left_join {
  my $self = shift;

  return (['regulatory_attribute', 'rf.regulatory_feature_id = ra.regulatory_feature_id']);
}



=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates RegulatoryFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
  #So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.
  # This code is ugly because caching is used to improve speed

  my ($sa, $reg_feat);#, $old_cs_id);
  $sa = ($dest_slice) ? $dest_slice->adaptor->db->get_SliceAdaptor() : $self->db->get_SliceAdaptor();

  #Some of this in now probably overkill as we'll always be using the DNADB as the slice DB
  #Hence it should always be on the same coord system
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
  my $fset_adaptor = $self->db->get_FeatureSetAdaptor();
  my (@features, $seq_region_id);
  my (%fset_hash, %slice_hash, %sr_name_hash, %sr_cs_hash, %ftype_hash);
  my $skip_feature = 0;

  my %feature_adaptors =
    (
     'annotated' => $self->db->get_AnnotatedFeatureAdaptor,
     'motif'     => $self->db->get_MotifFeatureAdaptor,
     #external
    );

  my $stable_id_prefix = $self->stable_id_prefix;

	my (
	    $dbID,                  
      $efg_seq_region_id,
	    $seq_region_start,      
      $seq_region_end,
	    $seq_region_strand,     
      $bound_start_length,
      $bound_end_length,  
      $display_label,
      $ftype_id,              
      $fset_id,
      $stable_id,             
      $attr_id,
      $attr_type,             
      $bin_string,
      $projected,             
      $has_evidence,
      $cell_type_count
     );

	$sth->bind_columns
    (
     \$dbID,              
     \$efg_seq_region_id,
     \$seq_region_start,  
     \$seq_region_end,
     \$seq_region_strand, 
     \$bound_start_length,
     \$bound_end_length,  
     \$display_label,
     \$ftype_id,          
     \$fset_id,
     \$stable_id,         
     \$bin_string,
     \$projected,         
     \$has_evidence,
     \$cell_type_count,
     \$attr_id,
     \$attr_type,
    );

	my ($asm_cs, $cmp_cs, $asm_cs_name);
	my ($asm_cs_vers, $cmp_cs_name, $cmp_cs_vers);

	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
  }

	my ($dest_slice_start, $dest_slice_end);
	my ($dest_slice_strand, $dest_slice_length, $dest_slice_sr_name);

	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}

  my $slice;
  my %reg_attrs = (
				   annotated => {},
				   motif     => {},
				   #external
				  );

  #Set 'unique' set of feature_set_ids
  my @fset_ids;

   if($self->{params_hash}{unique}){
	#FeatureSet have been pre-validated in the fetch method
	@fset_ids = map $_->dbID, @{$self->{params_hash}{feature_sets}};
  }


  my $skip_stable_id    = 0;#stable IDs are never 0
  my $no_skip_stable_id = 0;
  my @other_rf_ids;

 FEATURE: while ( $sth->fetch() ) {

	  #Handle non-unique skipping first

    if ( $stable_id &&
         ($skip_stable_id == $stable_id) ) {
      #Faster for queries which need to skip if we have this first
      next;
	  } elsif (@fset_ids) {
      #have no_skip_stable_id too
      #so we don't keep doing _fetch_other_feature_set_ids_by_stable_feature_set_ids
      #for ID s we have already checked

      if ($no_skip_stable_id != $stable_id) {
        @other_rf_ids = @{$self->_fetch_other_dbIDs_by_stable_feature_set_ids
                            ($stable_id,
                             \@fset_ids,
                             {
                              include_projected => $self->{params_hash}{include_projected}})};

        if (@other_rf_ids) {
          $skip_stable_id = $stable_id;
          #warn "skipping\n";
          next;
        } else {
          $no_skip_stable_id = $stable_id;
        }
      }
      #else don't skip this stable ID
	  }



	  if (! $reg_feat || ($reg_feat->dbID != $dbID)) {

      if ($skip_feature) {
        undef $reg_feat;        #so we don't duplicate the push for the feature previous to the skip feature
        $skip_feature = 0;
      }

      if ($reg_feat) {          #Set the previous attr cache and reset
        $reg_feat->attribute_cache(\%reg_attrs);
        push @features, $reg_feat;


        %reg_attrs = (
                      annotated => {},
                      motif     => {},
                      #external
                     );
      }

	    #Would need to build a slice adaptor cache here to enable mapping between assemblies
	    #Or if mapping between cs systems for a given schema_build
      #which would have to be handled by the core api

      #get core seq_region_id
      $seq_region_id = $self->get_core_seq_region_id($efg_seq_region_id);

      if (! $seq_region_id) {
        warn "Cannot get slice for eFG seq_region_id $efg_seq_region_id\n".
          "The region you are using is not present in the current dna DB";
        next;
      }

	    #if($old_cs_id && ($old_cs_id+ != $cs_id)){
	    #  throw("More than one coord_system for feature query, need to implement SliceAdaptor hash?");
	    #}
	    #$old_cs_id = $cs_id;
	    #Need to make sure we are restricting calls to Experiment and channel(i.e. the same coord_system_id)

      #Get the FeatureSet object
      $fset_hash{$fset_id} = $fset_adaptor->fetch_by_dbID($fset_id) if(! exists $fset_hash{$fset_id});

      #Get the FeatureType object
      $ftype_hash{$ftype_id} = $ft_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});

	    # Get the slice object
	    $slice = $slice_hash{'ID:'.$seq_region_id};

	    if (!$slice) {
	      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
	      $slice_hash{'ID:'.$seq_region_id} = $slice;
	      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
	      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
	    }

	    my $sr_name = $sr_name_hash{$seq_region_id};
	    my $sr_cs   = $sr_cs_hash{$seq_region_id};

	    # Remap the feature coordinates to another coord system if a mapper was provided
	    if ($mapper) {

	      throw("Not yet implmented mapper, check equals are Funcgen calls too!");

	      ($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand)
          = $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);

	      # Skip features that map to gaps or coord system boundaries
        if (! defined $sr_name) {
          $skip_feature = 1;
          next FEATURE;
        }

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
          #removed bound adjusts as this is now done dynamically

          if ($dest_slice_strand == 1) {
            $seq_region_start       = $seq_region_start - $dest_slice_start + 1;
            $seq_region_end         = $seq_region_end   - $dest_slice_start + 1;
          }
          else {
            my $tmp_seq_region_start       = $seq_region_start;
            $seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
            $seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
            $seq_region_strand      *= -1;
          }
	      }

	      # Throw away features off the end of the requested slice
        #Could account for bounds here. Currently this means we never
        #get just bounds in the region in detail
        #This would reintroduce calc here

	      if ($seq_region_end < 1 || $seq_region_start > $dest_slice_length
            || ( $dest_slice_sr_name ne $sr_name )) {
          $skip_feature = 1;
          next FEATURE;
        }

	      $slice = $dest_slice;
	    }


      #This stops an init warning when sid is absent for sid mapping
      my $sid = (defined $stable_id) ? sprintf($stable_id_prefix."%011d", $stable_id) : undef;

      $reg_feat = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new_fast
        ({
          'start'          => $seq_region_start,
          'end'            => $seq_region_end,
          '_bound_lengths' => [$bound_start_length, $bound_end_length],
          'strand'         => $seq_region_strand,
          'slice'          => $slice,
          'analysis'       => $fset_hash{$fset_id}->analysis(),
          'adaptor'        => $self,
          'dbID'           => $dbID,
          'display_label'  => $display_label,
          'binary_string'  => $bin_string,
          'projected'      => $projected,
          'set'            => $fset_hash{$fset_id},
          'feature_type'   => $ftype_hash{$ftype_id},
          'stable_id'      => $sid,
          'has_evidence'   => $has_evidence,
          'cell_type_count'=> $cell_type_count,
         });

	  }


	  #populate attributes cache
	  if (defined $attr_id  && ! $skip_feature) {

      $reg_attrs{$attr_type}->{$attr_id} = undef;
	  }
	}

  #handle last record
  if ($reg_feat) {
    $reg_feat->attribute_cache(\%reg_attrs);
    push @features, $reg_feat;
  }


  #reset params hash
  $self->{params_hash} = undef;

  return \@features;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::RegulatoryFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given RegulatoryFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Listref of stored RegulatoryFeatures
  Exceptions : Throws if a list of RegulatoryFeature objects is not provided or if
               the Analysis, CellType and FeatureType objects are not attached or stored.
               Throws if analysis of set and feature do not match
               Warns if RegulatoryFeature already stored in DB and skips store.
  Caller     : General
  Status     : At Risk

=cut

sub store{
  my ($self, @rfs) = @_;

  if (scalar(@rfs) == 0) {
	throw('Must call store with a list of RegulatoryFeature objects');
  }

  my $sth = $self->prepare("
		INSERT INTO regulatory_feature (
			seq_region_id,         
      seq_region_start,
			seq_region_end,        
      bound_start_length,
			bound_end_length,      
      seq_region_strand,
      display_label,         
      feature_type_id,
      feature_set_id,        
      stable_id,
      binary_string,         
      projected,
      has_evidence,          
      cell_type_count
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

 #bound_seq_region_start,	bound_seq_region_end


  my $sth2 = $self->prepare("
		INSERT INTO regulatory_attribute (
              regulatory_feature_id, attribute_feature_id, attribute_feature_table
		) VALUES (?, ?, ?)");

  my $db = $self->db();

  foreach my $rf (@rfs) {

	if( ! ref $rf || ! $rf->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature') ) {
	  throw('Feature must be an RegulatoryFeature object');
	}

	if ( $rf->is_stored($db) ) {
	  #does not accomodate adding Feature to >1 feature_set
	  warning('Skipping RegulatoryFeature [' . $rf->dbID() . '] as it is already stored in the database');
	  next;
	}


	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $rf->feature_set);


	my ($sid, $seq_region_id);
	($rf, $seq_region_id) = $self->_pre_store($rf);
	$rf->adaptor($self);#Set adaptor first to allow attr feature retreival for bounds
	#This is only required when storing


	#Actually never happens, as we always assign stable_ids after storing
	($sid = $rf->stable_id) =~ s/ENS[A-Z]*R0*// if defined $rf->stable_id;

	$sth->bind_param(1,  $seq_region_id,           SQL_INTEGER);
	$sth->bind_param(2,  $rf->start,               SQL_INTEGER);
	$sth->bind_param(3,  $rf->end,                 SQL_INTEGER);
	$sth->bind_param(4,  $rf->bound_start_length,  SQL_INTEGER);
	$sth->bind_param(5,  $rf->bound_end_length,    SQL_INTEGER);
	$sth->bind_param(6,  $rf->strand,              SQL_TINYINT);
	$sth->bind_param(7,  $rf->{display_label},     SQL_VARCHAR);#Deref so we don't store API default value
	$sth->bind_param(8,  $rf->feature_type->dbID,  SQL_INTEGER);
	$sth->bind_param(9,  $rf->feature_set->dbID,   SQL_INTEGER);
	$sth->bind_param(10, $sid,                     SQL_INTEGER);
	$sth->bind_param(11, $rf->binary_string,       SQL_VARCHAR);
  $sth->bind_param(12, $rf->is_projected,        SQL_BOOLEAN);
  $sth->bind_param(13, $rf->has_evidence,        SQL_BOOLEAN);
	$sth->bind_param(14, $rf->cell_type_count,     SQL_INTEGER);

	#Store and set dbID
	$sth->execute;
	$rf->dbID( $self->last_insert_id );


	#Store regulatory_attributes
	#Attr cache now only holds dbids not objs

	my %attrs = %{$rf->attribute_cache};

	foreach my $fclass(keys %attrs){

	  foreach my $attr_id(keys %{$attrs{$fclass}}){
		$sth2->bind_param(1, $rf->dbID, SQL_INTEGER);
		$sth2->bind_param(2, $attr_id,  SQL_INTEGER);
		$sth2->bind_param(3, $fclass,   SQL_VARCHAR);
		$sth2->execute();
	  }
	}
  }

  return \@rfs;
}





=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $regf_adaptor->fetch_all_by_Slice($slice);
  Description: Retrieves a list of features on a given slice, specific for the current
               default RegulatoryFeature set.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $fset) = @_;

  $fset ||= $self->_get_current_FeatureSet; #This get the MultiCell sets

  #Ternary operator essential here, as we don't want to die if there is no data!
  return (defined $fset) ? $self->fetch_all_by_Slice_FeatureSets($slice, [$fset]) : undef;
}


=head2 fetch_all_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Arg [3]    : optional HASHREF - params:
                   {
                    unique            => 0|1, #Get RegulatoryFeatures unique to these FeatureSets
                    include_projected => 0|1, #Consider projected features
                   }
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureSets($slice, \@fsets);
  Description: Simple wrapper to set unique flag.
               Retrieves a list of features on a given slice, specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : Throws if params hash is not valid or keys are not recognised
  Caller     : General
  Status     : At Risk

=cut

#These FeatureSets are not optional
#otherwise could have just re-implemented fetch_all_by_Slice
#with extra optional fsets args

#To implement unique flag, would either need to self join on stable ID to compare counts
#or set flag(would have to be array of fset IDs if we want more than one fset)
#to do look up in objs_from_sth method
#Would need to add wrapper for fetch_all_by_Slice_FeatureSets here

#This will slow down all regfeat fetches
#how can we do this without impacting obj_from_sth performance?
#is this worth worrying about? It's a fairly fast track anyway?
#could have method vars/code refs in obj_from_sth, but still adding method call instead of var test.
#impact likely negligable

#Could pre_store this in DB? 'Impossible' for 1 fset without doing self join comparison?

#We really need to expose the constraint here to enable more complex combined queries
#i.e. projected, FeatureType filter etc

#Need to do full re-implementation here rather than wrapper as we can't
#pass a preformed constraint to a method which might be exposed directly by the website
#as this may enable someone to inject SQL via a URL.
#Also can't handle unique/include_projected in generic SetFeatureAdaptor method

sub fetch_all_by_Slice_FeatureSets {
  my ($self, $slice, $fsets, $params_hash) = @_;

  if($params_hash){

	if(ref($params_hash) eq 'HASH'){
	  #Really only need feature_sets as unique is implicit at present
	  #define unique for clarify
	  $self->{params_hash}{unique}            = $params_hash->{unique};
	  $self->{params_hash}{feature_sets}      = $fsets;
	  $self->{params_hash}{include_projected} = $params_hash->{include_projected};
	  #include_projected as we also need to constrain
	  #_fetch_other_dbIDs_by_stable_feature_set_ids
	  #Waht happens when we have include_projected on it's own and set to 0
	}
	else{
	  throw("The params_hash argument must be a valid HASHREF, not:\t".ref($params_hash));
	}
  }


  my $constraint = 'rf.feature_set_id '.$self->_generate_feature_set_id_clause($fsets);
  $params_hash ||= {}; #To avoid deref fail below
  my $inc_proj = $params_hash->{include_projected};

  if( ($params_hash->{unique} && ( ! $inc_proj))  ||
	(defined $inc_proj && ($inc_proj == 0) )){
	$constraint .= ' AND '. ' rf.projected=0 ';
  }

  #explicit super call, just in case we ever re-implement in here
  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}




sub _fetch_other_dbIDs_by_stable_feature_set_ids{
  my ($self, $stable_id_int, $fset_ids, $params_hash) = @_;
  #Args and originating objects have been prevalidated

  #($stable_id = $stable_id) =~ s/^[A-Z0]+//;
  #Do this in caller, as we already have the stripped ID from _objs_from_sth

  my @fset_ids = @$fset_ids; #Deref here as we are pushing, and don't want modify in the caller
  my $projected_constraint = '';

  #This is internal, so we can assume $param_hash is a valid
  #HASHREF if defined.
  $params_hash ||= {};#quick way to prevent deref fail below

  if(! $params_hash->{include_projected}){
	$projected_constraint = ' AND projected=0  ';
  }

  #Handle MultiCell set, as this will always be present
  if(! $params_hash->{include_multicell}){
	push @fset_ids, $self->_get_current_FeatureSet->dbID;
  }


  my @other_rf_ids = @{$self->db->dbc->db_handle->selectcol_arrayref
						 ('SELECT regulatory_feature_id from regulatory_feature '.
						  "WHERE stable_id=${stable_id_int} ".$projected_constraint.
						  ' AND feature_set_id not in('.join(',', @fset_ids).')')};

  return \@other_rf_ids;
}




=head2 fetch_all_by_stable_ID

  Arg [1]    : string - stable ID e.g. ENSR00000000002
  Example    : my @cell_type_regfs = @{$regf_adaptor->fetch_all_by_stable_ID('ENSR00000000001');
  Description: Retrieves a list of RegulatoryFeatures with associated stable ID. One for each CellType or
               'core' RegulatoryFeature set which contains the specified stable ID.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_stable_ID {
  my ($self, $stable_id) = @_;

  #Add fsets here?

  throw('You must provide a stable_id argument') if ! $stable_id;
  $stable_id =~ s/[A-Z0]+//;

  $self->bind_param_generic_fetch($stable_id, SQL_INTEGER);
  return $self->generic_fetch('rf.stable_id=?');
}

=head2 fetch_all_by_attribute_feature

  Arg [1]    : Bio::Ensembl::Funcgen::AnnotatedFeature or MotifFeature
  Example    : my @regfs = @{$regf_adaptor->fetch_all_by_attribute_feature($motif_feature)};
  Description: Retrieves a list of RegulatoryFeatures which contain the givven attribute feature.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : Throws is argument not valid
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_attribute_feature {
  my ($self, $attr_feat) = @_;

  #add fsets here as optional arg

  my $attr_class = ref($attr_feat);

  if(! $valid_attribute_features{$attr_class}){
	#This assigns null to the hash value but we throw straight away
	throw("Attribute feature must be one of:\n\t".join("\n\t", keys(%valid_attribute_features)));
  }

  $self->db->is_stored_and_valid($attr_class, $attr_feat);
  my $attr_feat_table = $valid_attribute_features{$attr_class};


  #Don't retrict via existing left join as we want to get all
  #the reg_attrs not just those define by this query

  #Was originally doing a subselect, but this was doing a filesort on ALL rf with no key!
  #Separating the queries makes this a range query and uses the primary key
  #still files sort, but just on exact number of rows rather than ALL( I guess because it can't do it in the buffer for some reason)


  my ($rf_ids) = $self->db->dbc->db_handle->selectrow_array("SELECT group_concat(regulatory_feature_id) from regulatory_attribute ".
															 "WHERE attribute_feature_table='${attr_feat_table}' and attribute_feature_id=".$attr_feat->dbID);

  return (defined $rf_ids) ? $self->generic_fetch("rf.regulatory_feature_id in(${rf_ids})") : [];
}



=head2 fetch_type_config_by_RegulatoryFeatures

  Arg [1]    : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Example    : my $config = $regf_adaptor->fetch_type_config_by_RegualtoryFeature($rf);
  Description: Retrieves a config hash of CellType and FeatureType names and dbIDs supporting
               the given RegualtoryFeature
  Returntype : HASHREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub fetch_type_config_by_RegulatoryFeature{
  my ($self, $rf) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::RegulatoryFeature', $rf);


  my $sql = 'SELECT ft.name, ft.feature_type_id from '.
	'feature_type ft, feature_set fs, regulatory_attribute ra, annotated_feature af '.
	  'WHERE ft.feature_type_id=fs.feature_type_id AND fs.feature_set_id=af.feature_set_id AND '.
		'af.annotated_feature_id=ra.attribute_feature_id and ra.attribute_feature_table="annotated" AND '.
		  'ra.regulatory_feature_id=? group by ft.name order by ft.name';

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $rf->dbID, SQL_INTEGER);
  $sth->execute;
  my @ftype_config = @{$sth->fetchall_arrayref};
  $sth->finish;


  #Don't need cell type query here if we have a cell type sepcific set
  #What is quicker here?

  $sql = 'SELECT ct.name, ct.cell_type_id from '.
	'cell_type ct, feature_set fs, regulatory_attribute ra, annotated_feature af '.
	  'WHERE ct.cell_type_id=fs.cell_type_id AND fs.feature_set_id=af.feature_set_id AND '.
		'af.annotated_feature_id=ra.attribute_feature_id and ra.attribute_feature_table="annotated" AND '.
		  'ra.regulatory_feature_id=? group by ct.name order by ct.name';

  $sth = $self->prepare($sql);
  $sth->bind_param(1, $rf->dbID, SQL_INTEGER);
  $sth->execute;
  my @ctype_config = @{$sth->fetchall_arrayref};
  $sth->finish;

  return {
		  feature_types => \@ftype_config,
		  cell_types    => \@ctype_config,
		 };
}


1;
