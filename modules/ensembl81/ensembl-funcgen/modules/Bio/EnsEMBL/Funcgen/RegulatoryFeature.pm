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

Bio::EnsEMBL::Funcgen::RegulatoryFeature

=head1 SYNOPSIS

 use Bio::EnsEMBL::Registry;
 use Bio::EnsEMBL::Funcgen::RegulatoryFeature;

 my $reg = Bio::EnsEMBL::Registry->load_adaptors_from_db
     (
      -host    => 'ensembldb.ensembl.org',
      -user    => 'anonymous'
     );

 my $regfeat_adaptor = $reg->get_adaptor($species, 'funcgen', 'RegulatoryFeature');


 ### Creating/storing a RegulatoryFeature Set ###
 my $feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
     (
      -SLICE         => $chr_1_slice,
      -START         => 1000000,
	    -END           => 1000024,
      -STRAND        => 0,
      -DISPLAY_LABEL => $text,
      -FEATURE_SET   => $fset,
      -FEATURE_TYPE  => $reg_ftype,
      -ATTRIBUTE_CACHE => \%attr_cache,
     );

 my ($stored_feat) = @{$regfeat_adaptor->store([$feature])};


 ### Fetching some RegulatoryFeatures
 my @regfeats = @{$regfeat_adaptor->fetch_all_by_Slice_FeatureSets($slice, \@reg_feature_sets)};


 ### Print the bound and core loci
 print join(' - ', ($reg_feat->bound_start,
                    $reg_feat->start,
                    $reg_feat->end,
                    $reg_feat->bound_end)."\n";


 ### Getting some supporting evidence for a RegualtoryFeatures
 my @reg_attrs = @{$reg_feat->regulatory_attributes('annotated')};


=head1 DESCRIPTION

A RegulatoryFeature object represents the output of the Ensembl RegulatoryBuild:
    http://www.ensembl.org/info/docs/funcgen/regulatory_build.html

It may comprise many histone modification, transcription factor, polymerase and open
chromatin features, which have been combined to provide a summary view and
classification of the regulatory status at a given loci.


=head1 SEE ALSO

Bio::EnsEMBL:Funcgen::DBSQL::RegulatoryFeatureAdaptor
Bio::EnsEMBL::Funcgen::SetFeature

=cut


package Bio::EnsEMBL::Funcgen::RegulatoryFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw( Bio::EnsEMBL::Funcgen::SetFeature );


=head2 new

  Arg [-SLICE]             : Bio::EnsEMBL::Slice - The slice on which this feature is located.
  Arg [-START]             : int - The start coordinate of this feature relative to the start of the slice
                             it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]               : int -The end coordinate of this feature relative to the start of the slice
                    	     it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-FEATURE_SET]       : Bio::EnsEMBL::Funcgen::FeatureSet - Regulatory Feature set
  Arg [-FEATURE_TYPE]      : Bio::EnsEMBL::Funcgen::FeatureType - Regulatory Feature sub type
  Arg [-BINARY_STRING]     : (optional) string - Regulatory Build binary string
  Arg [-STABLE_ID]         : (optional) string - Stable ID for this RegualtoryFeature e.g. ENSR00000000001
  Arg [-DISPLAY_LABEL]     : (optional) string - Display label for this feature
  Arg [-ATTRIBUTE_CACHE]   : (optional) HASHREF of feature class dbID|Object lists
  Arg [-PROJECTED]         : (optional) boolean - Flag to specify whether this feature has been projected or not
  Arg [-dbID]              : (optional) int - Internal database ID.
  Arg [-ADAPTOR]           : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Example    : my $feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
                 (
									-SLICE         => $chr_1_slice,
									-START         => 1000000,
									-END           => 1000024,
                  -DISPLAY_LABEL => $text,
									-FEATURE_SET   => $fset,
                  -FEATURE_TYPE  => $reg_ftype,
                  -ATTRIBUTE_CACHE => \%attr_cache,
                 );


  Description: Constructor for RegulatoryFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($stable_id, $attr_cache, $bin_string, $projected, $has_evidence, $cell_type_count)
    = rearrange(['STABLE_ID', 'ATTRIBUTE_CACHE', 'BINARY_STRING', 'PROJECTED', 'HAS_EVIDENCE', 'CELL_TYPE_COUNT'], @_);

  #None of these are mandatory at creation
  #under different use cases
  $self->{binary_string} = $bin_string    if defined $bin_string;
  $self->{stable_id}     = $stable_id     if defined $stable_id;
  $self->{projected}     = $projected     if defined $projected;
  $self->{has_evidence}  = $has_evidence  if defined $has_evidence;
  $self->{cell_type_count}  = $cell_type_count  if defined $cell_type_count;
  $self->attribute_cache($attr_cache)     if $attr_cache;

  return $self;
}


=head2 display_label

  Example    : my $label = $feature->display_label;
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_label {
  my $self = shift;

  if(! defined $self->{display_label}){
    $self->{'display_label'}  = $self->feature_type->name.' Regulatory Feature';

    if( defined $self->cell_type ){
      $self->{display_label} .= ' - '.$self->cell_type->name;
    }
  }

  return $self->{display_label};
}


=head2 display_id

  Example    : print $feature->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. In this case the stable Id is
               preferred
  Returntype : String
  Exceptions : none
  Caller     : web drawing code, Region Report tool
  Status     : Stable

=cut

sub display_id {  return shift->{stable_id}; }


=head2 binary_string

  Arg [1]    : optional string - binary string from regulatory build
  Example    : my $bin_string = $feature->binary_string();
  Description: Getter for the binary_string for this feature.
  Returntype : String
  Exceptions : None
  Caller     : Regulatory build analyses
  Status     : At Risk - May change to BLOB

=cut

sub binary_string{ return shift->{binary_string}; }


=head2 stable_id

  Arg [1]    : (optional) string - stable_id e.g ENSR00000000001
  Example    : my $stable_id = $feature->stable_id();
  Description: Getter and setter for the stable_id attribute for this feature. 
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At risk - setter functionality to be removed

=cut

sub stable_id {
  my $self = shift;

  if (@_){
    #added v67
    warn "RegualtoryFeature::stable_id setter functionality is being removed\n";
    #Still used in stable_id_mapper.pl
    $self->{stable_id} = shift;
  }

  return $self->{stable_id};
}


=head2 regulatory_attributes

  Arg [1]    : String (optional) - Class of feature e.g. annotated or motif
  Example    : print "Regulatory Attributes:\n\t".join("\n\t", (map $_->feature_type->name, @{$feature->regulatory_attributes()}))."\n";
  Description: Getter for the regulatory_attributes for this feature.
  Returntype : ARRAYREF
  Exceptions : Throws if feature class not valid
  Caller     : General
  Status     : At Risk

=cut


sub regulatory_attributes{
  my ($self, $feature_class) = @_;
  my @fclasses;
  my %adaptors = ('annotated' => $self->adaptor->db->get_AnnotatedFeatureAdaptor,
                  'motif'     => $self->adaptor->db->get_MotifFeatureAdaptor     );

  if (defined $feature_class) {

    if (exists $adaptors{lc($feature_class)}) {
      @fclasses = (lc($feature_class));
    } 
    else {
      throw("The feature class you specified is not valid:\t$feature_class\n".
            "Please use one of:\t".join(', ', keys %adaptors));
    }
  } 
  else {
    @fclasses = keys %adaptors;
  }

  foreach my $fclass (@fclasses) {
    # Now structured as hash to facilitate faster has_attribute method
    # Very little difference to array based cache
    my @attr_dbIDs = keys %{$self->{attribute_cache}{$fclass}};
	
    if (scalar(@attr_dbIDs) > 0) {
	  
      if ( ! ( ref($self->{regulatory_attributes}{$fclass}->[0])  &&
               ref($self->{regulatory_attributes}{$fclass}->[0])->isa('Bio::EnsEMBL::Feature') )) {
      
        $adaptors{$fclass}->force_reslice(1); #So we don't lose attrs which aren't on the slice
        # fetch_all_by_Slice_constraint does relevant normalised Slice projection i.e. PAR mappingg
        $self->{'regulatory_attributes'}{$fclass} = 
          $adaptors{$fclass}->fetch_all_by_Slice_constraint
            ($self->slice,
             lc($fclass).'_feature_id in('.join(',', @attr_dbIDs).')' );

        # Forces reslice and inclusion for attributes not contained within slice 
        $adaptors{$fclass}->force_reslice(0);
      }
    } else {
      $self->{regulatory_attributes}{$fclass} = [];
    }
  }

  return [ map { @{$self->{regulatory_attributes}{$_}} } @fclasses ];
}

=head2 has_attribute

  Arg [1]    : Attribute Feature dbID
  Arg [2]    : Attribute Feature class e.g. motif or annotated
  Example    : if($regf->has_attribute($af->dbID, 'annotated'){ #do something here }
  Description: Identifies whether this RegualtoryFeature has a given attribute
  Returntype : Boolean
  Exceptions : Throws if args are not defined
  Caller     : General
  Status     : Stable

=cut

sub has_attribute{
  my ($self, $dbID, $fclass) = @_;

  throw('Must provide a dbID and a Feature class argument') if ! $dbID && $fclass;

  return exists ${$self->attribute_cache}{$fclass}{$dbID};
}


=head2 has_evidence

  Arg [1]     : None
  Returntype  : Boolean 
  Exceptions  : None
  Description : Returns 1 if has_evidence = 1, nothing otherwise

=cut

sub has_evidence { return shift->{has_evidence}; }


=head2 cell_type_count

  Arg [1]     : None 
  Returntype  : SCALAR 
  Exceptions  : None
  Description : Returns the amount of cells where this RegFeat is active

=cut

sub cell_type_count { shift->{cell_type_count}; }


=head2 get_focus_attributes

  Arg [1]    : None
  Example    : my @focus_attrs = @{$regf->get_focus_attributes};
  Description: Getter for the focus features of this RegualtoryFeature, used to defined the core region
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_focus_attributes{
  my $self = shift;

  if(! exists $self->{'focus_attributes'} ||
	 ! @{$self->{'focus_attributes'}}){
   $self->_sort_attributes;
  }

  return $self->{'focus_attributes'};
}


=head2 get_nonfocus_attributes

  Arg [1]    : None
  Example    : my @non_focus_attrs = @{$regf->get_nonfocus_attributes};
  Description: Getter for the non-focus features of this RegulatoryFeature, used to defined 
               the non core region i.e. the whiskers.
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_nonfocus_attributes{
  my $self = shift;

  #Test focus here as we may not have any nonfocus
  #But focus will show that we have sorted already
  if(! exists $self->{'focus_attributes'} ||
	 ! @{$self->{'focus_attributes'}}){
	$self->_sort_attributes;
  }

  return $self->{'nonfocus_attributes'};
}

#Add pod here

sub _sort_attributes{
  my $self = shift;

  $self->{'focus_attributes'} = [];
  $self->{'nonfocus_attributes'} = [];

  foreach my $attrf(@{$self->regulatory_attributes}){

    if($attrf->isa('Bio::EnsEMBL::Funcgen::MotifFeature') ||
       $attrf->feature_set->is_focus_set){
      push @{$self->{'focus_attributes'}}, $attrf;
    }
    else{
      push @{$self->{'nonfocus_attributes'}}, $attrf;
    }
  }

  return;
}


=head2 attribute_cache

  Arg [1]    : optional - HASHREF of attribute table keys with values as either a list of attribute 
               feature dbIDs or objects. If passing object, any MotifFeature objects should be in position
               order with respect to the slice.
  Example    : $feature->attribute_cache(\%attribute_feature_info);
  Description: Setter for the regulatory_attribute cache for this feature. This is a short cut method used by the 
               regulatory build and the webcode to avoid unnecessary fetching and enable enable lazy loading 
  Returntype : HASHREF
  Exceptions : Throws if trying to overwrite existing cache
  Caller     : RegulatoryFeatureAdaptor.pm and build_regulatory_features.pl
  Status     : At Risk

=cut


sub attribute_cache{
  my ($self, $attr_hash) = @_;

  if(defined $attr_hash){

	foreach my $fclass(keys %{$attr_hash}){

	  if(exists $self->{'attribute_cache'}{$fclass}){
		throw("You are trying to overwrite a pre-existing regulatory attribute cache entry for feature class:\t$fclass");
	  }
	  else{
		$self->{'attribute_cache'}{$fclass} = $attr_hash->{$fclass};
	  }
	}
  }

  return $self->{'attribute_cache'} || {};
}



# The bound_seq_region_start/end methods are dynamic to support
# HAP/PAR projection. Reverting this (not advised), would require
# changes in _objs_from_sth to remap the loci to a dest slice
# These are now 'lazy calculated' rather than pre-defined. Hence no
# change in performace unless they are called multiple times
# for the same slice


=head2 bound_seq_region_start

  Example    : my $bound_sr_start = $feature->bound_seq_region_start;
  Description: Getter for the seq_region bound_start attribute for this feature.
               Gives the 5' most start value of the underlying attribute
               features.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_seq_region_start { return $_[0]->seq_region_start - $_[0]->_bound_lengths->[0]; }


=head2 bound_seq_region_end

  Example    : my $bound_sr_end = $feature->bound_seq_region_end;
  Description: Getter for the seq_region bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_seq_region_end { return $_[0]->seq_region_end + $_[0]->_bound_lengths->[1]; }

# As this 'private' method is not exposed or required to be poylymorphic,
# it would theoretically, be quicker to have this as a sub.

sub _bound_lengths {
  my $self = shift;

  if(! defined  $self->{_bound_lengths}){

    my @af_attrs = @{$self->regulatory_attributes('annotated')};
    
    if (! @af_attrs) {
      throw('Unable to set bound length, no AnnotatedFeature attributes available for RegulatoryFeature: '
            .$self->dbID);
    }
    
    #Adding self here accounts for core region i.e.
    #features extending beyond the core may be absent on this cell type.
    my @start_ends;

    foreach my $feat (@af_attrs, $self) {
      push @start_ends, ($feat->seq_region_start, $feat->seq_region_end);
    }
        
    @start_ends = sort { $a <=> $b } @start_ends;
        
    $self->{_bound_lengths} = [ ($self->seq_region_start - $start_ends[0]),
                                ($start_ends[$#start_ends] - $self->seq_region_end) ];
  }

  return $self->{_bound_lengths};
}

#This appears to be slower than foreach
#my @start_ends = map { $_->start, $_->end } (@af_attrs, $self);




# The following methods are all dynamic wrt slice strand and projection/transfer
# and return local values

=head2 bound_start_length

  Example    : my $bound_start_length = $reg_feat->bound_start_length;
  Description: Getter for the bound_start_length attribute for this feature,
               with respect to the host slice strand
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_start_length {
  my $self = shift;
  return ($self->slice->strand == 1) ? $self->_bound_lengths->[0] : $self->_bound_lengths->[1];
}


=head2 bound_end_length

  Example    : my $bound_end_length = $reg_feat->bound_end_length;
  Description: Getter for the bound_end length attribute for this feature,
               with respect to the host slice strand.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_end_length {
  my $self = shift;
  return ($self->slice->strand == 1) ? $self->_bound_lengths->[1] : $self->_bound_lengths->[0];
}


=head2 bound_start

  Example    : my $bound_start = $feature->bound_start;
  Description: Getter for the bound_start attribute for this feature.
               Gives the 5' most start value of the underlying attribute
               features in local coordinates.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_start { return $_[0]->start - $_[0]->bound_start_length; }


=head2 bound_end

  Example    : my $bound_end = $feature->bound_start();
  Description: Getter for the bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features in local coordinates.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_end { return $_[0]->end + $_[0]->bound_end_length; }


=head2 is_projected

  Arg [1]    : optional - boolean
  Example    : if($regf->is_projected){ #do something different here }
  Description: Getter/Setter for the projected attribute.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : At risk - remove setter functionality

=cut

sub is_projected {
  my $self = shift;
  
  if(@_){
	#added v67
    warn "RegulatoryFeature::is_projected setter functionality is being removed\n";
    $self->{'projected'} = shift;
  }
  
  return $self->{'projected'};
}


=head2 get_underlying_structure

  Example    : my @web_image_structure = @{$regf->get_underlying_structure};
  Description: Getter for the bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : Arrayref
  Exceptions : None
  Caller     : Webcode
  Status     : At Risk

=cut

#Could precompute these as core region loci
#and store in the DB to avoid the MF attr fetch?

#This is also sensitive to projection/transfer after we have called it.
#Would have to do one of the following
#1 Projecting all motif_features. This could be done by extending/overwriting
#  Feature::project/transfer, and make all feature projection code use that e.g. BaseFeatureAdaptor
#2 Cache the start, end and strand of slice, and update when changed by transforming motif_feature_loci

# This is only ever used for web which will never call until any projection is complete.
# Hence no real need for this to be sensitive to pre & post projection calling
# Leave for now with above useage caveat

sub get_underlying_structure{
  my $self = shift;

  if (! defined $self->{underlying_structure}) {
    my @mf_loci;
    
    foreach my $mf (@{$self->regulatory_attributes('motif')}) {
      push @mf_loci, ($mf->start, $mf->end);
    }

    $self->{underlying_structure} = [ 
                                     $self->bound_start, $self->start,
                                     @mf_loci,
                                     $self->end, $self->bound_end
                                    ];
  }

  $self->{underlying_structure};
}


=head2 is_unique_to_FeatureSets

  Arg[1]     : optional - ARRAYREF of regualtory Bio::EnsEMBL::Funcgen::FeatureSet objects
                          Default is FeatureSet of given RegulatoryFeature, else need to be 
                          defined explicitly.
  Arg[2]     : optional - HASHREF Params hash: 
                                    {
                                     include_projected => 0|1, # Boolean, include 'projected' features
                                    }
  Example    : if($reg_feat->is_unique_to_FeatureSets($fsets)}{  
                   #then do some analysis here
               }
  Description: Identifies whether this RegualtoryFeature is unique to a set of FeatureSets.
  Returntype : boolean
  Exceptions : Throw is arguments not stored or valid.
  Caller     : General
  Status     : At risk

=cut

#Probably want to add in an FeatureType constraint here
#e.g. so we can compare active vs inactive or poised promoters

#omit include_multi doesn't make sense here

sub is_unique_to_FeatureSets{
  my ($self, $fsets, $params_hash) = @_;

  $fsets ||= [$self->feature_set];
  my @fset_ids;
  
  
  #define to avoid deref fails below.
  $params_hash ||= {};
  if(ref($params_hash) ne 'HASH'){
	throw("The params hash argument must be a valid HASHREF:\t".ref($params_hash));
  }


  foreach my $fset(@$fsets){
	#assume we have an adaptor set
	$self->adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	
	if($fset->feature_class ne 'regulatory'){
	  throw('Found non-regulatory FeatureSet');
	}

	push @fset_ids, $fset->dbID;
  }
  
  my $stable_id;
  ($stable_id = $self->stable_id) =~ s/^[A-Z0]+//;

  
  my @other_rf_ids = @{$self->adaptor->_fetch_other_dbIDs_by_stable_feature_set_ids
						 ($stable_id, 
							 \@fset_ids, 
						  { include_projected => $params_hash->{include_projected}} )};
  
  return (@other_rf_ids) ? 0 : 1;
}



=head2 get_other_RegulatoryFeatures

  Arg[1]     : optional - ARRAYREF of regualtory Bio::EnsEMBL::Funcgen::FeatureSet objects
                          Default is FeatureSet of given RegulatoryFeature, else need to be 
                          defined explicitly.
  Arg[2]     : optional - HASHREF Params hash: 
                                    {
                                     include_projected => 0|1, # Boolean, include 'projected' features
                                     include_multicell => 0|1, # Boolean, include MultiCell features
                                    }
  Example    : my @other_fsets = @{$reg_feat->get_other_FeatureSets($fsets)};
  Description: Gets other RegualtoryFeatures (linked via the stable ID) which are present in the 
               specified list of FeatureSets.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::RegulatoryFeature objects
  Exceptions : Throw is arguments not stored or valid.
  Caller     : General
  Status     : At risk

=cut

sub get_other_RegulatoryFeatures{
  my ($self, $fsets, $params_hash) = @_;
  
  #define to avoid deref fails below.
  $params_hash ||= {};
  if(ref($params_hash) ne 'HASH'){
	throw("The params hash argument must be a valid HASHREF:\t".ref($params_hash));
  }

  $fsets ||= [$self->feature_set];
  my @fset_ids;
  
  foreach my $fset(@$fsets){
	#assume we have an adaptor set
	$self->adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	
	if($fset->feature_class ne 'regulatory'){
	  throw('Found non-regulatory FeatureSet');
	}

	push @fset_ids, $fset->dbID;
  }
    
  my $stable_id;
  ($stable_id = $self->stable_id) =~ s/^[A-Z0]+//;

  my @other_fsets_ids = @{$self->adaptor->_fetch_other_dbIDs_by_stable_feature_set_ids
							($stable_id, \@fset_ids, 
							 {
							  include_projected => $params_hash->{include_projected},
							  include_multicell => $params_hash->{include_multicell},
							 })};

  return $self->adaptor->fetch_all_by_dbID_list(\@other_fsets_ids);
}


=head2 summary_as_hash

  Example       : $regf_summary = $regf->summary_as_hash;
  Description   : Retrieves a textual summary of this RegulatoryFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

#why has has_evidence been transposed to activity_evidence?

sub summary_as_hash {
  my $self   = shift;

  #Use cell_type count here to detect new build
  #but this is set to 0 for cell type specific features
  #my %build_specific_params;
  #my $ctype_cnt = $self->cell_type_count;

  #if($ctype_cnt){
  #  $build_specific_params{cell_type_count}    = $ctype_cnt;
  #  $build_specific_params{activity_envidence} = $self->has_evidence;
  #}
  #else{
  #  $build_specific_params{projected} = $self->projected;
  #}

  return
   {ID                      => $self->stable_id,
    cell_type               => $self->cell_type->name,
    bound_start             => $self->bound_seq_region_start,
    bound_end               => $self->bound_seq_region_end,
    start                   => $self->seq_region_start,
    end                     => $self->seq_region_end,
    strand                  => $self->strand,
    seq_region_name         => $self->seq_region_name,
    activity_evidence       => $self->has_evidence,
    description             => $self->feature_type->description,
    feature_type            => "regulatory",
    #%build_specific_params, 
   };
}

1;


