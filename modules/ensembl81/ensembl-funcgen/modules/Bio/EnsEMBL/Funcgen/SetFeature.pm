#
# Ensembl module for Bio::EnsEMBL::Funcgen::SetFeature
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

Bio::EnsEMBL::Funcgen::SetFeature - Base class for features of a Set.

=head1 SYNOPSIS

  # Would normally be created from an in inheriting class e.g. AnnotatedFeature.pm

  use base qw(Bio::Ensembl::Funcgen::SetFeature);

  sub new {
    my $caller = shift;
    my $class  = ref($caller) || $caller;
    my $self   = $class->SUPER::new(@_);
    # More construction here
  }

  # Alternative direct contruction

  my $feat = Bio::EnsEMBL::Funcgen::SetFeature->
    (
     -start         => 100,
     -end           => 220,
     -strand        => -1,
     -slice         => $slice,
     -set           => $fset,
     -feature_type  => $ftype,
     -display_label => $label,
    );

  # Accessing some attributes

 
  my $start     = $feat->start;
  my $end       = $feat->end;
  my $strand    = $feat->strand;
  my $fset      = $feat->set;
  my $cell_type = $feat->cell_type;

  # Printing some information

  print $feature->display_label.' has the FeatureType '.$feat->feature_type->name."\n";

=head1 DESCRIPTION

This is a base class for features which are contained within a Funcgen FeatureSet or ResultSet.
It provides generic methods for attributes which are common across all inheriting classes.

=cut


package Bio::EnsEMBL::Funcgen::SetFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw( Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable );


=head2 new


  Arg [-SET]          : Bio::EnsEMBL::Funcgen::ResultSet or FeatureSet.
  Arg [-DISPLAY_LABEL]: (optional) String - Display label for this feature
  Arg [-FEATURE_TYPE] : (optional) Bio::EnsEMBL::Funcgen::FeatureType.
                        Defaults to Feature/ResultSet FeatureType.
  
  #Bio::EnsEMBL::Feature arguments
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-STRAND]       : (optional) Int - The orientation of this feature relative to the 
                        strand it is on. Valid values are 1, -1 and 0.
  Arg [-START]        : Int - The start coordinate of this feature relative to the start of the slice
		                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : Int -The end coordinate of this feature relative to the start of the slice
	                      it is sitting on. Coordinates start at 1 and are inclusive. 
  Arg [-ANALYSIS]     : (optional) Bio::EnsEMBL::Analysis. Defaults to Feature/ResultSet Analysis.
  Arg [-dbID]         : (optional) Int - Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor



  Example             : my $feature = Bio::EnsEMBL::Funcgen::SetFeature->new
                                        (
                                         -SLICE         => $chr_1_slice,
                                         -START         => 1000000,
                                         -END           => 1000024,
                                         -STRAND        => -1,
							                           -DISPLAY_LABEL => $text,
							                           -SET           => $fset,
                                        );

  Description: Constructor for SetFeature objects. Should never be called directly, only by its children.
  Returntype : Bio::EnsEMBL::Funcgen::SetFeature
  Exceptions : Throws if no valid ResultSet or FeatureSet passed
               Throws if FeatureType is passed but not valid
  Caller     : General
  Status     : At Risk - FEATURE_SET arg to be removed, superceded by SET in v67

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($display_label, $fset, $ftype, $set)
    = rearrange(['DISPLAY_LABEL', 'FEATURE_SET', 'FEATURE_TYPE', 'SET'], @_);

  $set ||= $fset;

  if( ( ref($set) ne 'Bio::EnsEMBL::Funcgen::FeatureSet') &&
	  ( ref($set) ne 'Bio::EnsEMBL::Funcgen::ResultSet') ){
    throw("Must pass valid Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet object\n@_");
  }

  #Grab FeatureSet first so we can pass analysis to base Feature class
  #Funcgen analysis is currently always at the Set level
  #if this ever changes the SetFeature->analysis method will also need changing
  my $self = $class->SUPER::new(@_, -analysis => $set->analysis);
 
  if($ftype){
	
    if (ref($ftype) ne 'Bio::EnsEMBL::Funcgen::FeatureType') {
      throw("feature_type param must be a valid Bio::EnsEMBL::Funcgen::FeatureType\n@_");
    }
  
    $self->{feature_type} = $ftype;
  }
 
  #Setting attrs directly removes the need for setter code in methods
  $self->{set}           = $set;
  $self->{display_label} = $display_label if defined $display_label;
 	
  return $self;
}



=head2 feature_set

  Example    : my $set = $efeature->feature_set();
  Description: WARNING: Can now also return ResultSet aswell as FeatureSet attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet
  Exceptions : None
  Caller     : General
  Status     : At Risk - marked as to be removed in v67

=cut

sub feature_set {
  #??? deprecate
  #check webcode?

  return shift->{set};
}


=head2 set

  Example    : my $set = $set_feature->set();
  Description: Getter for the set attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub set {
  return shift->{set};
}

=head2 cell_type

  Example    : my $cell_name = $set_feature->cell_type->name;
  Description: Getter for the CellType attribute for the Set of this Feature.
  May not always be for some Set types e.g. ExternalFeatures.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : stable

=cut

sub cell_type{
	return shift->set->cell_type;
}


=head2 feature_type

  Example    : my $ft_name = $set_feature->feature_type->name;
  Description: Getter for the FeatureType attribute for this feature.
               If not explicitly set, defaults to the Set FeatureType
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : stable

=cut

sub feature_type{
  my $self = shift;

  if(! defined $self->{feature_type}){
    $self->{feature_type} = $self->set->feature_type;
  }
  
  return $self->{feature_type};
}


=head2 analysis

  Example    : my $analysis = $setfeature->analysis;
  Description: Getter for the Analysis attribute for this feature.
               Re-implementation of Bio::EnsEMBL::Feature->analysis.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : stable

=cut

#what about MFs? add as feature_set as MOODS/PWM analysis not represented
#This is a mandatory requirement for Bio::EnsEMBL::Feature
#Do we ever actually have analysis at the feature level?

sub analysis{
  return shift->set->analysis;
}


=head2 display_label

  Example    : my $label = $feature->display_label;
  Description: Getter for the display label of this feature.
  This will most likely be over-ridden by inheriting class
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_label{
  return shift->{display_label};
}


1;
