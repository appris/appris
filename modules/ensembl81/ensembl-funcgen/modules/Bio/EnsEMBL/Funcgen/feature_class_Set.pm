#
# Ensembl module for Bio::EnsEMBL::Funcgen::feature_class_Set
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

Bio::EnsEMBL::Funcgen::feature_class_Set - A module to add the feature_class role to a Set object.

=head1 SYNOPSIS

  use base qw( Bio::EnsEMBL::Funcgen::Set Bio::EnsEMBL::Funcgen::feature_class_Set )

  sub new {
    my $caller = shift;

    my $class = ref($caller) || $caller;

    my $self = $class->SUPER::new(@_);

    ...
    
    $self->_validate_feature_class(\@_);

    return $self;
  }

=head1 DESCRIPTION

This is an inheritance/mixin style 'role' class to handle the optional Set 
feature_class attribute.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::FeatureSet
Bio::EnsEMBL::Funcgen::ResultSet

=cut


package Bio::EnsEMBL::Funcgen::feature_class_Set;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );


=head2 _validate_feature_class

  Arg[0]     : Arrayref - Contructor arg
  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor for Set objects.
  Returntype : String - feature_class value
  Exceptions : Throws if feature_class now defined
  Caller     : Bio::Ensembl::Funcgen::Set sub class contructors
  Status     : At risk

=cut

sub _validate_feature_class {
  my $self = shift;
  my $args = shift;
  my ($fclass) = rearrange(['FEATURE_CLASS'], @$args);
  throw('Need to specify a -feature_class') if ! defined $fclass;
  $fclass = lc($fclass); 
  
    
  if(! $self->can('_valid_feature_classes')){
    throw(ucfirst($self->set_type).'Set inherits from feature_class_Set but'.
      ' does not have the mandatory _valid_feature_classes method');  
  }
  elsif(! grep(/^${fclass}$/, $self->_valid_feature_classes)){
    throw( "$fclass is not a valid ".ucfirst($self->set_type).
      "Set feature_class, valid feature classes are:\n\t".
      join("\t", $self->_valid_feature_classes ) );
  }
  
  $self->{feature_class} = $fclass;
  return $fclass;
}


=head2 feature_class

  Arg[0]     : String - feature class e.g. result, annotated, regulatory, external, dna_methylation or segmentation
  Example    : my $fclass = $set->feature_class;
  Description: Getter for the feature_type for this Set.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_class { return shift->{feature_class}; }


=head2 feature_class_name

  Example    : my $fclass_adaptor_method = 'get_'.$set->feature_class.'Adaptor';
  Description: Getter for the full feature class name for this Set e.g. AnnotatedFeature, RegulatoryFeature
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

#Test for adaptor
#Can't set this in new as we won't pass the adaptor
#unless creating from _obj_from_sth

sub feature_class_name{
  my $self = shift;

  if(! defined $self->{feature_class_name} ){
    $self->{feature_class_name} = $self->adaptor->build_feature_class_name($self->feature_class);
  }

  return $self->{feature_class_name};
}


1;

