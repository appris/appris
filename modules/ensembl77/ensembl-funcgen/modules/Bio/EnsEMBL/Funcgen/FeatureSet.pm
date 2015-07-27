#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureSet
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

Bio::EnsEMBL::Funcgen::FeatureSet - A container for a discrete set of experimental features.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureSet;

my $feature_set = Bio::EnsEMBL::Funcgen::FeatureSet->new
  (
   -dbid          => $dbid,
   -analysis      => $analysis,
   -feature_type  => $ftype,
   -cell_type     => $ctype,
   -name          => 'HeLa-S3_H3K4me3_ExperimentalGroup_Analysis',
   -feature_class => 'annotated',
   -description   => 'HeLa-S3 H3K4me3 SWEmbl peaks replicate 1',
   -display_label => 'HeLa-S3 H3K4me3',
   -input_set     => $iset,
  );


my @features = @{$feature_set->get_Features_by_Slice($slice)};

print $feature_set->display_label.' has '.
  scalar(@features).' on slice '.$slice->name."\n";


$feature_set        = $feature_set_adaptor->fetch_by_name('RegulatoryFeatures:HepG2');
my $feature_type    = $feature_type_adaptor->fetch_by_name('Promoter Associated');

#This may take a while and a lot of memory.
my @prom_like_regfs = @{$feature_set->get_all_by_FeatureType($feature_type)};

#This may take longer and more memory.
my @all_regfs       = @{$feature_set->get_all_Features()};

=head1 DESCRIPTION

A FeatureSet object defines a discrete set of features from a single analysis.
It is a container object providing access to meta data and convenience wrapper
methods to fetch the related features.

A FeatureSet is normally associated to a DataSet either as a product FeatureSet,
with the excpetion of FeatureSets that define external data. As these generally
have no supporting data, there is no need for a DataSet. FeatureSets can also be
associated with DataSets as supporting sets i.e. the input to a further analysis
such as the Ensembl Regulatory Build:
  http://www.ensembl.org/info/docs/funcgen/regulatory_build.html

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Set;
Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor;
Bio::EnsEMBL::Funcgen::DataSet;

=cut

package Bio::EnsEMBL::Funcgen::FeatureSet;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Funcgen::Set Bio::EnsEMBL::Funcgen::feature_class_Set);

=head2 new

  Args [1] : Hash of key value parameter as follows:

             MANDATORY:
               -NAME          => String - name for this Set.
               -FEATURE_TYPE  => Bio::EnsEMBL::Funcgen::FeatureType
               -FEATURE_CLASS => String - Class of feature e.g. result, annotated,
                                 regulatory, segmentation, external, mirna, or dna_methylation.
             OPTIONAL:
               -EXPERIMENT    => Bio::EnsEMBL::Funcgen::Experiment
               -EXPERIMENT_ID => Int
               -DESCRIPTION   => String
               -DISPLAY_NAME  => String
               -CELL_TYPE     => Bio::EnsEMBL::Funcgen::CellType
               -ANALYSIS      => Bio::EnsEMBL::Analysis
               -DBID          => Int
               -ADAPTOR       => Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor

  Example    : my $feature = Bio::EnsEMBL::Funcgen::FeatureSet->new
                 (
                  -dbid          => $dbid,
                  -analysis      => $analysis,
                  -feature_type  => $ftype,
                  -cell_type     => $ctype,
                  -name          => 'HeLa-S3_H3K4me3_ExperimentalGroup_Analysis',
                  -feature_class => 'annotated',
                  -description   => 'HeLa-S3 H3K4me3 SWEmbl peaks replicate 1',
                  -display_label => 'HeLa-S3 H3K4me3',
                  -experiment    => $exp,
			           );
  Description: Constructor for FeatureSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws if: FeatureType is not defined, Feature class is not valid
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = $class->SUPER::new(@_);

  my ($desc, $dlabel, $iset, $iset_id) = rearrange
   ( ['DESCRIPTION', 'DISPLAY_LABEL', 'INPUT_SET', 'INPUT_SET_ID'], @_ );


  #Mandatory params checks here (setting done in Set.pm)
  #explicit type check here to avoid invalid types being imported as NULL
  #subsequently throwing errors on retrieval
  my $type = $self->_validate_feature_class(\@_);

  if(! defined $self->cell_type){
    if( ($type eq 'annotated') ||  ($type eq 'regulatory') ||  ($type eq 'segmentation')){
      throw("FeatureSets with feature_class '$type' require a defined CellType");
    }
  }

  #Direct assignment to prevent need for set arg test in method
  $self->{description}   = $desc    if defined $desc;
  $self->{display_label} = $dlabel  if defined $dlabel;
 

  if(defined $iset || defined $iset_id){
    throw('FeatureSet not longer takes and InputSet or it\'s dbID in the constructor'.
          ', please upass an Experiment/dbID instead');
  }

  return $self;
}                               ## end sub new


sub _valid_feature_classes{
  return qw( annotated regulatory external segmentation mirna_target);
}


=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub new_fast { return bless ($_[1], $_[0]); }


=head2 description

  Example    : print "Feature set description is:\t".$fset->description."\n";
  Description: Getter for the description of this FeatureSet. e.g. Release 3.1
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub description { return shift->{description}; }



=head2 display_label

  Example    : print $rset->display_label;
  Description: Getter for the display_label attribute for this FeatureSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my $self = shift;

  if ( ! defined $self->{display_label} ) {

    if ( $self->feature_class eq 'regulatory' ) {
      $self->{display_label} = $self->name;
    } else {
      #This still fails here if we don't have a class or a cell_type set

      $self->{display_label} =
        $self->feature_type->name.' - '. $self->cell_type->name.' enriched sites';
    }
  }

  return $self->{display_label};
}




=head2 get_FeatureAdaptor

  Example    :
  Description: Retrieves and caches FeatureAdaptor of feature_set type
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_FeatureAdaptor{
  my $self = shift;

  if (! exists $self->{'adaptor_refs'}) {

    foreach my $valid_class ($self->_valid_feature_classes) {
      my $method = 'get_'.$self->adaptor->build_feature_class_name($valid_class).'Adaptor';
      $self->{'adaptor_refs'}{$valid_class} = $self->adaptor->db->$method;
    }
  }

  return $self->{'adaptor_refs'}->{$self->feature_class()};
}



=head2 get_Features_by_Slice

  Example    : my @features = @{$FeatureSet->get_Features_by_Slice($slice)};
  Description: Retrieves all Features for this FeatureSet for a given Slice
  Returntype : ARRAYREF containing Features of the feature_set type i.e. Annotated, Regulatory or Supporting;
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Features_by_Slice{
  my $self  = shift;
  my $slice = shift;

  return $self->get_FeatureAdaptor->fetch_all_by_Slice_FeatureSets($slice, [$self]);
}


=head2 get_Features_by_FeatureType

  Arg[0]     : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my @features = @{$FeatureSet->get_Features_by_FeatureType($ftype)};
  Description: Retrieves all Features for this FeatureSet for a given FeatureType
               or associated FeatureType. This is mainly used by external FeatureSets
               which can sometimes have more than one associated FeatureType.
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_Features_by_FeatureType{
  my $self = shift;
  my $type = shift;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureType_FeatureSets($type, [$self]);
}


=head2 get_all_Features

  Example    : my @features = @{$FeatureSet->get_all_Features};
  Description: Retrieves all Features for this FeatureSet
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_all_Features{
  my $self = shift;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureSets([$self]);
}


=head2 is_focus_set

  Args       : None
  Example    : if($fset->is_focus_set){ ... }
  Description: Returns true if FeatureSet is a focus set used to define
               the core region of a RegulatoryFeature
  Returntype : Boolean
  Exceptions : Throws if meta entry not present
  Caller     : General
  Status     : At Risk

=cut

#'focus' and 'core' have historically been interchageable


sub is_focus_set{
  my $self = shift;

  if (! defined $self->{focus_set}) {

    if (! defined $self->cell_type) {
      warn "FeatureSet without an associated CellType cannot be a focus set:\t".$self->name;
      $self->{focus_set} = 0;
    } else {
      $self->{focus_set} = $self->adaptor->fetch_focus_set_config_by_FeatureSet($self);
    }
  }

  return $self->{focus_set};
}


=head2 is_attribute_set

  Args       : None
  Example    : if($fset->is_attribute_set){ ... }
  Description: Returns true if FeatureSet is a supporting/attribute(focus or not) set used in the RegulatoryBuild
  Returntype : Boolean
  Exceptions : Throws if meta entry not present
  Caller     : General
  Status     : At Risk

=cut

sub is_attribute_set{
  my $self = shift;

  if (! defined $self->{attribute_set}) {

    if (! defined $self->cell_type) {
      warn "FeatureSet without an associated CellType cannot be a attribute set:\t".$self->name;
      $self->{attribute_set} = 0;
    } else {
      $self->{attribute_set} = $self->adaptor->fetch_attribute_set_config_by_FeatureSet($self);
    }
  }

  return $self->{attribute_set};
}


=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of FeatureSet method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name description display_label feature_class get_all_states
  Args[4]    : Arrayref - Optional list of FeatureSet method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: feature_type cell_type analysis experiment
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this FeatureSet to another based on the defined scalar
               and storable methods.
  Returntype : Hashref of key attribute/method name keys and values which differ.
               Keys will always be the method which has been compared.
               Values can either be a error string, a hashref of diffs from a
               nested object, or an arrayref of error strings or hashrefs where
               a particular method returns more than one object.
  Exceptions : None
  Caller     : Import/migration pipeline
  Status     : At Risk

=cut

sub compare_to {
  my ($self, $obj, $shallow, $scl_methods, $obj_methods) = @_;

  $scl_methods ||= [qw(name description display_label feature_class get_all_states)];
  $obj_methods ||= [qw(feature_type cell_type analysis experiment)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}



=head2 reset_relational_attributes

  Arg[1] : Hashref containing the following parameters.
           Mandatory:
            -analysis     => Bio::EnsEMBL::Analysis,
            -feature_type => Bio::EnsEMBL::Funcgen::FeatureType,

           Optional if not presently defined:
            -cell_type    => Bio::EnsEMBL::Funcgen::CellType,
            -experiment   => Bio::EnsEMBL::Funcgen::Experiment,

  Description: Resets all the relational attributes of a given FeatureSet.
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;
  my ($analysis, $feature_type, $cell_type, $exp) =
      rearrange(['ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE', 'EXPERIMENT'],
      %$params_hash);

  assert_ref($analysis,     'Bio::EnsEMBL::Analysis',             'Analysis');
  assert_ref($feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType', 'FeatureType');

  if( $self->cell_type && 
      ! check_ref($cell_type, 'Bio::EnsEMBL::Funcgen::CellType') ){
    throw("You must pass a valid Bio::EnsEMBL::Funcgen::CellType\n".
          "Passed:\t".ref($cell_type));
  }
  
  if( $self->experiment &&
      ! check_ref($exp, 'Bio::EnsEMBL::Funcgen::Experiment') ){
    throw("You must pass a valid Bio::EnsEMBL::Funcgen::Experiment\n".
          "Passed:\t".ref($exp));
  }
  
  if(defined $exp){
    #This will allow addition of an experiment
    #when it was not prevously defined
    $self->{experiment_id} = $exp->dbID;
    $self->{experiment}    = $exp;
  }

  $self->{cell_type}    = $cell_type;
  $self->{feature_type} = $feature_type;
  $self->{analysis}     = $analysis;


  #Finally undef the dbID and adaptor by default
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  #Reset cached dynamic attrs here, just in case there has been a change?
  return;
}




1;
