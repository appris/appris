#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSet
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

Bio::EnsEMBL::InputSet - A module to represent InputSet object.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::InputSet;

#Create an InputSet

my $inp_set = Bio::EnsEMBL::Funcgen::InputSet->new
                 (
                  -ANALYSIS     => $anal,
                  -CELL_TYPE    => $ctype,
                  -EXPERIMENT   => $exp,
                  -FEATURE_TYPE => $ftype,
                  -NAME         => 'SRR00000.fastq.gz',
                  -REPLICATE    => 1, # >0 for specific replicate or 0 for merged
                 );





=head1 DESCRIPTION

An InputSet object provides a generic container for any non-array based feature import,
allowing tracking of file import via the status table and integration into Data and FeatureSets to
provide traceability to the source experiment from a given FeatureSet.

=cut

package Bio::EnsEMBL::Funcgen::InputSet;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate);

use base qw(Bio::EnsEMBL::Funcgen::Set Bio::EnsEMBL::Funcgen::feature_class_Set);

=head2 new

  Example    : my $eset = Bio::EnsEMBL::Funcgen::InputSet->new
                 (
	                -DBID         => $dbID,
                  -ADAPTOR      => $self,
                  -ANALYSIS     => $anal,
                  -CELL_TYPE    => $ctype,
                  -EXPERIMENT   => $exp,
                  -FEATURE_TYPE => $ftype,
                  -NAME         => 'SRR00000.fastq.gz',
                  -FEATURE_CLASS => 'annotated',
                  -REPLICATE    => 1, # >0 for specific replicate or 0 for merged
                 );

  Description: Constructor for InputSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : Throws if no Experiment defined
               Throws if CellType or FeatureType are not valid or stored
  Caller     : General
  Status     : At risk

=cut

#### Methods which access subsets needs to test first if they exist, otherwise call _get_input_subsets
#### Support adding passing InputSubsets and storing them in iss and isiss
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new(@_);
  my ($exp, $rep) = rearrange(['EXPERIMENT', 'REPLICATE'], @_);

  $self->_validate_feature_class(\@_);	
  throw ('Must provide a  CellType') if ! defined $self->cell_type;

  if (! (ref $exp && $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') && $exp->dbID())){
    throw("Must specify a valid stored Bio::EnsEMBL::Funcgen::Experiment Passed: ".ref($exp));
  }

  #Set directly for speed
  $self->{replicate}  = $rep;
  $self->{experiment} = $exp;
  return $self;
}



sub _valid_feature_classes{
  return qw( annotated result segmentation dna_methylation );
}


=head2 get_Experiment

  Example    : my $exp = $inp_set->get_Experiment();
  Description: Getter for the Experiment of this InputSet.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Experiment{  return shift->{experiment}; }


=head2 get_InputSubsets

  Example    : my @subsets = @{$exp_set->get_InputSubsets()};
  Description: Getter for the InputSubsets for this InputSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_InputSubsets{
  my ($self)  = @_;

  if(! exists $self->{'subsets'}){
    my $iss_adaptor = $self->adaptor->db->get_InputSubsetAdaptor;

    my $subsets = $iss_adaptor->fetch_all_by_InputSet($self);

    for my $ss (@{$subsets}){
      if(exists $self->{'subsets'}->{$ss->name}){
        throw('Subsets with name ' . $ss->name . 'exists already');
      }

      $self->{'subsets'}->{$ss->name} = $ss;
    }
  }
  return [ values %{$self->{'subsets'}} ];
}

=head2 get_subset_by_name

  Example    : my $subsets = $exp_set->get_subset_by_name('subset1');
  Description: Getter for the subset of a given name for this InputSet.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_subset_by_name{
  my $self = shift;
  my $name = shift;
  return (exists $self->{'subsets'}{$name}) ? $self->{'subsets'}{$name} : undef;
}


=head2 get_subset_names

  Example    : my @subset_names = @{$exp_set->get_subset_names()};
  Description: Getter for the subset names for this InputSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_subset_names{
  my $self = shift;
  if(! exists $self->{'subsets'}){
    $self->get_InputSubsets;
  }
  return [ keys %{$self->{'subsets'}} ];
}



=head2 replicate

  Arg[1]     : Integer - replicate 0 = merged or NA, >0 refers to individual replicate
  Example    : if($iset->replicate){ #Do something replicate specific in here }
  Description: Getter for the replicate attribute of this InputSet.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub replicate {  return shift->{replicate}; }



=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of InputSet method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: vendor name feature_class replicate
  Args[4]    : Arrayref - Optional list of InputSet method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: analysis experiment feature_type cell_type
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this InputSet to another based on the defined scalar
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

  $scl_methods ||= [qw(vendor name feature_class replicate)];
  $obj_methods ||= [qw(analysis experiment feature_type cell_type)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}

=head2 reset_relational_attributes

  Arg[1]     : Hashref containing the following parameters:
                [Mandatory]
                -analysis     => Bio::EnsEMBL::Analysis,
                -feature_type => Bio::EnsEMBL::Funcgen::FeatureType,
                -cell_type    => Bio::EnsEMBL::Funcgen::CellType,
                -experiment   => Bio::EnsEMBL::Funcgen::Experiment
                [Optional]
                -support      => Arrayref of valid support objects (eg InputSet)
  Arg[2] ???
  Description: Resets all the relational attributes of a given InputSet.
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;

  my ($analysis, $cell_type, $experiment, $feature_type, $input_subsets) =
      rearrange(['ANALYSIS', 'CELL_TYPE', 'EXPERIMENT', 'FEATURE_TYPE', 'SUBSETS'],
      %$params_hash);

  if(! (defined $analysis &&
        ref($analysis) eq 'Bio::EnsEMBL::Analysis') ){
    throw('You must pass a valid Bio::EnsEMBL::Analysis. ' .
        'Passed: "' . ref($analysis) . '"');
  }

  if(! (defined $cell_type &&
        ref($cell_type) eq 'Bio::EnsEMBL::Funcgen::CellType') ){
    throw('You must pass a valid Bio::EnsEMBL::Funcgen::CellType. ' .
        'Passed: "' . ref($cell_type) . '"');
  }
  if(! (defined $experiment &&
        ref($experiment) eq 'Bio::EnsEMBL::Funcgen::Experiment') ){
    throw('You must pass a valid Bio::EnsEMBL::Funcgen::Experiment. ' .
        'Passed: "' . ref($experiment) . '"');
  }

  if(! (defined $feature_type &&
        ref($feature_type) eq 'Bio::EnsEMBL::Funcgen::FeatureType') ){
    throw('You must pass a valid Bio::EnsEMBL::Funcgen::FeatureType. ' .
        'Passed: "' . ref($feature_type) . '"');
  }

  if(! (defined $input_subsets &&
        ref($input_subsets) eq 'ARRAY') ){
    throw('You must pass an ARRAY_REF containing InputSubsets . ' .
        'Passed: "' . ref($input_subsets) . '"');
  }
  for my $iss(@{$input_subsets}){
    if(ref $iss ne 'Bio::EnsEMBL::Funcgen::InputSubset'){
      throw('You must pass a valid Bio::EnsEMBL::Funcgen::InputSubset. ' .
          'Passed: "' . ref($iss) . '"');
    }
    $self->{subsets}->{$iss->name} = $iss;
  }

  $self->{analysis}     = $analysis;
  $self->{cell_type}    = $cell_type;
  $self->{experiment}   = $experiment;
  $self->{feature_type} = $feature_type;



# Undef dbID and adaptor by default
    if(! $no_db_reset){
      $self->{dbID}    = undef;
      $self->{adaptor} = undef;
    }

    return;
}


### DEPRECATED/REMOVED ###

sub format {  deprecate('The InputSet::format method was removed in v74.');}

sub vendor {  deprecate('The InputSet::vendor method was removed in v74.'); }

sub _add_new_subset {  throw('The InputSet::_add_new_subset method was removed in v74'); }

1;
