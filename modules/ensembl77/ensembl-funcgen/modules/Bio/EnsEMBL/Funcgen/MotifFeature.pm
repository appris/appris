#
# Ensembl module for Bio::EnsEMBL::Funcgen::MotifFeature
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

=head1 NAME

Bio::EnsEMBL::MotifFeature - A module to represent a feature mapping as based
on a binding matrix e.g position weight matrix

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::MotifFeature;

my $feature = Bio::EnsEMBL::Funcgen::MotifFeature->new
 (
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
  -DISPLAY_LABEL => $text,
  -SCORE         => $score,
  -FEATURE_TYPE  => $ftype,
  -INTERDB_STABLE_ID    => 1,
 );

=head1 DESCRIPTION

A MotifFeature object represents the genomic placement of a sequence motif.
For example a transcription factor binding site motif associated with a
position weight matrix. These are generally associated with AnnotatedFeatures
of the corresponding FeatureType.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Funcgen::MotifFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable);


=head2 new

 
  Arg [-SCORE]          : (optional) int - Score assigned by analysis pipeline
  Arg [-SLICE]          : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-BINDING_MATRIX] : Bio::EnsEMBL::Funcgen::BindingMatrix - Binding Matrix associated to this feature.
  Arg [-START]          : int - The start coordinate of this feature relative to the start of the slice
		                it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]            : int -The end coordinate of this feature relative to the start of the slice
	                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-DISPLAY_LABEL]  : string - Display label for this feature
  Arg [-STRAND]         : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]           : (optional) int - Internal database ID.
  Arg [-ADAPTOR]        : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Example    : my $feature = Bio::EnsEMBL::Funcgen::MotifFeature->new(
                                									  -SLICE          => $chr_1_slice,
								                                	  -START          => 1_000_000,
                                									  -END            => 1_000_024,
								                                	  -STRAND         => -1,
                                									  -BINDING_MATRIX => $bm,
								                                	  -DISPLAY_LABEL  => $text,
                                  									  -SCORE          => $score,
                                                                      -INTERDB_STABLE_ID     => 1,
                                                                     );


  Description: Constructor for MotifFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Exceptions : Throws if BindingMatrix not valid
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  my $bmatrix;
  ($self->{score}, $bmatrix,   $self->{display_label}, $self->{interdb_stable_id}) 
   = rearrange(['SCORE', 'BINDING_MATRIX', 'DISPLAY_LABEL', 'INTERDB_STABLE_ID'], @_);
    

  if(! (ref($bmatrix) && $bmatrix->isa('Bio::EnsEMBL::Funcgen::BindingMatrix'))){
	throw('You must pass be a valid Bio::EnsEMBL::Funcgen::BindingMatrix');
  }

  $self->{binding_matrix} = $bmatrix;

  return $self;
}

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the calling code 
               is very disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  return bless ($_[1], $_[0]);
}


=head2 binding_matrix

  Example    : my $bmatrix_name = $mfeat->binding_matrix->name;
  Description: Getter for the BindingMatrix attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub binding_matrix{
  return shift->{binding_matrix};
}

=head2 feature_type

  Example    : my $TF_name = $motif_feature->feature_type->name;
  Description: Convenience method for accessing feature type of binding matrix
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut


sub feature_type{
  return shift->{binding_matrix}->feature_type;
}


=head2 score

  Example    : my $score = $feature->score();
  Description: Getter for the score attribute for this feature. 
  Returntype : double
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub score {
  return shift->{score};
}


=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_label {
  #If not set in new before store, a default is stored as:
  #$mf->binding_matrix->feature_type->name.':'.$mf->binding_matrix->name();
  return shift->{display_label};
}



=head2 associated_annotated_features

  Example    : my @associated_afs = @{$feature->associated_annotated_features()};
  Description: Getter/Setter for associated AnntoatedFeatures.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen:AnnotatedFeature objects
  Exceptions : None
  Caller     : General
  Status     : At risk - may change to associated_transcript_factor_features

=cut

sub associated_annotated_features{
  my $self = shift;
  my $afs  = shift;
  #Lazy load as we don't want to have to do a join on all features when most will not have any
 
  if (defined $afs) {

    if (ref($afs) eq 'ARRAY') {

      foreach my $af (@$afs) {
	
        if ( ! $af->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature') ) {
          throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::AnnotatedFeature objects');
        }
        #test is stored in adaptor
      }

      if (defined $self->{associated_annotated_features}) {
        warn('You are overwriting associated_annotated_features for the MotifFeature');
        #we could simply add the new ones and make them NR.
      }

      $self->{associated_annotated_features} = $afs;
    } 
    else {
      throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::AnnotatedFeature objects');
    }
  }


  if (! defined $self->{associated_annotated_features}) {

    if (defined $self->adaptor) {
      $self->{associated_annotated_features} = 
        $self->adaptor->db->get_AnnotatedFeatureAdaptor->fetch_all_by_associated_MotifFeature($self);
    }
  }
  
  #This has the potential to return undef, or an arrayref which may be empty.
  return $self->{associated_annotated_features};
}


=head2 is_position_informative

  Arg [1]    : int - 1-based position within the Motif
  Example    : $mf->is_position_informative($pos);
  Description: Indicates if a given position within the motif is highly informative
  Returntype : boolean
  Exceptions : throws if position out of bounds ( < 1 or > length of motif)
  Caller     : General
  Status     : At High risk

=cut

sub is_position_informative {
  my $self     = shift;
  my $position = shift || throw 'Need a position argument';

  if( ($position < 1) || 
      ($position > $self->binding_matrix->length) ){
    throw "Position outside of MotifFeature loci";
  }

  #if on the opposite strand, then need to reverse complement the position
  if($self->strand < 0){ 
    $position = $self->binding_matrix->length - $position + 1; 
  }

  return $self->binding_matrix->is_position_informative($position);
}


=head2 infer_variation_consequence

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature
  Arg [2]    : boolean - 1 if results in linear scale (default is log scale)
  Example    : my $vfs = $vf_adaptor->fetch_all_by_Slice($slice_adaptor->fetch_by_region('toplevel',$mf->seq_region_name,$mf->start,$mf->end,$mf->strand));
               foreach my $vf (@{$vfs}){
                   print $mf->infer_variation_consequence($vf)."\n";
               }

  Description: Calculates the potential influence of a given variation in a motif feature.
               Returns a value between -100% (lost) and +100% (gain) indicating the difference 
               in strength between the motif in the reference and after the variation.

               The variation feature slice needs to be the motif feature, including the strand
  Returntype : float
  Exceptions : throws if argument is not a  Bio::EnsEMBL::Variation::VariationFeature
               throws if the variation feature is not contained in the motif feature
  Caller     : General
  Status     : At High risk

=cut

sub infer_variation_consequence{
  my ($self, $vf, $linear) = @_;

  if(! $vf->isa('Bio::EnsEMBL::Variation::VariationFeature')){
    throw "We expect a Bio::EnsEMBL::Variation::VariationFeature object, not a ".$vf->class;
  }

  #See if these checks are required or if there are more efficient ways to do the checks...
  #if(($self->slice->seq_region_name ne $vf->slice->seq_region_name) ||
  #   ($self->slice->start != $vf->slice->start) || 
  #   ($self->slice->end != $vf->slice->end) ){
  #  throw "Variation and Motif are on distinct slices";
  #}
  #if(!(($vf->start >= $self->start) && ($vf->end <= $self->end ))){
  #  throw "Variation should be entirely contained in the Motif";
  #}

  if( ($vf->start < 1) || ($vf->end > $self->binding_matrix->length)){ throw "Variation not entirely contained in the motif feature"; }

  if(!($vf->allele_string =~ /^[ACTG]\/[ACTG]$/)){ throw "Currently only SNPs are supported"; }

  my $ref_seq = $self->seq;

  my $variant = $vf->allele_string;
  $variant =~ s/^.*\///;
  $variant =~ s/\s*$//;

  my ($vf_start,$vf_end) = ($vf->start, $vf->end);
  if($vf->strand == -1){
    #Needed for insertions
    $variant = reverse($variant);
    $variant =~ tr/ACGT/TGCA/;
  }
  my $var_seq = substr($ref_seq,0, $vf_start - 1).$variant.substr($ref_seq, $vf_start+length($variant)-1);

  my $bm = $self->{'binding_matrix'};
  return 100 * ($bm->relative_affinity($var_seq,$linear) - $bm->relative_affinity($ref_seq,$linear));
  
}

=head2 interdb_stable_id

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $idb_sid = $feature->interdb_stable_id();
  Description: Getter for the interdb_stable_id attribute for this feature.
               This is simply to avoid using internal db IDs for inter DB linking
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub interdb_stable_id {
  return shift->{interdb_stable_id};
}


=head2 summary_as_hash

  Example       : $motif_feature_summary = $motif_feature->summary_as_hash;
  Description   : Retrieves a textual summary of this MotifFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self = shift;
  my ($acc, $ftype);
  #split display_label as binding matrix may be lazy loaded and slow things down
  
  if ($self->display_label =~ /(.*[^:])(:)(.*)/o){ 
    $ftype = $1;
    $acc = $3;
  }
  else{
    warn "Failed to parse feature type and binding matric from display_id:\t".$self->display_id;
  }

  #Add bm.threshold in here?
  return
   {binding_matrix          => $acc,
    motif_feature_type      => $ftype,
    start                   => $self->seq_region_start,
    end                     => $self->seq_region_end,
    strand                  => $self->strand,
    seq_region_name         => $self->seq_region_name,
    score                   => $self->score            };
}


1;

