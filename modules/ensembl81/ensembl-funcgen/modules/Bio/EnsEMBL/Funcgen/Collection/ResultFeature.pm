#
# Ensembl module for Bio::EnsEMBL::Funcgen::Collection::ResultFeature
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

Bio::EnsEMBL::Funcgen::Collection::ResultFeature - A module to represent a lightweight ResultFeature collection

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Collection::ResultFeature;

my $rfeature = Bio::EnsEMBL::Funcgen::Collection::ResultFeature->new_fast({
   start  => $start, 
   end    => $end, 
   scores => [ $score ]
});

my @rfeatures = @{$rset->get_displayable_ResultFeatures_by_Slice($slice)};

foreach my $rfeature (@rfeatures){
    my $score = $rfeature->score();
    my $rf_start = $rfeature->start();
    my $rf_end = $rfeature->end();
}

=head1 DESCRIPTION

This is a Collection feature which is designed to store compressed/collected 
feature information for a defined window/bin size over a complete seq_region.
Or alternatively a single feature at the natural resolution i.e. window_size == 0.
The complete seq_region collections are cropped to provide a ResultFeature on any
given Slice. ResultFeatures are primarily stored in the result_feature table, 
but can also be generated on the fly from unprocessed data in the array result
tables.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::ResultFeatureAdaptor
Bio::EnsEMBL::Funcgen::Collector::ResultFeature

=cut


#This is distinct from a normal feature as the collection may have differen attributes and methods from the normal feature
#implementation.  For example a Bio::EnsEMBL::Collection::Gene would only have summary information over several genes.
#Altho', unlikely that we'll ever collect genes.

package Bio::EnsEMBL::Funcgen::Collection::ResultFeature;

use strict;
use warnings;
use base qw( Bio::EnsEMBL::Feature );

#This needs to inherit from Bio::EnsEMBL::Collection
#Which can host some of the below methods

#Reverted to hash implementation as we no longer deal with 
#huge amounts of features due to collections.
#This enables use of transform/seq_region_start/end methods
#and enable us to store on slices that do not begin at 1
#Altho need to remove code stipulating this



#To do
#Can we move any of these methods to a base Collection class?
#Should probably now use normal new method with validation?

=head2 new_fast

  Args       : Array with attributes start, end, strand, scores, probe, result_set_id, window_size, slice  IN THAT ORDER.
               WARNING: None of these are validated, hence can omit some where not needed
  Example    : none
  Description: Fast and list version of new. Only works if the code is very disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::ResultFeature
  Exceptions : None
  Caller     : ResultSetAdaptor
  Status     : At Risk

=cut

sub new_fast {
  #This is agnostic towards to type of reference
  return bless ($_[1], $_[0]);

}



=head2 scores

  Example    : my $score = $rf->score();
  Description: Getter of the scores attribute for ResultFeature
               objects
  Returntype : Arrayref.
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub scores {    return shift->{scores}; }



=head2 probe

  Example    : my $probe = $rf->probe();
  Description: Getter of the probe attribute for ResultFeature
               objects
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : None
  Caller     : General
  Status     : At Risk - This can only be used for Features with window 0.

=cut

#probe_id is currently not available in the result_feature table, so this would be a result/probe_feature query.

sub probe {    return shift->{probe}; }


sub result_set_id {    return shift->{result_set_id}; }

sub window_size {    return shift->{window_size}; }

sub get_min_max_scores{
  my $self = shift;

  if(! defined $self->{min_max_scores}){
    my @sorted_scores = sort { $a <=> $b } @{$self->{scores}};
    $self->{min_max_scores} = [$sorted_scores[0], $sorted_scores[$#sorted_scores]];
  }

  return $self->{min_max_scores};
}


1;

