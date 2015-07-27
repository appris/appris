#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor - A database aware adaptor responsible for
the creation of Slices in the context of eFG objects.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
 

  $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);

  $slice_adaptor = $db->get_SliceAdaptor();

  # get a slice on the entire chromosome X
  $gene_regulation_slice = $slice_adaptor->fetch_by_Gene_FeatureSets($gene, \@fsets);


=head1 DESCRIPTION

This module is simple wrapper class for the core SliceAdaptor, extending new
methods to generate Slices for eFG features associated with a given gene or transcript.

=cut


package Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::DBSQL::SliceAdaptor);

sub new {
  my ($caller, $efgdb) = @_;

  my $class = ref($caller) || $caller;

  #Now have to reset the efg db to the dnadb
  #Does work using just efg.seq_region_id
  #But this is undesirable as we want to perform 
  #other non seq_region table style queries on core DB
  my $self = $class->SUPER::new($efgdb->dnadb);
  $self->efgdb($efgdb);
  return $self;
}


sub efgdb{
  my ($self, $efgdb) = @_;
  
  if($efgdb){

	if(! (ref($efgdb) && $efgdb->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	  throw('Must provide a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
	}
	
	$self->{'efgdb'} = $efgdb;	
  }

  return $self->{'efgdb'};
}



=head2 fetch_by_Gene_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : my $gene_reg_slice = $efg_slice_adaptor->fetch_by_gene($gene);
  Description: Fetches a slice by finding the associated regulatory
               elements associated with a given gene(and it's transcripts and 
               translations), and extending the gene feature slice to 
               encapsulate the associated regulatory elements.
  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_Gene_FeatureSets {
  my ($self, $gene, $fsets) = @_;

  if(! ( ref($gene) && $gene->isa('Bio::EnsEMBL::Gene'))){
    throw("You must pass a valid Bio::EnsEMBL::Gene object");
  }
  
  #Now we want to get all the eFG feature associated with this gene and build the slice
  #resitrict to RegFeats for now.
  my $cs_name      = $gene->slice->coord_system_name;
  my $gene_chr     = $gene->seq_region_name;
  my $start        = $gene->seq_region_start;
  my $end          = $gene->seq_region_end;
  

  ($start, $end) = $self->_set_bounds_by_xref_FeatureSets
	($gene_chr, $start, $end, $gene, $fsets);

  #Now need to do this for all transcripts and proteins too?
  
  foreach my $trans(@{$gene->get_all_Transcripts}){
	($start, $end) = $self->_set_bounds_by_xref_FeatureSets
	  ($gene_chr, $start, $end, $trans, $fsets);

	my $translation = $trans->translation;
	($start, $end) = $self->_set_bounds_by_xref_FeatureSets
	  ($gene_chr, $start, $end, $translation, $fsets) if $translation;
  }

  return $self->fetch_by_region($cs_name, $gene_chr, $start, $end, $gene->strand);
}

=head2 fetch_by_Transcript_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : my $trans_reg_slice = $efg_slice_adaptor->fetch_by_transcript($transcript);
  Description: Fetches a slice by finding the associated regulatory
               elements associated with a given transcript(and it's translation), and 
               extending the gene feature slice to encapsulate the associated regulatory 
               elements.
  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : self
  Status     : at risk

=cut

sub fetch_by_Transcript_FeatureSets{
  my ($self, $transcript, $fsets) = @_;


  if(! ( ref($transcript) && $transcript->isa('Bio::EnsEMBL::Transcript'))){
    throw("You must pass a valid Bio::EnsEMBL::Transcript object");
  }
  
  #Now we want to get all the eFG feature associated with this gene and build the slice
  #resitrict to RegFeats for now.
  my $cs_name   = $transcript->slice->coord_system_name;
  my $trans_chr = $transcript->seq_region_name;
  my $start     = $transcript->seq_region_start;
  my $end       = $transcript->seq_region_end;
  

  ($start, $end) = $self->_set_bounds_by_xref_FeatureSets
	($trans_chr, $start, $end, $transcript, $fsets);

  #Now need to do this for the protein too
  my $translation = $transcript->translation;
  ($start, $end) = $self->_set_bounds_by_xref_FeatureSets
	($trans_chr, $start, $end, $translation, $fsets) if $translation;


  return $self->fetch_by_region($cs_name, $trans_chr, $start, $end, $transcript->strand);
}

#Do we need fetch_by_translation?

=head2 _set_bounds_by_xref_FeatureSets

  Arg [1]    : string - seq_region_name i.e. chromosome name.
  Arg [2]    : string - seq_region_start of current slice bound
  Arg [3]    : string - seq_region_end of current slice bound.
  Arg [4]    : Bio::EnsEMBL::Gene|Transcript|Translation
  Arg [5]    : arrayref - Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : ($start, $end) = $self->_set_bounds_by_regulatory_feature_xref
                              ($trans_chr, $start, $end, $transcript, $fsets);
  Description: Internal method to set an xref Slice bounds given a list of
               FeatureSets.
  Returntype : List - ($start, $end);
  Exceptions : throw if incorrent args provided
  Caller     : self
  Status     : at risk

=cut

sub _set_bounds_by_xref_FeatureSets{
  my ($self, $chr, $start, $end, $obj, $fsets) = @_;

  my ($extdb_name, $efg_feature);
  my $dbe_adaptor = $self->efgdb->get_DBEntryAdaptor;


  #do we need to test for start/end/chr here?
  #Is implicit if we test for obj and fsets, but 
  #does not check if defined.

  #Set ext_dbname and validate obj
  #Do we need a central store for ensembl db names?

  


  if($obj->isa('Bio::EnsEMBL::Gene')){
	$extdb_name = 'ensembl_core_Gene';
  }
  elsif($obj->isa('Bio::EnsEMBL::Transcript')){
	$extdb_name = 'ensembl_core_Transcript';
  }
  elsif($obj->isa('Bio::EnsEMBL::Translation')){
	$extdb_name = 'ensembl_core_Translation';
  }
  else{
	throw('Currently only handles Ensembl Gene, Transcript and Translation xrefs');
  }


  #warn "Setting bounds by $obj ".$obj->stable_id;

  #Set which eFG features we want to look at.

  if(ref($fsets) ne 'ARRAY' || scalar(@$fsets) == 0){
	throw('Must define an array of Bio::EnsEMBL::FeatureSets to extend xref Slice bound. You passed: '.$fsets);
  }

  my %feature_set_types;

  foreach my $fset(@$fsets){
	$self->efgdb->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	
	$feature_set_types{$fset->feature_class} ||= [];
	push @{$feature_set_types{$fset->feature_class}}, $fset;
  }


  #We can list the outer loop here and put in the BaseFeatureAdaptor, or possible storable as we do have FeatureType xrefs.
  #This would be useful for fetching all the efg features for a given xref and FeatureSets
  #Don't implement as a parent sub and call from here as this would mean looping through array twice.
  #Altho we could pass a code ref to do the filtering?
  #Just copy and paste for now to avoid obfuscation


  #Get xrefs for each eFG feature set type
  foreach my $fset_type(keys %feature_set_types){

	my $xref_method    = 'list_'.$fset_type.'_feature_ids_by_extid';
	#e.g. list_regulatory_feature_ids_by_extid

	my $adaptor_method = 'get_'.ucfirst($fset_type).'FeatureAdaptor';

	my $adaptor = $self->efgdb->$adaptor_method;

	my %feature_set_ids;
	map $feature_set_ids{$_->dbID} = 1, @{$feature_set_types{$fset_type}};

	#my $cnt = 0;


	#This should use fetch_all_by_external_name method

	foreach my $efg_feature(@{$adaptor->fetch_all_by_external_name($obj->stable_id, $extdb_name)}){
	  next if ! exists $feature_set_ids{$efg_feature->feature_set->dbID};
	  next if $efg_feature->seq_region_name ne $chr;
	  

	  #warn "found xref ".$efg_feature->display_label.' with start '.$efg_feature->seq_region_start;
	  #$cnt ++;
  	  
	  $start = $efg_feature->seq_region_start if $efg_feature->seq_region_start < $start;
	  $end   = $efg_feature->seq_region_end   if $efg_feature->seq_region_end   > $end;
	}
	
	#warn "Found $cnt $fset_type xrefs";

  }



  return ($start, $end);
}


1;
