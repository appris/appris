#
# Ensembl module for Bio::EnsEMBL::Funcgen::ProbeFeature
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

Bio::EnsEMBL::ProbeFeature - A module to represent an nucleotide probe
genomic mapping.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ProbeFeature;

my $feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new
 (
	-PROBE         => $probe,
	-MISMATCHCOUNT => 0,
	-SLICE         => $chr_1_slice,
	-START         => 1000000,
	-END           => 1000024,
	-STRAND        => -1,
  -ANALYSIS      => $analysis,
  -CIGAR_STRING  => '1U2M426D2M1m21M',
 );

=head1 DESCRIPTION

An ProbeFeature object represents the genomic placement of an Probe
object. The data are stored in the probe_feature table.

=cut

package Bio::EnsEMBL::Funcgen::ProbeFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( median );

use base qw( Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable );


=head2 new

  Arg [-PROBE]        : Bio::EnsEMBL::Funcgen::Probe - probe
        A ProbeFeature must have a probe. This probe must already be stored if
		you plan to store the feature.
  Arg [-MISMATCHCOUNT]: int
        Number of mismatches over the length of the probe. 
  Arg [-SLICE]        : Bio::EnsEMBL::Slice
        The slice on which this feature is.
  Arg [-START]        : int
        The start coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int
        The end coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-STRAND]       : int
        The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int
        Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
        Database adaptor.
  Example    : my $feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
				   -PROBE         => $probe,
				   -MISMATCHCOUNT => 0,
				   -SLICE         => $chr_1_slice,
				   -START         => 1_000_000,
				   -END           => 1_000_024,
				   -STRAND        => -1,
				   -ANALYSIS      => $analysis,
                   -CIGARLINE     => '15M2m3d4M', 
                   #Can represent transcript alignment as gapped genomic alignments
                   #D(eletions) representing introns
                   #Lowercase m's showing sequence mismatches
			   ); 
  Description: Constructor for ProbeFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::ProbeFeature
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  my ($probe, $mismatchcount, $pid, $cig_line)
    = rearrange(['PROBE', 'MISMATCHCOUNT', 'PROBE_ID', 'CIGAR_STRING'], @_);

  #remove mismatch?
  #mandatory args?
  
  #warn "creating probe feature with $pid";
  $self->{'probe_id'} = $pid if $pid;
  $self->probe($probe) if $probe;
  $self->mismatchcount($mismatchcount)  if defined $mismatchcount;#do not remove until probe mapping pipeline fixed
  $self->cigar_string($cig_line)          if defined $cig_line;
   
  #do we need to validate this against the db?  Grab from slice and create new if not present?  Will this be from the dnadb?
  
  #do we need this coordsys id if we're passing a slice?  We should have the method but not in here?

  return $self;
}



=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined. 
  Returntype : Bio::EnsEMBL::Funcgen::ProbeFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  bless ($_[1], $_[0]);
}

=head2 probeset

  Arg [1]    : (optional) string - probeset
  Example    : my $probeset = $feature->probeset();
  Description: Getter and setter for the probeset for this feature. Shortcut
               for $feature->probe->probeset(), which should be used instead.
			   Probeset is not persisted if set with this method.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk
             : Use $feature->probe->probeset() because this may be removed

=cut

sub probeset {
    my $self = shift;
	
    $self->{'probeset'} = shift if @_;
	
    if (! $self->{'probeset'}) {
	  $self->{'probeset'} = $self->probe()->probeset();
    }

	#We could bypass this entirely and call directly using proveset_id?

	
    return $self->{'probeset'};
}


#Only ever needs to be set in _objs_from_sth
#This is to allow linkage of probe_feature glyphs without retrieving the probeset.

sub probeset_id{
  my $self = shift;

  return $self->{'_probeset_id'};
}

=head2 mismatchcount

  Arg [1]    : int - number of mismatches
  Example    : my $mismatches = $feature->mismatchcount();
  Description: Getter and setter for number of mismatches for this feature.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub mismatchcount {
    my $self = shift;
	

    #replace with dynamic check of cigarline?

    $self->{'mismatchcount'} = shift if @_;
	
    return $self->{'mismatchcount'};
}


=head2 cigar_string

  Arg [1]    : str - Cigar line alignment annotation (M = Align & Seq match, m = Align matcht & Seq mismatch, D = Deletion in ProbeFeature wrt genome, U = Unknown at time of alignment)
  Example    : my $cg = $feature->cigar_string();
  Description: Getter and setter for number of the cigar line attribute for this feature.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub cigar_string {
  my $self = shift;
  
  $self->{'cigar_string'} = shift if @_;
	
  return $self->{'cigar_string'};
}


=head2 probe

  Arg [1]    : Bio::EnsEMBL::Funcgen::Probe - probe
  Example    : my $probe = $feature->probe();
  Description: Getter, setter and lazy loader of probe attribute for
               ProbeFeature objects. Features are retrieved from the database
			   without attached probes, so retrieving probe information for a
			   feature will involve another query.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub probe {
  my $self = shift;
  my $probe = shift;
  
  #can we not use _probe_id here?
  #why is probe_id not set sometimes?
  
  
  #warn "in pf and probe is ".$self->{'probe_id'};

  if ($probe) {

	#warn "Probe defined and is ".$probe. "and probe id is".$self->{'probe_id'};
    
    if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::Funcgen::Probe') ) {
      throw('Probe must be a Bio::EnsEMBL::Funcgen::Probe object');
    }
    $self->{'probe'} = $probe;
  }

  if ( ! defined $self->{'probe'}){
	# && $self->dbID() && $self->adaptor() ) {
    #$self->{'probe'} = $self->adaptor()->db()->get_ProbeAdaptor()->fetch_by_ProbeFeature($self);
	#warn "fetching probe with dbID ".$self->probe_id();
	$self->{'probe'} = $self->adaptor()->db()->get_ProbeAdaptor()->fetch_by_dbID($self->probe_id());
  }
  return $self->{'probe'};
}


=head2 probe_id

  Example    : my $probe_id = $pfeature->probe_id();
  Description: Getter for the probe db id of the ProbeFeature
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub probe_id{
  my $self = shift;

  return $self->{'probe_id'} || $self->probe->dbID();
}

=head2 get_results_by_channel_id

  Arg [1]    : int - channel_id (mandatory)
  Arg [2]    : string - Analysis name e.g. RawValue, VSN (optional)
  Example    : my @results = $feature->results();
  Description: Getter, setter and lazy loader of results attribute for
               ProbeFeature objects.
  Returntype : List ref to arrays containing ('score', 'Analysis logic_name');
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_results_by_channel_id {
    my $self = shift;
    my $channel_id = shift;
    my $anal_name = shift;

    warn("This method not fully implemented, remove/deprecate?");

    #$self->{'results'} ||= {};
    $self->{'results_complete'} ||= 0;
	
    if(! $self->{'results'} || ($anal_name && ! exists $self->{'results'}{$anal_name})){
      #fetch all, set complete set flag
      $self->{'results_complete'} ||= 1 	if(! $anal_name);
      
      foreach my $results_ref(@{$self->adaptor->fetch_results_by_channel_analysis($self->probe->dbID(), 
										  $channel_id, $anal_name)}){
	
	$self->{'results'}{$$results_ref[1]} = $$results_ref[0];
      }
    }
    
    return $self->{'results'}
}


#The experiment/al chip specificity has already been done by the ofa->fetch_all_by_Slice_Experiment
#This may be called with no preceding Experiment specificity
#this would return results for all experiments
#do we need to set a default Experiment?


#THis should return both Chip and Channel based results
#just Chip for now
#maybe retrieve and hash all if not Analysis object passed?  Then return what?  


=head2 get_result_by_Analysis_ExperimentalChips

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : listref - Bio::EnsEMBL::Funcgen::ExperimentalChip
  Example    : my $result = $feature->get_result_by_Analysis_ExperimentalChips($anal, \@echips);
  Description: Getter of results attribute for a given Analysis and set of ExperimentalChips
  Returntype : float
  Exceptions : Throws is no Analysis or ExperimentalChips are not passed?
  Caller     : General
  Status     : High Risk

=cut


#make ExperimentalChips optional?

#or have ResultSetAdaptor?  Do we need a ResultSet?
#may not have ExperimentalChip, so would need to return ec dbID aswell


######This will break/return anomalous if
#ECs are passed from different experiments
#ECs are passed from different Arrays


sub get_result_by_Analysis_ExperimentalChips{
    my ($self, $anal, $exp_chips) = @_;

    throw("Need to pass listref of ExperiemntalChips") if(scalar(@$exp_chips) == 0);
    throw("Need to pass a valid Bio::EnsEMBL::Analysis") if ! $anal->isa("Bio::EnsEMBL::Analysis");

    my (%query_ids, %all_ids, %ac_ids);
    my $anal_name = $anal->logic_name();
    
    foreach my $ec(@$exp_chips){
				
      throw("Need to pass a listref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects") 
	if ! $ec->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip");

		#my $tmp_id = $self->adaptor->db->get_ArrayAdaptor->fetch_by_array_chip_dbID($ec->array_chip_id())->dbID();
      
		#$array_id ||= $tmp_id;
      
      #throw("You have passed ExperimentalChips from different if($array_id != $tmp_id)
      
      #if(exists  $ac_ids{$ec->array_chip_id()}){
#	throw("Multiple chip query only works with contiguous chips within an array, rather than duplicates");
 #     }
      
      $ac_ids{$ec->array_chip_id()} = 1;
      $all_ids{$ec->dbID()} = 1;
      $query_ids{$ec->dbID()} = 1 if(! exists $self->{'results'}{$anal_name}{$ec->dbID()});
      
    }
    
    
    my @ec_ids = keys %query_ids;
    my @all_ids = keys %all_ids;
    
    
    #warn "ec ids @ec_ids\n";
    #warn "all ids @all_ids\n";
    
    #$self->{'results'} ||= {};
    #$self->{'results_complete'} ||= 0;#do we need this now?
    
    if((scalar(@all_ids) - scalar(@ec_ids))> 1){
      throw("DATA ERROR - There is more than one result stored for the following ExperimentalChip ids: @all_ids");
    }		
    elsif(! $self->{'results'} || (($anal_name && scalar(@ec_ids) > 0) && scalar(@all_ids) == scalar(@ec_ids))){
      #fetch all, set complete set flag
      #$self->{'results_complete'} ||= 1 	if(! $anal_name);
      #would need to look up chip and channel analyses here and call relevant fetch
      #or pass the chip and then build the query as = or IN dependent on context of logic name
      #if there are multiple results, last one will overwrite others
      #could do foreach here to deal with retrieving all i.e. no logic name
      #Can supply mutliple chips, but probe ids "should" be unique(in the DB at least) amongst contiguous array_chips
      #build the cache based on logic name and table_id
      #cahce key??  should we cat the ec_ids together?

      my @result_refs = @{$self->adaptor->fetch_results_by_Probe_Analysis_experimental_chip_ids($self->probe(), 
												$anal,
												\@ec_ids)};

      #Remove lines with no result
      while(@result_refs && (! $result_refs[0]->[0])){
	shift @result_refs;
      }

      my $num_results = scalar(@result_refs);
      my ($result, $mpos);
      #throw("Fetched more than one result for this ProbeFeature, Analysis and ExperimentalChips") if (scalar(@result_refs) >1);

      #No sort needed as we sort in the query

      if($num_results == 1){
	$result = $result_refs[0]->[0];
      }
      elsif($num_results == 2){#mean
	$result = ($result_refs[0]->[0] + $result_refs[1]->[0])/2;
    
      }
      elsif($num_results > 2){#median or mean of median flanks
	$mpos = $num_results/2;
    
	if($mpos =~ /\./){#true median
	  $mpos =~ s/\..*//;
	  $mpos ++;
	  $result =  $result_refs[$mpos]->[0];
	}else{
	  $result = ($result_refs[$mpos]->[0] + $result_refs[($mpos+1)]->[0])/2 ;
	}
      }
      
      $self->{'results'}{$anal_name}{":".join(":", @ec_ids).":"} = $result;
    }

	#do we return the ec ids here to, or do we trust that the user will know to only pass contiguous rather than duplicate chips

	#how are we going to retrieve the result for one of many possible ec id keys?
	#options, cat ec dbids as key, and grep them to find full key, then return result
	#this may hide the duplicate chip problem
	#If a query has already been made and cached,another query with one differing ID(duplicate result) may never be queried as we already have a cahced result
	#We shoulld pick up duplicates before this happens
	#If we try and mix ExperimentalChips from different experiments, then this would also cause multiple results, and hence hide some data
	
	my @keys;
	foreach my $id(@all_ids){
	  my @tmp = grep(/:${id}:/, keys %{$self->{'results'}{$anal_name}});
	  #Hacky needs sorting, quick fix for release!!

	  if(@tmp){
	    push @keys, grep(/:${id}:/, keys %{$self->{'results'}{$anal_name}});

	    last;
	  }

	}

    throw("Got more than one key for the results cache") if scalar(@keys) > 1;

    return $self->{'results'}{$anal_name}{$keys[0]};
}


#Will this be too slow, can we not do one query across all tables

sub get_result_by_ResultSet{
    my ($self, $rset) = @_;

    my $results = $rset->adaptor->fetch_results_by_probe_id_ResultSet($self->probe_id(), $rset);
   
    return median($results);
}

=head2 summary_as_hash

  Example       : $array_feature_summary = $array_feature->summary_as_hash;
  Description   : Retrieves a textual summary of this ArrayFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self = shift;
  my $probe = $self->probe;
 
  my %arrays = ();
  my $names_Arrays = $probe->get_names_Arrays; 
  foreach my $name (keys %$names_Arrays) {
    if (defined $names_Arrays->{$name}) {
      $arrays{$names_Arrays->{$name}->name} = $probe->get_probename($names_Arrays->{$name}->name);
    }
  }

  return
   {
    start                   => $self->seq_region_start,
    end                     => $self->seq_region_end,
    strand                  => $self->strand,
    seq_region_name         => $self->seq_region_name,
    feature_type            => 'array_probe',

    probe_length            => $probe->length,
    array_probe_names       => \%arrays
   };
}

1;
