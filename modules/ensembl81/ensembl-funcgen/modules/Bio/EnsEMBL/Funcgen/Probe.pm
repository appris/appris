#
# Ensembl module for Bio::EnsEMBL::Funcgen::Probe
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

Bio::EnsEMBL::Funcgen::Probe - A module to represent a nucleotide probe.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Probe;

my $probe = Bio::EnsEMBL::Funcgen::Probe->new
 (
  -PROBE_SET     => $probe_set,
	-NAME          => 'Probe-1',
  -ARRAY         => $array,
  -ARRAY_CHIP_ID => $ac_dbid,
	-CLASS         => "EXPERIMENTAL",
 );

=head1 DESCRIPTION

An Probe object represents an probe on a microarray. The data (currently the
name, probe_set_id, length, pair_index and class) are stored
in the oligo_probe table.

For Affy arrays, a probe can be part of more than one array, but only part of
one probeset. On each Affy array the probe has a slightly different name. For
example, two different complete names for the same probe might be
DrosGenome1:AFFX-LysX-5_at:535:35; and Drosophila_2:AFFX-LysX-5_at:460:51;. In
the database, these two probes will have the same oligo_probe_id. Thus the same
Affy probe can have a number of different names and complete names depending on
which array it is on.

=cut

package Bio::EnsEMBL::Funcgen::Probe;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use base qw( Bio::EnsEMBL::Funcgen::Storable );


=head2 new

  Arg [-NAME]          : string - probe name
        Used when the probe is on one array.
  Arg [-NAMES]         : Listref of strings - probe names
        Used when the probe is on multiple arrays.
  Arg [-ARRAY]          : Bio::EnsEMBL::Funcgen::Array
        Used when the probe is on one array.
  Arg [-ARRAYS]         : Listref of Bio::EnsEMBL::Funcgen::Array
        Used when the probe is on multiple arrays.
  Arg [-ARRAY_CHIP_ID]  : int - array_chip db ID
        Used when the probe is on one array.
  Arg [-ARRAY_CHIP_IDS]  : Listref of ints - array_chip dbIDs
        Used when the probe is on multiple array chips
  Arg [-NAMES]          : Listref of ints - arary_chip db IDs
        Used when the probe is on multiple arrays.
  Arg [-PROBE_SET]      : Bio::EnsEMBL::ProbeSet
        Each probe is part of one(and only one) probeset, if not probe set
        then probeset = probe i.e. probe_set size = 1
  Arg [-LENGTH]         : int - probe length
        Will obviously be the same for all probes if same probe
		is on multiple arrays.
  Arg [-CLASS]          : string - probe class e.g. CONTROL, EXPERIMENTAL
        Will be the same for all probes if same probe is on
		multiple arrays.
  Arg [-DESCRIPTION]    : (optional) string - description


  Example    : my $probe = Bio::EnsEMBL::Probe->new(
                   -NAME          => 'Probe-1',
				   -PROBE_SET     => $probe_set,
                   -ARRAY         => $array,
                   -ARRAY_CHIP_ID => $array_chip_id,
				   -LENGTH        => 25,
                   -CLASS         => 'EXPERIMENTAL',
                   -DESCRIPTION   => 'Some useful description',

               );
  Description: Creates a new Bio::EnsEMBL::Probe object.
  Returntype : Bio::EnsEMBL::Probe
  Exceptions : Throws if not supplied with probe name(s) and array(s)
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my (
      $names,          $name,
      $array_chip_ids, $array_chip_id,
      $arrays,         $array,
      $probeset,       $aclass,
      $length,         $desc
     ) = rearrange([
		    'NAMES',          'NAME',
		    'ARRAY_CHIP_IDS', 'ARRAY_CHIP_ID',
		    'ARRAYS',         'ARRAY',
		    'PROBE_SET',      'CLASS',
		    'LENGTH',         'DESCRIPTION'
		   ], @_);


  @$names = ($name) if(ref($names) ne "ARRAY");
  @$array_chip_ids = ($array_chip_id) if (ref($array_chip_ids) ne "ARRAY");
  @$arrays = ($array) if (ref($arrays) ne "ARRAY");

  #We need to record duplicates for each probe_set i.e. each array.
  #the relationship is really array_chip to name, as everything else stays the same
  #can't have same probe_set_id as this wouldn't maintain relationship
  #need unique ps id's or array_chip_id in probe table?
  #Then we can miss probeset id's out totally if required
  #or should we just duplicate everything with unique db IDs


  if (defined $$names[0]) {

    if(scalar(@$names) != scalar(@$array_chip_ids)){
      throw("You have not specified valid name:array_chip_id pairs\nYou need a probe name for each Array");
    }

    if(defined $$arrays[0]){
      if(scalar(@$names) != scalar(@$arrays)){
	throw("You have not specified valid name:Array pairs\nYou need a probe name for each Array\n");
      }
    }
    else{
      warn("You have not specified and Array objects, this will result in multiple/redundant queries based on the array_chip_id\nYou should pass Array objects to speed up this process");
	  #Is this true? We should cache this in the ArrayChip and make sure we're caching it in the caller.
    }

    # Probe(s) have been specified
    # Different names reflect different array

    for my $i(0..$#{$names}){
      $self->add_array_chip_probename($$array_chip_ids[$i], $$names[$i], $$arrays[$i]);
    }
  } else {
    throw('You need to provide a probe name (or names) to create an Probe');
  }

  $self->probeset($probeset) if defined $probeset;
  $self->class($aclass)      if defined $aclass;
  $self->length($length)     if defined $length;
  $self->description($desc)  if defined $desc;

  return $self;
}

#only takes single values for array and array_chip
#as we're shortcuting the constructor and simply blessing the hash
#therefore attr keys should not be lc and not prefix with '-'

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined. Cannot add array chip probe names unless we recreate
               the data structure in the caller.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  bless ($_[1], $_[0]);
}



=head2 add_array_chip_probename

  Arg [1]    : int - db ID of array_chip
  Arg [2]    : string - probe name
  Arg [3]    : Bio::EnsEMBL::Funcgen::Array
  Example    : $probe->add_array_chip_probename($ac_dbid, $probename, $array);
  Description: Adds a probe name / array pair to a probe, allowing incremental
               generation of a probe.
  Returntype : None
  Exceptions : None
  Caller     : General,
               Probe->new(),
               ProbeAdaptor->_obj_from_sth(),
			   AffyProbeAdaptor->_obj_from_sth()
  Status     : Medium Risk - Change to take ArrayChip object.

=cut

#mapping between probename and ac_dbid is conserved through array name between hashes
#only easily linked from arrays to probenames,as would have to do foreach on array name
  
#Can we change the implementation of this so we're only storing the array once, reverse
#the cache? But we want access to the array and using an object reference as a key is ????
#How would this impact on method functionality?
  
#We now handle multiple names per probe/array
#This will not capture the relationship between
#probe name and position on array!
#Not a problem for affy as name is position
#Currently not a problem for nimblegen as probes never have more than 1 name???
  
#This does not however accomodate multiple names on the same ArrayChip
#and there is currently no way of validating this with the current data model
#This is however, caught by ImportArrays


sub add_array_chip_probename {
  my $self = shift;
  my ($ac_dbid, $probename, $array) = @_;
  $self->{arrays}     ||= {};
  $self->{probenames} ||= {};
    
  if(! (ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array') && $array->dbID)){
	  throw('You must pass a valid Bio::EnsEMBL::Funcgen::Array. '.
          'Maybe you want to generate a cache in the caller?');
	}

  $self->{arrays}->{$ac_dbid}           = $array;
  $self->{probenames}->{$array->name} ||= [];
  push @{$self->{probenames}->{$array->name}}, $probename;
  
  return;
}


=head2 get_all_ProbeFeatures

  Args       : None
  Example    : my $features = $probe->get_all_ProbeFeatures();
  Description: Get all features produced by this probe. The probe needs to be
               database persistent.
  Returntype : Listref of Bio::EnsEMBL:Funcgen::ProbeFeature objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_ProbeFeatures {
	my $self = shift;
	if ( $self->adaptor() && $self->dbID() ) {
		return $self->adaptor()->db()->get_ProbeFeatureAdaptor()->fetch_all_by_Probe($self);
	} else {
		warning('Need database connection to retrieve Features');
		return [];
	}
}

=head2 get_all_Arrays

  Args       : None
  Example    : my $arrays = $probe->get_all_Arrays();
  Description: Returns all arrays that this probe is part of. Only works if the
               probe was retrieved from the database or created using
			   add_Array_probename (rather than add_arrayname_probename).
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_Arrays {
    my $self = shift;

	#Arrays are currently preloaded using a cache in _objs_from_sth
	return [ values %{$self->{'arrays'}} ];
}

=head2 get_names_Arrays

  Args       : None
  Example    : my %name_array_pairs = %{$probe->get_names_Arrays};
  Description: Returns Array name hash
  Returntype : hashref of probe name Bio::EnsEMBL::Funcgen::Array pairs
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_names_Arrays {
    my $self = shift;

	#Arrays are currently preloaded using a cache in _objs_from_sth
	return $self->{'arrays'};
}





=head2 get_all_probenames

  Arg [1]    : Optional - list of array names, defaults to all available
  Example    : my @probenames = @{$probe->get_all_probenames()};
  Description: Retrieves all names for this probe. Only makes sense for probes
               which share identical sequence for a given probeset and array.
               This can either be:
               1 A non-probeset array where probes with different names but identical 
                 sequence have beem merged, this is only possible if the probes in question 
                 share the same Array but are on seaprate ArrayChips.
               2 If they are part of a probeset (i.e. Affy probes), in which case 
                 get_all_complete_names() would be more appropriate.
  Returntype : Arrayref of strings
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_all_probenames {
  my ($self, @array_names) = @_;

	my @names;
	@array_names = keys %{$self->{'probenames'}} if ! @array_names;

	foreach my $array(@array_names){
	  push @names, @{$self->{'probenames'}->{$array}};
	}

  return \@names;
}



=head2 get_probename

  Arg [1]    : string - array name
  Example    : my $probename = $probe->get_probename('Array-1');
  Description: For a given array, retrieve the name for this probe.
  Returntype : string
  Exceptions : Throws if the array name is required but not specified
               Warns if probe has more than one name for the given array.
  Caller     : General
  Status     : Medium Risk

=cut


#we can have dulplicate probes on same array for Nimblegen
#what defines and unique probe?
#If we have a duplicate on the same array or even on the same array_chip, then we can still return the same name
#Needs more work

sub get_probename {
  my ($self, $arrayname) = @_;
	my $probename;

  if (! $arrayname){
    #Sanity check that there is only one non-AFFY array
    my @ac_ids = keys %{$self->{'arrays'}};

    if((scalar @ac_ids == 1) && ($self->get_all_Arrays()->[0]->vendor() ne "AFFY")){
      $arrayname = $self->get_all_Arrays()->[0]->name();
    }
    else{
      throw('Cannot retrieve name for Probe('.$self->dbID.") without arrayname if more than 1 array chip(@ac_ids) and not NIMBELGEN(".$self->get_all_Arrays()->[0]->vendor().")\n");
    }
  }

	return if ! exists ${$self->{'probenames'}}{$arrayname};

	my @names = @{$self->{'probenames'}->{$arrayname}};

	if(scalar(@names) > 1){
    my $p_info = '';

    if($self->probeset){
      $p_info = " probeset ".$self->probeset->name;
    }


	  #warn("Found replicate probes with different names for array ${arrayname}${p_info}.Returning comma separated string list:\t".join(',', @names)."\n");
	  return join(',', @names);	  
	}
	else{
	  ($probename) = @{$self->{'probenames'}->{$arrayname}};
	}

  return $probename;
}



=head2 get_all_complete_names

  Args       : None
  Example    : my @compnames = @{$probe->get_all_complete_names()};
  Description: Retrieves all complete names for this probe. The complete name
               is a concatenation of the array name, the probeset name and the
			   probe name.
  Returntype : Arrayref of strings
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_all_complete_names {
  my $self = shift;

  my ($probeset, @result);
	my $pset = $self->probeset;

	if ($pset) {
	  $probeset = $pset->name;
	}

  if(defined $probeset){
    $probeset = ':'.$probeset.':';
  }
  else{
    $probeset = ':';
  }


  #warn "For Nimblegen this need to be Container:Seqid::probeid?";

  while ( my (undef, $array) = each %{$self->{'arrays'}} ) {
    #would have to put test in here for $self->arrays()->vendor()
    #if($array->vendor() eq "AFFY"){

	  foreach my $name ( @{$self->{'probenames'}{$array->name()}} ) {
      push @result, $array->name . $probeset . $name;
	  }
  }

  return \@result;
}



#For affy this matters as name will be different, but not for Nimblegen
#Need to consolidate this
#have get name method which throws if there is more than one array
#detects array vendor and does appropriate method

=head2 get_complete_name

  Arg [1]    : string - array name
  Example    : my $compname = $probe->get_complete_name('Array-1');
  Description: For a given array, retrieve the complete name for this probe.
  Returntype : string
  Exceptions : Throws if the array name not specified or not known for this probe
  Caller     : General
  Status     : Medium Risk

=cut

sub get_complete_name {
    my $self = shift;
    my $arrayname = shift;


	throw('Must provide and array name argument to retreive the complete name') if ! defined $arrayname;

    my $probename = $self->get_probename($arrayname);

    if (!defined $probename) {
		throw('Unknown array name');
    }

	my $probeset = $self->probeset()->name();
	$probeset .= ':' if $probeset;

	return "$arrayname:$probeset$probename";
}

=head2 probeset

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::ProbeSet
  Example    : my $probe_set = $probe->probeset();
  Description: Getter and setter of probe_set attribute for Probe objects.
  Returntype : Bio::EnsEMBL::Funcgen::ProbeSet
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probeset {
    my $self = shift;

    $self->{'probe_set'} = shift if @_;
    return $self->{'probe_set'};
}

=head2 class

  Arg [1]    : (optional) string - class
  Example    : my $class = $probe->class();
  Description: Getter and setter of class attribute for Probe
               objects e.g. CONTROL, EXPERIMENTAL
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub class {
    my $self = shift;
    $self->{'class'} = shift if @_;
    return $self->{'class'};
}

=head2 length

  Arg [1]    : (optional) int - probe length
  Example    : my $probelength = $probe->length();
  Description: Getter and setter of length attribute for Probe
               objects.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub length {
    my $self = shift;
    $self->{'length'} = shift if @_;
    return $self->{'length'};
}

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $pdesc = $probe->description();
  Description: Getter and setter of description attribute for Probe
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub description {
  my $self = shift;
  $self->{'description'} = shift if @_;
  return $self->{'description'};
}


=head2 feature_count

  Arg[0]     : recount flag
  Example    : my $num_features = $probe->feature_count();
  Description: Counts the number of ProbeFeatures associated with this Probe
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut


sub feature_count{
  my ($self, $recount) = @_;

  if($recount ||
	 (! $self->{feature_count})){
	$self->{feature_count} = $self->adaptor->db->get_ProbeFeatureAdaptor->count_probe_features_by_probe_id($self->dbID);
  }

  return $self->{feature_count};
}


### DEPRECATED ###

### ARRAY DESIGN SPECIFIC METHODS

=head2 add_Analysis_score

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : string - analysis score (as string a precision may differ between analyses)??
  Example    : $probe->add_Analysis_score($analysis, $score);
  Description: Setter for probe analysis attributes from an array design
  Returntype : None
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub add_Analysis_score{
    my ($self, $anal, $score) = @_;
    
    throw('This functionality was removed in release 74');

    if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
    }

    throw("Must provide a score to add to the probe") if ! defined $score;

    $self->{'analysis'}{$anal->dbID()} = $score;

    return;
}

=head2 add_Analysis_CoordSystem_score

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : Bio::EnsEMBL::CoordSystem
  Arg [3]    : string - analysis score (as string a precision may differ between analyses)??
  Example    : $probe->add_Analysis_CoordSystem_score($analysis, $coord_sys, $score);
  Description: Setter for coord system dependant probe analysis attributes from an array design
  Returntype : None
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub add_Analysis_CoordSystem_score{
    my ($self, $anal, $cs, $score) = @_;

  throw('This functionality was removed in release 74'); 

    if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
    }

    if(! ($cs && $cs->dbID() && $cs->isa("Bio::EnsEMBL::Funcgen::CoordSystem"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::CoordSystem");
    }

    throw("Must provide a score to add to the probe") if ! defined $score;

    $self->{'analysis_coord_system'}{$anal->dbID()}{$cs->dbID()} = $score;

    return;
}

=head2 get_score_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Analysis
  Example    : my $anal_score = $probe->get_analysis_score($analysis);
  Description: Setter for probe analysis attributes from an array design
  Returntype : string
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub get_score_by_Analysis{
  my ($self, $anal) = @_;

   throw('This functionality was removed in release 74');

  $self->get_all_design_scores() if ! defined $self->{'analysis'};

  if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
    throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
  }


  return (exists $self->{'analysis'}{$anal->dbID()}) ? $self->{'analysis'}{$anal->dbID()} : undef;
}

=head2 get_score_by_Analysis_CoordSystem

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : Bio::EnsEMBL::CoordSystem
  Arg [3]    : string - analysis score (as string a precision may differ between analyses)??
  Example    : $probe->add_analysis($analysis, $coord_sys, $score);
  Description: Setter for coord system dependant probe analysis attributes from an array design
  Returntype : None
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub get_score_by_Analysis_CoordSystem{
    my ($self, $anal, $cs) = @_;

 throw('This functionality was removed in release 74');
    $self->get_all_design_scores() if ! defined $self->{'analysis_coord_system'};

    if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
    }

    if(! ($cs && $cs->dbID() && $cs->isa("Bio::EnsEMBL::Funcgen::CoordSystem"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::CoordSystem");
    }

    my $score = undef;

    if(exists $self->{'analysis_coord_system'}{$anal->dbID()} &&
       exists $self->{'analysis_coord_system'}{$anal->dbID()}{$cs->dbID()}){
      $score = $self->{'analysis_coord_system'}{$anal->dbID()}{$cs->dbID()};
    }

    return $score;
}


=head2 get_all_design_scores

  Arg [1]    : Boolean - No fetch flag, to fetch design scores from DB, used in adaptor
               To avoid testing DB for each probe when no design scores have been added.
  Example    : my @anal_score_coordsets = @{$probe->get_all_design_scores()};
  Description: Gets all design scores as analysis_id, score and optionally coord_system_id
  Returntype : ARRAYREF
  Exceptions : throws if no fetch flag is not defined and adaptor or probe is not defined and or stored.
  Caller     : General
  Status     : at risk

=cut

#not named get_all_Analysis_scores as this would imply only non-cs dependent scores
#hence named after table, as this returns simple table records

sub get_all_design_scores{
  my ($self, $no_fetch) = @_;

 throw('This functionality was removed in release 74');

  my ($analysis_id, $cs_id, $score, @design_scores);

  if(! $no_fetch){#can assume we have none stored already due to implementation of add methods

    throw("Probe must have and adaptor to fetch design scores from the DB") if(! $self->adaptor());

    foreach my $probe_analysis(@{$self->adaptor->fetch_all_design_records($self)}){
      #we can't use the add methods here as would be cyclical
      #nor do we need extra validation

      ($analysis_id, $cs_id, $score) = @$probe_analysis;

      if($cs_id){
	$self->{'analysis_coord_system'}{$analysis_id}{$cs_id} = $score;
      }else{
	$self->{'analysis'}{$analysis_id} = $score;
      }
    }
  }

  #populate array from attrs
  if(exists $self->{'analysis_coord_system'}){

    foreach $analysis_id(keys %{$self->{'analysis_coord_system'}}){

      foreach $cs_id(keys %{$self->{'analysis_coord_system'}{$analysis_id}}){
	push @design_scores, [$analysis_id, $self->{'analysis_coord_system'}{$analysis_id}{$cs_id}, $cs_id];
      }
    }
  }

  if(exists $self->{'analysis'}){

    foreach $analysis_id(keys %{$self->{'analysis'}}){

      push @design_scores, [$analysis_id, $self->{'analysis'}{$analysis_id}];
    }
  }


  return \@design_scores;

}


#do we need get_all methods for Analysis and Analysis_CoordSystem?
#maybe if we split into another Class and Adaptor

1;

