#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ProbeFeatureAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::ProbeFeatureAdaptor - A database adaptor for fetching and
storing ProbeFeature objects.

=head1 SYNOPSIS

my $ofa = $db->get_ProbeFeatureAdaptor();

my $features = $ofa->fetch_all_by_Probe($probe);
$features = $ofa->fetch_all_by_Slice_arrayname($slice, 'Array-1', 'Array-2');

=head1 DESCRIPTION

The ProbeFeatureAdaptor is a database adaptor for storing and retrieving
ProbeFeature objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::ProbeFeature

=cut


package Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;#DBI sql_types import
#use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);


my $true_final_clause = ' ORDER BY pf.seq_region_id, pf.seq_region_start, pf.probe_feature_id';
#Could drop pf.probe_feature_id from the ORDER as is implicit from the group?
#still uses filesort for ac clause
my $final_clause = $true_final_clause;


=head2 fetch_all_by_Probe

  Arg [1]    : Bio::EnsEMBL::Funcgen::Probe
  Example    : my $features = $ofa->fetch_all_by_Probe($probe);
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::PasteFeature objects
  Exceptions : Throws if argument is not a stored Probe object
  Caller     : Probe->get_all_ProbeFeatures()
  Status     : At Risk

=cut

sub fetch_all_by_Probe {
  my ($self, $probe, $coord_systems) = @_;

  if (! (ref($probe) && $probe->isa('Bio::EnsEMBL::Funcgen::Probe'))) {
    throw('fetch_all_by_Probe requires a Bio::EnsEMBL::Funcgen::Probe object');
  }

  if ( !defined $probe->dbID() ) {
    throw('fetch_all_by_Probe requires a stored Bio::EnsEMBL::Funcgen::Probe object');
  }

  return $self->fetch_all_by_probe_id($probe->dbID, $coord_systems);
}

=head2 fetch_all_by_probe_id

  Arg [1]    : int - Probe dbID
  Example    : my @features = @{$ofa->fetch_all_by_Probe_id($pid)};
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::PasteFeature objects
  Exceptions : Throws if argument not defined
  Caller     : Probe->get_all_ProbeFeatures()
  Status     : At Risk

=cut

sub fetch_all_by_probe_id {
  my ($self, $pid, $coord_systems) = @_;

  if ( ! defined $pid ) {
    throw('Need to specify a probe _id');
  }

  my @cs_ids = @{$self->_get_coord_system_ids($coord_systems)};
  $self->_tables([['seq_region', 'sr']]);

  my $cs_ids = join(', ', @cs_ids);
  my $constraint = " pf.probe_id=$pid AND pf.seq_region_id=sr.seq_region_id and sr.coord_system_id IN ($cs_ids)";
  $final_clause = ' GROUP by pf.probe_feature_id '.$final_clause;


  my $features = $self->generic_fetch($constraint);
  $self->reset_true_tables;
  $final_clause = $true_final_clause;


  return $features;
}



=head2 fetch_all_by_probeset_name

  Arg [1]    : String - probeset name
  Arg [2]    : ARRAYREF (optional) - Bio::EnsEMBL::CoordSystem objects
  Example    : my $features = $ofa->fetch_all_by_probeset('Set-1');
  Description: Fetchs all features that a given probeset creates.
  Returntype : Listref of Bio::EnsEMBL::ProbeFeature objects
  Exceptions : Throws if no probeset argument
  Caller     : General
  Status     : At Risk - add vendor/class to this?

=cut


sub fetch_all_by_probeset_name {
	my ($self, $probeset, $coord_systems) = @_;

	if (! $probeset) {
	  throw('fetch_all_by_probeset requires a probeset name argument');
	}

	#Restrict to default coord_systems
	#Can we remove the need for this by restricting the sr cache to default entries?
	my @cs_ids = @{$self->_get_coord_system_ids($coord_systems)};
    $self->_tables([['probe_set', 'ps'], ['seq_region', 'sr']]);

	#Need to protect against SQL injection here due to text params
	my $cs_ids = join(', ', @cs_ids);
	my $constraint = " ps.name=? AND ps.probe_set_id=p.probe_set_id AND pf.seq_region_id=sr.seq_region_id and sr.coord_system_id IN ($cs_ids)";
	$final_clause = ' GROUP by pf.probe_feature_id '.$final_clause;

	$self->bind_param_generic_fetch($probeset,  SQL_VARCHAR);
	my $features = $self->generic_fetch($constraint);
    $self->reset_true_tables;
	$final_clause = $true_final_clause;

	return $features;
}


=head2 fetch_all_by_ProbeSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::ProbeSet
  Arg [2]    : ARRAYREF (optional) - Bio::EnsEMBL::CoordSystem objects
  Example    : my @features = @{$probe_feature_adaptor->fetch_all_by_ProbeSet($pset)};
  Description: Fetches all ProbeFeatures from a given ProbeSet.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if no probeset argument
  Caller     : General
  Status     : At Risk - add vendor/class to this?

=cut

sub fetch_all_by_ProbeSet {
	my ($self, $pset, $coord_systems) = @_;
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ProbeSet', $pset);
	#Restrict to default coord_systems
	#Can we remove the need for this by restricting the sr cache to default entries?

	my @cs_ids = @{$self->_get_coord_system_ids($coord_systems)};
    $self->_tables([['seq_region', 'sr']]);

	my $cs_ids = join(', ', @cs_ids);
	my $constraint = ' p.probe_set_id='.$pset->dbID." AND pf.seq_region_id=sr.seq_region_id and sr.coord_system_id IN ($cs_ids)";
	$final_clause = ' GROUP by pf.probe_feature_id '.$final_clause;

	my $features = $self->generic_fetch($constraint);
    $self->reset_true_tables;
	$final_clause = $true_final_clause;

	return $features;
}


=head2 fetch_all_by_Slice_ExperimentalChips

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : ARRAY ref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Example    : my $features = $pfa->fetch_all_by_Slice_ExperimentalChips($slice, \@echips);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the given ExperimentalChips.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if args not valid
  Caller     :
  Status     : At Risk

=cut

sub fetch_all_by_Slice_ExperimentalChips {
  my ($self, $slice, $exp_chips) = @_;

  my %nr;

  foreach my $ec(@$exp_chips){

    throw("Need pass listref of valid Bio::EnsEMBL::Funcgen::ExperimentalChip objects")
      if ! $ec->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip");

    $nr{$ec->array_chip_id()} = 1;
  }

  my $constraint = " p.array_chip_id IN (".join(", ", keys %nr).") AND p.probe_id = pf.probe_id ";

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}



#Need to Group in the following methods as we may get array_chip
#to probe product if probe is presenton >1 array_chip.
#This will be slowing as GROUP implies order
#Does _objects_from_sth handle this without assuming order?



=head2 fetch_all_by_Slice_array_vendor

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : string - array name e.g. HG-U133A
  Arg [3]    : string - vendor e.g. AFFY
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_array_vendor($slice, $array_name, $vendor_name);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified array.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if no array name is provided
  #Caller     : Slice->get_all_ProbesFeatures()
  Status     : At Risk

=cut

sub fetch_all_by_Slice_array_vendor {
	my ($self, $slice, $array, $vendor) = @_;

	if(! ($array && $vendor)){
	  throw('You must provide and array name and a vendor name');
	}

    $self->_tables([['array', 'a'], ['array_chip', 'ac']]);

	#Need to protect against SQL injection here due to text params
	my $constraint = ' a.name=? and a.vendor=? and a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id';
	$final_clause  = ' GROUP by pf.probe_feature_id '.$final_clause;
	$self->bind_param_generic_fetch($array,  SQL_VARCHAR);
	$self->bind_param_generic_fetch($vendor, SQL_VARCHAR);

	my $features  = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
    $self->reset_true_tables;
	$final_clause = $true_final_clause;

	return $features;
}


=head2 fetch_all_by_Slice_Array

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Funcgen::Array
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $pfa->fetch_all_by_Slice_Array($slice, $array);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the given Array.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if no array name is provided
  Caller     :
  Status     : At Risk

=cut

sub fetch_all_by_Slice_Array {
  my ($self, $slice, $array) = @_;

  throw("Need pass a valid stored Bio::EnsEMBL::Funcgen::Array object")
	if (! (ref($array) && $array->isa("Bio::EnsEMBL::Funcgen::Array") && $array->dbID));

  $self->_tables([['array_chip', 'ac']]);
  my $constraint = ' ac.array_id='.$array->dbID.' and ac.array_chip_id=p.array_chip_id ';
  $final_clause = ' GROUP by pf.probe_feature_id '.$final_clause;

  my $features  = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
  $self->reset_true_tables;
  $final_clause = $true_final_clause;

  return $features;
}


=head2 fetch_all_by_Slice_Arrays

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects
  Arg [3]    : HASHREF - optional params hash e.g. {logic_name => 'AFFY_ProbeTranscriptAlign'}
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $pfa->fetch_all_by_Slice_Arrays($slice, \@arrays);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the given Arrays.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if ARRAYREF of arrays is not provided
  Caller     :
  Status     : At Risk

=cut

sub fetch_all_by_Slice_Arrays{
  my ($self, $slice, $arrays, $params) = @_;

  my $logic_name;
  $logic_name = $params->{'logic_name'} if exists ${$params}{'logic_names'};


  if(!(ref($arrays) eq 'ARRAY' &&  @$arrays)){
	throw('Must pass an ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects');
  }

  my $array_ids = join(',', (map $_->dbID, @$arrays));

  $self->_tables([['array_chip', 'ac']]);
  my $constraint = " ac.array_id IN ($array_ids) and ac.array_chip_id=p.array_chip_id ";

  $final_clause = ' GROUP by pf.probe_feature_id '.$final_clause;
  my $features  = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);
  $self->reset_true_tables;
  $final_clause = $true_final_clause;

  return $features;
}


=head2 fetch_Iterator_by_Slice_Arrays

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects
  Arg [3]    : HASHREF - optional params hash e.g. {logic_name => 'AFFY_ProbeTranscriptAlign'}
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $pfa->fetch_Iterator_by_Slice_Arrays($slice, \@arrays);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the given Array.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : Throws if ARRAYREF of arrays is not provided
  Caller     :
  Status     : At Risk

=cut

sub fetch_Iterator_by_Slice_Arrays{
  my ($self, $slice, $arrays, $params) = @_;


  return $self->fetch_Iterator_by_Slice_method
	($self->can('fetch_all_by_Slice_Arrays'),
	 [$slice, $arrays, $params],
	 0,#Slice idx
	 #500 #chunk length
	);
}

=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  return ([ 'probe_feature', 'pf' ], [ 'probe',   'p' ]);
}



=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  return qw(
            pf.probe_feature_id  pf.seq_region_id
			pf.seq_region_start  pf.seq_region_end
			pf.seq_region_strand pf.probe_id
			pf.analysis_id	   pf.mismatches
			pf.cigar_line        p.name
			p.probe_set_id
		    );
}

=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut
sub _default_where_clause {
	my $self = shift;

	return 'pf.probe_id = p.probe_id';
}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an ORDER BY clause. Sorting by probe_feature_id would be
			   enough to eliminate duplicates, but sorting by location might
			   make fetching features on a slice faster.
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : At Risk

=cut


sub _final_clause {
	return $final_clause;
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates ProbeFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::ProbeFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;

	#For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	#So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.



	# This code is ugly because caching is used to improve speed
	my ($seq_region_id);
	my $sa = $self->db->get_SliceAdaptor();
	$sa = $dest_slice->adaptor->db->get_SliceAdaptor() if($dest_slice);#don't really need this if we're using DNADBSliceAdaptor?

	#Some of this in now probably overkill as we'll always be using the DNADB as the slice DB
	#Hence it should always be on the same coord system, unless we're projecting

	my $aa = $self->db->get_AnalysisAdaptor();
	my @features;
	my (%analysis_hash, %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
	    $probe_feature_id,  $efg_seq_region_id,
	    $seq_region_start,  $seq_region_end,
	    $seq_region_strand, $mismatches,
		$probe_id,    	    $analysis_id,
		$probe_name,	    $cigar_line,
		$probeset_id
	);
	$sth->bind_columns(
					   \$probe_feature_id,  \$efg_seq_region_id,
					   \$seq_region_start,  \$seq_region_end,
					   \$seq_region_strand, \$probe_id,
					   \$analysis_id,       \$mismatches,
					   \$cigar_line,        \$probe_name,
					   \$probeset_id
	);

	my ($asm_cs, $cmp_cs, $asm_cs_name, $asm_cs_vers ,$cmp_cs_name, $cmp_cs_vers);
	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my ($dest_slice_start, $dest_slice_end, $dest_slice_strand);
	my ($dest_slice_length, $dest_slice_sr_name);
	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}

	#This has already been done by
	#build seq_region_cache based on slice
	#$self->build_seq_region_cache_by_Slice($slice);


	my $last_pfid;

	FEATURE: while ( $sth->fetch() ) {
		  #Need to build a slice adaptor cache here?
		  #Would only ever want to do this if we enable mapping between assemblies??
		  #Or if we supported the mapping between cs systems for a given schema_build, which would have to be handled by the core api

		#This is only required due to multiple records being returned
		#From nr probe entries due to being present on multiple ArrayChips
		#Group instead?
		next if($last_pfid && ($last_pfid == $probe_feature_id));
		$last_pfid = $probe_feature_id;

		#get core seq_region_id
		$seq_region_id = $self->get_core_seq_region_id($efg_seq_region_id);

		if(! $seq_region_id){
		  #warn "Cannot get slice for eFG seq_region_id $efg_seq_region_id for probe_feature $probe_feature_id\n".
		  #"The region you are using is not present in the current dna DB";
		  #This can happen as non slice fetches only restrict on cs_id
		  #Hence for the non-versioned cs's there may be seq_regions which have
		  #disappeared in the current assembly and hence won't be in the sr cache.
		  #We could get around this by adding an sr_id IN(all the sr_ids from this DB)
		  #but this will most likely just slow things down for data which is not present on
		  #just one assembly
		  #So preferable to clear old data!
		  next;
		}


		# Get the analysis object
		my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);

		# Get the slice object
		my $slice = $slice_hash{'ID:'.$seq_region_id};

		if (!$slice) {
			$slice                            = $sa->fetch_by_seq_region_id($seq_region_id);


			if(! $slice){
			  warn "Cannot get slice for seq_region_id $seq_region_id for probe_feature $probe_feature_id";
			}



			$slice_hash{'ID:'.$seq_region_id} = $slice;
			$sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
			$sr_cs_hash{$seq_region_id}       = $slice->coord_system();
		}

		#need to check once more here as it may not be in the DB,
		#i.e. a supercontig(non-versioned) may have been deleted between releases


		my $sr_name = $sr_name_hash{$seq_region_id};
		my $sr_cs   = $sr_cs_hash{$seq_region_id};

		# Remap the feature coordinates to another coord system if a mapper was provided
		if ($mapper) {

			throw("Not yet implmented mapper, check equals are Funcgen calls too!");

			($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand)
				= $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);

			# Skip features that map to gaps or coord system boundaries
			next FEATURE if !defined $sr_name;

			# Get a slice in the coord system we just mapped to
			if ( $asm_cs == $sr_cs || ( $cmp_cs != $sr_cs && $asm_cs->equals($sr_cs) ) ) {
				$slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"}
					||= $sa->fetch_by_region($cmp_cs_name, $sr_name, undef, undef, undef, $cmp_cs_vers);
			} else {
				$slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"}
					||= $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef, $asm_cs_vers);
			}
		}

		# If a destination slice was provided convert the coords
		# If the destination slice starts at 1 and is forward strand, nothing needs doing
		if ($dest_slice) {
			unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {
				if ($dest_slice_strand == 1) {
					$seq_region_start = $seq_region_start - $dest_slice_start + 1;
					$seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
				} else {
					my $tmp_seq_region_start = $seq_region_start;
					$seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
					$seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
					$seq_region_strand      *= -1;
				}
			}

			# Throw away features off the end of the requested slice
			next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
				|| ( $dest_slice_sr_name ne $sr_name );

			$slice = $dest_slice;
		}

		push @features, Bio::EnsEMBL::Funcgen::ProbeFeature->new_fast
		  ({
			'start'         => $seq_region_start,
			'end'           => $seq_region_end,
			'strand'        => $seq_region_strand,
			'slice'         => $slice,
			'analysis'      => $analysis,#we should lazy load this from analysis adaptor cache?
			'adaptor'       => $self,
			'dbID'          => $probe_feature_id,
			'mismatchcount' => $mismatches,
			'cigar_string'    => $cigar_line,
			'probe_id'      => $probe_id,
			#Do these need to be private?
			'_probeset_id'  => $probeset_id,#Used for linking feature glyphs
			'_probe_name'   => $probe_name,#?? There can be >1. Is this for array design purposes?
		   } );



	  }

	return \@features;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given ProbeFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if a list of ProbeFeature objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @ofs) = @_;

	if (scalar(@ofs) == 0) {
		throw('Must call store with a list of ProbeFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO probe_feature (
			seq_region_id,  seq_region_start,
			seq_region_end, seq_region_strand,
          	probe_id,  analysis_id,
			mismatches, cigar_line
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();
	my $analysis_adaptor = $db->get_AnalysisAdaptor();

	FEATURE: foreach my $of (@ofs) {

		if( ! ref $of || ! $of->isa('Bio::EnsEMBL::Funcgen::ProbeFeature') ) {
			throw('Feature must be an ProbeFeature object');
		}

		if ( $of->is_stored($db) ) {
			warn('ProbeFeature [' . $of->dbID() . '] is already stored in the database');
			next FEATURE;
		}

		if ( !defined $of->analysis() ) {
			throw('An analysis must be attached to the ProbeFeature objects to be stored.');
		}

		# Store the analysis if it has not been stored yet
		if ( !$of->analysis->is_stored($db) ) {
			$analysis_adaptor->store( $of->analysis() );
		}

		my $seq_region_id;
		($of, $seq_region_id) = $self->_pre_store($of);

		$sth->bind_param(1, $seq_region_id,        SQL_INTEGER);
		$sth->bind_param(2, $of->start(),          SQL_INTEGER);
		$sth->bind_param(3, $of->end(),            SQL_INTEGER);
		$sth->bind_param(4, $of->strand(),         SQL_TINYINT);
		$sth->bind_param(5, $of->probe_id(),       SQL_INTEGER);
		$sth->bind_param(6, $of->analysis->dbID(), SQL_INTEGER);
		$sth->bind_param(7, $of->mismatchcount(),  SQL_TINYINT);
		$sth->bind_param(8, $of->cigar_string(),   SQL_VARCHAR);

		$sth->execute();
		$of->dbID( $self->last_insert_id );
		$of->adaptor($self);


	}

	#No need to return this really as the dbID and adaptor has been
	#updated in the passed arrays of features via the object
	#reference
	return \@ofs
}



#Probe cache methods?

=head2 reassign_feature_to_probe

  Arg[0]     : ARRAYREF - feature dbIDs to reassign
  Arg[1]     : int - probe dbID to reassign to
  Example    : $ofa->reassign_feature_to_probe(\@fids, $pid);
  Description: Update features to link to given probe dbID
  Returntype : None
  Exceptions : Throws is args not met
  Caller     : Importer
  Status     : At Risk

=cut

sub reassign_feature_to_probe{
	my ($self, $fids_ref, $pid) = @_;

	if(! @$fids_ref || ! $pid){
	  throw('Need to pass a ref to an array of feature ids and a probe id to reassign to');
	}

	my $cmd = 'UPDATE probe_feature SET probe_id='.$pid.' WHERE probe_feature_id IN ('.join(',', @$fids_ref).')';
	$self->db->dbc->do($cmd);

	return;
}

=head2 delete_features

  Arg[0]     : ARRAYREF - feature dbIDs to reassign
  Example    : $pfa->delete_feature(\@fids);
  Description: Deletes feature with given probe_feature_ids
  Returntype : None
  Exceptions : Throws if not arg defines
  Caller     : Importer
  Status     : At Risk

=cut


#This does not rollback associated xrefs!

sub delete_features{
	my ($self, $fids_ref) = @_;

	if(! @$fids_ref){
	  throw('Need to pass a ref to an array of feature ids');
	}

	my $cmd = 'DELETE from probe_feature WHERE probe_feature_id IN ('.join(',', @$fids_ref).')';
	$self->db->dbc->do($cmd);

	return;
}


=head2 count_probe_features_by_probe_id

  Arg [1]    : string/int - id to count
  Example    : my $probe_feature_count = $pfa->count_features_by_probe_id($probe_id);
  Description: Returns a count of ProbeFeatures for a given probe id
  Returntype : string/int - count of features
  Exceptions : None
  Caller     : FeatureAdaptors
  Status     : At risk

=cut


sub count_probe_features_by_probe_id {
  my ($self, $probe_id) = @_;

  return $self->count_features_by_field_id('probe_id', $probe_id);
}

### DEPRECATED METHODS ###

sub fetch_all_by_probeset { #deprecated in v68
  my ($self, @args) = @_;

  deprecate('This method is deprecated, please use fetch_all_by_probeset_name or fetch_all_by_ProbeSet');

  return $self->fetch_all_by_probeset_name(@args);
}


1;

