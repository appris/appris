#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor - A database adaptor for fetching and
storing Probe objects.

=head1 SYNOPSIS

my $probe_adaptor = $db->get_ProbeAdaptor();

my $probe = $probe_adaptor->fetch_by_array_probe_probeset('Array-1', 'Probe-1');

=head1 DESCRIPTION

The ProbeAdaptor is a database adaptor for storing and retrieving
Probe objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Probe
Bio::EnsEMBL::Funcgen::ProbeFeature
Bio::EnsEMBL::Funcgen::ProbeSet
Bio::EnsEMBL::Funcgen::ArrayChip
Bio::EnsEMBL::Funcgen::Array

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_by_array_probe_probeset_name

  Arg [1]    : string - name of array
  Arg [2]    : string - name of probe
  Arg [3]    : (optional) string - name of probeset
  Example    : my $probe = $opa->fetch_by_array_probeset_probe('Array-1', 'Probe-1', 'ProbeSet-1');
  Description: Returns a probe given a combination of array name, probeset and
               probe name. This will uniquely define an Affy probe. Only one
			   probe is ever returned.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : throws if array or probe name not defined
  Caller     : General
  Status     : At Risk - rename to fetch_by_probe_array_probeset_name?

=cut

#This does not currently capture on plate replicate probes with different names
#Only returns the record corresponding to the given name and not the other replicate name

sub fetch_by_array_probe_probeset_name {
	my ($self, $array_name, $probe_name, $probeset_name) = @_;

	if(! (defined $array_name && defined $probe_name)){
	  throw('You must provide at least and array and probe name');
	}

	my $tables = 'probe p, array_chip ac, array a';
	$tables .= ', probe_set ps' if defined $probeset_name;

	my $sql = "SELECT distinct(p.probe_id) FROM $tables WHERE a.name=? and a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and p.name=?";
	$sql .= ' AND p.probe_set_id=ps.probe_set_id and ps.name=?' if defined $probeset_name;
	my $sth = $self->db->dbc->prepare($sql);
	$sth->bind_param(1, $array_name,    SQL_VARCHAR);
	$sth->bind_param(2, $probe_name,    SQL_VARCHAR);
	$sth->bind_param(3, $probeset_name, SQL_VARCHAR) if defined $probeset_name;
	$sth->execute;

	#This should only return one result
	#The only possible way this would not return one result
	#is if an identically named array(:probeset):probe which had a different sequence
	#As Import array would separate these based on the sequence hash
	my ($dbid) = $sth->fetchrow_array;

	return (defined $dbid) ? $self->fetch_by_dbID($dbid) : undef;
}




=head2 fetch_all_by_name

  Arg [1]    : string - probe name
  Example    : my @probes = @{$opa->fetch_all_by_name('Probe1')};
  Description: Convenience method to re-instate the functionality of
               $core_dbentry_adpator->fetch_all_by_External_name('probe_name');
               WARNING: This may not be the probe you are expecting as
               probe names are not unqiue across arrays and vendors.
               These should ideally be validated using the attached array
               information or alternatively use fetch_by_array_probe_probeset_name
               Returns a probe with the given name.
  Returntype : Arrayref
  Exceptions : Throws if name not passed
  Caller     : General
  Status     : At Risk

=cut


sub fetch_all_by_name{
  my ($self, $name) = @_;


  throw('Must provide a probe name argument') if ! defined $name;

  $self->bind_param_generic_fetch($name, SQL_VARCHAR);

  return $self->generic_fetch('p.name=?');
}




=head2 fetch_all_by_ProbeSet

  Arg [1]    : Bio::EnsEMBL::ProbeSet
  Example    : my @probes = @{$opa->fetch_all_by_ProbeSet($pset)};
  Description: Fetch all probes in a particular ProbeSet.
  Returntype : Listref of Bio::EnsEMBL::Probe objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ProbeSet {
	my ($self, $probeset) = @_;

	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ProbeSet', $probeset);
	return $self->generic_fetch('p.probe_set_id = '.$probeset->dbID);
}



=head2 fetch_all_by_Array

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array
  Example    : my @probes = @{$opa->fetch_all_by_Array($array)};
  Description: Fetch all probes on a particular array.
  Returntype : Listref of Bio::EnsEMBL::Probe objects.
  Exceptions : throws if arg is not valid or stored
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Array {
  my $self  = shift;
  my $array = shift;

   if(! (ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array') && $array->dbID())){
     throw('Need to pass a valid stored Bio::EnsEMBL::Funcgen::Array');
   }

  #get all array_chip_ids, for array and do a multiple OR statement with generic fetch

  return $self->generic_fetch('p.array_chip_id IN ('.join(',', @{$array->get_array_chip_ids()}).')');
}

=head2 fetch_all_by_ArrayChip

  Arg [1]    : Bio::EnsEMBL::Funcgen::ArrayChip
  Example    : my @probes = @{$opa->fetch_all_by_ArrayChip($array_chip)};
  Description: Fetch all probes on a particular ArrayChip.
  Returntype : Listref of Bio::EnsEMBL::Probe objects.
  Exceptions : throw is arg is not valid or stored
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ArrayChip {
  my $self  = shift;
  my $array_chip = shift;

  if(! (ref($array_chip) && $array_chip->isa('Bio::EnsEMBL::Funcgen::ArrayChip') && $array_chip->dbID())){
    throw('Need to pass a valid stored Bio::EnsEMBL::Funcgen::ArrayChip');
  }

  return $self->generic_fetch("p.array_chip_id =".$array_chip->dbID());
}


=head2 fetch_by_ProbeFeature

  Arg [1]    : Bio::EnsEMBL::Funcgen::ProbeFeature
  Example    : my $probe = $opa->fetch_by_ProbeFeature($feature);
  Description: Returns the probe that created a particular feature.
  Returntype : Bio::EnsEMBL::Probe
  Exceptions : Throws if argument is not a Bio::EnsEMBL::Funcgen::ProbeFeature object
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_ProbeFeature {
  my $self    = shift;
  my $feature = shift;

  if (! ref($feature) ||
      ! $feature->isa('Bio::EnsEMBL::Funcgen::ProbeFeature') ||
      ! $feature->{'probe_id'}
     ) {
    throw('fetch_by_ProbeFeature requires a stored Bio::EnsEMBL::Funcgen::ProbeFeature object');
  }

  return $self->fetch_by_dbID($feature->{'probe_id'});
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
  return (['probe', 'p']);
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
  return qw( p.probe_id p.probe_set_id p.name p.length p.array_chip_id p.class p.description);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Probe objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Probe objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $current_dbid, $arraychip_id, $probe_id, $probe_set_id, $name, $class, $probelength, $desc);
	my ($array, %array_cache, %probe_set_cache);

	$sth->bind_columns(\$probe_id, \$probe_set_id, \$name, \$probelength, \$arraychip_id, \$class, \$desc);


	#Complex query extension
	#We want the arrays, array_chip and probeset information setting
	#So the probe feature zmenu can just do one query to populate the zmenu unstead of 4
	#Let's just do this with probset to start with as this is more simple
	#The object creation for the linked adaptors need to be commodotised
	#So we can say something like $probeset_adaptor->create_obj_from_sth_args
	#Caches will need to be built in the calling adaptor rather than the true object adaptor

	#No group required as we will always want intermediate data
	#Therefore cannot use this in combined with simple extension???
	#Unless we explicitly state group by primary keys.

	#Need to build array of adaptor column hashes(order important)
	#Need to build and array of bound columns dependant
	#Then we need to parse output dependant on primary keys of each table
	#So we would need bolumns bound to hash values

	#We need a way to define the extended tables
	#Pass param hash to caller which uses BaseAdaptor method to add tables and columns
	#This would have to take into account anything added by the caller

	#Maybe the better way of doing this test would be to include all probes in a probeset?
	#Or maybe all probe_features for a probe?




	my $probe;

	while ( $sth->fetch() ) {

		#warn("Need to sort array cacheing, have redundant cache!!");
		#This is nesting array and probeset objects in probe!


		#Is this required? or should we lazy load this?
		#Should we also do the same for probe i.e. nest or lazy load probeset
		#Setting here prevents, multiple queries, but if we store the array cache in the adaptor we can overcome this
		#danger of eating memory here, but it's onld the same as would be used for generating all the probesets
		#what about clearing the cache?
		#also as multiple array_chips map to same array, cache would be redundant
		#need to store only once and reference.
		#have array_cache and arraychip_map
		#arraychip_map would give array_id which would be key in array cache
		#This is kinda reinventing the wheel, but reducing queries and redundancy of global cache
		#cache would never be populated if method not called
		#there for reducing calls and memory, increasing speed of generation/initation
		#if method were called
		#would slightly slow down processing, and would slightly increase memory as cache(small as non-redundant)
		#and map hashes would persist


	  ####MAKE THIS LAZY LOADED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  #Can we even do this given we won't then have the array context?
	  #We should just force this for efficiency and make people keep the array if they ever want to use that info?
	  #Will this affect any other methods?



	  #Can we not change this to an ArrayChip cache and just reimplement the array method?

		$array = $array_cache{$arraychip_id} || $self->db->get_ArrayAdaptor()->fetch_by_array_chip_dbID($arraychip_id);


		#I don't think we need this?  Certainly not for storing

		#$probe_set = $probe_set_cache{$probe_set_id} || $self->db->get_ArrayAdaptor()->fetch_by_array_chip_dbID($arraychip_id);
		#probe_set cache would be substantially bigger!!
		#potentially as many as the probes

		#Just build cache and nest for now,may want to just return ID and lazy load

		#This is a prime target for compound query extension
		#Either extend query by default and nest probe_set
		#Or lazy load probeset using cache somehow?
		#Use persistant probeset cache in ProbeSetAdaptor for dbID/Probe style queries

		my ($probeset);

		if($probe_set_id){
			$probeset = $probe_set_cache{$probe_set_id} || $self->db->get_ProbeSetAdaptor()->fetch_by_dbID($probe_set_id);
		}

		if (!$current_dbid || $current_dbid != $probe_id) {
			# New probe

			#UC??? or does rearrange handle this?

			$probe = Bio::EnsEMBL::Funcgen::Probe->new
			  (
			   -dbID          => $probe_id,
			   -name          => $name,
			   -array_chip_id => $arraychip_id,
			   -array         => $array,
			   -probe_set     => $probeset,
			   -length        => $probelength,
			   -class         => $class,
			   -description   => $desc,
			   -adaptor       => $self,
			);
			push @result, $probe;
			$current_dbid = $probe_id;
		} else {
		  # Extend existing probe
		  # Probe methods depend on preloading of Array objects
		  $probe->add_array_chip_probename($arraychip_id, $name, $array);
		}
	}
	return \@result;
}

=head2 store

  Arg [1]    : List of Bio::EnsEMBL::Funcgen::Probe objects
  Example    : $opa->store($probe1, $probe2, $probe3);
  Description: Stores given Probe objects in the database. Should only be
               called once per probe because no checks are made for duplicates
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : ARRAYREF
  Exceptions : Throws if arguments are not Probe objects
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my ($self, @probes) = @_;
  my $db = $self->db();
  throw('Must call store with a list of Probe objects') if (scalar @probes == 0);

  my $new_sth = $self->prepare
    (
     "INSERT INTO probe( probe_set_id, name, length, array_chip_id, class, description)".
     "VALUES (?, ?, ?, ?, ?, ?)"
    );
  
  my $existing_sth = $self->prepare
    (
     "INSERT INTO probe( probe_id, probe_set_id, name, length, array_chip_id, class, description )".
     "VALUES (?, ?, ?, ?, ?, ?, ?)"
    );


 PROBE: foreach my $probe (@probes) {
     
    if ( !ref $probe || ! $probe->isa('Bio::EnsEMBL::Funcgen::Probe') ) {
      throw("Probe must be an Probe object ($probe)");
    }
    
    if ( $probe->is_stored($db) ) {
      warning('Probe [' . $probe->dbID() . '] is already stored in the database');
      next PROBE;
    }
    
    # Get all the arrays this probe is on and check they're all in the database
    my %array_hashes;
    
    foreach my $ac_id (keys %{$probe->{'arrays'}}) {
      
      if (defined ${$probe->{'arrays'}}{$ac_id}->dbID()) {
        #Will this ever work as generally we're creating from scratch 
        #and direct access to keys above by passes DB fetch
        $array_hashes{$ac_id} = $probe->{'arrays'}{$ac_id};
      }
    }

    throw('Probes need attached arrays to be stored in the database') if ( ! %array_hashes );

    # Insert separate entry (with same probe_id) for each array/array_chip the probe is on

    foreach my $ac_id (keys %array_hashes) {
      my $ps_id = (defined $probe->probeset()) ? $probe->probeset()->dbID() : undef;

      foreach my $name (@{$probe->get_all_probenames($array_hashes{$ac_id}->name)}) {
      
        if ($probe->dbID) {    # Already stored
          $existing_sth->bind_param(1, $probe->dbID,        SQL_INTEGER);
          $existing_sth->bind_param(2, $ps_id,              SQL_INTEGER);
          $existing_sth->bind_param(3, $name,               SQL_VARCHAR);
          $existing_sth->bind_param(4, $probe->length(),    SQL_INTEGER);
          $existing_sth->bind_param(5, $ac_id,              SQL_INTEGER);
          $existing_sth->bind_param(6, $probe->class(),     SQL_VARCHAR);
          $existing_sth->bind_param(7, $probe->description, SQL_VARCHAR);
          $existing_sth->execute();
        } else {
          # New probe
          $new_sth->bind_param(1, $ps_id,              SQL_INTEGER);
          $new_sth->bind_param(2, $name,               SQL_VARCHAR);
          $new_sth->bind_param(3, $probe->length(),    SQL_INTEGER);
          $new_sth->bind_param(4, $ac_id,              SQL_INTEGER);
          $new_sth->bind_param(5, $probe->class(),     SQL_VARCHAR);
          $new_sth->bind_param(6, $probe->description, SQL_VARCHAR);
          $new_sth->execute();
          $probe->dbID($self->last_insert_id);
          $probe->adaptor($self);
        }
      }
    }
  }
  
  return \@probes;
}



=head2 fetch_all_design_scores

  Arg [1]    : Bio::EnsEMBL::Funcgen::Probe
  Example    : my @probe_analyses = @{$pa->fetch_all_design_scores($probe)};
  Description: Fetchs all probe design analysis records as analysis_id, score and coord_system_id
  Returntype : ARRAYREF
  Exceptions : throws if not passed a valid stored Probe
  Caller     : General
  Status     : at risk

=cut

sub fetch_all_design_scores{
  my ($self, $probe) = @_;

  if(! ($probe && $probe->isa('Bio::EnsEMBL::Funcgen::Probe') && $probe->dbID())){
    throw('Must pass a valid stored Bio::EnsEMBL::Funcgen::Probe');
  }

  my $sql = 'SELECT analysis_id, score, coord_system_id from probe_design WHERE probe_id='.$probe->dbID.';';
  return @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
}


1;

