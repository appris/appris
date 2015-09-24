#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);

  my $csa = $db->get_CoordSystemAdaptor();

  #
  # Fetch by name, schema_build and version(opt).
  #

  $cs = $csa->fetch_by_name_schema_build_version('chromosome', '39_36a', 'NCBI36');

  #As this is a multi-assembly DB, we have to accomodate the idea of schema versions, which will
  #enable a mapping from the feature table back to assembly/core DB of origin.


  #Old core methods, some may not work as they assume that there will only be one default version
  #where are there maybe multiple default versions, one for each assembly/schema_build

  #
  # Get all coord systems in the database:
  #
  foreach my $cs (@{$csa->fetch_all()}) {
    print $cs->name, ' ',  $cs->version, "\n";
  }

  #
  # Fetching by name:
  #

  #use the default version of coord_system 'chromosome' (e.g. NCBI33):
  $cs = $csa->fetch_by_name('chromosome');

  #get an explicit version of coord_system 'chromosome':
  $cs = $csa->fetch_by_name('chromsome', 'NCBI34');

  #get all coord_systems of name 'chromosome':
  foreach $cs (@{$csa->fetch_all_by_name('chromosome')}) {
     print $cs->name, ' ', $cs->version, "\n";
  }

  #
  # Fetching by rank:
  #
  $cs = $csa->fetch_by_rank(2);

  #
  # Fetching the pseudo coord system 'toplevel'
  #

  #Get the default top_level coord system:
  $cs = $csa->fetch_top_level();

  #can also use an alias in fetch_by_name:
  $cs = $csa->fetch_by_name('toplevel');

  #can also request toplevel using rank=0
  $cs = $csa->fetch_by_rank(0);

  #
  # Fetching by sequence level:
  #

  #Get the coord system which is used to store sequence:
  $cs = $csa->fetch_sequence_level();

  #can also use an alias in fetch_by_name:
  $cs = $csa->fetch_by_name('seqlevel');

  #
  # Fetching by id
  #
  $cs = $csa->fetch_by_dbID(1);


=head1 DESCRIPTION

The Funcgen CoordSystemAdaptor works slighty different to the core version. As
the Funcgen DB stores features mapped to multiple core/dna DBs the schema and 
data versions(i.e. the last bit of the DB name) have to be stored. This maintains
a link between the seq_region_id stored in the Funcgen DB and the seq_region and assembly
tables stored in the core DB on which the features were originally built.

Default versions or ranking has not yet been tested.

This adaptor allows the querying of information from the coordinate system
adaptor.

Note that many coordinate systems do not have a concept of a version
for the entire coordinate system (though they may have a per-sequence version).
The 'chromosome' coordinate system usually has a version (i.e. the
assembly version) but the clonal coordinate system does not (despite having
individual sequence versions).  In the case where a coordinate system does
not have a version an empty string ('') is used instead.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Funcgen::CoordSystem;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my %cs_warnings;

=head2 new

  Arg [1]    : See BaseAdaptor for arguments (none specific to this
               subclass)
  Example    : $cs = $db->get_CoordSystemAdaptor(); #better than new()
  Description: Creates a new CoordSystem adaptor and caches the contents
               of the coord_system table in memory.
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor
  Exceptions : none
  Caller     :
  Status     : At risk

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #
  # Cache the entire contents of the coord_system table cross-referenced
  # by dbID and name
  #

  #Funcgen specific
  #Added extra key on schema_build for all
  #keyed on name, list of coord_system value
  $self->{'_name_cache'} = {};

  #keyed on id, coord_system value
  $self->{'_dbID_cache'} = {};

  #keyed on rank
  #$self->{'_rank_cache'} = {};

  #keyed on id, 1/undef values
  $self->{'_is_sequence_level'} = {};
  $self->{'_is_default_version'} = {};
  
  my $sql = 'SELECT coord_system_id, name, rank, version, attrib, schema_build, core_coord_system_id FROM coord_system';
  my @args;
  if($self->is_multispecies()) {
  	$sql.=' where species_id =?';
  	push(@args, $self->species_id());
  }
	$sql.=' order by coord_system_id';
  my $sth = $self->prepare($sql);
  $sth->execute(@args);

  my ($dbID, $name, $rank, $version, $attrib, $sbuild, $ccs_id, $cs);
  $sth->bind_columns(\$dbID, \$name, \$rank, \$version, \$attrib, \$sbuild, \$ccs_id);

  while($sth->fetch()) {
    my $seq_lvl = 0;
    my $default = 0;


	#what we need is an add schema_build, seq_level, default, rank method
	#name and version shoudl be same for one CS

	if($attrib) {
      foreach my $attrib (split(',', $attrib)) {
        $self->{"_is_$attrib"}->{$dbID} = 1;
        if($attrib eq 'sequence_level') {
          $seq_lvl = 1;
        } elsif($attrib eq 'default_version') {
          $default = 1;
        }
      }
    }


	#Found new name, version pair
	if(! $cs || ($dbID != $cs->dbID())){

	  if($cs){

		#handle caching here
		
		#the get methods which utilise these caches need to sort the results based on the latest schema build.
		#or maybe instead of just having one name, where cat the schema_build, but point to the same cs
		#so loop through all the schema build for one CS?


		$self->{'_dbID_cache'}->{$cs->dbID()} = $cs;



		#Right then
		#Unless we're querying by cs_id from the eFG DB then we will always need
		#schema_build&level||rank or name&version

		#No point in having NR rank cache, need to resolve with schema_build?
		#Name 
		#have schema_build as optional arg in all methods, get from BDAdaptor if not defined?
		#This will just match the schema build to the current eFG DB

		$self->{'_name_cache'}->{lc($cs->name())} ||= [];
		#$self->{'_rank_cache'}->{$rank} ||= [];
		#push @{$self->{'_rank_cache'}->{$rank}}, $cs;
		push @{$self->{'_name_cache'}->{lc($cs->name())}}, $cs;
	  }


	  $cs = Bio::EnsEMBL::Funcgen::CoordSystem->new
		(-DBID           => $dbID,
		 -ADAPTOR        => $self,
		 -NAME           => $name,
		 -VERSION        => $version,
		 #-IS_CURRENT     => $is_current, #Currently only used by mart direct SQL
	  );
	}


	#could we fetch the actual core CS here, and add it to the eFG coord sys?
	#or should we just handle the individual args?
	#do we need to write generic method in DBAdaptor for this, then we can use the 
	#CSAdaptor as a cache for all DBAdaptor(CSs) should we not use reg for this?

	#we could populate objects from new rather than from db, then create adaptor as required?
	#still need to store is stored in CD?  and also we need to test everytime to see if we have an adaptor

	$cs->add_core_coord_system_info(
									-RANK                 => $rank,
									-SEQUENCE_LEVEL       => $seq_lvl,
									-DEFAULT              => $default,
									-SCHEMA_BUILD         => $sbuild,
									-CORE_COORD_SYSTEM_ID => $ccs_id,
									-IS_STORED            => 1,
									);


	#can we change these caches to use just the name and version rather than schema_build?

	#$self->{'_sb_name_cache'}->{$sbuild.":".lc($name)} ||= [];
    #$self->{'_dbID_cache'}->{$dbID} = $cs;
    #$self->{'_sb_rank_cache'}->{$sbuild.":".$rank} = $cs;
    #push @{$self->{'_sb_name_cache'}->{$sbuild.":".lc($name)}}, $cs;
  }


  #handle last cs
 
  if($cs){
    $self->{'_dbID_cache'}->{$cs->dbID()} = $cs;
    #push @{$self->{'_rank_cache'}->{$rank}}, $cs;
    push @{$self->{'_name_cache'}->{lc($cs->name())}}, $cs;
  }

  $sth->finish();


  #Removed mapping path support as this is handled by core API
  #This also handled toplevel

  return $self;
}


=head2 fetch_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  Alternatively
               this may be an alias for a real coordinate system.  Valid
               aliases are 'toplevel' and 'seqlevel'.
  Arg [2]    : optional - string $version
               The version of the coordinate system to retrieve.  If not
               specified the default version for the appropriate schema_build
               will be used.
  Example    : $coord_sys = $csa->fetch_by_name('chromosome', 'NCBI36');
               # toplevel is an pseudo coord system representing the highest
               # coord system in a given region
               # such as the chromosome coordinate system
               $coord_sys = $csa->fetch_by_name('toplevel');
               #seqlevel is an alias for the sequence level coordinate system
               #such as the clone or contig coordinate system
               $coord_sys = $csa->fetch_by_name('seqlevel');
  Description: Retrieves a coordinate system by its name
  Returntype : Bio::EnsEMBL::Funcgen::CoordSystem
  Exceptions : throw if no name argument provided
               warning if no version provided and default does not exist
  Caller     : general
  Status     : At risk

=cut


#we need the schema_build for the top/sequence_level!!!!!!!!!!!!!!!!!!!!!!!!!
#if schema_build not defined them we need to use ->db->dnadb schema_build
#careful, this could be using the default dnadb already
#but this is the desired behaviour is it not?

#need a generic method to fetch the best cs based on dnadb or latest schema_build
#also need generic method in DBAdaptor to set dnadb by Experiment
#need to populate schema_build in Experiment?
#how can we do this dynamically?  all Experiment(ec, channel, rset) based methods should set dnadb appropriately?
#could this potentially mean this is called too many times for one query?
#or we could just let the user manage it?
#we need to check whether different/non-comparable schema_builds are added to the same result set
#use latest schema_build i.e. gene set or original schema_build.

sub fetch_by_name{
  my $self = shift;
  my $name = lc(shift);
  my $version = shift;
  my $sbuild = $self->db->_get_schema_build($self->db->dnadb);
  my ($cs, $found_cs);
  
  throw(q(Mandatory argument 'name') ) if ! defined $name;

  #Set default_version if not specified
  if(! defined $version){
    my $core_cs = $self->db->dnadb->get_CoordSystemAdaptor->fetch_by_name($name);

    if(! defined $core_cs){
      return;
    }

    $version = $core_cs->version;
  }

  $version = lc($version);


  #can we not just use
  #if(($name eq 'toplevel' || $name eq 'seqlevel') && ! $schema_build){
  #	throw('To access toplevel or seqlevel you must provide a the third schema_build argument');
  # }

  warn "Using dnadb(${sbuild}) to acquire $name" if $name =~ /level/;

  if($name eq 'seqlevel') {
    return $self->fetch_sequence_level_by_schema_build($sbuild);
  } elsif($name eq 'toplevel') {
    return $self->fetch_top_level_by_schema_build($sbuild);
  }

  if(! exists($self->{_name_cache}->{$name})) {
    if($name =~ /top/) {
      warn( q(Did you mean 'toplevel' coord system instead of '$name'?) );
    } elsif($name =~ /seq/) {
      warn( q(Did you mean 'seqlevel' coord system instead of '$name'?) );
    }
    return undef;
  }


  my @coord_systems = @{$self->{_name_cache}->{$name}};
 
  #Filter versions if or get the default for the schema_build or comparable

  #This will only get non-versioned CSs if there are already loaded on a given schema_build
  #Hence we can never retrieve a 'comparable' supercontig if it has not been loaded onto the current schema_build
  #Hence we end up loading a new CS for each non-versioned level.
  

  foreach $cs (@coord_systems) {
	#Versions are only relevant to assembled levels e.g. chromosome & scaffold?
	#No some DBs have version for contig, supercontig etc if they have no assembled level
	#Issues around unassembled levels with version i.e. not being able to access data from 
	#old version, even tho' 'assembly' of contig should be identical.  This should not be a big 
	#problem as the only species we are likely to do this with will have mappings between levels
	#or we can just re import for new version of unassembled level? This will give redundant data in 
	#the DB for those contigs which appear in both versions and are identical.

    if($version) {#Assembled level
      
      if( lc($cs->version()) eq lc($version) ){
        #This will pick the right CS even if the dnadb schema_build is not present
        $found_cs = $cs;
        last;
      }
    }
    else{#We have an unassembled/non-versioned level and can use any as there should only be 1
      $found_cs = $cs;
      last;
    }
  }
  
  if(! $found_cs){
  
    if($version) {
      warn "No coord system found for $sbuild version '$version'";
      return undef;
    }
    else{
      warn "Could not find $name CoordSystem.";
      return undef;
    }
  }

  return $found_cs;
}



=head2 fetch_all

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [2]    : optional string - attribute e.g. DEFAULT, SEQUENCE_LEVEL, RANK
  Example    : foreach my $cs (@{$csa->fetch_all()}) {
                 print $cs->name(), ' ', $cs->version(), "\n";
               }
  Description: Retrieves every coordinate system defined in the DB. Will 
               restrict to those from a particular dnadb if one is passed.
               #These will be returned in ascending order of rank. I.e.
               #The highest coordinate system with rank=1 would be first in the
               #array.
  Returntype : listref of Bio::EnsEMBL::Funcgen::CoordSystems
  Exceptions : none
  Caller     : general
  Status     : at risk - make arg optional so we can truly retrieve all

=cut

sub fetch_all {
  my ($self, $dnadb, $attribute) = @_;

  if(! $dnadb || ! (ref($dnadb) && $dnadb->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
	throw('Not yet implement full fetch_all, please pass a dnadb');
  }

  my @coord_systems;
  my $schema_build = $self->db->_get_schema_build($dnadb);

  foreach my $cs(values %{$self->{'_dbID_cache'}}){

	if ($cs->contains_schema_build($schema_build)){
	  next if($attribute && ! $cs->get_coord_system_attribute(uc($attribute), $dnadb));
	  push @coord_systems, $cs;
	}
  }

  ##order the array by rank in ascending order
  #foreach my $rank (sort {$a <=> $b} keys %{$self->{'_rank_cache'}}) {
  #  push @coord_systems, $self->{'_rank_cache'}->{$rank};
  #}

  return \@coord_systems;
}



=head2 fetch_by_rank

  Arg [1]    : int $rank
  Example    : my $cs = $coord_sys_adaptor->fetch_by_rank(1);
  Description: Retrieves a CoordinateSystem via its rank. 0 is a special
               rank reserved for the pseudo coordinate system 'toplevel'.
               undef is returned if no coordinate system of the specified rank
               exists.
  Returntype : Bio::EnsEMBL::Funcgen::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub fetch_by_rank {
  my $self = shift;
  my $rank = shift;

  thrw('not implemented rank cache yet');

  throw("Rank argument must be defined.") if(!defined($rank));
  throw("Rank argument must be a non-negative integer.") if($rank !~ /^\d+$/);

  if($rank == 0) {
    return $self->fetch_top_level();
  }

  return $self->{'_rank_cache'}->{$rank};
}



=head2 fetch_all_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  This can be
               the name of an actual coordinate system or an alias for a
               coordinate system.  Valid aliases are 'toplevel' and 'seqlevel'.
  Example    : foreach my $cs (@{$csa->fetch_all_by_name('chromosome')}){
                 print $cs->name(), ' ', $cs->version();
               }
  Description: Retrieves all coordinate systems of a particular name
  Returntype : listref of Bio::EnsEMBL::Funcgen::CoordSystem objects
  Exceptions : throw if no name argument provided
  Caller     : general
  Status     : Medium

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = lc(shift); #case insensitive matching

  throw('Name argument is required') if(!$name);

  if($name eq 'seqlevel') {
    return [$self->fetch_sequence_level()];
  } elsif($name eq 'toplevel') {
    return [$self->fetch_top_level()];
  }

  return $self->{'_name_cache'}->{$name} || [];
}


=head2 fetch_by_dbID

  Arg [1]    : int dbID
  Example    : $cs = $csa->fetch_by_dbID(4);
  Description: Retrieves a coord_system via its internal
               identifier, or undef if no coordinate system with the provided
               id exists.
  Returntype : Bio::EnsEMBL::Funcgen::CoordSystem or undef
  Exceptions : thrown if no coord_system exists for specified dbID
  Caller     : General
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;
  throw('dbID argument is required') if ! defined $dbID;
  return $self->{'_dbID_cache'}->{$dbID};
}



=head2 fetch_top_level

  Arg [1]    : none
  Example    : $cs = $csa->fetch_top_level();
  Description: Retrieves the toplevel pseudo coordinate system.
  Returntype : a Bio::EnsEMBL::Funcgen::CoordSystem object
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub fetch_top_level {
  my $self = shift;

  throw("Not yet implemented with schema_build");

  return $self->{'_top_level'};
}


=head2 fetch_sequence_level

  Arg [1]    : none
  Example    : ($id, $name, $version) = $csa->fetch_sequence_level();
  Description: Retrieves the coordinate system at which sequence
               is stored at.
  Returntype : Bio::EnsEMBL::Funcgen::CoordSystem
  Exceptions : throw if no sequence_level coord system exists at all
               throw if multiple sequence_level coord systems exists
  Caller     : general
  Status     : At risk

=cut

sub fetch_sequence_level {
  my $self = shift;

  throw("Not yet implemented with schema_build");

  my @dbIDs = keys %{$self->{'_is_sequence_level'}};

  throw('No sequence_level coord_system is defined') if(!@dbIDs);

  if(@dbIDs > 1) {
    throw('Multiple sequence_level coord_systems are defined.' .
          'Only one is currently supported');
  }

  return $self->{'_dbID_cache'}->{$dbIDs[0]};
}




=head2 get_mapping_path

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs1
  Arg [2]    : Bio::EnsEMBL::CoordSystem $cs2
  Example    : foreach my $cs @{$csa->get_mapping_path($cs1,$cs2);
  Description: Given two coordinate systems this will return a mapping path
               between them if one has been defined.  Allowed Mapping paths are
               explicitly defined in the meta table.  The following is an
               example:

          mysql> select * from meta where meta_key = 'assembly.mapping';
          +---------+------------------+--------------------------------------+
          | meta_id | meta_key         | meta_value                           |
          +---------+------------------+--------------------------------------+
          |      20 | assembly.mapping | chromosome:NCBI34|contig             |
          |      21 | assembly.mapping | clone|contig                         |
          |      22 | assembly.mapping | supercontig|contig                   |
          |      23 | assembly.mapping | chromosome:NCBI34|contig|clone       |
          |      24 | assembly.mapping | chromosome:NCBI34|contig|supercontig |
          |      25 | assembly.mapping | supercontig|contig|clone             |
          +---------+------------------+--------------------------------------+

               For a one-step mapping path to be valid there needs to be
               a relationship between the two coordinate systems defined in
               the assembly table.  Two step mapping paths work by building
               on the one-step mapping paths which are already defined.

               The first coordinate system in a one step mapping path must
               be the assembled coordinate system and the second must be
               the component.

               Example of use:
               my $cs1 = $cs_adaptor->fetch_by_name('contig');
               my $cs2 = $cs_adaptor->fetch_by_name('chromosome');

               my @path = @{$cs_adaptor->get_mapping_path($cs1,$cs2)};

               if(!@path) {
                 print "No mapping path.";
               }
               elsif(@path == 2) {
                 print "2 step mapping path.";
                 print "Assembled = " . $path[0]->name() . "\n";
                 print "Component = " . $path[1]->name() . "\n";
               } else {
                 print "Multi step mapping path\n";
               }

  Returntype : reference to a list of Bio::EnsEMBL::CoordSystem objects

  Exceptions : none
  Caller     : general
  Status     : At risk

=cut


#Need to be redirected to the core/dnadb of interest

sub get_mapping_path {
  my $self = shift;
  my $cs1 = shift;
  my $cs2 = shift;

  if(!ref($cs1) || !ref($cs2) ||
     !$cs1->isa('Bio::EnsEMBL::CoordSystem') ||
     !$cs2->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('Two Bio::EnsEMBL::CoordSystem arguments expected.');
  }

  my $key1 = $cs1->name() . ":" . $cs1->version();
  my $key2 = $cs2->name() . ":" . $cs2->version();

  my $path = $self->{'_mapping_paths'}->{"$key1|$key2"};

  return $path if($path);

  $path = $self->{'_mapping_paths'}->{"$key2|$key1"};

  if(!$path) {
    # No path was explicitly defined, but we might be able to guess a
    # suitable path.  We only guess for missing 2 step paths.

    my %mid1;
    my %mid2;

    foreach my $path (values(%{$self->{'_mapping_paths'}})) {
      next if(@$path != 2);

      my $match = undef;

      if($path->[0]->equals($cs1)) {
        $match = 1;
      } elsif($path->[1]->equals($cs1)) {
        $match = 0;
      }

      if(defined($match)) {
        my $mid = $path->[$match];
        my $midkey = $mid->name() . ':' . $mid->version();

        # is the same cs mapped to by other cs?
        if($mid2{$midkey}) {
          my $path = [$cs1,$mid,$cs2];
          $self->{'_mapping_paths'}->{"$key1|$key2"} = $path;
          $key1 =~ s/\:$//;
          $key2 =~ s/\:$//;
          $midkey =~ s/\:$//;
          warning("Using implicit mapping path between '$key1' and '$key2' " .
                  "coord systems.\n" .
                  "An explicit 'assembly.mapping' entry should be added " .
                  "to the meta table.\nExample: " .
                  "'$key1|$midkey|$key2'\n");
          return $path;
        } else {
          $mid1{$midkey} = $mid;
        }
      }

      $match = undef;

      if($path->[0]->equals($cs2)) {
        $match = 1;
      } elsif($path->[1]->equals($cs2)) {
        $match = 0;
      }


      if(defined($match)) {
        my $mid = $path->[$match];
        my $midkey = $mid->name() . ':' . $mid->version();

        # is the same cs mapped to by other cs?
        if($mid1{$midkey}) {
          my $path = [$cs2,$mid,$cs1];
          $self->{'_mapping_paths'}->{"$key2|$key1"} = $path;

          $key1 =~ s/\:$//;
          $key2 =~ s/\:$//;
          $midkey =~ s/\:$//;
          warning("Using implicit mapping path between '$key1' and '$key2' " .
                  "coord systems.\n" .
                  "An explicit 'assembly.mapping' entry should be added " .
                  "to the meta table.\nExample: " .
                  "'$key1|$midkey|$key2'\n");

          return $path;
        } else {
          $mid2{$midkey} = $mid;
        }
      }
    }
  }

  return $path || [];
}


=head2 _fetch_by_attribute

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 
  Status     : At risk

=cut

sub _fetch_by_attrib {
  my $self = shift;
  my $attrib = shift;
  my $version = shift;

  $version = lc($version) if($version);

  my @dbIDs = keys %{$self->{"_is_$attrib"}};

  throw("No $attrib coordinate system defined") if(!@dbIDs);

  foreach my $dbID (@dbIDs) {
    my $cs = $self->{'_dbID_cache'}->{$dbID};
    if($version) {
      return $cs if(lc($version) eq $cs->version());
    } elsif($self->{'_is_default_version'}->{$dbID}) {
      return $cs;
    }
  }

  #specifically requested attrib system was not found
  if($version) {
    throw("$attrib coord_system with version [$version] does not exist");
  }

  #coordsystem with attrib exists but no default is defined:
  my $dbID = shift @dbIDs;
  my $cs = $self->{'_dbID_cache'}->{$dbID};
  my $v = $cs->version();
  warning("No default version for $attrib coord_system exists. " .
          "Using version [$v] arbitrarily");

  return $cs;
}

=head2 _fetch_all_by_attribute

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 
  Status     : At risk

=cut

sub _fetch_all_by_attrib {
  my $self = shift;
  my $attrib = shift;

  my @coord_systems = ();
  foreach my $dbID (keys %{$self->{"_is_$attrib"}}) {
    push @coord_systems, $self->{"_dbID_cache"}->{$dbID};
  }

  return \@coord_systems;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Funcgen::CoordSystem
  Example    : $csa->store($coord_system);
  Description: Stores a CoordSystem object in the database.
  Returntype : none
  Exceptions : Warning if CoordSystem is already stored in this database.
  Caller     : none
  Status     : At risk

=cut

sub store {
  my $self = shift;
  my $cs = shift;

  if(!$cs || !ref($cs) || !$cs->isa('Bio::EnsEMBL::Funcgen::CoordSystem')) {
    throw('CoordSystem argument expected.');
  }

  my $sth;
  my $db = $self->db();
  my $name = $cs->name();
  my $version = $cs->version();

  if($name eq 'toplevel' || $name eq 'seqlevel' || !$name) {
    throw("[$name] is not a valid name for a storable CoordSystem.");
  } 
  
  foreach my $sbuild(keys %{$cs->{'core_cache'}}){
	my $rank    = $cs->{'core_cache'}->{$sbuild}->{'RANK'};
	my $seqlevel = $cs->{'core_cache'}->{$sbuild}->{'SEQUENCE_LEVEL'};
	my $default  = $cs->{'core_cache'}->{$sbuild}->{'DEFAULT'};
	my $ccs_id = $cs->{'core_cache'}->{$sbuild}->{'CORE_COORD_SYSTEM_ID'};

	#
	# Do lots of sanity checking to prevent bad data from being entered
	#
	
	if($cs->{'core_cache'}->{$sbuild}->{'IS_STORED'}) {
	  #Doesn't this only check on dbID?
	  next;
	}

 
	#if($seqlevel && keys(%{$self->{'_is_sequence_level'}})) {
	#  throw("There can only be one sequence level CoordSystem.");
	#}

	#if(exists $self->{'_name_cache'}->{lc($name)}) {
	#  my @coord_systems = @{$self->{'_name_cache'}->{lc($name)}};

	#  foreach my $c (@coord_systems) {
	#	if(lc($c->version()) eq lc($version)) {
	#	  warning("CoordSystem $name $version is already in db.\n");
	#	  return;
	#	}
	#	if($default && $self->{'_is_default_version'}->{$c->dbID()}) {
	#	  throw("There can only be one default version of CoordSystem $name");
	#	}
	#  }
	#}

	if($rank !~ /^\d+$/) {
	  throw("Rank attribute must be a positive integer not [$rank]");
	}

	if($rank == 0) {
	  throw("Only toplevel CoordSystem may have rank of 0.");
	}

	#if(defined($self->{'_rank_cache'}->{$rank})) {
	#  throw("CoordSystem with rank [$rank] already exists.");
	#}

	my @attrib;

	push @attrib, 'default_version' if($default);
	push @attrib, 'sequence_level' if($seqlevel);

	my $attrib_str = (@attrib) ? join(',', @attrib) : undef;

	#
	# store the coordinate system in the database
	#
	
	if(! $cs->dbID()){
		
	  $sth = $self->prepare('insert into coord_system (name, version, attrib, rank, schema_build, core_coord_system_id, species_id) values (?,?,?,?,?,?,?)');  

	  $sth->bind_param(1, $name,               SQL_VARCHAR);
	  $sth->bind_param(2, $version,            SQL_VARCHAR);
	  $sth->bind_param(3, $attrib_str,         SQL_VARCHAR);
	  $sth->bind_param(4, $rank,               SQL_INTEGER);
	  $sth->bind_param(5, $sbuild,             SQL_VARCHAR);
	  $sth->bind_param(6, $ccs_id,             SQL_INTEGER);
	  $sth->bind_param(7, $self->species_id(), SQL_INTEGER);


		#Here we are getting failures due to concurrent processes storing the same CS.
		#There is no abolsolute way of protecting again this unless we lock the tables
		#before we query.
		#This could happen with seq_region also, but it is higly unlikley that we will write at the same time.
		
	  $sth->execute();
	  
	  #eval { 	$sth->execute() };
	  
	  #if($@){
	  
	  #}
	  
	  
	  my $dbID = $sth->{'mysql_insertid'};
	  $sth->finish();
	  
	  if(!$dbID) {
		throw("Did not get dbID from store of CoordSystem.");
	  }
	  
	  $cs->dbID($dbID);
	  $cs->adaptor($self);
	}else{
	  #can we prep this out of the loop
	  #we don't know until we're in it
	  my $sql = 'insert into coord_system (coord_system_id, name, version, attrib, rank, schema_build, core_coord_system_id, species_id) values (?,?,?,?,?,?,?,?)';
	  $sth = $db->dbc->prepare($sql);

	  $sth->bind_param(1, $cs->dbID(),         SQL_INTEGER);
	  $sth->bind_param(2, $name,               SQL_VARCHAR);
	  $sth->bind_param(3, $version,            SQL_VARCHAR);
	  $sth->bind_param(4, $attrib_str,         SQL_VARCHAR);
	  $sth->bind_param(5, $rank,               SQL_INTEGER);
	  $sth->bind_param(6, $sbuild,             SQL_VARCHAR);
	  $sth->bind_param(7, $ccs_id,             SQL_INTEGER);
	  $sth->bind_param(8, $self->species_id(), SQL_INTEGER);
	  $sth->execute();
	  $sth->finish();
	}
	
	$cs->{'core_cache'}{$sbuild}{'IS_STORED'} = 1;


  }

 
  #
  # update the internal caches that are used for fetching
  #
  #$self->{'_is_default_version'}->{$dbID} = 1 if($default);
  #$self->{'_is_sequence_level'}->{$dbID} = 1 if($seqlevel);

  $self->{'_name_cache'}->{lc($name)} ||= [];
  #$self->{'_rank_cache'}->{$rank} ||= [];
  $self->{'_dbID_cache'}->{$cs->dbID()} = $cs;
  #this will duplicate CS in cache if we add a core cs and then store
  #same with rank cache, need to replace
  my $push = 1;

  foreach my $name_cs(@{$self->{'_name_cache'}->{lc($name)}}){
	
	if($name_cs->version() eq $cs->version()){
	  $push = 0;
	  $name_cs = $cs;
	}
  }

  push @{$self->{'_name_cache'}->{lc($name)}}, $cs if $push;

  #$push = 1;

  

  #this could result in mixed rank cs in the same rank cache
  #push @{$self->{'_rank_cache'}->{$rank}}, $cs;
  #need to rethink rank cache?  make it schema_rank cache
  

  return $cs;
}


=head2 validate_and_store_coord_system

  Arg [1]    : Bio::EnsEMBL::CoordSystem (could also be Funcgen::CoordSystem)
  Example    : my $funcgen_cs = $csa->validate_coord_system($core_cs);
  Description: Given a CoordSystem retrieves the corresponding Funcgen CoordSystem 
               or generates new one
  Returntype : Bio::EnsEMBL::Funcgen::CoordSystem
  Exceptions : throw if arg not valid and stored
  Caller     : general
  Status     : At risk - just have validate and let DBAdaptor store totally new CSs?

=cut

#currently get cs from slice, and need to validate for dnadb too
#can take FGCoordSystem or CoordSystem
 
sub validate_and_store_coord_system{
  my ($self, $cs) = @_;
	
  if (! (ref($cs) && $cs->isa('Bio::EnsEMBL::CoordSystem') && $cs->dbID())) {
    throw('Must provide a valid stored Bio::EnsEMBL::CoordSystem');
  }
  

  #Need to add to Funcgen coord_system here
  #check if name and version are present and reset coord_system_id to that one, else get last ID and create a new one
  #coord_system_ids will not match those in core DBs, so we need ot be mindful about this.
  #can't use is_stored as this simply checks the dbID
  #seq_region_ids may change between schemas with the same assembly version
  #Store schema_version in coord_system and create seq_region translation
  #table to maintain the seq_region_id mapping back to each core DB

  
  #Do we need to check the the dnadb and the slice db match?
  #Do we have to have specified a dnadb at this point?  No.
  #But need to put checks in place for dnadb methods i.e. seq/slice retrieval



  my $sbuild = $self->db->_get_schema_build($cs->adaptor->db());

  #this should implicitly use the current schema_build
  #hence providing specificty for non-version CS's e.g. supercontig etc...
  my $fg_cs = $self->fetch_by_name($cs->name(), $cs->version());

  #this needs to satify both schema_build and version
  #retrieving by name version should retunr the lastest schema_build unless the it is not the toplevel or highest expected rank?
  
  my $version;

  if (! $fg_cs) {
	
		#if($cs->name ne 'clone' && (! $cs->version)){
    #NO VERSION for assembled level !!
    #Assume the default version
    #we could get this from meta, but is unreliable
    #get from default chromosome version
    #my $tmp_cs = $cs->adaptor->fetch_by_name('chromosome');
    #$version = $tmp_cs->version;
		#}

	
		$fg_cs = Bio::EnsEMBL::Funcgen::CoordSystem->new(
                                                     -NAME    => $cs->name(),
                                                     -VERSION =>  $cs->version(),
                                                    );

		warn "Created new CoordSystem:\t".$fg_cs->name().":".$fg_cs->version()."\n";
  }


  #This is done in BaseFeatureAdaptor->_pre_store
  #to avoid users without write permission trying 
  #to store olf assemblies on new shema_builds 
  #on old schema_build which are already present
  #If the CS can't be found then you're probably 
  #importing for the first time and have write permissions

  #re-instated as we don't want any extra calls in _pre_store
  #as this iterates over every features stored
  #increasing import time.

  if (! $fg_cs->contains_schema_build($sbuild)) {
	
    #Need to set all attribs here.

    $fg_cs->add_core_coord_system_info(
                                       -RANK                 => $cs->rank(), 
                                       -SEQUENCE_LEVEL       => $cs->is_sequence_level(), 
                                       -DEFAULT              => $cs->is_default(), 
                                       -SCHEMA_BUILD         => $sbuild,
                                       -CORE_COORD_SYSTEM_ID => $cs->dbID(),
                                       -IS_STORED            => 0,
                                      );
	
    eval {  $fg_cs = $self->store($fg_cs) };
	
    if ($@ && (! exists $cs_warnings{$fg_cs->name.':'.$fg_cs->version})) {
      $cs_warnings{$fg_cs->name.':'.$fg_cs->version} = 1;
      warning("$@\nYou do not have permisson to store the CoordSystem for schema_build $sbuild\n".
              "Using comparable CoordSystem:\t".$fg_cs->name.':'.$fg_cs->version."\n");
    }
  }
  
  

  return $fg_cs;
}


1;
