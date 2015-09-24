#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::CoordSystem
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

Bio::EnsEMBL::Funcgen::CoordSystem

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);

  my $csa = $db->get_CoordSystemAdaptor();

  #
  # Get default chromosome coord system for the 39_36a DB:
  #
  my $cs = $csa->fetch_by_name_schema_build_version('chromosome', '39_36a');
  my $str = join ':', $cs->name(),$cs->version(),$cs->dbID();
  print "$str\n";


=head1 DESCRIPTION

This has been adapted from the core CoordSystem object to accomodate the multi-assembly
aspects of the eFG schema, namely hadnling the schema_build of the referenced core DB.

This is a simple object which contains a few coordinate system attributes:
name, internal identifier, version and schema_build.  A coordinate system is 
uniquely defined by its name and version and which DB it came from i.e. schema_build.  
A version of a coordinate system applies to all sequences within a coordinate system.  
This should not be confused with individual sequence versions.

Take for example the Human assembly.  The version 'NCBI33' applies to
to all chromosomes in the NCBI33 assembly (that is the entire 'chromosome'
coordinate system).  The 'clone' coordinate system in the same database would
have no version however.  Although the clone sequences have their own sequence
versions, there is no version which applies to the entire set of clones.

Coordinate system objects are immutable. Their name and version, and other
attributes may not be altered after they are created.

=cut

package Bio::EnsEMBL::Funcgen::CoordSystem;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Storable);

my %warnings;


=head2 new

  Arg [..]   : List of named arguments:
               -NAME      - The name of the coordinate system
               -VERSION   - (optional) The version of the coordinate system.
                            Note that if the version passed in is undefined,
                            it will be set to the empty string in the
                            resulting CoordSystem object.
               -RANK      - The rank of the coordinate system. The highest
                            level coordinate system should have rank 1, the
                            second highest rank 2 and so on.  An example of
                            a high level coordinate system is 'chromosome' an
                            example of a lower level coordinate system is
                            'clone'.
               -SCHEMA_BUILD - The schema and data build version of the DB of
                               origin.
               -TOP_LEVEL - (optional) Sets whether this is a top-level coord
                            system. Default = 0. This should only be set to
                            true if you are creating an artificial toplevel
                            coordsystem by the name of 'toplevel'
               -SEQUENCE_LEVEL - (optional) Sets whether this is a sequence
                                 level coordinate system. Default = 0
               -DEFAULT      - (optional)
                               Whether this is the default version of the 
                               coordinate systems of this name. Default = 0
               -DBID         - (optional) The internal identifier of this
                                coordinate system
               -ADAPTOR      - (optional) The adaptor which provides database
                               interaction for this object
  Example    : $cs = Bio::EnsEMBL::CoordSystem->new(-NAME    => 'chromosome',
                                                    -VERSION => 'NCBI33',
                                                    -RANK    => 1,
                                                    -DBID    => 1,
                                                    -SCHEMA_BUILD => '39_36a',
                                                    -ADAPTOR => adaptor,
                                                    -DEFAULT => 1,
                                                    -SEQUENCE_LEVEL => 0);
  Description: Creates a new CoordSystem object representing a coordinate
               system.
  Returntype : Bio::EnsEMBL::Funcgen::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);


  #Can we just hadnle schema_build here and call super->new for the rest.
  #We will also have to handle the top/default levels issues with multiple DBs


  #my ($name, $version, $sbuild, $top_level, $sequence_level, $default, $rank) =
  #  rearrange(['NAME','VERSION', 'SCHEMA_BUILD','TOP_LEVEL', 'SEQUENCE_LEVEL',
  #             'DEFAULT', 'RANK'], @_);

  my ($name, $version) =  rearrange(['NAME','VERSION'], @_);


  throw('A name argument is required') if(! $name);
  

  $version = '' if(!defined($version));
  

  #$top_level       = ($top_level)      ? 1 : 0;
  #$sequence_level  = ($sequence_level) ? 1 : 0;
  #$default         = ($default)        ? 1 : 0;
  #$rank ||= 0;

  #if($top_level) {
  #  if($rank) {
  #    throw('RANK argument must be 0 if TOP_LEVEL is 1');
  #  }

  #  if($name) {
  #    if($name ne 'toplevel') {
  #      throw('The NAME argument must be "toplevel" if TOP_LEVEL is 1')
  #    }
  #  } else {
  #    $name = 'toplevel';
  #  }

  #  if($sequence_level) {
  #    throw("SEQUENCE_LEVEL argument must be 0 if TOP_LEVEL is 1");
  #  }

  #  $default = 0;

#  } else {
#    if(!$rank) {
#      throw("RANK argument must be non-zero if not toplevel CoordSystem");
#    }
#    if($name eq 'toplevel') {
#      throw("Cannot name coord system 'toplevel' unless TOP_LEVEL is 1");
#    }
#  }

#  if($rank !~ /^\d+$/) {
#    throw('The RANK argument must be a positive integer');
#  }


  $self->{'core_cache'} = {};
  $self->{'version'} = $version;
  $self->{'name'} = $name;
  #$self->{'schema_build'} = $sbuild;
  #$self->{'top_level'} = $top_level;
  #$self->{'sequence_level'} = $sequence_level;
  #$self->{'default'} = $default;
  #$self->{'rank'}    = $rank;




  return $self;
}


=head2 add_core_coord_system_info

  Arg [1]    : mandatory hash:
                              
							 	-RANK                 => $rank,
								-SEQUENCE_LEVEL              => $seq_lvl,
								-DEFAULT          => $default,
								-SCHEMA_BUILD         => $sbuild,
								-CORE_COORD_SYSTEM_ID => $ccs_id,
								-IS_STORED            => $stored_status,
							  
  Example    : $cs->add_core_coord_system_info(
									-RANK                 => $rank,
									-SEQUENCE_LEVEL              => $seq_lvl,
									-DEFAULT          => $default,
									-SCHEMA_BUILD         => $sbuild,
									-CORE_COORD_SYSTEM_ID => $ccs_id,
									-IS_STORED            => 1,
									);

  Description: Setter for core coord system information
  Returntype : none
  Exceptions : throws if:
               rank not 0 when toplevel
               name not 'TOPLEVEL" when toplevel
               sequence level and top level
               no schema_build defined
               no rank
               rank 0 when not toplevel
               name 'TOPLEVEL' when not toplevel
               
  Caller     : Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor and ?
  Status     : at risk - replace with add_core_CoordSystem? implement top level?

#this does not check name and version!


=cut

sub add_core_coord_system_info {
  my ($self) = shift;

  my ($sbuild, $top_level, $sequence_level, $default, $rank, $stored, $ccs_id) =
    rearrange(['SCHEMA_BUILD','TOP_LEVEL', 'SEQUENCE_LEVEL',
               'DEFAULT', 'RANK', 'IS_STORED', 'CORE_COORD_SYSTEM_ID'], @_);


  throw('Must provide a schema_build') if ! $sbuild;
  throw('Must provide a core_coord_system_id') if ! $ccs_id;

  
 #$top_level       = ($top_level)      ? 1 : 0;
  $sequence_level  = ($sequence_level) ? 1 : 0;
  $default         = ($default)        ? 1 : 0;
  $stored ||=0;

  $rank ||= 0;

  if($top_level) {
    if($rank) {
      throw('RANK argument must be 0 if TOP_LEVEL is 1');
    }

    if($self->name()) {
      if($self->name() ne 'toplevel') {
        throw('The NAME argument must be "toplevel" if TOP_LEVEL is 1')
      }
    } else {
	  throw('toplevel not yet implemented');
      #$name = 'toplevel';
    }

    if($sequence_level) {
      throw("SEQUENCE_LEVEL argument must be 0 if TOP_LEVEL is 1");
    }

    $default = 0;

  } else {
    if(!$rank) {
      throw("RANK argument must be non-zero if not toplevel CoordSystem");
    }
    if($self->name() eq 'toplevel') {
      throw("Cannot name coord system 'toplevel' unless TOP_LEVEL is 1");
    }
  }

  if($rank !~ /^\d+$/) {
    throw('The RANK argument must be a positive integer');
  }


  #We can add unstored coord systems here
  #But will these ever have valid entries in the seq_region cache
  #Initialising this cache key turning off the warning in equals about
  #Using the nearest coord_system

  $self->{'core_cache'}{$sbuild} = {(
									 RANK                 => $rank,
									 SEQUENCE_LEVEL       => $sequence_level,
									 DEFAULT              => $default,
									 CORE_COORD_SYSTEM_ID => $ccs_id,
									 IS_STORED            => $stored,
									)};




  return;
}


#remove all but schema_buil and equals?
#depends on how we handle levels

=head2 name

  Arg [1]    : (optional) string $name
  Example    : print $coord_system->name();
  Description: Getter for the name of this coordinate system
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  return $self->{'name'};
}


=head2 get_latest_schema_build

  Example    : my $db_schema_build = $coord_system->get_latest_schema_build();
  Description: Getter for the most recent schema_build of this coordinate system
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub get_latest_schema_build {
  my $self = shift;

  return (sort (keys %{$self->{'core_cache'}}))[0];
}


=head2 contains_schema_build

  Example    : if ($coord_system->contains_schema_build('43_36e')){..do some coord system things ..};
  Description: Returns true is the CoordSystem maps to the corresponding core CoordSystem
  Returntype : Boolean
  Exceptions : throws if schema_build not defined
  Caller     : general
  Status     : at risk

=cut

sub contains_schema_build {
  my ($self, $schema_build) = @_;
  
  throw('Must pass a schema_build') if ! $schema_build;

  return (exists $self->{'core_cache'}{$schema_build}) ? 1 : 0;
}

=head2 version

  Arg [1]    : none
  Example    : print $coord->version();
  Description: Getter/Setter for the version of this coordinate system.  This
               will return an empty string if no version is defined for this
               coordinate system.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
  my $self = shift;

  return $self->{'version'};
}




=head2 equals

  Arg [1]    : Bio::EnsEMBL::Funcgen::CoordSystem $cs
               The coord system to compare to for equality.
  Example    : if($coord_sys->equals($other_coord_sys)) { ... }
  Description: Compares 2 coordinate systems and returns true if they are
               equivalent.  The definition of equivalent is sharing the same
               name and version.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub equals {
  my $self = shift;
  my $cs = shift;

  if(!$cs || !ref($cs) || 
	 (! $cs->isa('Bio::EnsEMBL::Funcgen::CoordSystem')  && 
	  ! $cs->isa('Bio::EnsEMBL::CoordSystem'))){
	throw('Argument must be a Bio::EnsEMBL::Funcgen::CoordSystem');
  }
  
  
  #need to add check on schema_build here
  #all schema_builds should have been added by BaseFeatureAdaptor during import
  #fails if we are using two different versions with the same cs's

  if((lc($self->version()) eq lc($cs->version())) &&
     (lc($self->name()) eq lc($cs->name()))){

	#we need to make sure these are default CS, otherwise we can get into trouble with
	#re-used or mismatched seq_region_ids between DB wih different default assemblies

	if (! $self->contains_schema_build($self->adaptor->db->_get_schema_build($cs->adaptor->db))) {

	  #Only warn first time this is seen
	  my $warning_key = $self->adaptor->db->_get_schema_build($cs->adaptor->db).':'.$self->name().':'.$self->version;

	  if(! exists $warnings{$warning_key}){
      warn 'You are using a schema_build('.$self->adaptor->db->_get_schema_build($cs->adaptor->db).
        ') which has no CoordSystem stored for '.$cs->version.". Defaulting to closest name version match.\n";
      $warnings{$warning_key} = 1;
	  }
	}
    return 1;
  }

  return 0;
}




=head2 is_top_level

  Arg [1]    : none
  Example    : if($coord_sys->is_top_level()) { ... }
  Description: Returns true if this is the toplevel pseudo coordinate system.
               The toplevel coordinate system is not a real coordinate system
               which is stored in the database, but it is a placeholder that
               can be used to request transformations or retrievals to/from
               the highest defined coordinate system in a given region.
  Returntype : 0 or 1
  Exceptions : none
  Caller     : general
  Status     : at risk - not implemented yet

=cut

sub is_top_level {
  my $self = shift;

  throw('Not yet implmented, need to test against the core cache using dnadb/schema_build');

  return $self->{'top_level'};
}


#These attribute methods are largely redundant
#is_default is used by Feature Adaptors to restrict features to
#current default assembly for non slice based methods
#Especially redundant now we have implemented this in fetch_all

=head2 is_sequence_level

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : if($coord_sys->is_sequence_level($dnadb)) { ... }
  Description: Returns true if this is a sequence level coordinate system
               for a given dnadb
  Returntype : 0 or 1
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub is_sequence_level {
  my ($self, $dnadb) = @_;

  return $self->get_coord_system_attribute('sequence_level', $dnadb);
}


=head2 is_default

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : if($coord_sys->is_default($dnadb)) { ... }
  Description: Returns true if this coordinate system is the default
               version of the coordinate system of this name for a given dnadb.
  Returntype : 0 or 1
  Exceptions : none
  Caller     : general - Used 
  Status     : at risk

=cut

sub is_default {
  my ($self, $dnadb) = @_;

  return $self->get_coord_system_attribute('default', $dnadb);
}

sub get_coord_system_attribute{
  my($self, $attr_name, $dnadb) = @_;
  
  if(! ($dnadb && ref($dnadb) && $dnadb->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
	throw("You must pass a dnadb to access the CoordSystem attribute:\t $attr_name");
  }
  
  my $schema_build = $self->adaptor->db->_get_schema_build($dnadb);
  
  if(! $self->contains_schema_build($schema_build)){
	throw("CoordSystem does not contain the schema_build:\t$schema_build");
  }

  return $self->{'core_cache'}{$schema_build}{uc($attr_name)};

}


=head2 rank

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : if($cs1->rank($dnadb) < $cs2->rank($dnadb)) {
                 print $cs1->name(), " is a higher level coord system than",
                       $cs2->name(), "\n";
               }
  Description: Returns the rank of this coordinate system for a given dnadb. 
               A lower number is a higher coordinate system.  The highest level coordinate
               system has a rank of 1 (e.g. 'chromosome').  The toplevel
               pseudo coordinate system has a rank of 0.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : at risk - not yet implemented

=cut

sub rank {
  my ($self, $dnadb) = @_;
  return $self->get_coord_system_attribute('rank', $dnadb);

}

1;
