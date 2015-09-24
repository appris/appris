#
# Ensembl module for Bio::EnsEMBL::Funcgen::Storable
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

Bio::EnsEMBL::Funcgen::Storable

=head1 SYNOPSIS

  my $dbID = $storable_object->dbID();
  my $adaptor = $storable_object->adaptor();
  if($storable_object->is_stored($db_adaptor))) {

  }
=head1 DESCRIPTION

This is a simple wrapper class to provide convenience methods for the StorableAdaptor.
Only get type methods have been implemented here to avoid obfuscating DB writes which
should only be done by the specific 'Storable'Adaptors.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor

=cut

package Bio::EnsEMBL::Funcgen::Storable;

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref );

use base qw( Bio::EnsEMBL::Storable );

=head2 new

  Arg [-STATES]  : Arrayref of states
  Arg [-dbID]    : database internal id
  Example        : none
  Caller         : internal calls
  Description    : create a new Storable object
  Returntype     : Bio::EnsEMBL::Storable
  Exceptions     : Adaptor not a Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor
  Status         : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($states, $assoc_ftypes) = rearrange(['STATES', 'ASSOCIATED_FEATURE_TYPES'] ,@_);

  if ($self->adaptor() && (! $self->adaptor->isa("Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor"))){
    throw("Adaptor muct be a valid Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor");
  }

  if($states){
    assert_ref($states, 'ARRAY');
    #Deref so we don't get unwanted update behaviour
    @{$self->{'states'}} = @$states if $states;
  }
  
  $self->associated_feature_types($assoc_ftypes) if(defined $assoc_ftypes);

  return $self;
}




=head2 has_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Example    : if($experimental_chip->has_status('IMPORTED'){ ... skip import ... };
  Description: Tests whether storable has a given status
  Returntype : BOOLEAN
  Exceptions : Throws if not status is provided
  Caller     : general
  Status     : At risk

=cut

sub has_status{
   my ($self, $status) = @_;

   throw("Must provide a status to check") if ! $status;

   my @state = grep(/^$status$/, @{$self->get_all_states()});
   my $boolean = scalar(@state);#will be 0 or 1 due to table contraints

   return $boolean;
}

#There is a potential to create an obj from scratch which may already exist in the db
#If we add a state to this (obj has not dbID so will not retrieve stored states)
# and then try and store it, this will result in adding the state to the previously stored obj.
#The behaviour is silent and could cause problems.

#To resolve this the adaptor implementations must throw if we find a matching object
#We must force the user to generate the obj from the db(use recover) rather than from scratch
#to make them aware of the situation.  This is useful to protect objects where we do not want to overwrite previous data
#e.g. experiment, experimental_chip, channel
#For objects which are routinely resued, we must make sure we always try the db first(not just when recover is set)
#Then warn/throw if there are differing attributes

#This is not possible for set objects, but is not a problem as it will just create another set entry rather than overwriting
#All update/store_states methods should be okay so long as we have a dbID first.


=head2 get_all_states

  Example    : my @ec_states = @{$experimental_chip->get_all_states()};
  Description: Retrieves of all states
  Returntype : Arrayref
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut

sub get_all_states{
   my ($self) = @_;

   #This could miss states in the DB for storables which have been created and had states added
   #but already exist with states in the DB
   #The way to get around this is to throw if we try and store an object without a dbID which matches
   #something in the DB.
   #Remove func in adaptors(ec and channel only?) to automatically use prestored objects, throw instead if no dbID and matches.
   #force use of recover to retrieve object from DB and then skip to relevant step based on states.
   #Have states => next method hash in Importer/ArrayDefs?

  if(! defined $self->{states}){

    if(! $self->adaptor){
      #throw('Cannot get_all_states for a Storable without an associated adaptor');
      $self->{states} = []; #Need to set this for compare_to
    }
    else{ #Populate cache

      if( $self->is_stored($self->adaptor->db) ){
        $self->{states} = $self->adaptor->fetch_all_states($self);
      }
    }
   }

   return $self->{states};
}


=head2 add_status

  Example    : $ec->add_state('DISPLAYABLE');
  Description: Adds a state to a new or previously stored Storable
  Returntype : None
  Exceptions : Throws if no status supplied
  Caller     : general
  Status     : At risk

=cut



sub add_status{
   my ($self, $status) = @_;

   throw("Must pass a status to add e.g. 'DISPLAYABLE'") if ! $status;



   #this does not resolve the problem!!???
   #can add a status to an unstored object which

   if($self->adaptor && $self->is_stored($self->adaptor->db()) && ! $self->{'states'}){
     @{$self->{'states'}} = @{$self->adaptor->fetch_all_states($self)};
   }

   push @{$self->{'states'}}, $status;

   return;
}

sub is_displayable{
  my $self = shift;

  return $self->has_status('DISPLAYABLE');
}

#These DBEntry methods are here to enable xrefs to FeatureType, Probe, ProbeSet & ProbeFeature
#They only work as SetFeature has Storable as the second element of @ISA and Bio::EnsEMBL::Feature
#get_all_DBEntries be incorporated into Bio::EnsEMBL::Storable as generic method
#With the other wrapper methods in the Storables of the non-core APIs?
#Or can these be moved to core as a supra core class?
#i.e. Is common to all non-core APIs but not relevant for core?
#Can we bring these together to stop code propogation?


=head2 get_all_Gene_DBEntries

  Example    : my @gene_dbentries = @{ $storable->get_all_Gene_DBEntries };
  Description: Retrieves Ensembl Gene DBEntries (xrefs) for this Storable.
               This does _not_ include the corresponding translations
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

#Need to add optional Transcript/Gene param so we can filter
#Filter here or would be better to restrict in sql query ni DBEntryAdaptor?

sub get_all_Gene_DBEntries {
  my $self = shift;


  #We wouldn't need this if we made the xref schema multi species
  #my $species = $self->adaptor->db->species;
  my $species = Bio::EnsEMBL::Registry->get_alias($self->adaptor->db->species);

  if(!$species){
	throw('You must specify a DBAdaptor -species to retrieve DBEntries based on the external_db.db_name');
  }

  #safety in case we get Homo sapiens
  ($species = lc($species)) =~ s/ /_/;



  return $self->get_all_DBEntries($species.'_core_Gene');
}

=head2 get_all_Transcript_DBEntries

  Arg[0]     : optional - Bio::EnsEMBL::Transcript to filter DBEntries on.
  Example    : my @transc_dbentries = @{ $set_feature->get_all_Transcript_DBEntries };
  Description: Retrieves ensembl Transcript DBEntries (xrefs) for this Storable.
               This does _not_ include the corresponding translations
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the Storable (i.e. they have not already been added or
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

#change the to ensembl_core when we implement Gene/Transcript/Protein|Translation links on the same external_db

sub get_all_Transcript_DBEntries {
  my ($self, $transcript) = @_;

  #We wouldn't need this if we made the xref schema multi species
  my $species = Bio::EnsEMBL::Registry->get_alias($self->adaptor->db->species);
  #Need to make sure this is latin name


  if(!$species){
	throw('You must specify a DBAdaptor -species to retrieve DBEntries based on the external_db.db_name');
  }

  #safety in case we get Homo sapiens
  #($species = lc($species)) =~ s/ /_/;

  my $dbes = $self->get_all_DBEntries($species.'_core_Transcript');

  #This needs to be moved to the DBEntryAdaptor and restrict the query using the
  #dbprimary_acc

  if($transcript){
	my $sid = $transcript->stable_id;

	#Test for sid here?

	if(ref($transcript) && $transcript->isa('Bio::EnsEMBL::Transcript')){
	  my @dbes;

	  foreach my $dbe(@$dbes){
		if($dbe->primary_id eq $sid){
		  push @dbes, $dbe;
		}
	  }
	  $dbes = \@dbes;
	}
	else{
	  throw('Transcript argument must be a valid Bio::EnsEMBL::Transcript');
	}
  }


  return $dbes;
}


=head2 get_all_UnmappedObjects

  Example    : my @uos = @{$storable->get_all_UnmappedObjects };
  Description: Retrieves UnamappedObjects for this Storable.
  Returntype : arrayref of Bio::EnsEMBL::UnmappedObject objects
  Exceptions : none
  Caller     : general
  Status     : At risk - move to Bio::Ensembl::Storable?

=cut

sub get_all_UnmappedObjects{
  my $self = shift;
  #Do we want to add external_db or analysis param here?

  my $class = ref($self);
  $class =~ s/.*:://;

  return $self->adaptor->db->get_UnmappedObjectAdaptor->fetch_all_by_object_type_id($class, $self->dbID);
}


=head2 get_all_DBEntries

  Arg[1]     : string - External DB name e.g. ensembl_core_Gene
  Arg[2]     : string - External DB type
  Example    : my @dbentries = @{ $set_feature->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this SetFeature.
               This does _not_ include the corresponding translations
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the SetFeature (i.e. they have not already been added or
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general, get_all_DBLinks
  Status     : Stable - at risk move to storable

=cut


#We could add 3rd arg here which would be xref(info_)type e.g. Gene/Transcript etc.
#Move info_type to ox.linkage_type to sit along side linkage_annotated


sub get_all_DBEntries {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my $cache_name = "dbentries";

  if(defined($ex_db_exp)){
    $cache_name .= $ex_db_exp;
  }
  if(defined($ex_db_type)){
    $cache_name .= $ex_db_type;
  }

  #Need to add tests for valid objects for xrefs

  # if not cached, retrieve all of the xrefs for this gene

  #This is not using the caching optimally
  #It seems for naive(ex_db_exp,ex_db_type) queries we create a naive cache
  #This means that further more specific queries will make another query and not use the cache


  if( (! defined $self->{$cache_name}) && $self->adaptor() ){

	my @tables = $self->adaptor->_tables;
	@tables = split/_/, $tables[0]->[0];#split annotated_feature
	my $object_type = join('', (map ucfirst($_), @tables));#change to AnnotatedFeature

    $self->{$cache_name} =
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $object_type, $ex_db_exp, $ex_db_type);
  }
  elsif( ! defined $self->{$cache_name} ){
	throw('You must have set and adaptor to be able to get_all_DBEntries');
  }


  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntry->new(...);
               $transcript->add_DBEntry($dbe);
  Description: Associates a DBEntry with this object. Note that adding
               DBEntries will prevent future lazy-loading of DBEntries for this
               storable (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

=cut

sub add_DBEntry {
  my $self = shift;
  my $dbe = shift;

  unless($dbe && ref($dbe) && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}


=head2 associated_feature_types

  Example    : my @associated_ftypes = @{$feature->associated_feature_types()};
  Description: Getter/Setter for other associated FeatureTypes.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen:FeatureType objects
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub associated_feature_types{
  my ($self, $ftypes) = @_;

  #Lazy load as we don't want to have to do a join on all features when most will not have any


  if(defined $ftypes){

	if(ref($ftypes) eq 'ARRAY'){

	  foreach my $ftype(@$ftypes){

		if( ! $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType') ){
		  throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::FeatureType objects');
		}
		#test is stored in adaptor
	  }

	  if(defined $self->{'associated_feature_types'}){
		warn('You are overwriting associated feature types');
		#we could simply add the new ones and make them NR.
	  }

	  $self->{'associated_feature_types'} = $ftypes;
	}
	else{
	  throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::FeatureType objects');
	}
  }


  if(! defined $self->{'associated_feature_types'}){
	#This will fail if we have not stored yet

	if(defined $self->adaptor){
	  $self->{'associated_feature_types'} = $self->adaptor->db->get_FeatureTypeAdaptor->fetch_all_by_association($self);
	}

  }


  #This has the potential to return undef, or an arrayref which may be empty.
  return $self->{'associated_feature_types'};
}


=head2 compare_scalar_methods

  Arg [1]    : Object to compare to, not necessarily stored.
  Arg [2]    : Arrayref - Method names which return a Scalar, or an Array/Arrayref
               of Scalars.
  Example    : my %diffs = %{$self->compare_to($other_obj,
                                               [qw(name, get_scalars_method)]);
  Description: Compares the returned string values of the specified methods
               between this object and the object passed.
  Returntype : Hashref of method name keys error string values which show
               the differences between self and the other object (in that order)
  Exceptions : Throws if arguments not defined and valid.
               Throws if this Storable or the object passed cannot call the
               specified methods.
               Throws if methods do not return scalars.
  Caller     : Storable::compare_to
  Status     : At risk

=cut

sub compare_scalar_methods {
  my ($self, $obj, $methods) = @_;

  if(! (defined $methods &&
        ref($methods) &&
        ref($methods) eq 'ARRAY')){
    throw('You must pass an Arrayref of methods to compare');
  }

  my $diffs = {};

  foreach my $method(@$methods){

    if ( ! $self->can($method) ) {
      throw(ref($self)." cannot call method '".$method."' to compare");
    }

    my ($these_scls, $other_scls, $diff, $multi) =
      @{$self->_compare_method_return_types($obj, $method)};

    if($diff){
      $diffs->{$method} = $diff;
      next;
    }

    foreach my $i(0..$#{$these_scls}){
      #Check all are scalar
      if (ref(\($these_scls->[$i])) ne 'SCALAR'){
        throw("$method does not return a SCALAR value or an ARRAY or ARRAYREF ".
              "of SCALAR values:\t".$these_scls->[$i]);
      }

      #equating strings with == results in all strings evaling as 0, therefore
      #would match any string (unless prefixed with numbers in which case it uses those).
      #ne/eq converts numbers into strings accurately, so this is safe.
      if($these_scls->[$i] ne $other_scls->[$i]){
        $diffs->{$method} = $these_scls->[$i].' - '.$other_scls->[$i];
      }
    }
  }

  return $diffs;
}

=head2 compare_storable_methods

  Arg [1]    : Storable to compare to. (Mandatory)
  Arg [2]    : Arrayref - Method names which return a Storable or an Array
               or Arrayref of Storables. (Mandatory)
  Example    : my %diffs = %{$self->compare_to($other_obj,
                                               [qw(get_storable, get_storables)]);
  Description: Compares the returned Storable(s) of the specified methods. This is
               done checking whether the are stored in the same DB as this Storable,
               and have matching dbIDs.
  Returntype : Hashref of method name keys and error string values which show
               the differences between self and the other Storable (in that order)
  Exceptions : Throws if arguments not defined and valid.
               Throws if this Storable or the object passed cannot call the
               specified methods.
               Throws if methods do not return stored Bio::EnsEMBL::Storables with an adptor
  Caller     : Storable::compare_to
  Status     : At risk

=cut

#removed deep mode from here as this was only valid for storing between DBs
#As it is not possible to equivalently sort multiple storables from
#different DBs based on their dbIDs, hence this comparison would fail.
#Would need to be done with unique key values, which are not generically accessibly here
#Deep inter-DB comparisons need to be done by iterative shallow comparisons
#through an object hierarchy in a wrapper method i.e. the wrapper method will know
#the appropriate unique keys to sort on.
#This also removes the need to track previous comparisons to prevent deep recursion
#Does this mean we can now add missing object methods which would have caused circular refs?
#Not need to do this in either inter/intraDB comparison

#Could remove require ment for them being storables if we can call compare_to on the
#returned objects. This would require some 'Role' definition as we can't guarantee what
#another object compare method might be name or in fact do.

#TODO Make test names (value string prefixes) accessable for validation?

sub compare_storable_methods {
  my ($self, $obj, $methods) = @_;
  my $diffs = {};

  if(! (defined $methods &&
        ref($methods) &&
        ref($methods) eq 'ARRAY')){
    throw('You must pass an Arrayref of methods to compare');
  }

  foreach my $method(@$methods){
    my $obj_diffs = {};

    my ($these_objs, $other_objs, $diff, $multi) =
      @{$self->_compare_method_return_types($obj, $method, 1)};
    #1 is 'storable' flag

    if(! (@$these_objs && @$other_objs)){
      next;
    }
    elsif($diff){
      $diffs->{$method} = $diff;
      next;
    }
    #else we have same number of storables
    
    #Need to account for a single undef value, here?
    

    for my $i (0..$#{$these_objs}){
      #We can't alter $method key here for test discrimination
      #else we won't be able to easily identify what method gave the error
      #TP ID specific test result, have to pattern match
      #the start of the error string

      #These always have to be stored, so can use DB from either object
      if(! $other_objs->[$i]->adaptor){
        throw('Could not access DBAdaptor from self('.ref($other_objs->[$i]).
          ") for $method is_stored check");
      }


      if(! $these_objs->[$i]->is_stored($other_objs->[$i]->adaptor->db)){
        $diffs->{$method} = 'DBs are not the same for '.ref($other_objs);
        last if $multi;
      }
      elsif(ref($these_objs->[$i]) ne ref($other_objs->[$i])){ #Different return types
        #As there is not requirement for $self and $obj to be the same class
        #there is no guarantee they will return the same object
        $diffs->{$method} = "Namespace mismatch:\n\t".
            ref($these_objs->[$i]).' - '.ref($other_objs->[$i]);
        last if $multi;
      }
      elsif($these_objs->[$i]->dbID != $other_objs->[$i]->dbID){
          $diffs->{$method} = "dbID mismatch:\t".
            join(', ', (map $_->dbID, @$these_objs))."\t-\t".
            join(', ', (map $_->dbID, @$other_objs));
          last if $multi;
      }
    }
  }

  return $diffs;
}


#We don't enforce that $self and $obj are the same class, just that they have the
#same method

sub _compare_method_return_types{
  my ($self, $obj, $method, $storable) = @_;

  if ( ! $self->can($method) ) {
      throw(ref($self).' cannot call method '.$method);
  }
    elsif ( ! $obj->can($method) ) {
      throw(ref($obj).' cannot call method '.$method);
  }

  my ($diff, $multi, @these_values, @other_values);

  #|| () avoids handing a array of a single undef value
  @these_values = $self->$method || ();
  @other_values = $obj->$method  || ();

  #Handle return types
  if( (scalar(@these_values) == 1) &&
      (scalar(@other_values) == 1) ){

    if( ref($these_values[0]) ){
      #This is allowed to be an Arrayref for scalar or Arrayref and Storable for storable

      #Storable test done in compare_storable_methods

      if(ref($these_values[0]) eq 'ARRAY'){
        @these_values = @{$these_values[0]};
        @other_values = @{$other_values[0]};
      }
    }
  }
  
  #Compare sizes and sort
  if(scalar(@these_values) != scalar(@other_values) ) {
    $diff = "Return size mismatch:\t".
              scalar(@these_values).', '.scalar(@other_values);
  }
  elsif(scalar(@these_values)){
    $multi = 1;

    if($storable){
        # Do storable check here as this is required for dbID sort. adds iteration
        # Move to compare_storable_methods?

        if(! eval { map { die() if ! $_->isa('Bio::EnsEMBL::Storable') } 
                     (@these_values, @other_values); 
                    1;}){
          throw($method.' method does not return Storable(s)');
        }

        @these_values = sort {$a->dbID <=> $b->dbID} @these_values;
        @other_values = sort {$a->dbID <=> $b->dbID} @other_values;
    }
    else{
      #scalar check done in compare_scalar_methods (for speed)
      @these_values = sort @these_values;
      @other_values = sort @other_values;
    }
  }

  return [\@these_values, \@other_values, $diff, $multi];
}


=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Mandatory list of Storable method names each
               returning a Scalar or an Array or Arrayref of Scalars.
  Args[4]    : Arrayref - Mandatory ist of Storable method names each
               returning a Storable or an Array or Arrayref of Storables.
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset,
                                                       1,
                                                       [qw(get_scalar get_scalars)],
                                                       [qw(get_storable get_storables)]
                                                       )};
  Description: Compare this Storable to another based on the defined scalar
               and storable methods.
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this Storable to another.
  Returntype : Hashref of key attribute/method name keys and values which differ.
               Keys will always be the method which has been compared.
               Values can either be a error string, a hashref of diffs from a
               nested object, or an arrayref of error strings or hashrefs where
               a particular method returns more than one object.
  Exceptions : Throws if object class does not match self.
               Throws if depth level invalid.
  Caller     : Storables with compare_to wrapper defining method arguments
  Status     : At Risk

=cut

#We could potentially allow compare_to between different classes
#But this would require some Role definition i.e. we can guarantee
#that methods with the same name between classes will return the same types

sub compare_to {
  my ($self, $obj, $shallow, $scl_methods, $obj_methods) = @_;

  if(! (defined $obj &&
        ref($obj)    &&
        $obj->isa(ref($self))) ){
      throw('You must pass a valid '.ref($obj).' to compare_to '.ref($self));
  }

  my $diffs = $self->compare_scalar_methods($obj, $scl_methods);

  if(! $shallow){
    %$diffs = (
        %$diffs,
        %{$self->compare_storable_methods (
          $obj,
          $obj_methods
          )
        }
      );
  }

  return $diffs;
}
1;
