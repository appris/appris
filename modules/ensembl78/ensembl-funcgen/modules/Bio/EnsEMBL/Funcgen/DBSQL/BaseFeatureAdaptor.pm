#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor - Funcgen Feature Adaptor base class

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for Funcgen feature adaptors. This base class is simply a way
to redefine some methods to use with the Funcgen DB.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

use vars qw(@EXPORT); #require Exporter is done in BaseAdaptor
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

### To do ###
#Add Translation method for fetch by FeatureSets methods?
#Test and implement feature mapping between coord systems
#Correct/Document methods!!!
#Implement externaldb db_name version registry

my %warnings;

=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Wrapper method for core BaseAdaptor, build seq_region cache for features
  Returntype : ARRAYREF of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : FeatureAdaptor classes
  Status     : at risk

=cut

sub generic_fetch {
  my $self = shift;

  #need to wrap _generic_fetch to always generate the
  #seq_region_cache other wise non-slice based fetch methods fail

  #build seq_region cache here once for entire query
  #Using default schema_build here
  #So would need to set dnadb appropriately
  #This is cleaning tmp cache values such that
  #nested feature fetches cause failure e.g. when regfeats retrieve their reg_attrs
  #We need a way of always generating the tmp cache, or having it persist?
  #This is because we haven't built the tmp cache in non_slice based methods i.e. we haven't run get_seq_region_id_by_Slice
  #This is the same for all non-Slice based methods
  #And the is no way around it as we are not providing that info about the new DB by passing a slice!
  #The only way to get around this is to make the tmp_cache persistant

  $self->build_seq_region_cache();

  return $self->SUPER::generic_fetch(@_);
}


=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : String (optional) - An SQL query constraint (i.e. part of the WHERE clause)
               OR
               Hashref - Params hashref containing constraints and/or optional contraints
               See BaseAdaptor::compose_constraint for more information.
  Arg [3]    : String (optional) - The logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which
               are on the Slice defined by $slice and fulfill the SQL
               constraint defined by $constraint. If logic name is defined,
               only features with an analysis of type $logic_name will be
               returned.
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
  my($self, $slice, $params, $logic_name) = @_;

  my @result;

  if (! (ref($slice) && $slice->isa('Bio::EnsEMBL::Slice'))) {
    throw('Bio::EnsEMBL::Slice argument expected.');
  }

  #Maintain $params, as this may contain more than the constraint config
  my $constraint = '';

  if((defined $params) &&
     (ref($params) eq 'HASH')){ #This is a scalar status arg
    #overwriting logic_names here

    if(defined $logic_name){
      $params->{constraints}{logic_names} = [$logic_name];
    }

    $constraint = $self->compose_constraint_query($params);

    #This does not yet return undef for unknown logic_names
    #So the query will still be performed
  }
  else{#still supporting passing a constraint string 
    $constraint = $self->_logic_name_to_constraint($params, $logic_name);
    #if the logic name was invalid, undef was returned
    return [] if ! defined $constraint; 
  }

  my $fg_cs = $self->db->get_FGCoordSystemAdaptor->fetch_by_name
   ($slice->coord_system->name,
    $slice->coord_system->version);

  if (! defined $fg_cs) {
    warn "No CoordSystem present for ".$slice->coord_system->name().":".$slice->coord_system->version();
    return \@result;
  }

  #build seq_region cache here once for entire query
  $self->build_seq_region_cache($slice);
  my $syn = $self->_main_table->[1];


  #check the cache and return if we have already done this query
  #Added schema_build to feature cache for EFG
  my $key;
  my $cache;
  if ( !$self->db->no_cache() ) {
    my $bind_params = $self->bind_param_generic_fetch();
    $key = uc(join(':', $slice->name, $constraint, $self->db->_get_schema_build($slice->adaptor->db())));

    if ( defined($bind_params) ) {
      $key .= ':'
        . join( ':', map { $_->[0] . '/' . $_->[1] } @{$bind_params} );
    }
    $cache = $self->_slice_feature_cache();
    if (exists($cache->{$key})) {
      $self->{'_bind_param_generic_fetch'} = ();
      return $cache->{$key};
    }
  }

  my $sa = $slice->adaptor();

  # Hap/PAR support: retrieve normalized 'non-symlinked' slices
  my @proj = @{$sa->fetch_normalized_slice_projection($slice)};

  if (@proj == 0) {
    throw('Could not retrieve normalized Slices. Database contains ' .
          'incorrect assembly_exception information.');
  }

  # Want to get features on the FULL original slice
  # as well as any symlinked slices

  # Filter out partial slices from projection that are on
  # same seq_region as original slice

  my $sr_id = $slice->get_seq_region_id();
  @proj = grep { $_->to_Slice->get_seq_region_id() != $sr_id } @proj;
  my $segment = bless([1,$slice->length(),$slice ],
                      'Bio::EnsEMBL::ProjectionSegment');
  push( @proj, $segment );

  # construct list of Hap/PAR boundaries for entire seq region
  my @bounds;
  my $ent_slice = $sa->fetch_by_seq_region_id($sr_id);
  $ent_slice    = $ent_slice->invert() if($slice->strand == -1);
  my @ent_proj  = @{$sa->fetch_normalized_slice_projection($ent_slice)};

  shift @ent_proj;              # skip first
  @bounds = map {$_->from_start - $slice->start() + 1} @ent_proj;


  # fetch features for the primary slice AND all symlinked slices
  foreach my $seg (@proj) {

    my $offset = $seg->from_start();
    my $seg_slice  = $seg->to_Slice();
    my $features = $self->_slice_fetch($seg_slice, $constraint); ## NO RESULTS? This is a problem with the cs->equals method?

    # if this was a symlinked slice offset the feature coordinates as needed
    if ($seg_slice->name() ne $slice->name()) {

    FEATURE:
      foreach my $f (@$features) {

        if ($offset != 1) {
          $f->{start} += $offset-1;
          $f->{end}   += $offset-1;
        }

        # discard boundary crossing features from symlinked regions
        foreach my $bound (@bounds) {
          if ($f->{'start'} < $bound && $f->{'end'} >= $bound) {

            next FEATURE;
          }
        }

        $f->{'slice'} = $slice;
        push @result, $f;
      }
    } else {
      push @result, @$features;
    }
  }

  if (defined $key) {
    $cache->{$key} = \@result;
  }

  $self->reset_true_tables; #for compose_constraint support

  return \@result;
}

=head2 build_seq_region_cache

  Arg [1]    : optional - Bio::EnsEMBL::Slice
               the slice from which to obtain features
  Example    : $self->build_seq_region_cache();
  Description: Builds the seq_region_id caches based on the species and schema_build
  Returntype : None
  Exceptions : thrown if optional Slice argument is not valid
  Caller     : self
  Status     : At risk - should be private _build_seq_region_cache? Change arg to DBAdaptor? Or remove if we are building the full cache?

=cut

#do we even need to have the coord system?
#so long as we are only using one schema build i.e. one dnadb defualt = current
#slice and fg_cs are optional
#need to look at this

sub build_seq_region_cache{
  my ($self, $slice) = @_;

  #instead of building this for each query
  #we need to store the current cached coord_system
  #and test whether it contains the query coord_system
  #No as this then restricts us to one coord_system
  #where as we want to cover all from one DB/schema_build

  if (defined $slice) {
    throw('Optional argument must be a Bio::EnsEMBL::Slice') if(! ( ref($slice) && $slice->isa('Bio::EnsEMBL::Slice')));
  }

  my $dnadb = (defined $slice) ? $slice->adaptor->db() : $self->db->dnadb();
  my $schema_build = $self->db->_get_schema_build($dnadb);
  my $sql = 'select sr.core_seq_region_id, sr.seq_region_id from seq_region sr';
  my @args = ($schema_build);

  if ($self->is_multispecies()) {
  	$sql.= ', coord_system cs where sr.coord_system_id = cs.coord_system_id and cs.species_id=? and';
  	unshift(@args, $self->species_id());
  } else {
  	$sql.= ' where';
  }
  $sql.=' sr.schema_build =?';

  #Check we have not already got the right cache
  my $cache_key = join(':', @args);

  #Do we already have a valid cache?
  return if($self->cache_key eq $cache_key);


  #Can't maintain these caches as we may be adding to them when storing
  #we should just maintain the cache with the cache_key

  #$self->{'seq_region_cache'}{$cache_key} ||= {};
  #$self->{'core_seq_region_cache'} = {};

  my $sth = $self->prepare($sql);
  $sth->execute(@args);
  while (my $ref = $sth->fetchrow_arrayref()) {
    $self->{seq_region_cache}->{$cache_key}->{$ref->[0]} = $ref->[1];
    $self->{core_seq_region_cache}->{$cache_key}->{$ref->[1]} = $ref->[0];
  }
  $sth->finish();

  $self->cache_key($cache_key);

  return;
}


=head2 cache_key

  Arg [1]    : optional string - species_id.schema_build e.g. 1.55_37
               the slice from which to obtain features
  Example    : $self->build_seq_region_cache();
  Description: Getter/Setter for the seq_region cache_key
  Returntype : string
  Exceptions : None
  Caller     : self
  Status     : At risk

=cut

sub cache_key{
  my ($self, $key) = @_;

  $self->{'cache_key'} = $key if $key;
  return $self->{'cache_key'} || '';
}



#Need to separate these into by slice and by cs?
#We could just use an old cs/schema_build to grab the correct slice based on the name
#however we want to at least warn that the db needs updating
#The problem is that we want to detect whether the seq_region_id is present for _pre_store
#but automatically use a comparable seq_region for the normal fetch_methods.

#So we need to split or add $?

sub get_seq_region_id_by_Slice{
  my ($self, $slice, $fg_cs, $test_present) = @_;

  if (! ($slice && ref($slice) && $slice->isa("Bio::EnsEMBL::Slice"))) {
    throw('You must provide a valid Bio::EnsEMBL::Slice');
  }

  #We really need to validate the schema_build of the slice to make sure it
  #present in the current default coord_system i.e. the one which was used to
  #generate the seq_region_cache
  #This is set with the dnadb or with a slice query
  #This may not always have been done.

  #Now all we have to do is test the cache_key
  #Or can we just build_seq_region_cache as this checks and rebuilds if not correct
  #This may generate a mistmach between the dnadb and the schema_build used to generate the cache
  #This will be reset if required for subesquent queries using the cache key
  $self->build_seq_region_cache($slice);

  my ($core_sr_id,  $fg_sr_id);


  #Slice should always have an adaptor, no?
  if ( $slice->adaptor() ) {
    $core_sr_id = $slice->adaptor()->get_seq_region_id($slice);
  } else {
    $core_sr_id = $self->db->dnadb->get_SliceAdaptor()->get_seq_region_id($slice);
  }


  #This does not work!! When updating for a new schema_build we get the first
  #seq_region stored, than for each subsequent one, it arbitrarily assigns a value from the hash even tho the
  #the exists condition isn't met!
  #my $fg_sr_id = $self->{'seq_region_cache'}{$core_sr_id} if exists $self->{'seq_region_cache'}{$core_sr_id};
  #Can't replicate this using a normal hash

  #This cache has been built based on the schema_build
  if (exists $self->{'seq_region_cache'}{$self->cache_key}{$core_sr_id}) {
    $fg_sr_id = $self->{'seq_region_cache'}{$self->cache_key}{$core_sr_id};
  }


  if (! $fg_sr_id && ref($fg_cs)) {
    #This is used to store new seq_region info along side previous stored seq_regions of the same version

    if ( ! $fg_cs->isa('Bio::EnsEMBL::Funcgen::CoordSystem')) {
      throw('Must pass as valid Bio::EnsEMBL::Funcgen:CoordSystem to retrieve seq_region_ids for forwards compatibility, passed '.$fg_cs);
    }

    my $sql = 'select seq_region_id from seq_region where coord_system_id =? and name =?';
    my $sth = $self->prepare($sql);
    $sth->execute($fg_cs->dbID(), $slice->seq_region_name());

    #This may not exist, so we need to catch it here?
    ($fg_sr_id) = $sth->fetchrow_array();
    $sth->finish();

    #if we are providing forward compatability
    #Then we know the eFG DB doesn't have the core seq_region_ids in the DB
    #Hence retrieving the slice will fail in _obj_from_sth
    #So we need to set it internally here
    #Then pick it up when get_core_seq_region_id is called for the first time(from where?)
    #and populate the cache with the value

    #This only works if there is a comparable slice
    #If we are dealing with a new assembly, then no $fg_sr_id will be returned
    #So need to catch this in the caller

    #Can we remove this now the cache is persistant?
    #i.e. Cache is not regenerated everytime, hence we don't lose the data?

    if ($fg_sr_id) {
      $self->{'_tmp_core_seq_region_cache'}{$self->cache_key} = {(
                                                                  $fg_sr_id => $core_sr_id
                                                                 )};
    }
  } elsif (! $fg_sr_id && ! $test_present) {
		#This generally happens when using a new core db with a efg db that hasn't been updated
		#Default to name match or throw if not present in DB
		my $schema_build = $self->db->_get_schema_build($slice->adaptor->db());
		my $core_cs = $slice->coord_system;

		#Avoids mapping of core to efg seq_region_ids
		#via schema_build(of the new core db) as we are matching directly to the seq_name

		my $sql = 'select distinct(seq_region_id) from seq_region sr, coord_system cs where sr.coord_system_id=cs.coord_system_id and sr.name=? and cs.name =?';
		my @args = ($slice->seq_region_name(), $core_cs->name());

		if ($core_cs->version()) {
		  $sql.= ' and cs.version =?';
		  push(@args, $core_cs->version());
		}
		if ($self->is_multispecies()) {
			$sql.=' and cs.species_id=?';
			push(@args, $self->species_id());
		}

		my $sth = $self->prepare($sql);
		$sth->execute(@args);
		($fg_sr_id) = $sth->fetchrow_array();
		$sth->finish();

		if (! $fg_sr_id) {
		  #Warn instead of throw so we can catch absent seq_region without eval
		  #warn('Cannot find previously stored seq_region for: '.$core_cs->name.':'.$core_cs->version.':'.$slice->seq_region_name.
      #		   "\nYou need to update your eFG seq_regions to match your core DB using: update_DB_for_release.pl\n");
		}

		#Only warn first time this is seen
		my $warning_key = $core_cs->name.':'.$core_cs->version.':'.$slice->seq_region_name;

		if (! exists $warnings{$warning_key}) {
		  warn 'Defaulting to previously store seq_region for: '.$warning_key.
			  "\nYou need to update your eFG seq_regions to match your core DB using: update_DB_for_release.pl\n";
		  $warnings{$warning_key} = 1;
	  }
  }

  return $fg_sr_id;
}

sub get_core_seq_region_id{
  my ($self, $fg_sr_id) = @_;

  #This is based on what the current schema_build is
  #to enable multi-schema/assmelby look up
  #we will have to nest these in schema_build caches
  #and use the schema_build of the slice which is passed to acquire the core_seq_region_id
  #likewise for reverse, i.e. when we store.

  my $core_sr_id = $self->{'core_seq_region_cache'}{$self->cache_key}{$fg_sr_id};

  if (! defined $core_sr_id && exists $self->{'_tmp_core_seq_region_cache'}{$self->cache_key}{$fg_sr_id}) {
    #Do we need to test the cache_key here, might it have changed since get_seq_region_id_by_Slice?
    #Or will build_seq_region_cache handle this?

    #Why do we reset this in the main cache here if we are returning the value
    #Is this not corrupting the main cache? But we always want this value?
    $self->{'core_seq_region_cache'}{$self->cache_key}{$fg_sr_id} = $self->{'_tmp_core_seq_region_cache'}{$self->cache_key}{$fg_sr_id};

    #Delete here so we don't have schema_build specific sr_ids hanging around for another query?
    #delete  $self->{'_tmp_core_seq_region_cache'}{$self->cache_key}{$fg_sr_id};

    #These are valid values and would most likely be required for subsequent queries.
    #Cache key is now also schema_build specific so this is no longer a problem
    #Removed this as this causes RegFeat retrieval to fail
    #As all non slice based methods would fail due to the lack of info about unstored DB i.e. we need a slice


    $core_sr_id = $self->{'core_seq_region_cache'}{$self->cache_key}{$fg_sr_id};

  } elsif (! defined $core_sr_id) { #No cache entry, so get from dnadb

    my $sql = "select distinct(name) from seq_region where seq_region_id=$fg_sr_id";
    my $slice_name = $self->db->dbc->db_handle->selectall_arrayref($sql);
    ($slice_name) = @{$slice_name->[0]};

    #if(scalar(@names) != 1){#This should never happen

    #We should really grab the coord sys name above too.
    my $slice_adaptor = $self->db->dnadb->get_SliceAdaptor;
    $core_sr_id = $slice_adaptor->get_seq_region_id($slice_adaptor->fetch_by_region(undef, $slice_name));

    #Set this in the cache so we don't have to retrieve it again
    $self->{'core_seq_region_cache'}{$self->cache_key}{$fg_sr_id} = $core_sr_id;
  }

  return $core_sr_id;
}

=head2 _pre_store

  Arg [1]    : Bio::EnsEMBL::Feature
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Helper function containing some common feature storing functionality
               Given a Feature this will return a copy (or the same feature if no changes
	           to the feature are needed) of the feature which is relative to the start
               of the seq_region it is on. The seq_region_id of the seq_region it is on
               is also returned.  This method will also ensure that the database knows which coordinate
               systems that this feature is stored in.  This supercedes the core method, to trust the
               slice the feature has been generated on i.e. from the dnadb.  Also handles multi-coordsys
               aspect, generating new coord_system_ids as appropriate and assembly projection.
  Returntype : Bio::EnsEMBL::Feature and the seq_region_id it is mapped to
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::"Type"FeatureAdaptors
  Status     : At risk

=cut


sub _pre_store {
  my ($self, $feature, $new_assembly) = @_;
  #May want to add cs_level arg?
  #What about ignore length flag?

  if (!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand());

  my $db = $self->db();
  my $slice = $feature->slice();

  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region
  if ($slice->start != 1 || $slice->strand != 1) {

    #throw("You must generate your feature on a slice starting at 1 with strand 1");
    #We did remove this transfer it uses direct hash access which
    #did not work with old array based ResultFeatures

    #move feature onto a slice of the entire seq_region
    $slice = $slice->adaptor->fetch_by_region($slice->coord_system->name(),
                                              $slice->seq_region_name(),
                                              undef, #start
                                              undef, #end
                                              undef, #strand
                                              $slice->coord_system->version());
    $feature = $feature->transfer($slice);

    if (!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }



  #Project here before we start building sr caches and storing CSs
  if ($new_assembly) {
     #warn "Projecting ".$feature->start.'-'.$feature->end." on "..$feature->slice->name." to $new_assembly";


    my @segments = @{$feature->feature_Slice->project($slice->coord_system->name, $new_assembly)};
    # do some sanity checks on the projection results:
    # discard the projected feature if
    #   1. it doesn't project at all (no segments returned)
    #   2. the projection is fragmented (more than one segment)
    #   3. the projection doesn't have the same length as the original
    #      feature

    # this tests for (1) and (2)
    if (scalar(@segments) == 0) {
      warn "Feature doesn't project to $new_assembly\n";
      return;
    } elsif (scalar(@segments) > 1) {
      warn "Feature projection is fragmented in $new_assembly\n";
      return;
    }

    # test (3)
    my $proj_slice = $segments[0]->to_Slice;

    if ($feature->length != $proj_slice->length) {

      #if(! $conf->param('ignore_length')){
      warn "Feature projection is wrong length in $new_assembly\n";
      return;
      #  }
    }

    #warn "proj ".$proj_slice->name;

    # everything looks fine, so adjust the coords of your feature
    #Have to generate new_slice here as we are not sure it is going to be
    #on the same slice as the old assembly
    $slice = $proj_slice->adaptor->fetch_by_region($proj_slice->coord_system->name, $proj_slice->seq_region_name);

    #These are just callers for ResultFeature!
    #For speed.

    $feature->start($proj_slice->start);
    $feature->end($proj_slice->end);
    $feature->slice($slice);

    #warn  "new feature ".$feature->feature_Slice->name;
  }


  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system; #from core/dnadb

  #retrieve corresponding Funcgen coord_system and set id in feature
  my $csa = $self->db->get_FGCoordSystemAdaptor(); #had to call it FG as we were getting the core adaptor
  my $fg_cs = $csa->validate_and_store_coord_system($cs);
  $fg_cs = $csa->fetch_by_name($cs->name(), $cs->version()); #Why are we refetching this?
  my $tabname = $self->_main_table->[0];


  #Need to do this for Funcgen DB
  my $mcc = $db->get_MetaCoordContainer();

  $mcc->add_feature_type($fg_cs, $tabname, $feature->length);


  #build seq_region cache here once for entire query
  $self->build_seq_region_cache($slice);

  #Now need to check whether seq_region is already stored
  #1 is test present flag
  my $seq_region_id = $self->get_seq_region_id_by_Slice($slice, undef, 1);


  if (! $seq_region_id) {
    #check whether we have an equivalent seq_region_id
    $seq_region_id = $self->get_seq_region_id_by_Slice($slice, $fg_cs);
    my $schema_build = $self->db->_get_schema_build($slice->adaptor->db());
    my $sql;
    my $core_sr_id = $slice->get_seq_region_id;
    my @args = ($slice->seq_region_name(), $fg_cs->dbID(), $core_sr_id, $schema_build);

    #Add to comparable seq_region
    if ($seq_region_id) {
      $sql = 'insert into seq_region(seq_region_id, name, coord_system_id, core_seq_region_id, schema_build) values (?,?,?,?,?)';
      unshift(@args, $seq_region_id);
    }
    #No compararble seq_region
    else {
      $sql = 'insert into seq_region(name, coord_system_id, core_seq_region_id, schema_build) values (?,?,?,?)';
    }

    my $sth = $self->prepare($sql);

    #Need to eval this
    eval{$sth->execute(@args);};

    if (!$@) {
      $seq_region_id =  $self->last_insert_id;
    }

    #Now we need to add this to the seq_region caches
    #As we are not regenerating them every time we query.
    $self->{seq_region_cache}{$self->cache_key}{$core_sr_id} = $seq_region_id;
    $self->{core_seq_region_cache}{$self->cache_key}{$seq_region_id} = $core_sr_id;
  }

  #Need to return seq_region_id as they  are not stored
  #in the slice retrieved from slice adaptor
  return ($feature, $seq_region_id);
}


#This is causing problems with remapping
#Need to merge this back into the core code, subbing out the relevant parts
#required for funcgen support

sub _slice_fetch {
  my $self = shift;
  my $slice = shift;
  my $orig_constraint = shift;


  my $slice_start  = $slice->start();
  my $slice_end    = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs     = $slice->coord_system();
  #my $slice_seq_region = $slice->seq_region_name();

  #We Need to translate the seq_regions IDs to efg seq_region_ids
  #we need to fetch the seq_region ID based on the coord_system id and the name
  #we don't want to poulate with the eFG seq_region_id, jsut the core one, as we need to maintain core info in the slice.

  #we need to cache the seq_region_id mappings for the coord_system and schema_build,
  #so we don't do it for every feature in a slice
  #should we just reset it in temporarily in fetch above, other wise we lose cache between projections
  #should we add seq_region cache to coord_system/adaptor?
  #Then we can access it from here and never change the slice seq_region_id
  #how are we accessing the eFG CS?  As we only have the core CS from the slice here?
  #No point in having in eFG CS if we're just having to recreate it for everyquery
  #CSA already has all CS's cached, so query overhead would be mininal
  #But cache overhead for all CS seq_region_id caches could be quite large.
  #we're going to need it in the _pre_store, so having it in the CS would be natural there
  #can we just cache dynamically instead?
  #we're still calling for csa rather than accessing from cache
  #could we cache the csa in the base feature adaptor?




  #my $slice_seq_region_id = $slice->get_seq_region_id();

  #get the synonym and name of the primary_table
  my ($tab_name, $tab_syn) = @{$self->_main_table};

  #find out what coordinate systems the features are in
  my $mcc = $self->db->get_MetaCoordContainer();
  my @feat_css=();


  my $mca = $self->db->get_MetaContainer();
  my $value_list = $mca->list_value_by_key( $tab_name."build.level" );
  if ( @$value_list and $slice->is_toplevel()) {
    push @feat_css, $slice_cs;
  } else {
    @feat_css = @{$mcc->fetch_all_CoordSystems_by_feature_type($tab_name)};
  }

  my $asma = $self->db->get_AssemblyMapperAdaptor();
  my @features;


  # fetch the features from each coordinate system they are stored in
 COORD_SYSTEM: foreach my $feat_cs (@feat_css) {

    my $mapper;
    my @coords;
    my @ids;

    #This is not returning true if the feat_cs is no default
    #



    if ($feat_cs->equals($slice_cs)) {
      # no mapping is required if this is the same coord system

      #eFG change
      #We want to set this to undef if we are dealing with result_features which
      #are not wsize==0!
      #This is only required if we load array data at the natural resolution!
      #This test will effect ever single feature query, can we omit this somehow?
      my $max_len;


      if (! (($self->can('_window_size') &&
              $self->_window_size))) {
        $max_len = $self->_max_feature_length() ||
          $mcc->fetch_max_length_by_CoordSystem_feature_type($feat_cs,$tab_name);
      }

      my $constraint = $orig_constraint;
      my $sr_id = $self->get_seq_region_id_by_Slice($slice, $feat_cs);


      #If we have no slice id here, we know we don't have anything in the DB for this slice
      return [] if ! $sr_id;

      #Only add the slice clause if we are not forcing reslice
      #force_reslice means we have dbID in the query for reg attrs
      #which maybe of the end of the query slice, but we still need them for
      #building reg features which overlap the end of the query slice

      if (! $self->force_reslice) {


        $constraint .= " AND " if($constraint);
        $constraint .=
          "${tab_syn}.seq_region_id = $sr_id AND " .
            "${tab_syn}.seq_region_start <= $slice_end AND " .
              "${tab_syn}.seq_region_end >= $slice_start";


        if ($max_len) {

          my $min_start = $slice_start - $max_len;
          $constraint .=
            " AND ${tab_syn}.seq_region_start >= $min_start";
        }
      }

      my $fs = $self->generic_fetch($constraint,undef,$slice);
      # features may still have to have coordinates made relative to slice
      # start




      $fs = $self->_remap($fs, $mapper, $slice);
      push @features, @$fs;
    } else {

      #Table contains some feature on a CS that differs from the Slice CS
      #can't do CS remapping yet as AssemblyMapper expects a core CS
      #change AssemblyMapper?
      #or do we just create a core CS just for the remap and convert back when done?

      #  $mapper = $asma->fetch_by_CoordSystems($slice_cs, $feat_cs);

      #   next unless defined $mapper;

      # Get list of coordinates and corresponding internal ids for
      # regions the slice spans
      #   @coords = $mapper->map($slice_seq_region, $slice_start, $slice_end,
      #                         $slice_strand, $slice_cs);

      #  @coords = grep {!$_->isa('Bio::EnsEMBL::Mapper::Gap')} @coords;

      #  next COORD_SYSTEM if(!@coords);

      # @ids = map {$_->id()} @coords;
      ##coords are now id rather than name
      ##      @ids = @{$asma->seq_regions_to_ids($feat_cs, \@ids)};

      # When regions are large and only partially spanned by slice
      # it is faster to to limit the query with start and end constraints.
      # Take simple approach: use regional constraints if there are less
      # than a specific number of regions covered.

      #  if(@coords > $MAX_SPLIT_QUERY_SEQ_REGIONS) {
      #    my $constraint = $orig_constraint;
      #    my $id_str = join(',', @ids);
      #    $constraint .= " AND " if($constraint);
      #    $constraint .= "${tab_syn}.seq_region_id IN ($id_str)";

      #   my $fs = $self->generic_fetch($constraint, $mapper, $slice);

      #       $fs = _remap($fs, $mapper, $slice);

      #      push @features, @$fs;

      #   } else {
      # do multiple split queries using start / end constraints

      #    my $max_len = $self->_max_feature_length() ||
      #    $mcc->fetch_max_length_by_CoordSystem_feature_type($feat_cs,
      #                                                         $tab_name);
      #   my $len = @coords;
      # for(my $i = 0; $i < $len; $i++) {
      #   my $constraint = $orig_constraint;
      #   $constraint .= " AND " if($constraint);
      #   $constraint .=
      #       "${tab_syn}.seq_region_id = "     . $ids[$i] . " AND " .
      #       "${tab_syn}.seq_region_start <= " . $coords[$i]->end() . " AND ".
      #       "${tab_syn}.seq_region_end >= "   . $coords[$i]->start();

      #   if($max_len) {
      #     my $min_start = $coords[$i]->start() - $max_len;
      #     $constraint .=
      #       " AND ${tab_syn}.seq_region_start >= $min_start";
      #   }

      #  my $fs = $self->generic_fetch($constraint,$mapper,$slice);

      #         $fs = _remap($fs, $mapper, $slice);

      #          push @features, @$fs;
      #    }
      #  }
    }
  }                             #COORD system loop

  return \@features;
}




#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#

#We have to have this here as this is a sub not a method, hence the _remap
#in the core BaseFeatureAdaptor is not available here

#Need to over-ride this for Iterator and remap per feature?
#Can we skip this altogether if the check test is done before hand?
#Would have to account for each CS set of features
#but we are currently only handling the query slice CS, so we can skip this completely for now


sub _remap {
  my ($self, $features, $mapper, $slice) = @_;

  #check if any remapping is actually needed
  if (@$features && (!$features->[0]->isa('Bio::EnsEMBL::Feature') ||
                     $features->[0]->slice == $slice)) {
    return $features;
  }

  #remapping has not been done, we have to do our own conversion from
  #to slice coords

  my @out;

  my $slice_start = $slice->start();
  my $slice_end   = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs    = $slice->coord_system();

  my ($seq_region, $start, $end, $strand);

  #my $slice_seq_region_id = $slice->get_seq_region_id();
  my $slice_seq_region = $slice->seq_region_name();

  foreach my $f (@$features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();

    if (!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region = $fslice->seq_region_name();
    my $fseq_region_id = $fslice->get_seq_region_id();
    my $fcs = $fslice->coord_system();

    if (!$slice_cs->equals($fcs)) {
      #slice of feature in different coord system, mapping required

      ($seq_region, $start, $end, $strand) =
        $mapper->fastmap($fseq_region_id,$f->start(),$f->end(),$f->strand(),$fcs);

      # undefined start means gap
      next if(!defined $start);
    } else {
      $start      = $f->start();
      $end        = $f->end();
      $strand     = $f->strand();
      $seq_region = $f->slice->seq_region_name();
    }

    # maps to region outside desired area
    next if ($start > $slice_end) || ($end < $slice_start) ||
      ($slice_seq_region ne $seq_region);

    #shift the feature start, end and strand in one call
    if ($slice_strand == -1) {
      $f->move( $slice_end - $end + 1, $slice_end - $start + 1, $strand * -1 );
    } else {
      $f->move( $start - $slice_start + 1, $end - $slice_start + 1, $strand );
    }

    $f->slice($slice);

    push @out,$f;
  }

  return \@out;
}


=head2 fetch_all_by_stable_Storable_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Storable
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

#We really only want this method to work on storables
#Namely all SetFeatures, could extend to FeatureTypes
#This is a reimplementation of what has been used in the
#Funcgen::SliceAdaptor::_set_bounds_by_xref_FeatureSet

#we could have a fourth arg here which is a coderef to filter the feature returned?
#This could be used by the SliceAdaptor to only return features which exceed a current slice
#We would still have to test the start and end in the caller, so not a massive win
#for too much code obfucation.  Altho we are maintaining two versions of this code
#which could get out of sync, epsecially with respect to external_db names

#This is currently feature generic
#So can call from one adaptor but will return features from all adaptors?
#Where can we put this?
#DBEntry adaptor?
#The only logical place for this is really a convinience method in the core BaseFeatureAdaptor
#This removes the problem of the generic nature of the returned data, as the focus is with the
#caller which is a singular feature, therefore, not mixing of data types

#restrict to one feature type for now. but this will be slower as we will have to make the list call
#for each type of feature_set passed, so if these are mixed they will be redundant
#Let's just put a warning in
#Keep this structure so we can port it to a convinience method somewhere
#and slim down later?


#This does not descend Gene>Transcript>Translation

#Now we have conflicting standards
#FeatureSets was kep as a list to prevent having to pass [$fset] for one FeatureSet
#This also elimated the need to test ref == ARRAY
#But now we want an extra optional descend flag
#Should we change back?
#Or should we wrap this method in one which will decend?
#Write wrappers for now, as we will have to write the descend loop anyway
#But this does not resolve the @fsets vs [$fset] issue

#Will have to do this in the wrappers anyway, so change back
#and live with the dichotomy of FeatureSet implementations for now

sub fetch_all_by_stable_Storable_FeatureSets{
  my ($self, $obj, $fsets) = @_;

  my ($extdb_name);
  my $dbe_adaptor = $self->db->get_DBEntryAdaptor;


  #Do we need a central registry for ensembl db names, what about versions?


  if (ref($obj) && $obj->isa('Bio::EnsEMBL::Storable') && $obj->can('stable_id')) {
    my @tmp = split/::/, ref($obj);
    my $obj_type = pop @tmp;
    my $group = $obj->adaptor->db->group;
    #Could sanity check for funcgen here?


    if (! defined $group) {
      throw('You must pass a stable Bio::EnsEMBL::Feature with an attached DBAdaptor with the group attribute set');
    }

    $extdb_name = 'ensembl_'.$group.'_'.$obj_type;

    #warn "ext db name is $extdb_name";

  } else {
    throw('Must pass a stable Bio::EnsEMBL::Feature, you passed a '.$obj);
  }


  #Set which eFG features we want to look at.

  if (ref($fsets) ne 'ARRAY' || scalar(@$fsets) == 0) {
    throw('Must define an array of Bio::EnsEMBL::FeatureSets to extend xref Slice bound. You passed: '.$fsets);
  }

  my %feature_set_types;

  foreach my $fset (@$fsets) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);

    $feature_set_types{$fset->feature_class} ||= [];
    push @{$feature_set_types{$fset->feature_class}}, $fset;
  }


  #We can list the outer loop here and put in the BaseFeatureAdaptor, or possible storable as we do have FeatureType xrefs.
  #This would be useful for fetching all the efg features for a given xref and FeatureSets
  #Don't implement as a parent sub and call from here as this would mean looping through array twice.
  #Altho we could pass a code ref to do the filtering?
  #Just copy and paste for now to avoid obfuscation

  my @features;

  #Get xrefs for each eFG feature set type
  foreach my $fset_type (keys %feature_set_types) {

    #my $xref_method    = 'list_'.$fset_type.'_feature_ids_by_extid';
    #e.g. list_regulatory_feature_ids_by_extid


    #Do type test here
    my $adaptor_type = ucfirst($fset_type).'FeatureAdaptor';

    #remove this for generic method
    next if ref($self) !~ /$adaptor_type/;

    #This is for generic method
    #my $adaptor_method = 'get_'.ucfirst($fset_type).'FeatureAdaptor';
    #my $adaptor = $self->db->$adaptor_method;

    my %feature_set_ids;
    map $feature_set_ids{$_->dbID} = 1, @{$feature_set_types{$fset_type}};

    my $cnt = 0;

    #Change this self to adaptor for generic method
    foreach my $efg_feature (@{$self->fetch_all_by_external_name($obj->stable_id, $extdb_name)}) {

      #Skip if it is not in one of our FeatureSets
      next if ! exists $feature_set_ids{$efg_feature->feature_set->dbID};
      push @features, $efg_feature;
    }
  }


  return \@features;
}


sub fetch_all_by_Gene_FeatureSets{
  my ($self, $gene, $fsets, $dblinks) = @_;

  if (! ( ref($gene) && $gene->isa('Bio::EnsEMBL::Gene'))) {
    throw("You must pass a valid Bio::EnsEMBL::Gene object");
  }

  my @features = @{$self->fetch_all_by_stable_Storable_FeatureSets($gene, $fsets)};

  if ($dblinks) {

    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      push @features, @{$self->fetch_all_by_Transcript_FeatureSets($transcript, $fsets, $dblinks)};
    }
  }

  return \@features;
}

sub fetch_all_by_Transcript_FeatureSets{
  my ($self, $transc, $fsets, $dblinks) = @_;

  if (! ( ref($transc) && $transc->isa('Bio::EnsEMBL::Transcript'))) {
    throw("You must pass a valid Bio::EnsEMBL::Transcript object");
  }


  my @features = @{$self->fetch_all_by_stable_Storable_FeatureSets($transc, $fsets)};

  if ($dblinks) {
    my $translation = $transc->translation;
    push @features, @{$self->fetch_all_by_stable_Storable_FeatureSets($translation, $fsets)} if $translation;
  }

  return \@features;
}



=head2 fetch_all_by_display_label

  Arg [1]    : String $label - display label of feature to fetch
  Example    : my $feat = $adaptor->fetch_all_by_display_label('abd-A_dpp:REDFLY:TF000092');
  Description: Returns the features with the given display label.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Feature objects
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

#There is no result_feature.display_label attribute.
#Move this to individual feature adaptors to avoid over-riding?


sub fetch_all_by_display_label {
  my ($self, $label) = @_;

  throw('You must defined a display_label argument') if ! defined $label;

  my $table_syn = $self->_main_table->[1];
  my $constraint = "${table_syn}.display_label = ?";
  $self->bind_param_generic_fetch($label, SQL_VARCHAR);

  return $self->generic_fetch($constraint);
}




=head2 _get_coord_system_ids

  Arg [1]    : optional - ARRAYREF of Bio::EnsEMBL::Funcgen::CoordSystem objects
  Example    : my $cs_ids = $self->_get_coord_system_ids([$cs1, $cs2], ...);
  Description: Returns the IDs of the CoordSystems specified or all default CoordSystems.
  Returntype : ARRAYREF
  Exceptions : Throws if CoordSystem object are not stored and valid
  Caller     : BaseFeatureAdaptor.pm
  Status     : At risk

=cut

#Can we remove the need for this by restricting the sr cache to default entries?

sub _get_coord_system_ids{
  my ($self, $coord_systems) = @_;

  my @cs_ids;

  if ($coord_systems) {

    if (ref($coord_systems) eq 'ARRAY' &&
        (scalar(@$coord_systems) >0)) {

      foreach my $cs (@$coord_systems) {
        $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CoordSystem', $cs);
        push @cs_ids, $cs->dbID;
      }
    } else {
      throw('CoordSystems parameter must be an arrayref of one or more Bio::EnsEMBL::Funcgen::CoordSystems');
    }
  } else {

    #Get current default cs's
    foreach my $cs (@{$self->db->get_FGCoordSystemAdaptor->fetch_all($self->db->dnadb, 'DEFAULT')}) {
      push @cs_ids, $cs->dbID;
    }
  }

  #This should never happen
  if (scalar(@cs_ids) == 0) {
    throw('Could not find any default CoordSystems');
  }

  return \@cs_ids;
}


=head2 count_features_by_field_id

  Arg [1]    : string     - table field to count
  Arg [2]    : string/int - id to count
  Example    : my $probe_feature_count = $pfa->count_feature_by_field_id('probe_id', $probe_id);
  Description: Returns a count of the features with the acciated field id
  Returntype : string/int - count of features
  Exceptions : Throws args are not set
  Caller     : FeatureAdaptors
  Status     : At risk

=cut

#This does not assume one record per feature
#But does assume primary key if ${table_name}_id
#Can't move to core due to cs issues, but could mirror implementation.

sub count_features_by_field_id{
  my ($self, $count_field, $count_id) = @_;
  #Any other params here?

  if (! ($count_field && $count_id)) {
    throw('You must provide a count name and a count id to count by');
  }

  my ($table_name, $table_syn) = @{$self->_main_table};
  my $table_id   = "${table_name}_id";
  my @cs_ids     = @{$self->_get_coord_system_ids};
  my $sql = "SELECT count(distinct($table_id)) from $table_name $table_syn, seq_region sr where ${table_syn}.${count_field}=? and ${table_syn}.seq_region_id in(select distinct(seq_region_id) from seq_region where coord_system_id in(".join(',', @cs_ids).'))';
  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $count_id, SQL_INTEGER);
  $sth->execute;

  return $sth->fetchrow_array;
}



=head2 force_reslice

  Arg [1]    : Optional - Boolean
  Example    : if($self->force_reslice){
                    # Reslice features past ends of destination Slice
               }
  Description: Sets/Returns force_reslice boolean
  Returntype : Boolean
  Exceptions : None
  Caller     : FeatureAdaptors::_objs_from_sth
  Status     : At risk

=cut

sub force_reslice{
  my $self  = shift;
  my $force = shift;

  if (defined $force) {
    $self->{force_reslice} = $force;
  }

  return $self->{force_reslice};
}


# Currently over-riding this method in SetFeatureAdaptor
# But might be nice to remove that and genericise this by getting 
# the analysis_id table via a method which can be over-ridden e.g. 
# SetFeatureAdaptor::_analysis_id_table_syn which would return the relevant set table name and syn
# BaseFeatureAdaptor::_analysis_id_table_syn would simply return the _main_table syn
# 

#Should this use the cache via AnalysisAdpator::fetch_by_logic_name?
#As _logic_name_to_constriant does?
#_logic_name_to_constraint returns undef if logic_name is not known
# allowing caller to return [] early.
# How would this work in compose_constraint?
# simply delete this first if exists
# How does this work wrt redundancy between constraints and optional_constraints

# For generic methods, there is a slight redundancy in calling
# the _main_table method
# Can we call this in compose_constraint and pass to contrain method

sub _constrain_logic_names {
  my $self        = shift;
  my $logic_names = shift;
  assert_ref($logic_names, 'ARRAY');
  
  if(! @$logic_names){
    throw('Must pass an Arrayref of logic_name strings to contrain by');  
  }

  for my $lname(@$logic_names){
    if(! defined $lname){
      throw('Found undefined logic_name value');  
    }  
  }

  my $syn = $self->_main_table->[1];
  $self->_tables([['analysis', 'a']]);

  my $constraint = $syn.'.analysis_id = a.analysis_id and a.logic_name in ("'.
    join('", "', @$logic_names).'")';
    
  return ($constraint, {});   #{} = not further constraint conf
}

1;


