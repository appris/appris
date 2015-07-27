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

Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor

=head1 SYNOPSIS

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => "ensembldb.ensembl.org",
   -dbname => "mus_musculus_funcgen_67_37",
   -species => "Mus_musculus",
   -user => "anonymous",
   -dnadb => $mouse_core_db,
   -port => '3307',
  );

my $experiment_adaptor = $db->get_ExperimentAdaptor();

=back

=head1 DESCRIPTION

An adaptor to access the funcgen database and expose other available adaptors.

=cut

################################################################################


package Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use strict;
use warnings;
use DBI; #for resolving core DB
use Bio::EnsEMBL::Utils::Exception         qw( throw deprecate warning );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( assert_ref_do );
use Bio::EnsEMBL::Registry;

use base qw(Bio::EnsEMBL::DBSQL::DBAdaptor);
my $reg = "Bio::EnsEMBL::Registry";


=head2 new

  Arg [-DNADB_HOST]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_USER]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_PASS]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_PORT]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_NAME]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_ASSEMBLY] : String - Overrides defaults (ensembldb or registry)

  Arg [...]         : Other args are passed to superclass Bio::EnsEMBL::DBSQL::DBAdaptor
  Example1          : $db = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
                              (
                               -user   => 'readonly', #No password
                               -dbname => 'pog',
                               -host   => 'caldy',
                              );

                      #If dnadb is not defined in registry, this will automatically
                      #set it from the default dnadb_host (e.g. ensembldb)

  Exmaple2          : $db = new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
                               -user           => 'write',
                               -pass           => 'password',
                               -dbname         => 'pog',
                               -host           => 'caldy',
                               -dnadb_assmebly => '36',
						                  );
                       #This will specifically look for a dnadb with assembly version 36
                       #on the default dnadb_host

   Exmaple2          : $db = new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
                               -user           => 'write',
                               -pass           => 'password',
                               -dbname         => 'pog',
                               -host           => 'caldy',
                               #The following will over-ride the default dnadb setting
                               -dnadb_name     => 'my_homo_sapiens_funcgen_67_37',
                               -dnadb_host     => 'my_host',
                               #Can add more dnadb params if required
                               #-dnadb_host
                               #-dnadb_user
                               #-dnadb_port
						                  );

  Example3            : $db new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
                               -user           => 'write',
                               -pass           => 'password',
                               -dbname         => 'pog',
                               -host           => 'caldy',
                               -dnadb          => $dnadb_object,
                              );

  Description: Constructor for DBAdaptor. Will automatically set the dnadb based on dnadb params.
               This makes some assumptions about how the DBs name are defined i.e. the last part must
               conform to the _RELEASE_ASSEMBLY format e.g. 67_37 in homo_sapiens_funcgen_67_37
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : Throws if conflicting dnadb params found
  Caller     : general
  Status     : Stable

=cut

#a main problem here is that we try to autoset the species
#this should really be done in Config Registry or possible SUPER new via reg->find_and_add_alias
#
#In the interim do we need to force the specification of species?
#This is required for automatic dnadb setting
#we could try and parse it from the dbname (trinomial names would be an issue)

#Trick is to do dnadb setting stuff first?
#still won't help as this doesn't set the species properly

#Woo woo! AY fixed this in the registry
#Species will never be DEFAULT
#Hence species.production_name will be mandatory
#We can change this if required
#We don't need to set species if we have a dnadb (params) passed
#but we can't validate. We could do some soft validation regex vs dbname? Too much?

sub new {
  my ($class, @args) = @_;

  #Force group to be funcgen as this is the only valid group.
  #This is now also in ConfigRegistry::gen_load
  my $self = $class->SUPER::new(@args, '-group', 'funcgen');

  #DEFAULT species handling now done in ConfigRegistry.pm
  #via gen_load which is called from SUPER new
  #gen_load throws if no species param is defined
  #and species.production_name is also absent

  #Should we need to do this if we define a valid dnadb?
  #reg aliases still not loaded


  my ( $dnadb_host, $dnadb_user, $dnadb_port, $dnadb_pass, $dnadb_assm, $dnadb_name, $dnadb)
    = rearrange([
                 'DNADB_HOST', 'DNADB_USER',
                 'DNADB_PORT', 'DNADB_PASS',
                 'DNADB_ASSEMBLY', 'DNADB_NAME',
                 'DNADB'
                ], @args );

  my $default_dnadb = $self->SUPER::dnadb;
  my ($default_host, $default_port, $default_user, $default_pass, $default_assm, $default_name, $efg_assm);
  my ($dnadb_predefined, $dnadb_params);

  if ( $dnadb_host || $dnadb_user || $dnadb_port || $dnadb_pass || $dnadb_assm || $dnadb_name) {
    $dnadb_params = 1;

    if ($dnadb) {
      throw('You cannot specific -dnadb and other dnadb params');
    }
  }


  #This is currently searching for the DB one ensembldb despite being predefined as -dnadb!

  if ($default_dnadb->group eq 'core') {
    #is 'funcgen' if not already by passing dnadb or loading registry
    $default_host = $default_dnadb->dbc->host;
    $default_port = $default_dnadb->dbc->port;
    $default_user = $default_dnadb->dbc->username;
    $default_pass = $default_dnadb->dbc->password;
    $default_name = $default_dnadb->dbc->dbname;
    ($default_assm = (split/_/, $self->_get_schema_build($default_dnadb))[1]) =~ s/[a-z]//;
    $dnadb_predefined = 1;
  }


  #We need to test dnadb_name vs dnadb_assm here before we over-ride we set defaults
  #We can expect a mistmatch if we have only defined one or the other
  if ( ($dnadb_assm && $dnadb_name) &&
       ($dnadb_name !~ /_${dnadb_assm}_/) ) {
    throw("You have specified conflicting -dnadb_name(${dnadb_name}) and -dnadb_assembly(${dnadb_assm}) parameters");
  } elsif ($dnadb_name) {       #Get dnadb_assm from name
    #This is not strictly required, but means we don't set is incorrectly below

    if ($dnadb_name =~ /.*_([0-9]+[a-z]*)$/) {
      $dnadb_assm = $1;
    }
  }

  #Defaults now set in dnadb as we want to test registry first;
  #These will pick default_dnadb values if params not explicitly set
  $self->{'dnadb_host'} = $dnadb_host || $default_host || 'ensembldb.ensembl.org';
  $self->{'dnadb_port'} = $dnadb_port || $default_port || 3306;
  $self->{'dnadb_user'} = $dnadb_user || $default_user || 'anonymous';
  $self->{'dnadb_pass'} = $dnadb_pass || $default_pass || undef;
  $self->{'dnadb_name'} = $dnadb_name || $default_name || undef;
  ($efg_assm = (split/_/, $self->_get_schema_build($self))[1]) =~ s/[a-z]//;


  $self->{'dnadb_assm'} = $dnadb_assm || $default_assm || $efg_assm;
  $dnadb_assm = $self->{'dnadb_assm'}; #reset here as we use below.

  #This only tries to _set_dnadb if we set some dnadb_params
  #or the dnadb_assm doesn't match the default/predefined dnadb

  if ($dnadb_params ||
      ($self->_get_schema_build($self->dnadb()) !~ /[0-9]+_${dnadb_assm}[a-z]*$/) ) {

    if (! $dnadb_params) {
      warn ':: WARNING: Unable to match assembly version between the dnadb name ('.
        $self->dnadb->dbc->dbname.') and the specified -dnadb_assm '.$self->dnadb_assembly.
          "\nMaybe you need to rename your DBs according to the Ensembl naming convention e.g. myprefix_homo_sapiens_55_37";
    } elsif ($dnadb_predefined) {
      #No can't be dnadb as we throw if we have conflicting -dnadb and dnadb params
      warn ":: Over-riding pre-defined dnadb regsitry values with dnadb params:\t".
        $self->dnadb_user.'@'.$self->dnadb_host.':'.$self->dnadb_port;
    }

    $self->_set_dnadb;
  }

  return $self;
}


#This should be in Storable to mirror core is_stored method?
#These do not fit in Storable, move these stored methods to BaseAdaptor?

=head2 is_stored_and_valid

  Arg [1]    : String - class namespace
  Arg [2]    : Bio::EnsEMBL::Funcgen::Storable e.g. ResultSet etc.
  Arg [3]    : String (optional) - Name of variable to use in error output (for use with assert_ref)
  Example    : $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  DESCRIPTION: Validates object class and stored status
  Returntype : None
  Exceptions : Throws if Storable is not stored
  Caller     : General
  Status     : At risk

=cut

sub is_stored_and_valid{
  my ($self, $class, $obj, $name) = @_;
  assert_ref($obj, $class, $name);
  throw("$obj is not stored") if ! $obj->is_stored($self);
  return;
}


=head2 are_stored_and_valid

  Arg [1]    : String - Namespace of class
  Arg [2]    : ARRAYREF os Bio::EnsEMBL::Funcgen::Storable objects e.g. ResultSet
  Arg [3]    : String : return value method name
  Example    : $db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', \@rsets);
  DESCRIPTION: Wrapper for is_stored_and_valid. Will optionally return array of values
               defined by calling method name arg on each object passed
  Returntype : ARRAYREF - contents defined by optional method name arg
  Exceptions : Throws if object list is not an ARRAY with at least one element
  Caller     : general
  Status     : At risk

=cut

sub are_stored_and_valid{
  my ($self, $class, $obj_list, $method_name) = @_;
  assert_ref($obj_list, 'ARRAY', 'object list');

  if(scalar(@$obj_list) == 0){
   throw('Objects Arrayref is empty');
  }

  my @return_vals;

  foreach my $obj (@$obj_list) {
    $self->is_stored_and_valid($class, $obj);

    if(! $method_name){
      assert_ref($obj, $class, 'object');
    }
    else{
      push @return_vals, assert_ref_do($obj, $class, $method_name, 'object');
    }
  }

  return \@return_vals;


}



#Move these to Helper.pm! Check method dependencies first!

=head2 load_table_data

  Arg [1]    : string - table name
  Arg [1]    : string - file path for file to load
  Example    : $db->load_table_data("result",  $self->get_dir($results_dir)."/result.txt");
  DESCRIPTION: Generic method to load a file into a specified table
  Returntype : none
  Exceptions : Throws if argument not supplied
  Caller     : general
  Status     : At risk - only used by for results at present, to be removed

=cut

sub load_table_data{
  my ($self, $table, $file, $ssh) = @_;
  my $cmd = 'mysqlimport -L '.$self->connect_string().' '.$file;
  system($cmd) && throw("Failed to load data from $file\nExit code:\t".($?>>8)."\n$!");
  return;
}



=head2 get_available_adaptors

  Example    : my %pairs = %{$dba->get_available_adaptors()};
  Description: gets a hash of the available adaptors
  ReturnType : reference to a hash
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::ConfigRegistry
  Status     : Stable

=cut


#will adding SliceAdaptor here use the dna DB? i.e. the core DB rather than the efg DB?

sub get_available_adaptors{
  my ($self) = shift;

  my %pairs = (
               'AnnotatedFeature'       => 'Bio::EnsEMBL::Funcgen::DBSQL::AnnotatedFeatureAdaptor',
               'Array'                  => 'Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor',
               'ArrayChip'              => 'Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor',
               'BindingMatrix'          => 'Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor',
               'CellType'               => 'Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor',
               'CoordSystem'            => 'Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor',
               'DataSet'                => 'Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor',
               'DataSet'                => 'Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor',
               'DBEntry'                => 'Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor',
               'DNAMethylationFeature'  => 'Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor',
               'Experiment'             => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor',
               'ExperimentalChip'       => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor',
               'ExperimentalGroup'      => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalGroupAdaptor',
               'ExternalFeature'        => 'Bio::EnsEMBL::Funcgen::DBSQL::ExternalFeatureAdaptor',
               'FeatureSet'             => 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor',
               'FeatureType'            => 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor',
               'FGCoordSystem'          => 'Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor',
               'InputSet'               => 'Bio::EnsEMBL::Funcgen::DBSQL::InputSetAdaptor',
               'InputSubset'            => 'Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor',
               'MirnaTargetFeature'     => 'Bio::EnsEMBL::Funcgen::DBSQL::MirnaTargetFeatureAdaptor',
               'MetaCoordContainer'     => 'Bio::EnsEMBL::Funcgen::DBSQL::MetaCoordContainer',
               'MotifFeature'           => 'Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor',
               'Probe'                  => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor',
               'ProbeFeature'           => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureAdaptor',
               'ProbeSet'               => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor',
               'RegulatoryFeature'      => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor',
               'ResultFeature'          => 'Bio::EnsEMBL::Funcgen::DBSQL::ResultFeatureAdaptor',
               'ResultSet'              => 'Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor',
               'SegmentationFeature'    => 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationFeatureAdaptor',
               'Slice'                  => 'Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor',
               #prepended FG to override core adaptor. Now fixed in registry. Maintained for backwards compatibility
               'Channel'                => 'Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor',

               #add required EnsEMBL(core) adaptors here
               #Should write/retrieve from efg not dna db
               'UnmappedObject'     => 'Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor',
               'Analysis'           => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
               "MetaContainer"      => 'Bio::EnsEMBL::DBSQL::MetaContainer',
              );

  return (\%pairs);
}


=head2 _get_schema_build

  Arg [1]    : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor or Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : my $schema_build = $db->_get_schema_build($slice->adaptor->db());
  DESCRIPTION: Given a DBAdaptor, parses the the 'SCHEMA_BUILD' string from the DB name.
  Returntype : String or undef - 'SCHEMA_BUILD' as per the standard DB naming.
  Exceptions : Throws if DBAdaptor not passed 
               Warns if unable to match expected schema_build string format from the DB name
               and returns undef.
  Caller     : General
  Status     : At risk - replace with MetaContainer method

=cut

# This is actually more accurately the DB suffix, as 
# 1 The build is now just normally the assembly number (with some exceptions in EG)
# 2 It may not be the current default schema_build used in the API if the coord_system.is_current flag
#   has been moved to a different schema_build i.e. after the update_DB_for_release.pl script has been run

sub _get_schema_build{
  my $self = shift;
  my $db   = shift;

  assert_ref($db, 'Bio::EnsEMBL::DBSQL::DBAdaptor');
  my $schema_build;
  my $name = $db->dbc->dbname;

  if ($name =~ /.*_([0-9]+_[0-9]+[a-z]*)$/o) {
    $schema_build = $1;
  
    if (length($schema_build) > 10) {
      $schema_build = substr($schema_build, -10);
    }
  }
  else {
    warning("Wrong format: '$name' Release & Assembly expected at the end of dbname, e.g.: *_core_75_37");
  }

  return $schema_build;
}


#Remove all this dnadb setter functionality

sub dnadb_host{
  my ($self, $host) = @_;
  $self->{'dnadb_host'} = $host if $host;
  return $self->{'dnadb_host'};
}

sub dnadb_name{
  return $_[0]->{dnadb_name};
}

sub dnadb_port{
  my ($self, $port) = @_;
  $self->{'dnadb_port'} = $port if $port;
  return $self->{'dnadb_port'};
}

sub dnadb_pass{
  my ($self, $pass) = @_;
  $self->{'dnadb_pass'} = $pass if $pass;
  return $self->{'dnadb_pass'};
}

sub dnadb_user{
  my ($self, $user) = @_;
  $self->{'dnadb_user'} = $user if $user;
  return $self->{'dnadb_user'};
}

sub dnadb_assembly{
  my ($self, $assm) = @_;
  $self->{'dnadb_assm'} = $assm if $assm;
  return $self->{'dnadb_assm'};
}

#Redefine dbadb here to add coordsystem

=head2 dnadb

  Arg [1]:     Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [2]:     string - coord_system name e.g. chromosome
  Usage :      my $dnadb = $db->dnadb();
  Description: returns the database adaptor where the dna lives i.e. the core db for a given species
               There are at least 2 cases where you need to set this explicitly
               1.  If you want to retrieve features on an assembly which is not the default in
               the correspeonding core DB with matching schema_build
               2.  If the corresponding core DB is not available on the default ensembl DB
               server(ensembldb/ens-livemirror) i.e. before a new release.
  Status :     At risk. - Might remove validation of CS

=cut

#why was the SUPER::dnadb->group ever ne core?

sub dnadb {
  my ($self, $dnadb, $cs_name) = @_;

  #warn "dnadb x $dnadb x || ".$self->SUPER::dnadb->group();
  #dnadb is passed directly to this method from super new

  if ($dnadb || $self->SUPER::dnadb->group() ne 'core') {
    #what

    if (! $dnadb) {             #Guess/set dnadb by assembly or dnadb params
      #warn "calling _set_dnadb";
      return $self->_set_dnadb;
    }

    #warn "setting via SUPER";

    $self->SUPER::dnadb($dnadb);

    #set default coordsystem here
    #there might not be a chromosome level if we just have a scaffold assembly
    #supercontig is already loaded as we use toplevel?
    #This should really get all the default CSs from the core CSAdaptor?
    #We should also do this during the update?
    #This will also enable people to query using clone/contig/supercontig level slices
    #How will this work? Where will the mapping between CSs be done?

    my @cs_names;
    @cs_names = ($cs_name) if $cs_name;
    #$cs_name ||= 'chromosome';

    if (! $cs_name) {

      foreach my $cs (@{$dnadb->get_CoordSystemAdaptor->fetch_all_by_attrib('default_version')}) {
        push @cs_names, $cs->name;
      }
    }


    foreach my $cs_name (@cs_names) {

      my $cs;
      eval { $cs = $dnadb->get_CoordSystemAdaptor->fetch_by_name($cs_name)};
      my $error = $@;

      if ($error) {
        my ($schema, $build) = split/_/, $self->_get_schema_build($dnadb);
        $build =~ s/[a-z]//;
        throw("It appears that the schema of ".$dnadb->dbc->dbname.
              ' is incompatible with your current core API version('.$reg->software_version.
              ").  You could try using the $schema version of the core API, or alternatively try specifying ".
              "different -dnadb/registry_host parameters to point to a make recent version containing build $build\n");

      }


      #this will only add the default assembly for this DB, if we're generating on another we need to add it separately.

      #!!! This is a non-obvious store behaviour !!!
      #This can result in coord_system entries being written unknowingly if you are using the efg DB with a write user
      $self->get_FGCoordSystemAdaptor->validate_and_store_coord_system($cs);
    }
  }

  return $self->SUPER::dnadb(); #never pass @_ here!
}


=head2 _set_dnadb

  Usage :      $self->_set_dnadb;
  Description: Sets the dnadb to the latest version given the assembly version
  Exceptions:  Throws if no assembly version provided or cannot for appropriate dnadb on ensembldb
  Caller:      Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor::new or dnadb
  Status :     At risk

=cut


sub _set_dnadb{
  my $self       = shift;
  my $assm_ver   = $self->dnadb_assembly;
  my $dnadb_name = $self->dnadb_name;

  #Sanity check, but this should always be the case when called from new/dnadb
  if (! ($assm_ver || $dnadb_name)) {
    throw("You need to have define at least a dnadb_assembly($assm_ver) or dnadb_name($dnadb_name) before calling this method");
  }

  throw('Must provide and assembly version to set the dnadb') if ! defined $assm_ver;

  my $reg_lspecies = $reg->get_alias($self->species());
  #The registry has incremented the species as we have recreated the efg DB
  #possibly using a different schema_build
  #This set true lspecies to allow dnadb detection
  #in multi DB environments e.g. DAS server
  my $lspecies = $reg_lspecies;
  $lspecies =~ s/[0-9]+$// if($lspecies =~ /[0-9]$/);

  if( $lspecies eq 'default'){
    throw('Either specify a species parameter or set species.production_name'.
      " in the meta table (DB: $dnadb_name) to set dnadb automatically, alternatively pass a dnadb parameter");
  }

  #Start with lastest MySQL instances
  #We are over-riding specified port here, only for known hosts
  #we should really account for this and make it nr
  my @ports;

  if($self->dnadb_port){
    @ports = ($self->dnadb_port);
  }
  elsif ($self->dnadb_host eq 'ensdb-archive') { #
    @ports = (5304, 3304);
  } elsif ($self->dnadb_host eq 'ensembldb.ensembl.org') {
    @ports = (5306, 3306);
  }

  my $match_name = $dnadb_name || $lspecies.'_core_%_'.$assm_ver.'%';


  my $sql = 'show databases like "'.$match_name.'"';
  my ($dbh, @dbnames, $port, $host_port);


  foreach $port(@ports) {
    #This is probably duplicating connections and over-riding any connection
    #pooling going on in the base DBConnection if we are using the same host port
    #as the registry connection

    $dbh = DBI->connect('DBI:mysql:host='.$self->dnadb_host.";port=${port}",
                        $self->dnadb_user,
                        $self->dnadb_pass,
                        {
                         'RaiseError' => 1});

    eval { @dbnames = map {$_ = "@$_"} @{$dbh->selectall_arrayref($sql)}; };

    if ($@) {
      throw('Failed to fetch dna DB names from '.$self->dnadb_host.":${port}"."\n$@");
    }


    #Will always take the latest release, not the latest genebuild version
    #Which is probably what we want anyway

    if(! $dnadb_name){
      @dbnames = grep(/core_[0-9]/, sort @dbnames);
    }


    if (scalar(@dbnames)==0) {
      warn(':: Failed to find dnadb like '.$match_name.', using '
           .$self->dnadb_user.'@'.$self->dnadb_host.':'.$port);
    } else {
      $host_port = $port;
      last;
    }
  }

  throw("Failed to find dnadb like $match_name.") if(scalar(@dbnames)==0);
  warn ":: Auto-selecting build $assm_ver core DB as:\t".
    $self->dnadb_user.'@'.$dbnames[$#dbnames].':'.$self->dnadb_host.':'.$host_port."\n";

  my $db = $reg->reset_DBAdaptor($reg_lspecies, 'core', $dbnames[$#dbnames],
                                 $self->dnadb_host, $host_port,
                                 $self->dnadb_user, $self->dnadb_pass);
  $self->dnadb($db);
  return $db;
}

=head2 connect_string

  Example    : my $import_cmd = 'mysqlimport '.$db->connect_string()." $table_file";
  Description: Retrieves the mysql cmdline connection string
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub connect_string{
  my $self = shift;

  return '-h'.$self->dbc->host().' -u'.$self->dbc->username().' -p'.$self->dbc->password()
    .' -P'.$self->dbc->port().' '.$self->dbc->dbname();
}


=head2 rollback_table

  Arg [1]    : String - SQL to execute
  Arg [2]    : String - Main table name
  Arg [3]    : String (optional) - Table ID for refresh_table
  Arg [4]    : Boolean (optional) - Do not refresh_table
  Arg [5]    : Boolean (optional) - Force clean up even if no records deleted
  Example    : $db->rollback_table($sql, 'result_set', 'result_set_id');
  Description: Executes the SQL and calls refresh_table method
  Returntype : None
  Exceptions : Throws if SQL or table name arguments are not set.
  Caller     : Helper rollback methods.
  Status     : At risk

=cut

#table is not really mandatory, and this can delete from multiple tables at a time

#TODO return $row_cnt?
#Use SQLHelper?

sub rollback_table {
  my ( $self, $sql, $table, $id_field, $no_clean_up, $force_clean_up ) = @_;
  my $row_cnt;

  if(! ($sql && $table)){
   throw('Must provide at least the SQL and table name arguments');
  }

  #warn $sql;
  if(! eval { $row_cnt = $self->dbc->do($sql); 1 }){
    throw("Failed to rollback table $table using sql:\t$sql\n$@");
  }

  $row_cnt = 0 if $row_cnt eq '0E0';
  #$self->log("Deleted $row_cnt $table records");

  if ( $force_clean_up || ( $row_cnt && !$no_clean_up ) ) {
    $self->refresh_table( $table, $id_field );
  }

  return;
}


=head2 refresh_table

  Arg [2]    : String - table name
  Arg [3]    : String (optional) - Table ID for refresh_table
  Example    : $db->refresh_table('result_set', 'result_set_id');
  Description: Calls reset_autoinc if table_id set, optimzes and analyses table
  Returntype : None
  Exceptions : Throws if table argument is not set
  Caller     : Helper rollback methods
  Status     : At risk

=cut

#Now separated so that we can do this once at the end of a rollback of many Sets

sub refresh_table {
  my ( $self, $table, $id_field ) = @_;

  if(! $table){
    throw('Must provide a table argument');
  }

  #This only works if the new calue is available
  #i.e. do not need lock for this to be safe
  $self->reset_table_autoinc( $table, $id_field ) if $id_field;

  #$self->log("Optimizing and Analyzing $table");

  $self->dbc->do("optimize table $table");  #defrag data, sorts indices, updates table stats
  $self->dbc->do("analyze  table $table");  #analyses key distribution
  return;
}


=head2 reset_table_autoinc

  Arg [1]    : String - table name
  Arg [2]    : String - Autoinc field
  Example    : $db->reset_table_autoinc('result_set', 'result_set_id');
  Description: Resets the autoinc field of a table to the next logical number.
               This is useful when importing and rolling back many entries.
               A tidyness thing really.
  Returntype : None
  Exceptions : Throws if table argument or auto_inc field arguments are not set
  Caller     : refresh_table and Helper rollback methods
  Status     : At risk

=cut

sub reset_table_autoinc {
  my ( $self, $table_name, $autoinc_field ) = @_;

  if ( !( $table_name && $autoinc_field ) ) {
    throw('You must pass a table_name and an autoinc_field'.
      ' to reset the autoinc value');
  }


  #Unsafe to do this in two queries as parallel jobs may add in between select and alter
  #in fact this needs a table lock to be totally safe
  #although current ALTER will just fail silently if this happens
  my $sql = "select $autoinc_field from $table_name order by $autoinc_field desc limit 1";
  my ($current_auto_inc) = $self->dbc->db_handle->selectrow_array($sql);
  my $new_autoinc = ($current_auto_inc) ? ( $current_auto_inc + 1 ) : 1;
  $sql = "ALTER TABLE $table_name AUTO_INCREMENT=$new_autoinc";
  $self->dbc->do($sql);
  return;
} ## end sub reset_table_autoinc



1;

