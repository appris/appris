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

Bio::EnsEMBL::Funcgen::Parsers::BaseImporter

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the start of merging the Importer with the  InputSet & BaseExternal Parsers.
This class holds all generic methods used for importing data into the funcgen schema.
Move all generic methods from Importer to here, and move format specific methods to new parsers.
Then remove Importer completely.

=head1 SEE ALSO


=cut

package Bio::EnsEMBL::Funcgen::Parsers::BaseImporter;

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( validate_path dump_data );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::FeatureType;

use base qw( Bio::EnsEMBL::Funcgen::Utils::Helper );

#new
#edb_release over-ride, to enable loading of old data.
#How do we handle set rollback? (with force for associated sets)
#

#Force was actually specific to store_window_bin_by_Slice_Parser
#

=head2 new

  Arg [1]    : HASH containing attributes:
                  -rollback    Performs full rollback of imported features.
                  -recover     Performs rollback of features/sets which do 
                               not have an associated Imported status.
           
                  -config_file Path to file containing config hash.
             
                  -
  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Bed class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Simple
  Exceptions : throws if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Parsers:Simple
  Status     : at risk

=cut




#Also implement rollback at the script level? Then we can re-use the methods to perform the rollback?


#Would just need to revoke states before calling validate_files
#Can do this manually before import


#Inheritance fix
#See ENSREGULATION-18 ticket in JIRA

#We really want to remove many of the mandatory params from new here and probably
#just enforce the generic params e.g. db etc but not import specific stuff i.e. ctyp/ftype
#write new method first before updating existing parser to use new init_method

sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new(@_);

  my $reg = "Bio::EnsEMBL::Registry"; 
  my ($config_file, $clobber, $rollback, $species, $fset_desc,
      $user, $host, $port, $pass, $dbname, $db, $ssh,
      $assm_version, $release, $reg_config, $verbose, $slices,
      $reg_db, $reg_host, $reg_port, $reg_user, $reg_pass,
      $ftype_name, $ctype_name, $feature_analysis, $no_disconnect,
      $input_files);
  

  #Set some directly here to allow faster accessor only methods
  #These will over-ride defaults on Parser config
  

  ($ftype_name,       $ctype_name,          $feature_analysis,    $self->{feature_set_desc},
   $species,          $db,                  $user,                $host,
   $port,             $pass,                $dbname,              $assm_version, $ssh,
   $release,          $reg_config,          $reg_db,              $reg_host,
   $reg_port,         $reg_user,            $reg_pass,            $verbose,
   $slices,           $self->{recover},     $clobber,             $rollback,
   $config_file,      $self->{ucsc_coords}, $self->{_dump_fasta}, $no_disconnect,
   $input_files,      $self->{input_dir},   $self->{data_dir},    $self->{output_dir},
   $self->{name}
  ) = rearrange
    (
     ['FEATURE_TYPE_NAME', 'CELL_TYPE_NAME', 'FEATURE_ANALYSIS', 'FEATURE_SET_DESCRIPTION',
      'species',           'db',             'user',             'host',
      'PORT',              'PASS',           'DBNAME',           'ASSEMBLY',      'SSH',
      'RELEASE',           'REG_CONFIG',     'REGISTRY_DB',      'REGISTRY_HOST',
      'REGISTRY_PORT',     'REGISTRY_USER',  'REGISTRY_PASS',    'VERBOSE',
      'slices',            'recover',        'clobber',          'rollback',
      'config_file',       'ucsc_coords',    'DUMP_FASTA',       'no_disconnect',
      'input_files',       'input_dir',      'data_dir',         'output_dir',
      'name'],
     @_);
  
  
 
  #Other stuff to bring in here:
  #  dirs
  #  prepared/batch/farm
  #  analysis/cell/feature_type/feature_set_description
  #    This is redundant wrt config leave for now and catch

  #Can we move all of this DB/reg handling stuff out
  #So it can be re-used by other modules/scripts?

  #Need to compare these to BaseExternalParser

  
  if($db){
    #Don't load the registry as this will 'increment' the db species
    #Currently assuming all is fine here, but ultimately. to avoid redundancy
    #migrate Base::_set_out_db and db to this module
    $self->{species} = $db->species;
  }
  else{

    $species || throw('Mandatory param -species not met');
 
    # Registry and DB handling - re/move to separate method?
    
    
    if(! $db){
    
    if ($reg_host && $self->{'reg_config'}) {
      warn "You have specified registry parameters and a config file:\t".$self->{'reg_config'}.
        "\nOver-riding config file with specified paramters:\t${reg_user}@${reg_host}:$reg_port";
    }
  
  
    #### Set up DBs and load and reconfig registry
    
    # Load Registry using assembly version
    # Then just redefine the efg DB
  
    #We have problems here if we try and load on a dev version, where no dev DBs are available on ensembldb
    #Get the latest API version for the assembly we want to use
    #Then load the registry from that version
    #Then we can remove some of the dnadb setting code below?
    #This may cause problems with API schema mismatches
    #Can we just test whether the current default dnadb contains the assembly?
    #Problem with this is that this will not have any other data e.g. genes etc 
    #which may be required for some parsers
  
    #How does the registry pick up the schema version??
  
    #We should really load the registry first given the dnadb assembly version
    #Then reset the eFG DB as appropriate
    $self->{'reg_config'} = $reg_config || ((-f "$ENV{'HOME'}/.ensembl_init") ? "$ENV{'HOME'}/.ensembl_init" : undef);
  
    if ($reg_host || ! defined $self->{'_reg_config'}) {
      #defaults to current ensembl DBs
      $reg_host ||= 'ensembldb.ensembl.org';
      $reg_user ||= 'anonymous';
  
      #Default to the most recent port for ensdb
      if ( (! $reg_port) && 
           ($reg_host eq 'ensdb-archive') ) {
        $reg_port = 5304;
      }
  
      #This will try and load the dev DBs if we are using v49 schema or API?
      #Need to be mindful about this when developing
      #we need to tip all this on it's head and load the reg from the dnadb version
  
      my $version = $release || $reg->software_version;
      $self->log("Loading v${version} registry from $reg_user".'@'.$reg_host);
  
      #Note this defaults API version, hence running with head code
      #And not specifying a release version will find not head version
      #DBs on ensembldb, resulting in an exception from reset_DBAdaptor below
     
      #This currently loads from reg_host
      #Which triggers the funcgen DB to try and _set_dnadb
      #even if we have specified as db/dnadb already,as this is reset after this
      #we probably just want to use the dnadb_host as the default reg_host
      #This also need cleaning up wrt DBAdaptor behaviours
      #Simply and/or remove
  
      my $num_dbs = $reg->load_registry_from_db
        (
         -host       => $reg_host,
         -user       => $reg_user,
         -port       => $reg_port,
         -pass       => $reg_pass,
         -db_version => $version,
         -verbose    => $verbose,
        );
  
      if (! $num_dbs) {
        #only throw if we don't have any other db params passed?
        throw("Failed to load any DBs from $reg_user".'@'.$reg_host." for release $version.\n".
             "This will result in a failure to validate the species.\n".
             "Please define a valid -release for the registry and/or registry params/config\n".
             "Or select a -registry_host which matches the API version:\t".$reg->software_version);
      }
  
      if ((! $dbname) && (! $db)) {
        throw('Not sensible to set the import DB as the default eFG DB from '
              .$reg_host.', please define db params');
      }
    } else {
      $self->log("Loading registry from:\t".$self->{'_reg_config'});
      $reg->load_all($self->{'_reg_config'}, 1);
    }
  
  
    #Need to test the DBs here, as we may not have loaded any!
    #get_alias will fail otherwise
    #This is a cyclical dependancy as we need alias to get species which we use to grab the DB
    #alias is dependant on core DB being loaded with relevant meta entries.
    #revise this when we split the Importer
  
    #Validate species
    my $alias = $reg->get_alias($species) || throw("Could not find valid species alias for $species");
    #You might want to clean up:\t".$self->get_dir('output'));
    $self->{'param_species'} = $species; #Only used for dir generation  
    $self->{species} = $alias;
    

    #define eFG DB from params or registry

    if ($reg_db) {              #load eFG DB from reg

      if ($dbname) {
        throw("You cannot specify DB params($dbname) and load from the registry at the same time.");
      }

      $self->log('WARNING: Loading eFG DB from Registry');
      $db = $reg->get_DBAdaptor($self->species(), 'funcgen');
      throw("Unable to retrieve ".$self->species." funcgen DB from the registry") if ! $db;
    } else {                    #resets the eFG DB in the custom or generic registry

      $dbname || throw('Must provide a -dbname if not using default custom registry config');
      $pass   || throw('Must provide a -pass parameter');
	 
      #remove this and throw?
      if (! defined $host) {
        $self->log('WARNING: Defaulting to localhost');
        $host = 'localhost';
      }
	  
      $port ||= 3306;
      my $host_ip = '127.0.0.1'; #is this valid for all localhosts?
	  
      if ($ssh) {
        $host = `host localhost`; #mac specific? nslookup localhost wont work on server/non-PC 
        #will this always be the same?
		
        if (! (exists $ENV{'EFG_HOST_IP'})) {
          warn "Environment variable EFG_HOST_IP not set for ssh mode, defaulting to $host_ip for $host";
        } else {
          $host_ip = $ENV{'EFG_HOST_IP'};
        }
		
        if ($self->host() ne 'localhost') {
          warn "Overriding host ".$self->host()." for ssh connection via localhost($host_ip)";
        }
      }
	

      #data version is only used if we don't want to define the dbname
      #This should never be guessed so don't need data_version here
      #$dbname ||= $self->species()."_funcgen_".$self->data_version();


      #Remove block below when we can
      my $dbhost = ($ssh) ? $host_ip : $host;

      #This isn't set yet!?
      #When we try to load e.g. 49, when we only have 48 on ensembldb
      #This fails because there is not DB set for v49, as it is not on ensembl DB
      #In this case we need to load from the previous version
      #Trap this and suggest using the -schema_version/release option
      #Can we autodetect this and reload the registry?
      #We want to reload the registry anyway with the right version corresponding to the dnadb
      #We could either test for the db in the registry or just pass the class.

      $db = $reg->reset_DBAdaptor($self->species, 'funcgen', $dbname, $dbhost, $port, $user, $pass,
                                  {
                                   -dnadb_host => $reg_host,
                                   -dnadb_port => $reg_port,
                                   -dnadb_assembly => $assm_version,
                                   -dnadb_user => $reg_user,
                                   -dnadb_pass => $reg_pass,
                                  });
      }
    }
  }


  #Test connections
  $self->{db} = $db; 
  $db->dbc->db_handle;
  $db->dnadb->dbc->db_handle;
  #Set re/disconnect options

  #These really need setting dependant on the import parser
  $db->dbc->disconnect_when_inactive(1)        if ! $no_disconnect;
  $db->dnadb->dbc->disconnect_when_inactive(1) if ! $no_disconnect;


  
  #Catch config clashes/redundancy
  
  if ($self->feature_set_description  &&
      ($config_file || exists ${$self}{static_config}) ) {
    throw('You have specified a -feature_set_description alongside user/static_config. Please define this in the config');
  }

  if ( ($feature_analysis || $ctype_name || $ftype_name) &&
       exists ${$self}{static_config} ) {
    throw('You have specified a analysis/cell/feature_type params alongside static_config. Please define this in the static_config');
  }

  #Catch no config here? Or can this be imported via Parser specific meta/config files


  ### Check analyses/feature_type/cell_type
  if ($feature_analysis) {
    my $fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($feature_analysis);
    throw("The Feature Analysis $feature_analysis does not exist in the database") if(!$fanal);
    $self->feature_analysis($fanal);

    #This currently fails before the config gets loaded!
    #Need to load config before this validation!
  }

  if ($ctype_name) {
    my $ctype = $self->db->get_CellTypeAdaptor->fetch_by_name($ctype_name);
    throw("The CellType $ctype_name does not exist in the database") if(!$ctype);
    $self->cell_type($ctype);
  }

  if ($ftype_name) {
    my $ftype = $self->db->get_FeatureTypeAdaptor->fetch_by_name($ftype_name);
    throw("The FeatureType $ftype_name does not exist in the database") if(!$ftype);
    $self->feature_type($ftype);
  }



  #Set some attrs to allow setter only methods
  $self->{slice_adaptor} = $db->dnadb->get_SliceAdaptor;
  $self->slices($slices) if defined $slices;
  $self->{rollback}         = $rollback || $clobber;
  $self->{counts}           = {};
  $self->{seq_region_names} = []; #Used for slice based import


  # USER CONFIG #
  #Here we need to read config based on external file
  #Should do something similar to set_feature_sets
  #and validate_and_store_feature_types in BaseExternalParser
  #but we are using define and validate sets instead

  #BaseExternalParser and BaseImporter really need to be merged
  #After we have stripped out all the array/experiment specific stuff


  if ($config_file) {
    my $config;

    $self->log("Reading config file:\t".$config_file);

    if (! ($config = do "$config_file")) {
      throw("Couldn't parse config file:\t$config_file:\n$@") if $@;
      throw("Couldn't do config:\t$config_file\n$!")          if ! defined $config;
      throw("Couldn't compile config_file:\t$config_file")    if ! $config;
    }

    #At least check it is hash
    if (ref($config) ne 'HASH') {
      throw("Config file does not define a valid HASH:\t$config_file");
    }
	
    $self->{user_config} = $config;	

    #Can call validate_and_store_config directly ehre once we have remove set_config stuff

  }


  $self->{'data_dir'} ||= $ENV{'EFG_DATA'};


  #if( (! $self->input_files($input_files)) &&
  #     (! defined $self->get_dir('input') ) ){ 
  #  #Set default input_dir if we have not specified files
  #   #This is dependant on name which is not mandatory yet!
  #    $self->{'input_dir'} = $self->get_dir("data").'/input/'.
  #    $self->{'param_species'}.'/'.$self->vendor().'/'.$self->name();   
  #}

  if(defined $self->get_dir('input')){
    validate_path($self->get_dir('input'), 1); #dir flag
  }
    
  


  #$self->debug(2, "BaseImporter class instance created.");
  #$self->debug_hash(3, \$self);
  
  return $self;
}







=head2 validate_and_store_config

  Args       : None
  Example    : $self->validate_and_store_config;
  Description: Imports feature types defined by import_sets.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk 

=cut

#Taken from InputSet
#Needs to support 'static' config from external parsers

#Need to get the config has format working with this and set_feature_sets?
#Or just migrate define_and_validate_sets here now?

#Now takes arrayref, need to change in other callers

#All analyses and ftypes now stored here(inc fset only defined), so don't need to do this in the define sets method

#define/set sets method should write to user config first, call validate_store, before defining sets?

#Should we allow empty hashes to fetch from DB?

sub validate_and_store_config{
  my ($self, $fset_names) = @_;

  if ( (ref($fset_names) ne 'ARRAY') ||
       (scalar(@$fset_names) == 0) ) {
    throw('Must pass FeatureSet names to validate_and_store_config');
  }

  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;


  #Do we want to enable no config? Just don't call this if we have no config!

  my ($static_config, $user_config, $fset_config);

  #Set here to avoid auto-vivifying in tests below
  $user_config   = $self->{user_config} if exists ${$self}{user_config};
  $static_config = $self->{static_config} if exists ${$self}{static_config};
  my $config = $user_config || $static_config;

  if (! $config) {
    #throw('No user or static config found');
    warn('No user or static config found');
  } elsif ($user_config && $static_config) {
    throw('BaseImporter does not yet support overriding static config with user config');
    #Would need to over-ride on a key by key basis, account for extra config in either static or user config?
  } else {
    #Store config for each feature set
    #inc associated feature_types and analyses
    #add cell_types in here?

    foreach my $import_set (@$fset_names) {
      #warn "validating config for $import_set";

	
      if (exists ${$config}{feature_sets}{$import_set}) {
        $fset_config = $config->{feature_sets}{$import_set};
      }
	
      if (! $fset_config) {
        throw("Could not find config for:\t$import_set");
      }

      $self->log("Validating and storing config for:\t$import_set");

	  
      #If we grab the ftype and analysis from the feature set first
      #Then we can remove the reundancy in the config
      #would need to handle case i.e. -ANALYSIS -analysis
      #use rearrange!

      #Set in analyses and feature_types config, should use same ref if already exist
      #else we have a duplicated entry using the same name, which may have different attrs!
	
      #This is going to be a problem as we are assigning a new value to the key?
      #Do we need to use references in feature_sets?

      if (exists ${$fset_config}{feature_set}) {
        #my ($fset_analysis, $fset_ftype) = rearrange(['ANALYSIS', 'FEATURE_TYPE'], %{$fset_config->{feature_set}});

        #Can't use rearrange for key we are setting as we need to now the case
        #This is grabbing top level config keys, not hashes
        #config names refer to top level keys and may not match logic_name or ftype name
        my $fset_params          = $fset_config->{feature_set};
        my $fset_analysis_key    = (exists ${$fset_params}{-analysis})            ? '-analysis'                        : '-ANALYSIS';
        my $analysis_config_name = (exists  ${$fset_params}{$fset_analysis_key} ) ? $fset_params->{$fset_analysis_key} : undef;
        my $fset_ftype_key       = (exists ${$fset_params}{-feature_type})        ? '-feature_type'                    : '-FEATURE_TYPE';
        my $ftype_config_name    = (exists ${$fset_params}{$fset_ftype_key})      ? $fset_params->{$fset_ftype_key}    : undef;
	  	  

        #Check we have these in feature_sets analyses/feature_types already?
        #Should that be a hash with empty {} values
        #Then we can use the top level analyses hash if required
        #Setting feature_set analysis/feature_type only if not present i.e. we don't over-write the top level
        #defining hash
        #Postponing checking here will lose context if we need to throw
        #i.e. we won't know where the error has come from if the name is not present in the top level hash
        #Can we call direct from here instead?

  

        #Leave to store to catch these mandatory attrs

        if ($analysis_config_name) {
		
          if (ref(\$analysis_config_name) ne 'SCALAR') {
            throw("You must set feature_set($import_set) -analysis config as a string value i.e. a key referencing the top level analyses config");
          }
		
          if (! exists ${$fset_config}{analyses}{$analysis_config_name}) { #Set in analyses to validate/store below
            $fset_config->{analyses}{$analysis_config_name}     = {}; #Don't have to set a ref to top level here, as these will all get set to the same obj ref in validate/store
          }

          #Only use refs in the feature_set
          #top level and feature_set analyses get set to same obj at same time		  
          $fset_config->{feature_set}{$fset_analysis_key} = \$fset_config->{analyses}{$analysis_config_name};
		  
        }

	  
        if ($ftype_config_name) {		
	
          if (ref(\$ftype_config_name) ne 'SCALAR') {
            throw("You must set feature_set($import_set) -feature_type config as a string value i.e. a key referencing the top level feature_types config");
          }

          if (! exists ${$fset_config}{feature_types}{$ftype_config_name}) { #Set in analyses to validate/store below
            $fset_config->{feature_types}{$ftype_config_name} = {};
          }
		
          $fset_config->{feature_set}{$fset_ftype_key}  = \$fset_config->{feature_types}{$ftype_config_name};

        }
      }

	
      #Can self ref user config if 'do' will work with %config, specified as the last line
      #Merge these two loops?
      #Need to account for additional config keys in user/static config

	
      if (exists ${$fset_config}{'analyses'}) {
	  
        foreach my $logic_name (keys %{$fset_config->{'analyses'}}) {
          $fset_config->{'analyses'}{$logic_name} =
            $self->validate_and_store_analysis($logic_name, $fset_config->{'analyses'}{$logic_name});
        }
      }

      if (exists ${$fset_config}{'feature_types'}) {
	  
        foreach my $ftype_key (keys %{$fset_config->{'feature_types'}}) {
          $fset_config->{'feature_types'}{$ftype_key} = 
            $self->validate_and_store_feature_type($ftype_key, $fset_config->{'feature_types'}{$ftype_key});
        }
      }
	
      #if(exists ${$fset_config}{'feature_set'}){
      #  #Just ignore this for now, as set_feature_set and define_and_validate_sets deal with this repectively
      #  #Solution is to extend define_and_validate_set to support external feature static config
      #  #Then remove set_feature_sets
      # Also need to consider InputSet::define_sets
      #}
    }
  }
  return;
}


#Need some method to get and add analyses and ftypes to reduce hardcoding of key strings
#Can we add Class->new to config to provide some of this?
#Would proliferate constructor calls into config
#May also prevent some logic imposed by validate/set methods i.e. defaults?
#Has to be one or other as we can't assume we can use methods if we are dealing with a hash

#Should we allow empty hashes to default to DB entry?
#Move HASH test into these methods

#Need to use rearrange in these methods for case safety
#rearrange is quite slow, fine for import

#Currently can't use empty hash to default to DB, as this gives us no key in the config
#Can we set this to the string required for the fetch method in the referenced hash(analyses/feature_types)?
# eq test would work, but would have to skip HASH test

#Maybe we just set the analysis in the feature_set by the logic_name string, rather than a ref
#Then we test for presence in the relevant hash
#Would have to set obj and ref when validating/storing
#Would also have to do this in the feature_set feature_types/analyses i.e. we would be able to ref the whole hash
#would have to be an arrayrefs

#Do we want to support no entries in top level definition if we can find it in the DB?


#Could merge these two using a config hash to define adaptor method, params etc?
#Would be easy to add more types to store e.g. cell_type

sub validate_and_store_analysis{
  my ($self, $logic_name, $analysis_params) = @_; 
  my $analysis;
  eval {$self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis_params)};

  if (! $@) {                   #Can assume we have already set this valid Analysis
    $analysis = $analysis_params
  } else {
    my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
	
    ### Validate config entries
    #Catches undef or empty {} at feature_sets and top levels, defaults to DB
    my $invalid_entry = 1;
    my $obj_logic_name;
    my $got_config = 0;


    if ( (! defined $analysis_params ) ||
         (ref($analysis_params) eq 'HASH') ) { #We have a valid entry
	
      if (  (ref($analysis_params) eq 'HASH') && 
            (%{$analysis_params}) ) { #number of keys
        $got_config = 1;
        $invalid_entry = 0;
      } else {                  #empty feature_set analyses hash -> check top level analyses first
        #$invalid_entry = 1;

        if (exists ${$self->{static_config}{analyses}}{$logic_name}) {
          $analysis_params = $self->{static_config}{analyses}{$logic_name};
		
          if ( (! defined $analysis_params ) ||
               (ref($analysis_params) eq 'HASH') ) {
            $invalid_entry = 0;
		  
            if ( (ref($analysis_params) eq 'HASH') &&
                 (%{$analysis_params}) ) { #number of keys
              $got_config = 1;
            }                   #else is empty top level {}

          }                     #else is invalid
        } else {                #No top level config, assume we want to use the DB
          $invalid_entry = 0;
        }
      }
    }
  
  
    if ($invalid_entry) {
      throw("You have defined a none HASH value in your config for analysis:\t$logic_name\n".
            "Please define config as HASH, or use empty HASH or undef to use existing default config or DB");
    }
	
    if ($got_config) {
      ($obj_logic_name) = rearrange(['LOGIC_NAME'], %{$analysis_params});
    } else {
      $obj_logic_name   = $logic_name;
    }
  	
  
    if ($logic_name ne $obj_logic_name) { #Not a show stopper as this is just the config key
      warn "Found analysis key name - logic_name mismatch in config:\t$logic_name vs $obj_logic_name\n";
    }
 
    $analysis         = $analysis_adaptor->fetch_by_logic_name($obj_logic_name);
  
  

    if ($got_config) {

      my $config_anal      = Bio::EnsEMBL::Analysis->new(%{$analysis_params});
	
      if (! defined $analysis) {	
        $self->log('Analysis '.$obj_logic_name." not found in DB, storing from config");		
        $analysis_adaptor->store($config_anal);
        $analysis = $analysis_adaptor->fetch_by_logic_name($obj_logic_name);	
      } else {
	  
        my $not_same = $analysis->compare($config_anal);
        #Analysis::compare returns the opposite of what you expect!
	  
        if ($not_same) {
          throw('There is a param mismatch between the '.$obj_logic_name.
                ' Analysis in the DB and config. Please rectify or define a new logic_name');
        }
      }
    } elsif (! defined $analysis) {
      throw("Cannot fetch $obj_logic_name analysis from DB, please check your config key or define new top level analyses config");
    }
  }

  return $analysis;
}

sub validate_and_store_feature_type{
  my ($self, $ftype_name, $ftype_params) = @_;
  my $ftype;
  eval {$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype_params)};
  
  if (! $@) {                   #Can assume we have already set this valid Analysis
    $ftype = $ftype_params
  } else {
    my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
	
    #Validate config entries
    #Catches undef or empty {} at feature_sets and top levels, defaults to DB
    my $invalid_entry = 1;
    my ($name, $class);
  
    my $got_config = 0;

    if ( (! defined $ftype_params ) ||
         (ref($ftype_params) eq 'HASH') ) { #We have a valid entry
	
      if (  (ref($ftype_params) eq 'HASH') && 
            (%{$ftype_params}) ) { #number of keys
        $got_config = 1;
        $invalid_entry = 0;
      } else {                  #empty feature_set ftype hash -> check top level ftype first
        #$invalid_entry = 1;
	  
        if (exists ${$self->{static_config}{feature_types}}{$ftype_name}) {
          $ftype_params = $self->{static_config}{feature_types}{$ftype_name};
		
          if ( (! defined $ftype_params ) ||
               (ref($ftype_params) eq 'HASH') ) {
            $invalid_entry = 0;
		  
            if ( (ref($ftype_params) eq 'HASH') &&
                 (%{$ftype_params}) ) { #number of keys
              $got_config = 1;
            }                   #else is empty top level {}

          }                     #else is invalid
        } else {                #No top level config, assume we want to use the DB
          $invalid_entry = 0;
        }
      }
    }
  
 
    if ($invalid_entry) {
      throw("You have defined a none HASH value in your config for feature_type:\t$ftype_name\n".
            "Please define config as HASH, or use empty HASH or undef to use existing default config or DB");
    }
	
    if ($got_config) {
      ($name, $class) = rearrange(['NAME', 'CLASS'], %{$ftype_params});
    } else {
      $name   = $ftype_name;
    }
  	

   
    #Can't use rearrange for key we are setting as we need to now the case
    my $analysis_key = (exists ${$ftype_params}{-analysis}) ? '-analysis' : '-ANALYSIS';
    my $analysis;

    if (exists ${$ftype_params}{$analysis_key}) {
      #This is slightly redundant as we may have already validated this analysis
      my ($lname) = rearrange(['LOGIC_NAME'], %{$ftype_params->{$analysis_key}});
      $ftype_params->{$analysis_key} = $self->validate_and_store_analysis($lname, $ftype_params->{$analysis_key});
      $analysis = $ftype_params->{$analysis_key};
    }
  
    my @ftypes = $ftype_adaptor->fetch_by_name($name, $class, $analysis);
   
    if (scalar(@ftypes) > 1) {
      throw("Unable to fetch unique feature_type $name. Please specify top level config to define class (and analysis)");
    } else {
      $ftype = $ftypes[0];      #can be undef
    }
  
  
    if ($got_config) {
      my $config_ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(%{$ftype_params});

      if ($ftype) {
	  
        if (! $ftype->compare($config_ftype)) {
          my $label = $name."($class";
          $label .= (defined $analysis) ? ' '.$analysis->logic_name.')' : ')';
		
          throw('There is a param mismatch between the '.$name.
                ' FeatureType in the DB and config. Please rectify in the config.');
        }
      } else {
        $self->log('FeatureType '.$name." not found in DB, storing from config");		
        ($ftype) = @{$ftype_adaptor->store($config_ftype)};
      }
    } elsif (! defined $ftype) {
      throw("Cannot fetch $name feature_type from DB, please check your config key or define new top level feature_types config");
    }
  }

  return $ftype;
}





sub counts{
  my ($self, $count_type) = @_;

  if ($count_type) {
    $self->{'_counts'}{$count_type} ||=0;
    return 	$self->{'_counts'}{$count_type};
  }
 
  return $self->{'_counts'}
}



sub slices{
  my ($self, $slices) = @_;

  if (defined $slices) {
    
    if (ref($slices) ne 'ARRAY') {
      throw("-slices parameter must be an ARRAYREF of Bio::EnsEMBL::Slices (i.e. not $slices)");
    }

    foreach my $slice (@$slices) {
      
      if (! ($slice && ref($slice) && $slice->isa('Bio::EnsEMBL::Slice'))) {
        throw("-slices parameter must be Bio::EnsEMBL::Slices (i.e. not $slice)");
      }
      
      #Removed cache_slice from here as this was
      #preventing us from identifying the seq_region in an input file
      
      my $full_slice = $self->slice_adaptor->fetch_by_name($slice->name);

      if (($slice->start != 1) ||
          ($slice->end != $full_slice->end)) {
        throw("InputSet Parser does not yet accomodate partial Slice based import i.e. slice start > 1 or slice end < slice length:\t".$slice->name);
		
      }

      push @{$self->{seq_region_names}}, $slice->seq_region_name;
    }
    $self->{slices} = $slices;
  }

  return $self->{slices} || [];
}


sub count{
  my $self       = shift;
  my $count_type = shift;;

  $self->{'_counts'}{$count_type} ||=0;
  $self->{'_counts'}{$count_type}++;
  return $self->{'_counts'}{$count_type};
}


sub rollback{ return shift->{rollback}; }

sub recovery{ return shift->{recover}; }


=head2 db
  
  Example    : my $db = $imp->db;
  Description: Getter for the db attribute
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut

sub db{ return shift->{db}; }


=head2 species
  
  Example    : my $species = $imp->species;
  Description: Getter for species attribute
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut

sub species{ return shift->{species}; }


=head2 ucsc_coords
  
  Example    : $start += 1 if $self->ucsc_coords;
  Description: Getter for UCSC coordinate usage flag
  Returntype : Boolean
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub ucsc_coords{ return shift->{ucsc_coords}; }


=head2 dump_fasta
  
  Example    : if($self->dump_fasta()){...do fasta dump...}
  Description: Getter for the dump_fasta flag
  Returntype : Boolean
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub dump_fasta{ return shift->{_dump_fasta}; }

sub slice_adaptor{ return shift->{slice_adaptor}; }


=head2 name
  
  Example    : $imp->name('Experiment1');
  Description: Getter for the name of this import
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name{ return shift->{name}; }


=head2 feature_set_description
  
  Example    : $imp->description("ExperimentalSet description");
  Description: Getter for the FeatureSet description
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub feature_set_description{ return shift->{feature_set_desc}; }


=head2 project_feature

  Args [0]   : Bio::EnsEMBL::Feature
  Args [1]   : string - Assembly e.g. NCBI37
  Example    : $self->project($feature, $new_assembly);
  Description: Projects a feature to a new assembly via the AssemblyMapper
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Throws is type is not valid.
  Caller     : General
  Status     : At risk 

=cut

#Was in BaseExternalParser

# --------------------------------------------------------------------------------
# Project a feature from one slice to another
sub project_feature {
  my $self         = shift;
  my $feat         = shift;
  my $new_assembly = shift;
  my $feat_slice   = $feat->feature_Slice;  # project feature to new assembly

  if (! $feat_slice) {
    throw('Cannot get Feature Slice for '.$feat->start.':'.$feat->end.':'.$feat->strand.' on seq_region '.$feat->slice->name);
  }

  my @segments = @{ $feat_slice->project('chromosome', $new_assembly) };
  
  if (! @segments) {
    $self->log("Failed to project feature:\t".$feat->display_label);
    return;
  } elsif (scalar(@segments) >1) {
    $self->log("Failed to project feature to distinct location:\t".$feat->display_label);
    return;
  }
  
  my $proj_slice = $segments[0]->to_Slice;
  
  if ($feat_slice->length != $proj_slice->length) {
    $self->log("Failed to project feature to comparable length region:\t".$feat->display_label);
    return;
  }
  
  
  # everything looks fine, so adjust the coords of the feature
  $feat->start($proj_slice->start);
  $feat->end($proj_slice->end);
  $feat->strand($proj_slice->strand);
  my $slice_new_asm = $self->slice_adaptor->fetch_by_region('chromosome', $proj_slice->seq_region_name, undef, undef, undef, $new_assembly);
  $feat->slice($slice_new_asm);
  
  return $feat;

}



=head2 registry_host

  Example    : my $reg_host = $imp->registry_host;
  Description: Accessor for registry host attribute
  Returntype : string e.g. ensembldb.ensembl.org
  Exceptions : None
  Caller     : general
  Status     : at risk 

=cut

sub registry_host{
  return $_[0]->{'reg_host'};
}

=head2 registry_user

  Example    : my $reg_user = $imp->registry_user;
  Description: Accessor for registry user attribute
  Returntype : string e.g. ensembldb.ensembl.org
  Exceptions : None
  Caller     : general
  Status     : at risk 

=cut

sub registry_user{
  return $_[0]->{'reg_user'};
}

=head2 registry_port

  Example    : my $reg_port = $imp->registry_port;
  Description: Accessor for registry port attribute
  Returntype : string e.g. ensembldb.ensembl.org
  Exceptions : None
  Caller     : general
  Status     : at risk 

=cut

sub registry_port{
  return $_[0]->{'reg_port'};
}

=head2 registry_pass

  Example    : my $reg_pass = $imp->registry_pass;
  Description: Accessor for registry pass attribute
  Returntype : string e.g. ensembldb.ensembl.org
  Exceptions : None
  Caller     : general
  Status     : at risk 

=cut

sub registry_pass{
  return $_[0]->{'reg_pass'};
}


#init_methods to go here
#each will check for mandatory params
#and set an init_flag? Such that any dependant methods
#can check it to see if mandatory params have been checked
#alterantively split this out into separate classes?

sub check_base_mandatory_params {
  my $self = $_[0];
}



=head2 input_files
  
  Example    : $imp->input_files(\@files);
  Description: Getter for the input file paths
  Returntype : Arrayref
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub input_files{
  my ($self, $input_files) = @_;
  #coudl also be really kind and handle being passed an Array rather than an Arrayref
  
  
  if(defined $input_files){
    my $ref = ref($input_files);
    
    if($ref ne 'ARRAY'){
      throw('-input_files parameter is not an Arrayref:\t'.$input_files);
    }
    elsif($ref eq ''){  #Be kind to scalars too
      $input_files = [$input_files] 
    }
       
    map { validate_path($_) } @$input_files;
    $self->log("Processing files:\n\t\t".join("\n\t\t", @$input_files)); 
    
    #We don't yet support multiple files
    if (scalar(@$input_files) >1) {
      warn("Found more than one input file:\n".join("\n\t", @{$self->input_files}).
        "\nThe InputSet parser does not yet handle multiple input files(e.g. replicates).\n".
        "We need to resolve how we are going handle replicates with random cluster IDs");
      #do we even need to?
    }
    
    $self->{input_files} = $input_files;
  }

  return $self->{input_files}; 
}
## !!!!
## Dirty Fix for release 75
## !!!!!! Method copied from Importer.pm

=head2 get_dir

  Example    : $imp->get_dir("import");
  Description: Retrieves full path for given directory
  Arg [1]    : mandatory - dir name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk - move to Helper?

=cut

sub get_dir{
  my ($self, $dirname) = @_;
  return $self->get_data("${dirname}_dir");
}


1;
