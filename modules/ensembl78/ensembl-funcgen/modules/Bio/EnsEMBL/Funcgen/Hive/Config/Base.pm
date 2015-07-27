
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

    Bio::EnsEMBL::Funcgen::Hive::Config::Base;

=head1 SYNOPSIS

   # TODO: this could be easily merged with the Peaks pipeline...
   # TODO: allow subfolders which will represent replicates...
   # Allow semaphores so jobs can be run truly in parallel (see SemaStart and SemaLongMult_conf)

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the base config file for all funcgen pipeline config

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    Please see the implementation details in 'Runnable/DB' modules themselves.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut

package Bio::EnsEMBL::Funcgen::Hive::Config::Base;

use strict;
use warnings;
use base qw(Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf);
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of 
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used
                  to initialize default options.
    Caller      : DependantOptions::process_options (init_pipeline.pl)

=cut

#These are used by init_pipeline.pl and according to the docs should the pipelinedb spec
#Anything set (e.g. test) in default_options will only be available where it is reused(substituted)
# into the rest of the config below
#cmdline parameters over-ride these defaults
#This is basically mechanism to remove the needs for doing soemthing like this f
#or defaults in pipeline_wide_params: test => $self->o('test') || default_value
#or: test => (defined $self->o('test')) ? $self->o('test').'/ing' : '/default/ing'
#$self->o also captures undefined cmdline params, if they have not been set in defaults


#It looks like it is possible to refer to a previously defined option within the same hash

#Never assign a command line opt to it's namesake param key in default_options e.g.
# pipeline_name => $self->o('pipeline_name')
#This will create a circular dependancy!

#pipeline_db is special and is the only thing that will actually be used in the init_pipeline
#would need to define in reset of config if we want to access it as a param


#pdb to mysql args string method in long mutl example
#$self->db_execute_command('pipeline_db', $SQL)
#only works in piepline_create_Commands
#also only available in conf, flag here is to 
#append dbname, then it spearately append dbname usign ->o?
# 'cmd'           => 'mysqldump '.$self->dbconn_2_mysql('pipeline_db', 0).' '
#.$self->o('pipeline_db','-dbname').' >#work_dir#/#filename#',


### WARNING: Want to add data to an existing HASH?
# You may be tempted to do this:
#   my_hash => {%{$self->o(my_hash)}, new_key => new_value},  ### DO NOT DO THIS!
# This would work with a normal hash, but the 'o' method is not a simple hash accessor
# and will result in a circular reference!
#
# It is possible to redefine just one value in a hash, but it has to exist already and 
# reference another config value e.g.
#   my_hash => {some_key = $self->o('some_value')}
# You can now redefine some_value, and this will be sibstituted post-hoc into my_hash e.g
#   some_value => new_value,
# So undef values could be added for later definition

#### Analysis top up config behaviour
# Overlapping pipeline wide paramters will be over-written by the 
# config which is topped up last.
# Conversly analysis parameters are maintained as per the original config.
# Hence parameters for overlapping analyses should be identical, or
# at least do the required thing dependant on the order in which the configs
# are topped up.


#To do
# 1 Move species & assembly to BaseDB?

sub default_options {
  my $self = $_[0];  
  
  return {
    %{$self->SUPER::default_options},    
    #ensembl_cvs_root_dir => $ENV{'SRC'},   #Now set in env as is mandatory  requirement for loading default config 
    
    
    bin_dir              => undef,#'/software/ensembl/funcgen',
    #This is preventing things being picked up from $PATH
    #bin_dir should really only be defined if all programs are located in the same dir
    #or an an analysis level, to override $PATH or specify something which isn't in $PATH
    #this should probably only be defined at the analysis level
    #There is no way the user will know which Runnables use bin_dir
    #Best thing to do is to probably expect this to be set in the environment
    #and set this to undef here by default, so this can be used explicitly to 
    #run with a custom bin_dir environment
    
    
    pdb_port             => undef,
    #no_write             => undef, #For use with runWorker.pl -no_write, so we can write some STDOUT in run
    
       
    #Will need to define this below if we want access to it as a param
    pipeline_db => 
	  {
     #  test => $self->o('GRR'),
	    #CR enable different db params from output db
	   -host   => $self->o('pdb_host'),
	   -port   => $self->o('pdb_port'),
	   -user   => $self->o('pdb_user'),
	   -pass   => $self->o('ENV', 'PDB_PASS'),#Via env for security
	   -dbname => $self->o('pipeline_name'),
	   -driver => 'mysql',
	   #todo deal with this in the env and don't corrupt pipeline_name add $ENV{USER}?
	  },
    
	use_tracking_db     => 1,     
	species             => undef,
	
	#This does not work?!!
	#Other confs still have uc species
	#do it in the env instead
	#Set to lc here, so we don't have to do this for all other confs
	#species             => ($self->o('species')) ? lc($self->o('species')) : undef,#This can be inferred/fetched from the db
    
	assembly            => undef,
	

	#Keep these separate so we don't have to mess around with tidy up when archiving
	#These are also defined in the modules as a fail safe
	#but keep here so they og into pipeline wide params, rather than data flowing
	#for each job
	work_root_dir     => $self->o('data_root_dir').'/output/'.$self->o('pipeline_name'),
	hive_output_dir   => $self->o('data_root_dir').'/output/'.$self->o('pipeline_name').'/hive_debug',
  alt_data_root_dir => undef,

    #data_root_dir is omited, so it is made mandatory when loading the config
    #and we don't want to store a specific path in here
	
	#CR investigate usage of FileSpec
	
	
	#Mandatory DB params will now be caught by constructors via _set_out_db, 
	#rather than in config, as having this optional for Base means we can't 
	#enforce it for BaseDB :/
    
    
	#Maybe we shoudl merge after all, what pipeline won't use the DB?
	#alignment may need it for dumping seq
	
	
	#undef_params =>
	

    #Other HiveGeneric_conf vars may be set to thing we don't expect
    #Redefine here so that we don't use them unknowingly?
    #password, host?
    #How is password not throwing and error as we don't have ENSADMIN_PSW set in the env?

    verbose => undef,
    #no_tidy => undef, #todo implement no_tidy
    #Set this to stop the pipeline removing intermediates
    #which maybe required for rerunning certain jobs
    #leaving this unset greatly reduces the running footprint of the pipeline
    #This is useful for debugging
    
    archive_root     => undef,
    allow_no_archive => 0,
    #Set this to 1 to allow analyses to skip
    #archiving if archive_root is not specified
    
   };
}

=head2 resource_classes

    Description : Implements resource_classes() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the LSF resource classes available

=cut


sub resource_classes {
  my $self = shift;
  return 
    {#todo add in lsf group spec here to top LSF warning output
     #todo pass DB_HOST_LSFNAME as param
     
     default                 => { 'LSF' => '' },    
     #urgent                  => { 'LSF' => '-q yesterday' },
     #Should never use this in the pipleine, best to bswitch after submission if required
     normal_2GB              => { 'LSF' => ' -M2000 -R"select[mem>2000] rusage[mem=2000]"' },
     normal_monitored        => { 'LSF' => " -R\"select[$ENV{DB_HOST_LSFNAME}<1000] ".
                                            "rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1]\"" },
     normal_high_mem         => { 'LSF' => ' -M5000 -R"select[mem>5000] rusage[mem=5000]"' },
     normal_high_mem_2cpu    => { 'LSF' => ' -n2 -M5000 -R"select[mem>5000] rusage[mem=5000] span[hosts=1]"' },
     normal_monitored_2GB    => {'LSF' => " -M2000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>2000]".
                                                " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=2000]\"" },
     normal_monitored_4GB    => {'LSF' => " -M4000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>4000]".
                                            " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=4000]\"" },  
     normal_monitored_8GB    => {'LSF' => " -M8000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>8000]".
                                            " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=8000]\"" },   
     normal_monitored_16GB    => {'LSF' => " -M16000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>16000]".
                                            " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=16000]\"" },    
     normal_10gb_monitored    => {'LSF' => " -M10000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>10000]".
                                            " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=10000]\"" },                                                                         
     normal_5GB_2cpu_monitored => {'LSF' => " -n2 -M5000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>5000]".
                                                " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=5000] span[hosts=1]\"" },
     normal_10gb             => { 'LSF' => ' -M10000 -R"select[mem>10000] rusage[mem=10000]"' },
     long_monitored          => { 'LSF' => "-q long -R\"select[$ENV{DB_HOST_LSFNAME}<1000] ".
                                            "rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1]\"" },
     long_high_mem           => { 'LSF' => '-q long -M4000 -R"select[mem>4000] rusage[mem=4000]"' },
     long_monitored_high_mem => { 'LSF' => "-q long -M4000 -R\"select[$ENV{DB_HOST_LSFNAME}<1000 && mem>4000]".
                                                " rusage[$ENV{DB_HOST_LSFNAME}=10:duration=10:decay=1,mem=4000]\"" },
    };
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.
    Caller      : HiveGeneric_conf::run (init_pipeline.pl)

=cut

#$self->o here will take the cmdline option else default to what was set in default_options, else with fail if not defined

#Will init_pipeline with analysis_topup, use previously stored opts so we don't have to define them again?

#All these are written to meta table and available via param unless redefined for a given
#analysis (pipeline_analyses parameters stored in analysis_base_parameters), or 
#for a given job e.g. input_id


sub pipeline_wide_parameters {
  my ($self) = @_;
                            
  return {
    %{$self->SUPER::pipeline_wide_parameters},  # inheriting pipeline_name(not any more)
    pipeline_name => $self->o('pipeline_name'),
    
    
    
    #Will need to catch as mandatory params in fetch input where required
    species        => $self->o('species'),
    default_gender => 'male', #Used for defining reference files when gender is unkown
    #assembly is pipeline wide for the most part, but may want to load on two assemblies?
    #unlikely, more likely that we project from an old assembly onto a single default assembly
    
    assembly     => $self->o('assembly'),
    #Used for non DB stuff, and loading old data(projection) in DB analyses
    #this should use -old_assembly and -assembly
    
    
	#These can't be modified by each conf, as they will be in the same DB things will get messy
	#todo These should be used as root dirs and specific confs should define new keys i.e. 'confname_data_dir'
    #They can however be redefined for a given analysis! or job!
    

    #Allow over-writing defaults set above
       
    #No default for data dir, as we don't want to check in internal paths
    data_root_dir     => $self->o('data_root_dir'), #'/lustre/scratch109/ensembl/funcgen',
    work_root_dir     => $self->o('work_root_dir'), #'/lustre/scratch109/ensembl/funcgen',
    alt_data_root_dir => $self->o('alt_data_root_dir'), 
    
    
    archive_root     => $self->o('archive_root'),
    allow_no_archive => $self->o('allow_no_archive'),
    
       
    bin_dir           => $self->o('bin_dir'),
  
    #root_output_dir => $self->o('root_output_dir'),
    #need root_output_dir for non-DB specific analyses e.g. alignments
    #don't need this as we can define outputdir local to an analysis
    #in the parameters config, this will not be set in the meta table?
    #Would be useful to have a view on all the output data which is not dbspecific
	#can do this by having the following structure 
	#output/species/assembly/alignments/
	#output/species/assembly/peaks
	#output/species/assembly/dbname/alignments
	#output/species/assembly/dbname/peaks
	#Then createing the necessary links in the top level alignments/peaks dirs
	#from the dbspecific dirs
	#would need to handle duplicates
      
    #db_output_dir   => $self->o('db_output_dir'),
    #This will be defined as the output_dir for a given analysis
    
    hive_output_dir => $self->o('hive_output_dir'), 
    #enable work_root_dir in here?
       
    use_tracking_db => $self->o('use_tracking_db'),	
        
    #These are params which should be defined when seeding a batch of inputs
    #they will be automatically dataflowed to the appropriate analyses for
    #that batch, but not effect other batches i.e. they are batch_wide
    #not pipeline wide
    
    #We can still have pipeline_wide params for these, this just
    #ensures that any over-ride when seeding will be passed onto subsequent 
    #analyses 
    
    #These should always be parsed by all RunnableDBs at least?
    #Only necessary if they are to be used
    #generic data flow method will catch them in write_output
    
    batch_param_names => ['no_write', #For use with runWorker.pl -no_write, so we can write some STDOUT in run
                                 #is this already available in the job, or is it just passed ot the worker?
                    ],
    
    #This should really only ever be defined
    #when running in debug mode, but here for clarity.
    #This will not be passed as a param, so would either need to reseed the job
    #or pass no_tidy as the debug level
    #we could translate that to a real debug level via a hash in Base
    #not_tidy, no_tidy_2, no_tidy_3
    #This would work really nicely, but the hive code currently barfs as it expects a number
    #Will have to reseed job instead!
    #Currently only implemented in:
    #MergeQCAlignements - to maintain fastq and bam chunk files
    #no_tidy => $self->o('no_tidy'),  
  };
}



=head2 pipeline_create_commands

    Description : Interface method that should return a list of command lines to be run in order to create and set up the pipeline database.
                  Please see existing PipeConfig modules for examples.

=cut


#Copied from HiveGeneric_conf 
#commented out CREATE DB as we always have a pre-existing DB for funcgen work
#commented out foreign keys to enable DropPipeline (fails otherwise due tofk constraints
#re-enabled this as we are now using a separate DB for the hive tables

#default pipeline_create_commands defined in HiveGeneric_conf.pm
#specify more in sub-class conf if required but always do @{$self->SUPER::pipeline_create_commands},

#CR consider adding hive tables to output DB, therefore associating pipeline data with output data
#Need to consider overwriting/cleaning existing hive tables
#This will remove the possibility of having two hives run on the same DB. Is this a good thing?

## WARNING!!
## Currently init_pipeline.pl doesn't run this method when a pipeline is created with the -analysis_topup option
# configure_hive.pl handles this



1;

