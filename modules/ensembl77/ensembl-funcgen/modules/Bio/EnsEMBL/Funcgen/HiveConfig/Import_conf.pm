
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf;

=head1 SYNOPSIS

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the Config file for the Import Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Import pipeline consists of several "analysis":
        * SetupImportPipeline is equivalent to the "prepare" in parse_and_import.pl
        * RunImport loads the reads per each slice...
        * WrapUpImport finalizes when all partial imports are done...

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut

sub default_options {
  my ($self) = @_;
  return {
	  'ensembl_cvs_root_dir' => $ENV{'SRC'},

	  'pipeline_db' => {
	  		    -host   => $self->o('host'),
	  		    -port   => $self->o('port'),
	  		    -user   => $self->o('user'),
	  		    -pass   => $self->o('pass'),
	  		    #-dbname => $ENV{USER}.'_peaks_'.$self->o('dbname'),
			    -dbname => $self->o('pipedb_name'),
	  		   },

	  'verbose'  => 0,

	  'format'   => "SEQUENCING",
	  'vendor'   => "SEQUENCING",
	  'parser'   => "Bed", #Bed
	  'group'    => "efg",
	  'location' => "Hinxton",
	  'contact'  => 'http://lists.ensembl.org/mailman/listinfo/dev',

	  'input_feature_class' => 'result',
	  'registry_host'       => 'ens-livemirror',
	  'registry_port'       => 3306,
	  'registry_user'       => 'ensro',
	  #'assembly' => 37,

	  'host'              => 'ens-genomics1',
	  'port'              => 3306,
	  'feature_analysis'  => 'bwa_samse',
	  'recover'           => 1,

	  'data_dir'   =>'/lustre/scratch109/ensembl/funcgen',
	  'output_dir' =>'/lustre/scratch109/ensembl/funcgen/output/'.$self->o('dbname'),

	  #you can only add a slice at a time... see if list of values can be passed as params...
	  #'slice' => undef,

	  #why is this needed??
	  'farm' => 1,


	 };
}

=head2 resource_classes

    Description : Implements resource_classes() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the LSF resource classes available

=cut

sub resource_classes {
    my ($self) = @_;
    return {
        'default'                    => { 'LSF' => '' },
        'urgent'                     => { 'LSF' => '-q yesterday' },
        'normal_monitored'           => { 'LSF' => "        -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
        'long_monitored'             => { 'LSF' => "-q long -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
        'long_high_memory'           => { 'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
        'long_monitored_high_memory' => { 'LSF' => "-q long -M4000000 -R\"select[$ENV{LSF_RESOURCE_HOST}<1000 && mem>4000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1,mem=4000]\"" },

#	    0 => { -desc => 'default',                         'LSF' => '' },
#	    1 => { -desc => 'urgent',                          'LSF' => '-q yesterday' },
#	    2 => { -desc => 'normal ens-genomics1',            'LSF' => '-R"select[myens_genomics1<1000] rusage[myens_genomics1=10:duration=10:decay=1]"' },
#	    3 => { -desc => 'long ens-genomics1',              'LSF' => '-q long -R"select[myens_genomics1<1000] rusage[myens_genomics1=10:duration=10:decay=1]"' },
#	    4 => { -desc => 'long high memory',                'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
#	    5 => { -desc => 'long ens-genomics1 high memory',  'LSF' => '-q long -M6000000 -R"select[myens_genomics1<600 && mem>6000] rusage[myens_genomics1=12:duration=5:decay=1,mem=6000]"' },
	   };
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {

	  'pipeline_name'   => $self->o('pipedb_name'),  # name used by the beekeeper to prefix job names on the farm
	  #Add a parameter $self->o('output_dir')
	  'output_dir'      => $self->o('output_dir')."/".$self->o('vendor'),
	  'hive_output_dir' => $self->o('output_dir')."/".$self->o('vendor')."/hive_debug",
	  "species"         => $self->o('species'),


	  "dbname"          => $self->o("dbname"),
	  "user"            => $self->o("user"),
	  "pass"            => $self->o("pass"),
	 'host' => $self->o('host'),
	  'port' => $self->o('port'),


#Added DNADB parameters
      "dnadb_host"   => $self->o('dnadb_host'),
	  "dnadb_port"   => $self->o('dnadb_port'),
	  "dnadb_user"   => $self->o('dnadb_user'),
	  "dnadb_name"   => $self->o('dnadb_name'),


	  "input_dir"       => $self->o('input_dir'),

	  'verbose'  => $self->o('verbose'),

	  'format'   => $self->o('format'),
	  'vendor'   => $self->o('vendor'),
	  'parser'   => $self->o('parser'),
	  'group'    => $self->o('group'),
	  'location' => $self->o('location'),
	  'contact'  => $self->o('contact'),

	  'input_feature_class' => $self->o('input_feature_class'),
	  'registry_host' => $self->o('registry_host'),
	  'registry_port' => $self->o('registry_port'),
	  'registry_user' => $self->o('registry_user'),
		  'registry_version' => $self->o('registry_version'),
	  'assembly' => $self->o('assembly'),


	  'feature_analysis' => $self->o('feature_analysis'),
	  'recover' => $self->o('recover'),

	  'data_dir' => $self->o('data_dir'),

	  #you can only add a slice at a time... see if list of values can be passed as params...
	  #'slice' => undef,

	  #why is this needed??
	  'farm' => $self->o('farm'),


	 };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands
      that will create and set up the Hive database.

=cut


sub pipeline_create_commands {
 my ($self) = @_;


  return [

	  #HiveGeneric assumes ensembl-hive folder while if you use the stable version its ensembl-hive_stable!
	  @{$self->SUPER::pipeline_create_commands},
	  # inheriting database and hive tables creation

	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 0)." -e 'CREATE DATABASE ".$self->o('pipeline_db', '-dbname')."'",

	  # standard eHive tables and procedures:
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/tables.sql',
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/procedures.sql',

	  #Create hive output folders as required
	  'mkdir -p '.$self->o('output_dir')."/".$self->o("vendor"),
	  'mkdir -p '.$self->o('output_dir')."/".$self->o('vendor')."/hive_debug",

	 ];
}


=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.


=cut

sub pipeline_analyses {
  my ($self) = @_;

  return [
	  {

	   -logic_name => 'setup_pipeline',
	   -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::Import',
	   -parameters => { 'batch_job' => 0, 'prepared' => 0 },
	   -input_ids => [
			  # No initial input_ids... these will be added as needed by init_pipeline -job_topup
			  { 'cell_type' => $self->o('cell_type'), 'feature_type' => $self->o('feature_type'), 'input_set' => $self->o('input_set'), 'result_file' => $self->o('result_file')  },
			  #allow slice(s) for partial import...
			  #{ 'cell_type' => $self->o('cell_type'), 'feature_type' => $self->o('feature_type'), 'input_set' => $self->o('input_set'), 'result_file' => $self->o('result_file'), 'slice' => $self->o('slice')  },
			 ],
	   -flow_into => {
			  1 => [ 'import_reads' ],
			  2 => [ 'wrap_up_pipeline' ],
			 },
	   -hive_capacity => 10,
           -rc_name => 'default',
#	   -rc_id => 0,
	   #this really need revising as this is sorting the bed files
	   #Need to change resource to reserve tmp space
	  },

	  {
	   -logic_name    => 'import_reads',
	   -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::Import',
	   -parameters    => { 'batch_job' => 1, 'prepared' => 1 },
	   -input_ids     => [
				 # (jobs for this analysis will be flown_into via branch-1 from 'setup_pipeline' jobs above)
			     ],
	   -hive_capacity => 50,
	   #Control files should be handled by setup_pipeline.
           -rc_name => 'long_monitored_high_memory',
            #-rc_id => 5, # Better safe than sorry... size of datasets tends to increase...
	   #use semaphores...
	   #-wait_for => [ 'setup_pipeline' ]
	  },

	  {
	   -logic_name => 'wrap_up_pipeline',
	   -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::Import',
	   -parameters => { 'wrap_up' => 1 },
	   -input_ids => [
	  		  # (jobs for this analysis will be flown_into via branch-2 from 'setup_pipeline' jobs above)
	  		 ],
	   -hive_capacity => 10,
           -rc_name => 'default',
#	   -rc_id => 0,
	   #Use semaphores...
	   #-wait_for => [ 'run_peaks_DNAse', 'run_peaks' ],
	  },
	 ];
}

1;

