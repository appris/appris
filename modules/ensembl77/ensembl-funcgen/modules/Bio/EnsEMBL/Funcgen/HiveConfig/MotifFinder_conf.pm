
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::MotifFinder_conf;

=head1 SYNOPSIS

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::*_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::*_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::*_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the Config file for the Motif Finder Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Motif Finder pipeline consists of several "analysis":
        * SetupMotifPipeline
        * InferSubMotifs
        * ClusterMotifs

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::HiveConfig::MotifFinder_conf;

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
	  		    -host   => $self->o('dbhost'),
	  		    -port   => $self->o('dbport'),
	  		    -user   => $self->o('pipeuser'),
	  		    -pass   => $self->o('pipepass'),
	  		    #-dbname => $ENV{USER}.'_peaks_'.$self->o('dbname'),
			    -dbname => $self->o('pipedb_name'),
	  		   },

	  'dnadb_host' => 'ens-livemirror',
	  'dnadb_port' => 3306,
	  "dnadb_user" => 'ensro',

	  'bin_dir' => "/software/ensembl/funcgen",

	  'bin_size'  => 500,
	  'window_size' => 50,

	 };
}

=head2 resource_classes

    Description : Implements resource_classes() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the LSF resource classes available

=cut

sub resource_classes {
    my ($self) = @_;
    return {
   'default'            => { 'LSF' => '' },
   'urgent'             => { 'LSF' => '-q yesterday' },
   'normal_monitored'   => { 'LSF' => "        -R\"select[$ENV{LSF_RESOURCE_HOST}<600 && $ENV{LSF_DNA_DB}<600] rusage[$ENV{LSF_DNA_DB}=10:$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
   'long_monitored'     => { 'LSF' => "-q long -R\"select[$ENV{LSF_RESOURCE_HOST}<1000]                        rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
   'long_high_memory'   => { 'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },

#	    0 => { -desc => 'default',          'LSF' => '' },
#	    1 => { -desc => 'urgent',           'LSF' => '-q yesterday' },
#	    2 => { -desc => 'normal ens-genomics1',  'LSF' => '-R"select[myens_genomics1<600 && myens_livemirror<600] rusage[myens_livemirror=10:myens_genomics1=10:duration=10:decay=1]"' },
#	    3 => { -desc => 'long ens-genomics1',    'LSF' => '-q long -R"select[myens_genomics1<1000] rusage[myens_genomics1=10:duration=10:decay=1]"' },
#	    4 => { -desc => 'long high memory',      'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },

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

	  'output_dir'      => $self->o('work_dir')."/motifs/results",
	  'hive_output_dir' => $self->o('work_dir')."/motifs/hive_output",

	  'dbhost' => $self->o('dbhost'),
	  'dbport' => $self->o('dbport'),
	  "dbuser" => $self->o("dbuser"),
	  "dbname" => $self->o("dbname"),

	  'dnadb_host' => $self->o('dnadb_host'),
	  'dnadb_port' => $self->o('dnadb_port'),
	  "dnadb_user" => $self->o("dnadb_user"),
	  "dnadb_name" => $self->o("dnadb_name"),

	  "species" => $self->o("species"),

	  'bin_dir' => $self->o('bin_dir'),

	  'bin_size' => $self->o('bin_size'),
	  'window_size' => $self->o('window_size'),
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

	  @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables creation
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 0)." -e 'CREATE DATABASE ".$self->o('pipeline_db', '-dbname')."'",
	  # standard eHive tables and procedures:
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_hive_root_dir').'/sql/tables.sql',
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_hive_root_dir').'/sql/procedures.sql',

	  #Create hive output folders as required
	  'mkdir -p '.$self->o('work_dir')."/motifs/results",
	  'mkdir -p '.$self->o('work_dir')."/motifs/hive_output",

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
	   -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::SetupMotifInference',
	   -parameters => { 'batch_job' => 0, 'prepared' => 0 },
	   -input_ids => [
			  # No initial input_ids... these will be added as needed by init_pipeline -job_topup
			  { 'feature_set' => $self->o('feature_set') },
			 ],
	   -flow_into => {
			  2 => [ 'infer_submotifs' ],
			  3 => [ 'cluster_motifs' ],
			 },
	   -hive_capacity => 10,
	   -rc_name => 'default',
	  },

	  {
	   #This basically consists on running a command...
	   -logic_name    => 'infer_submotifs',
	   -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::InferMotifs',
	   -parameters    => { },
	   -input_ids     => [
				 # (jobs for this analysis will be flown_into via branch-1 from 'setup_pipeline' jobs above)
			     ],
	   -hive_capacity => 100,
	   -rc_name => 'normal_monitored',
	   #use semaphores...
	   #-wait_for => [ 'setup_pipeline' ]
	  },

	  {
	   #This basically consists on running a command...
	   -logic_name => 'cluster_motifs',
	   -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ClusterMotifs',
	   -parameters => { },
	   -input_ids => [
	  		  # (jobs for this analysis will be flown_into via branch-2 from 'setup_pipeline' jobs above)
	  		 ],
	   -hive_capacity => 10,
	   -rc_name => 'default',
	   #Use semaphores...
	   #-wait_for => [ 'run_peaks_DNAse', 'run_peaks' ],
	  },
	 ];
}

1;

