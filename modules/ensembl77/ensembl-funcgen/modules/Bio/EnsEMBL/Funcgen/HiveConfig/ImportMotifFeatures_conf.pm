
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::ImportMotifFeatures_conf;

=head1 SYNOPSIS

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::*_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::*_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::*_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the Config file for the Import Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Import pipeline consists of several "analysis":
        * SetupPipeline is equivalent to the "prepare" in parse_and_import.pl
        * LoadMotifFeatures loads motif features per each slice...
        * WrapUpPipeline finalizes when all partial imports are done...

    Please see the implementation details in LoadMotifFeatures Runnable module

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::HiveConfig::ImportMotifFeatures_conf;

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
  return
    {
     %{$self->SUPER::default_options},    
    'hive_root_dir' => $ENV{'SRC'}.'/ensembl-hive',
     # some Compara developers might prefer $ENV{'HOME'}.'/ensembl_main'
    'pipeline_name' => $ENV{'PDB_NAME'},

    #'pipeline_db' =>
    # {
    #  -host   => $self->o('host'),
    #  -port   => $self->o('port'),
    #  -user   => $self->o('user'),
    #  -pass   => $self->o('pass'),
    #  -driver => 'mysql',
    #  -dbname => $ENV{'PDB_NAME'},
    # },

#Will need to define this below if we want access to it as a param
    pipeline_db => 
    {
     #  test => $self->o('GRR'),
      #CR enable different db params from output db
     -host   => $self->o('host'),
     -port   => $self->o('port'),
     -user   => $self->o('user'),
     -pass   => $self->o('pass'),
     -dbname => $self->o('ENV', 'PDB_NAME'),
     -driver => 'mysql',
     #todo deal with this in the env and don't corrupt pipeline_name add $ENV{USER}?
    },


     #Need to change this to use $ENV{OUT_ROOT} so we can switch scratch usage easily
     'output_dir' => undef,
     'slices'     => '',

	 };
}

=head2 resource_classes

    Description : Implements resource_classes() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the LSF resource classes available

=cut

sub resource_classes {
  my ($self) = @_;
  return
    {
     'default'                    => { 'LSF' => '' },
     'urgent'                     => { 'LSF' => '-q yesterday' },
     'normal_monitored'           => { 'LSF' => "-M1000 -R\"select[$ENV{LSF_RESOURCE_HOST}<1000 && mem>1000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1,mem=1000]\"" },
     'long_monitored'             => { 'LSF' => "-q long -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
     'long_high_memory'           => { 'LSF' => '-q long -M4000 -R"select[mem>4000] rusage[mem=4000]"' },
     'long_monitored_high_memory' => { 'LSF' => "-q long -M4000 -R\"select[$ENV{LSF_RESOURCE_HOST}<600 && mem>4000] rusage[$ENV{LSF_RESOURCE_HOST}=12:duration=5:decay=1,mem=4000]\"" },

#     0 => { -desc => 'default',          'LSF' => '' },
#     1 => { -desc => 'urgent',           'LSF' => '-q yesterday' },
#     2 => { -desc => 'normal ens-genomics1',  'LSF' => '-M1000000 -R"select[myens_genomics1<1000 && mem>1000] rusage[myens_genomics1=10:duration=10:decay=1:mem=1000]"' },
#     3 => { -desc => 'long ens-genomics1',    'LSF' => '-q long -R"select[myens_genomics1<1000] rusage[myens_genomics1=10:duration=10:decay=1]"' },
#     4 => { -desc => 'long high memory',      'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
#     5 => { -desc => 'long ens-genomics1 high memory',  'LSF' => '-q long -M4000000 -R"select[myens_genomics1<600 && mem>4000] rusage[myens_genomics1=12:duration=5:decay=1:mem=4000]"' },

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

	  'pipeline_name'   => $self->o('pipeline_name'),  # name used by the beekeeper to prefix job names on the farm
	  'hive_output_dir' => $self->o('output_dir')."/motif_features/hive_output",
	  'output_dir' => $self->o('output_dir')."/motif_features/results",

	  'host'   => $self->o('host'),
	  'port'   => $self->o('port'),
	  'user'   => $self->o('user'),
	  'pass'   => $self->o('pass'),
	  'dbname' => $self->o('dbname'),

	  'dnadb_host'  => $self->o('dnadb_host'),
	  'dnadb_port'  => $self->o('dnadb_port'),
	  'dnadb_user'  => $self->o('dnadb_user'),
	  'dnadb_name'  => $self->o('dnadb_name'),

	  'efg_src'    => $self->o('efg_src'),

	  'slices'     => $self->o('slices'),

	 };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands
      that will create and set up the Hive database.

=cut


sub pipeline_create_commands {
 my ($self) = @_;

  return
    [
     #HiveGeneric assumes ensembl-hive folder while if you use the stable version its ensembl-hive_stable!
     @{$self->SUPER::pipeline_create_commands},
     # inheriting database and hive tables creation rather than doing the following
     #'mysql '.$self->dbconn_2_mysql('pipeline_db', 0)." -e 'CREATE DATABASE ".$self->o('pipeline_db', '-dbname')."'",
     # standard eHive tables and procedures:
     #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/tables.sql',
     #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/procedures.mysql',


     #Create hive output folders as required
     'mkdir -p '.$self->o('output_dir')."/motif_features/hive_output",
     'mkdir -p '.$self->o('output_dir')."/motif_features/results",
    ];
}


=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.


=cut

sub pipeline_analyses {
  my ($self) = @_;

  return
    [
     {
      -logic_name    => 'run_import',
      -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::ImportMotifFeatures',
      -parameters    => { },
      -analysis_capacity => 50,   # allow several workers to perform identical tasks in parallel
      -batch_size    => 1,
      -input_ids     => [
                         #For the moment it only receives the matrix, and deduces feature_type(s) from there...
                         { 'matrix' => $self->o('matrix'), 'file' => $self->o('file') }
                        ],
      -rc_name => 'normal_monitored',
     },
    ];
}

1;

