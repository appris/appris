
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::Peaks_conf;

=head1 SYNOPSIS

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Peaks_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Peaks_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Peaks_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the Config file for the Peaks Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Peaks pipeline consists of several "analysis":
        * SetupPeaksPipeline verifies the existence of experiments etc...
        * RunSWEmbl make the peak calling and stores the annotated features...
        * WrapUpSWEmbl do some filtering when needed and QC

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::HiveConfig::Peaks_conf;

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
    %{$self->SUPER::default_options},  # inheriting database and hive tables creation

	  'pipeline_db' => {
	  		    -host   => $self->o('dbhost'),
	  		    -port   => $self->o('dbport'),
	  		    -user   => $self->o('dbuser'),
	  		    -pass   => $self->o('dbpass'),
	  		    #-dbname => $ENV{USER}.'_peaks_'.$self->o('dbname'),
			    -dbname => $self->o('pipedb_name'),
	  		   },

	  #'efgdb_host' => ...
	  #We can add default values for all these but tend to avoid since people quickly forget what those are...

	  'skip_control' => 0,
	  'control_feature' => 'WCE',
	  'bin_dir' => '/software/ensembl/funcgen',
	  'control_file' => ''

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

#	    1 => { -desc => 'urgent',                           'LSF' => '-q yesterday' },
#	    2 => { -desc => 'normal ens-genomics1',             'LSF' => '-R"select[myens_genomics1<1000] rusage[myens_genomics1=10:duration=10:decay=1]"' },
#	    3 => { -desc => 'long ens-genomics1',               'LSF' => '-q long -R"select[myens_genomics1<1000] rusage[myens_genomics1=10:duration=10:decay=1]"' },
#	    4 => { -desc => 'long high memory',                 'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
#	    5 => { -desc => 'long ens-genomics1 high memory',   'LSF' => '-q long -M4000000 -R"select[myens_genomics1<1000 && mem>4000] rusage[myens_genomics1=10:duration=10:decay=1:mem=4000]"' },
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
      %{$self->SUPER::pipeline_wide_parameters},  # inheriting database and hive tables creation
	    'pipeline_name' => $self->o('pipedb_name'),  # name used by the beekeeper to prefix job names on the farm


	    'work_dir'        => $self->o('work_dir'),            # data directories and filenames
	    #'output_dir'      => $self->o('work_dir').'/ehive/'.$self->o('efgdb_name').'/hive_results',
	    #'hive_output_dir' => $self->o('work_dir').'/ehive/'.$self->o('efgdb_name').'/hive_debug',
	    #'output_dir'      => $self->o('output_dir').'/ehive/'.$self->o('efgdb_name').'/hive_results',
	    #'hive_output_dir' => $self->o('output_dir').'/ehive/'.$self->o('efgdb_name').'/hive_debug',
	    'output_dir'      => $self->o('output_dir').'/peaks/results',
	    'hive_output_dir' => $self->o('output_dir').'/peaks/hive_debug',
	    'bin_dir'         => $self->o('bin_dir'),

	    #Maybe use parameters instead of ENV Variables... use ENV variables in the pipeline ENV
	    "dnadb"	   => {
			       "-host"   => $self->o('dnadb_host'),
			       "-port"   => $self->o('dnadb_port'),
			       "-user"   => $self->o('dnadb_user'),
			       "-dbname" => $self->o('dnadb_name'),
			      },
	    "efgdb"	   => {
			       "-host"   => $self->o('dbhost'),
			       "-port"   => $self->o('dbport'),
			       "-user"   => $self->o('dbuser'),
			       "-pass"   => $self->o('dbpass'),
			       "-dbname" => $self->o('dbname'),
			      },
	    #This could be inferred from the db, but it's probably safer(?) to pass as parameter...
	    "species"      => $self->o('species'),
	    #May pass this to input_id... to allow for files of different assemblies in the same pipeline run.
	    "assembly"     => $self->o('assembly'),

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
	  #HiveGeneric assumes ensembl-hive folder while if you use the stable version it's ensembl-hive_stable!
	  @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables creation

	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 0)." -e 'CREATE DATABASE ".$self->o('pipeline_db', '-dbname')."'",

	  # standard eHive tables and procedures:
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/tables.sql',
	  #'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/procedures.sql',

	  #Create hive output folders as required
	  'mkdir -p '.$self->o('output_dir').'/peaks/results',
	  'mkdir -p '.$self->o('output_dir').'/peaks/hive_debug',

	  #This command makes the init_pipeline script fail.
	  #Any value in meta can be added with pipeline_wide_parameters()
	  #'mysql '.$self->dbconn_2_mysql($pipeline_db, 1).' -e "INSERT INTO meta (meta_key, meta_value) VALUES (\'hive_output_dir\',\''.$self->o('work_dir').'/ehive/'.$self->o('efgdb_name').'/hive_debug\');"',

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
	   -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::SetupPeaksPipeline',
	   -parameters => {},
	   -input_ids => [
			  # No initial input_ids... these will be added as needed by init_pipeline -job_topup
			  { 'cell_type' => $self->o('cell_type'), 'feature_type' => $self->o('feature_type'), 'experiment_name' => $self->o('experiment_name'), 'file_type' => $self->o('file_type'), 'analysis' => $self->o('analysis_name'), 'skip_control' => $self->o('skip_control'),  'control_feature' => $self->o('control_feature'), 'control_file' => $self->o('control_file') },
# , 'control_file' =>$self->o('cell_type')."_".$self->o('control_feature')."_".$self->o('experiment_name').".samse.".$self->o('file_type').".gz" },
			 ],
	   -flow_into => {
			  3 => [ 'run_peaks' ],
			  4 => [ 'run_peaks_wide' ],
			  5 => [ 'run_macs' ],
			  #2 => [ 'wrap_up_pipeline' ],
			 },
	   #These jobs cannot run in parallel due to race conditions! Do NOT change this setting unless you know what you're doing
	   -hive_capacity => 1,
           -rc_name => 'default',
	  },

	  {
	   -logic_name    => 'run_peaks',
	   -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::RunSWEmbl',
	   -parameters    => { },
	   -input_ids     => [
				 # (jobs for this analysis will be flown_into via branch-3 from 'setup_pipeline' jobs above)
			     ],
	   -hive_capacity => 10,
	   #Control files should be handled by setup_pipeline.
	   -rc_name => 'long_monitored_high_memory', # Better safe than sorry... size of datasets tends to increase...
	   -wait_for => [ 'setup_pipeline' ]
	  },

	  {
	   -logic_name    => 'run_peaks_wide',
	   -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::RunCCAT',
	   -parameters    => { },
	   -input_ids     => [
			      # (jobs for this analysis will be flown_into via branch-3 from 'setup_pipeline' jobs above)
			     ],
	   -hive_capacity => 10,
	   #Control files should be handled by setup_pipeline.
	   -rc_name => 'normal_monitored', # CCAT does not need much
	   -wait_for => [ 'setup_pipeline' ]
	  },

    {
	   -logic_name    => 'run_macs',
	   -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::RunMACS',
	   -parameters    => { },
	   -input_ids     => [
				 # (jobs for this analysis will be flown_into via branch-5 from 'setup_pipeline' jobs above)
			     ],
	   -hive_capacity => 10,
	   #Control files should be handled by setup_pipeline.
	   -rc_id => 5, # Better safe than sorry... size of datasets tends to increase...
	   -wait_for => [ 'setup_pipeline' ]
	  },

	  #{
	  # -logic_name => 'wrap_up_pipeline',
	  # -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::WrapUpPeaksPipeline',
	  # -parameters => {},
	  # -input_ids => [
	  #		  # (jobs for this analysis will be flown_into via branch-1 from 'setup_pipeline' jobs above)
	  #		 ],
	  # #TODO see if it can run in paralell, usually we should be able to
	  # -hive_capacity => 10,
	  # -wait_for => [ 'run_peaks_DNAse', 'run_peaks' ],
	  #},
	 ];
}

1;

