
=pod

=head1 NAME

  Bio::EnsEMBL::Hive::PipeConfig::Dnase_profile_conf

=head1 SYNOPSIS


=head1 DESCRIPTION

    This is an example pipeline put together from basic building blocks:

    Analysis_1: JobFactory.pm is used to turn the list of tables of the given database into jobs

        these jobs are sent down the branch #2 into the second analysis

    Analysis_2: SystemCmd.pm is used to run these dumping+compression jobs in parallel.

=head1 CONTACT

Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut

package Bio::EnsEMBL::Funcgen::HiveConfig::Dnase_profile_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly

=head2 default_options

    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                   o('password')       your read-write password for creation and maintenance of the hive database
=cut

sub default_options {
    my ($self) = @_;
    return {
        'ensembl_cvs_root_dir' => $ENV{'SRC'},     # some Compara developers might prefer
	'pipeline_name' => 'dnase_profile',                 # name used by the beekeeper to prefix job names on the farm

        'pipeline_db' => {                                  # connection parameters
            -host   => 'ens-genomics1',
            -port   => 3306,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),                     # a rule where a previously undefined parameter is used (which makes either of them obligatory)
            -dbname => $ENV{USER}.'_'.$self->o('pipeline_name'),    # a rule where a previously defined parameter is used (which makes both of them optional)
	    },

	 'is_male'       => 0,                                          # include table creation statement before inserting the data

	 'work_dir'   => '/lustre/scratch101/ensembl/ds19/Dnase_Footprint_ENCODE',                                        #
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
	    'work_dir'        => $self->o('work_dir'),            # data directories and filenames
	    'hive_output_dir' => $self->o('work_dir').'/hive_debug',
    };
}

=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates a directory for storing the output.

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
	    'mysql '.$self->dbconn_2_mysql('pipeline_db', 0)." -e 'CREATE DATABASE ".$self->o('pipeline_db', '-dbname')."'",

	    # standard eHive tables and procedures:
	    'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/tables.sql',
	    'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).' <'.$self->o('ensembl_cvs_root_dir').'/ensembl-hive/sql/procedures.sql',

	    'mkdir -p '.$self->o('work_dir').'/hive_debug',
	   ];
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
      'normal_high_memory' => { 'LSF' => '-M8000000  -R"select[mem>8000]  rusage[mem=8000]"' },
      'normal_huge_memory' => { 'LSF' => '-M12000000 -R"select[mem>12000] rusage[mem=12000]"' },
#	    0 => { -desc => 'default',          'LSF' => '' },
#	    1 => { -desc => 'urgent',           'LSF' => '-q yesterday' },
#	    2 => { -desc => 'normal high memory',      'LSF' => '-M8000000 -R"select[mem>8000] rusage[mem=8000]"' },
#	    3 => { -desc => 'normal huge memory',      'LSF' => '-M12000000 -R"select[mem>12000] rusage[mem=12000]"' },
	   };
}

=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
 Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.

=cut

sub pipeline_analyses {
    my ($self) = @_;
    return [
	    {   -logic_name => 'make_profile',
		-module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MakeDnaseProfile',
		-parameters => {
			       },
		-input_ids => [
			       {
				'matrix'  => $self->o('matrix'),
				'dnase'   => $self->o('dnase'),
				'is_male' => $self->o('is_male'),
			       },
			      ],
		-hive_capacity => 100,       # allow several workers to perform identical tasks in parallel
		-rc_name => 'normal_high_memory',
		-flow_into => {
			       2 => [ 'run_centipede' ],   # will create a fan of jobs
			      },
	    },

	    {   -logic_name    => 'run_centipede',
		-module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::RunCentipede',
		-parameters    => {
				  },
		-hive_capacity => 100,       # allow several workers to perform identical tasks in parallel
		-rc_name => 'normal_high_memory',
		-input_ids     => [
				   # (jobs for this analysis will be flown_into via branch-2 from 'get_tables' jobs above)
				  ],
	    },
	   ];
  }

1;

