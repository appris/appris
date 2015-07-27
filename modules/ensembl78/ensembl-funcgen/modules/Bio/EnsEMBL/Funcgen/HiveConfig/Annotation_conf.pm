
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::Annotation_conf;


=head1 DESCRIPTION

    This is the Config file for the Annotation Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Annotation pipeline implements Damian;s scripts in /scripts/regulatory_annotation:

    Please see the implementation details in the Runnable modules.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut

package Bio::EnsEMBL::Funcgen::HiveConfig::Annotation_conf;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
  my ($self) = @_;
  return {
	  'ensembl_cvs_root_dir' => $ENV{'SRC'},                  # some Compara developers might prefer $ENV{'HOME'}.'/ensembl_main'

	  'pipeline_db' => {
	  		    -host   => $self->o('dbhost'),
	  		    -port   => $self->o('dbport'),
	  		    -user   => $self->o('dbuser'),
	  		    -pass   => $self->o('dbpass'),
	  		    -dbname => $ENV{USER}.'_regfeat_annotation_'.$self->o('dbname'),
			    #-dbname => $self->o('pipedb_name'),
	  		   },

	 };
}

sub resource_classes {
    my ($self) = @_;
    return {
      'default'                    => { 'LSF' => '' },
      'urgent'                     => { 'LSF' => '-q yesterday' },
      'normal_monitored'           => { 'LSF' => "        -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
      'long_monitored'             => { 'LSF' => "-q long -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
      'long_high_memory'           => { 'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
      'long_monitored_high_memory' => { 'LSF' => "-q long -M4000000 -R\"select[$ENV{LSF_RESOURCE_HOST}<1000 && mem>4000] rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1,mem=4000]\"" },

#	    0 => { -desc => 'default',          'LSF' => '' },
#	    1 => { -desc => 'urgent',           'LSF' => '-q yesterday' },
#	    2 => { -desc => 'normal ens-genomics2',  'LSF' => '-R"select[myens_genomics2<1000] rusage[myens_genomics2=10:duration=10:decay=1]"' },
#	    3 => { -desc => 'long ens-genomics2',    'LSF' => '-q long -R"select[myens_genomics2<1000] rusage[myens_genomics2=10:duration=10:decay=1]"' },
#	    4 => { -desc => 'long high memory',      'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
#	    5 => { -desc => 'long ens-genomics2 high memory',  'LSF' => '-q long -M4000000 -R"select[myens_genomics2<1000 && mem>4000] rusage[myens_genomics2=10:duration=10:decay=1:mem=4000]"' },
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

	    'pipeline_name' => $self->o('pipeline_db', '-dbname'),  # name used by the beekeeper to prefix job names on the farm

	    'work_dir'        => $self->o('work_dir'),   # data directories and filenames
	    #use this as scratch dir or create one specifically?
	    'output_dir'      => $self->o('output_dir').'/annotation/results',
	    'hive_output_dir' => $self->o('output_dir').'/annotation/hive_debug',

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

	    "workdb"	   => {
			       "-host"   => $self->o('workdb_host'),
			       "-port"   => $self->o('workdb_port'),
			       "-user"   => $self->o('workdb_user'),
			       #workdbpass?
			       "-pass"   => $self->o('dbpass'),
			      },
	    #This could be inferred from the db, but it's probably safer(?) to pass as parameter...
	    "species"      => $self->o('species'),

	    #"release"      => $self->o('release'),


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
	  'mkdir -p '.$self->o('output_dir').'/annotation/results',
	  'mkdir -p '.$self->o('output_dir').'/annotation/hive_debug',


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
	   -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::SetupAnnotationPipeline',
	   -parameters => {},
	   -input_ids => [
			  # No initial input_ids... these will be added as needed by init_pipeline -job_topup
			  { },

			 ],
	   -flow_into => {
			  2 => [ 'annotate_regulatory_features' ],
			  #3 => [ 'wrap_up_pipeline' ],
			 },
	   #These jobs cannot run in parallel due to race conditions! Do NOT change this setting unless you know what you're doing
	   -hive_capacity => 1,
	   -rc_name => 'default',
	  },

	  {
	   -logic_name    => 'annotate_regulatory_features',
	   -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::AnnotateRegulatoryFeatures',
	   -parameters    => { },
	   -input_ids     => [
				 # (jobs for this analysis will be flown_into via branch-2 from 'setup_pipeline' jobs above)
			     ],
	   #Since all the weight is in the database it is safer to run only one at a time... or a small number at least
	   -hive_capacity => 1,
	   #Control files should be handled by setup_pipeline.
	   -rc_name => 'long_monitored_high_memory', # Better safe than sorry... size of datasets tends to increase...
	   -wait_for => [ 'setup_pipeline' ]
	  },

	 ];
}

1;
