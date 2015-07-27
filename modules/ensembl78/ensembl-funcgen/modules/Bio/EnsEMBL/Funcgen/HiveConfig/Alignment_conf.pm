
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf;

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

    This is the Config file for the Alignment Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Alignment pipeline consists of several "analysis":
        * SetupAlignmentPipeline verifies the existence of the files and creates alignment jobs ...
        * RunAlignment makes the alignment...
        * WrapUpAlignment merges the alignments, some QC and fills in the data tracking db

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf;

use strict;
use warnings;
use Data::Dumper;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut

sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options},

    'pipeline_db' => {
          -host   => $self->o('dbhost'),
          -port   => $self->o('dbport'),
          -user   => $self->o('dbuser'),
          -pass   => $self->o('dbpass'),
          #The aligments are independent of the EFG DB but since we should call the collections and the peaks, we can keep it
          #-dbname => $ENV{USER}.'_alignments_'.$self->o('efgdb_name'),
          -dbname => $self->o('pipedb_name'),
         },

    'bin_dir' => '/software/ensembl/funcgen',

    #'efgdb_host' => ...
    #We can add default values for all these but tend to avoid since people quickly forget what those are...

    #We could add a dummy default dataset, just so we can create an "empty" pipeline and add sets as needed...
    #This dummy dataset would have to be specifically detected and ignored in a setup step...
    #'experiment_name' => 'Dummy',
    #'cell_type' => 'Dummy',
    #'feature_type' => 'Dummy',

   };
}

=head2 resource_classes

    Description : Implements resource_classes() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the LSF resource classes available

=cut

sub resource_classes {
    my ($self) = @_;
    return {

#Use this section when running on Sanger Farm
   'default'            => { 'LSF' => '' },
   'urgent'             => { 'LSF' => '-q yesterday' },
   'long_high_memory'   => { 'LSF' => '-q long -M5000000 -R"select[mem>5000] rusage[mem=5000]"' },
   'normal_high_memory' => { 'LSF' => '        -M5000000 -R"select[mem>5000] rusage[mem=5000]"' },

#Use this section when running on EBI cluster
#    0 => { -desc => 'default',          'LSF' => '' },
#    1 => { -desc => 'long high memory',      'LSF' => '-M5000 -R"select[mem>5000] rusage[mem=5000]"' },
#    2 => { -desc => 'normal high memory',    'LSF' => '-M5000 -R"select[mem>5000] rusage[mem=5000]"' },

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
      'pipeline_name' => 'alignments'.'_'.$self->o('dbname'),  # name used by the beekeeper to prefix job names on the farm


      'work_dir'        => $self->o('work_dir'),            # data directories and filenames
      #This will be used for temp files and debug
      #'output_dir'      => $self->o('output_dir').'/ehive/'.$self->o('dbname').'/hive_results',
      #'hive_output_dir' => $self->o('output_dir').'/ehive/'.$self->o('dbname').'/hive_debug',
      'output_dir'      => $self->o('output_dir').'/alignments/results',
      'hive_output_dir' => $self->o('output_dir').'/alignments/hive_debug',


      #Maybe use parameters instead of ENV Variables... use ENV variables in the pipeline ENV
      "dnadb"   => {
             "-host"   => $self->o('dnadb_host'),
             "-port"   => $self->o('dnadb_port'),
             "-user"   => $self->o('dnadb_user'),
             "-dbname" => $self->o('dnadb_name'),
            },
      "efgdb"  => {
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

      #Make sure the bwa_indexes were generated with the same version!!!
      #Just use the default bwa: should be in in /software/varinfo/bin
      "bwa_bin"      => $self->o('bin_dir')."/bwa",
      #"bwa_bin"      => "/nfs/users/nfs_d/ds19/src/bwa-0.5.8a/bwa",
      #get new versions of bwa in /software/ensembl/bin/bwa
      'bin_dir' => $self->o('bin_dir'),

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      "fastq_chunk_size"      => "16000000" #This should run in 30min-1h

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
    #'mkdir -p '.$self->o('work_dir').'/ehive/'.$self->o('efgdb_name').'/hive_debug',
    #'mkdir -p '.$self->o('work_dir').'/ehive/'.$self->o('efgdb_name').'/hive_results',
    'mkdir -p '.$self->o('output_dir').'/alignments/results',
    'mkdir -p '.$self->o('output_dir').'/alignments/hive_debug',

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
     -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::SetupAlignmentPipeline',
     -parameters => {},
     -input_ids => [
        # No initial input_ids... these will be added as needed by init_pipeline -job_topup
        { 'cell_type' => $self->o('cell_type'), 'feature_type' => $self->o('feature_type'), 'experiment_name' => $self->o('experiment_name')  },
       ],
     -flow_into => {
        '1->A' => [ 'run_alignments' ],
        'A->2' => [ 'wrap_up_pipeline' ],
       },
     # These jobs can be run in parallell... don't put too many since it may generate many jobs...jobs
     -limit => 1,
     -batch_size => 1,
     -hive_capacity => 10,
     -rc_name => 'default',
    },

    {
     -logic_name    => 'run_alignments',
     -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::RunBWA',
     -parameters    => { },
     -input_ids     => [
         # (jobs for this analysis will be flown_into via branch-3 from 'setup_pipeline' jobs above)
           ],
     -hive_capacity => 100,
     # Better safe than sorry... size of datasets tends to increase...
     -rc_name => 'normal_high_memory',
     #No need to wait once since it is independent from everything else
     #-wait_for => [ 'setup_pipeline' ]
    },

    {
     -logic_name => 'wrap_up_pipeline',
     -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::WrapUpAlignment',
     -parameters => { },
     -input_ids  => [
          # (jobs for this analysis will be flown_into via branch-1 from 'setup_pipeline' jobs above)
         ],
     -hive_capacity => 10,
     -rc_name => 'long_high_memory',
     # No need to wait, if we use semaphores...
     # -wait_for => [ 'run_alignments' ],
    },

    #Temporary thing, while replicates are not handled properly...
    #{
    # -logic_name => 'converge_replicates',
    # -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ConvergeReplicates',
    # -parameters => { },
    # -input_ids  => [
          # (jobs for this analysis will be flown_into via branch-1 from 'setup_pipeline' jobs above)
    #jobs ],
    # -hive_capacity => 10,
    # -rc_id => 1,
     # No need to wait, if we use semaphores...
     # -wait_for => [ 'run_alignments' ],
    #},

   ];
}

1;

