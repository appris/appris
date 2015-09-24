
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::ReadAlignment;

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


package Bio::EnsEMBL::Funcgen::Hive::Config::ReadAlignment;

use strict;
use warnings;
use Data::Dumper;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut

sub default_options {
  my $self = shift;
  return {
    %{$self->SUPER::default_options},

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      fastq_chunk_size      => 16000000, #This should run in 30min-1h
      alignment_analysis    => 'bwa_samse',
      bwa_samse_param_methods     => ['sam_ref_fai'],

      #Use the same set up as for the peak caller wrt conditionaly data flowing


      #do we really need this? May aswell do it now as it costs nothing
      #although this would mean much more added complexity in DefineResultSets?
      #no, we just need to build the config correctly and the branch descriptor
      #This will be used in DefineResultSet and also Run_ALIGNER_and_QC_control analyses
      #alignment_branches =>
      # {
      #  bwa_control   => 2,
      #  bwa_merged    => 3,
      #  bwa_idr       => 4,

        #add more aligner triplets in here
        #custom support will not support >1 branch unless we alter branch_output_id
        #custom => 100,
      # },

      fastq_root_dir => $self->o('data_root_dir').'/fastq',

   };
}




=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
    my $self = shift;
    return
     {
      %{$self->SUPER::pipeline_wide_parameters},  # inheriting database and hive tables creation

      #Make sure the bwa_indexes were generated with the same version!!!
      #Just use the default bwa: should be in in /software/varinfo/bin
      #"bwa_bin"      => $self->o('bin_dir')."/bwa",
      #set this in the PeakCaller?

      #get new versions of bwa in /software/ensembl/bin/bwa? Hardened bwa?

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      fastq_chunk_size      => $self->o('fastq_chunk_size'),   #Change to batch specific
      alignment_analysis    => $self->o('alignment_analysis'), #Nope we may want this to be batch specific!
      #aligner_param_methods => $self->o($self->o('alignment_analysis').'_param_methods'),

      # this is currently preaking the option processing
      aligner_param_methods => $self->o('bwa_samse_param_methods'),

      #This is stricly not required anymore as we use the local_url from the tracking tables
      fastq_root_dir      => $self->o('fastq_root_dir'),
      #This will should be set to one in downstream config
      can_PreprocessIDR              => 0,
      can_run_SWEmbl_R0005_replicate => 0,
      can_DefineMergedDataSet       => 0,
    };
}




=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.


=cut

sub pipeline_analyses {
  my $self = shift;
  return
   [
    @{$self->SUPER::pipeline_analyses}, #To pick up BaseSequenceAnalysis-DefineMergedOutputSet

    {-logic_name => 'IdentifyAlignInputSubsets',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',
     -meadow_type => 'LOCAL',#should always be uppercase
     -parameters =>
      {set_type                     => 'InputSubset',
       feature_set_analysis_type    => 'peak',
       default_feature_set_analyses => $self->o('default_peak_analyses'),

       #Define these here, rather than hardcoding in IdentifySetInputs as is it generic
       #need to dataflow these just incase they have been redefined for a batch
       #why not batch flow as they will be required or flowed through
       #IdentifyAlignInputSubsets, DownloadInputSubsetsFactory, DefineInputSets, Run_BWA_and_QC_control
       #also if we don't use alignment_analysis in the input_set, then we would also need to flow it
       #to DefineMergedOutputSet, which might be lost across config(similar) to the IDR analysis
       #Hence let's change InputSet analysis to the alignment
       #This is basically reciprocating the ResultSet and we should probably j
       #ust drop one of InputSubset/Set
       dataflow_param_names => ['no_idr'], 
       # 'alignment_analysis'], #now batch flown
       #, 'broad_peak_feature_types'], Removed this for now as we currently don't allow over
       #of any of the defaults hashes apart from at initialisation
       #this is probably a case for using the hive batch flow
       #but this is still risky with respect to segregation of recover/rollback
       #flags which may not be safe to maintain across the analyses.
       #we could get around this by either changing there names to be more specific
       #or unsetting them after they have been used
      },
     -flow_into =>
      {
       '3->A' => [ 'DownloadInputSubset' ],
       'A->2' => [ 'DefineResultSets' ],#This will not always have pre-req Download jobs
      },
     -analysis_capacity => 1, #For safety, and can only run 1 LOCAL?
     -rc_name => 'default',   #NA as LOCAL?
     -failed_job_tolerance => 100,
     #We don't care about these failing, as we expect them too
    },


   # {-logic_name  => 'DownloadInputSubsetFactory',
   #  -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
   #  -meadow_type => 'LOCAL',
   #  #This either:
   #  #1 Downloads the shared controls and signal for all relevant inputsets
   #  # and semaphores the DefineInputSets job
   #  #2 Downloads a single experiment subsets which does not have a controls
   #  # and semaphores the DefineInputSets job

#     #For the download we are going to have to wait for everything (shared ctrl and all associated singal)
     #so we can put everything through the alignment at the same time
     #mirroring the alignment set up here (where we run the control first, then the reps and
     #flow directly on to the peaks), would mean direct flow of individual replicate sets
     #to the alignment jobs, meaning we have parallel jobs with the same controls
     #which would cause problems i.e.we might get parallel jobs trying to align the same control

#     -parameters =>
#      {dataflow_param_names => ['no_idr', 'alignment_analysis']},
#     -flow_into =>
#      {
#       #Is this the right markup for a factory?
#       '2->A' => [ 'DownloadInputSubset' ],
#       'A->1' => [ 'DefineResultSets' ],
#      },

#      -analysis_capacity => 100,       #Will this just run 1 as LOCAL?
#      -rc_name           => 'default', #NA as LOCAL?
#    },


    {-logic_name => 'DownloadInputSubset',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
     -parameters => {
                    },
     -analysis_capacity => 100,       #Could this go higher?
     -rc_name           => 'default', #TODO change this to long?
    },


    {-logic_name => 'DefineResultSets',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineResultSets',
     -meadow     => 'LOCAL',
     #-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
     #Inputs are either:
     #1 A set of controls and a series of InputSubset batches(replicates) which share the control
     #2 A series of unrelated InputSubset batches(replicates) which have no controls
    -parameters => {dataflow_param_names => ['no_idr']},#this is now batch flowed as we need it in DefineMergedOutputSet????
    -flow_into =>
     {#Needs to know about -no_idr and default_broad_peak_feature_types so it can flow correctly
      #To get broad_peak feature types to run IDR you would have to unset broad_peak_feature_types
      #We could alternately specify force_idr? Although redefining broad_peak_feature_types
      #has more finegrained control

      #We currently need to run the control jobs first!!! This is non-optimal
      #Would be better to have both the control job/analysis and the replicate job/analysis
      #form a multi semaphore for the down stream peak calling
      #2 => ['Preprocess_bwa_samse_control'], #This merges and aligns controls and then...
      #flows to either to:
      #1 A factory which will submits the replicate batches which share the controls
      # each job will semaphore the IDR job for that batch. Each of the replicate jobs
      # will flow directly to DefineReplicateOutputSet.
      #2 Merged alignment jobs which share this controls. These can flow directly to
      # DefineMergedOutputSet

      #For InputSubsets with no control we also need to short cut this,
      #and flow direct to the downstream analyses of Run_BWQ_and_QC_control
      'A->2'  => ['PreprocessIDR'],#move to 2 'A->2' but make sure I flow this last for each control group
      '10'    => ['Preprocess_bwa_samse_control'],
      '11'    => ['Preprocess_bwa_samse_merged'],
      '12->A' => ['Preprocess_bwa_samse_replicate'],

      #Could add another alignment analysis here easily e.g.
      #'20'    => ['Preprocess_NEWALIGNER_control'],
      #'21'    => ['Preprocess_NEWALIGNER_merged'],
      #'22->A' => ['Preprocess_NEWALIGNER_replicate'],

      #Could support custom analyses like this
      #i.e. before we start the decimal blocks of
      #known resource optimised/branched analyses
      #Assuming we will never have more the 7 outflow analyses
      #allow_custom_branching would need to be an analysis level param
      #'3'    => ['Preprocess_custom_control'],
      #'4'    => ['Preprocess_custom_merged'],
      #'5->A' => ['Preprocess_custom_replicate'],


      #Also need to replicate this for custom aligner? sigh, leave for now
      #'100' => [ 'Run_QC_and_CustomAligner_merged' ],
     },
     -analysis_capacity => 100,       #Will this just run 1 as LOCAL?
     -rc_name           => 'default', #NA as LOCAL?
    },

    ### PREPROCESS ANALYSES ###
    #Currently averaging ~22mins
    #But this depends entirely on the number and size of the reps in the control

    {-logic_name => 'Preprocess_bwa_samse_control',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
     #-parameters => {run_controls => 1}, #This is now implicit with the flow of signal 'result_set_groups'
     #rename this to control_set_groups or something more obvious?
     -flow_into =>
      {
       '2->A' => ['Run_bwa_samse_control_chunk'],
       'A->3' => ['MergeControlAlignments_and_QC'],
       },
     -batch_size => 1, #max parallelisation???
     -analysis_capacity => 200,
     -rc_name => 'normal_high_mem'},


     {-logic_name => 'Preprocess_bwa_samse_merged',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
     -parameters => {merge => 1},
     -flow_into =>
      {
       '2->A' => ['Run_bwa_samse_merged_chunk'],
       'A->3' => [ 'MergeAlignments_and_QC' ], #This is in Peaks.pm conf
       },
     -batch_size => 1, #max parallelisation???
     -analysis_capacity => 200,
     -rc_name => 'normal_high_mem'},


     {-logic_name => 'Preprocess_bwa_samse_replicate',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
     #-parameters => { },
     -flow_into =>
      {
       '2->A' => ['Run_bwa_samse_replicate_chunk'],
       'A->3' => [ 'MergeReplicateAlignments_and_QC' ],
       },
     -batch_size => 1, #max parallelisation???
     -analysis_capacity => 200,
     -rc_name => 'normal_high_mem'},


    ### RUN BWA ANALYSES ###
    #These are currently averaging ~35mins, so can probably up to batch_size 2 to reduce submission time?
    #Although this will just throttle this step to 70mins, regardless of load
    #if we maintain it at 1, then all the analysis capacity will be used
    #meaning batching will effectively be done by the workers, which will probably take little time
    #to be submitted, but we suffer from this anyway.

    {-logic_name => 'Run_bwa_samse_control_chunk',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
     -batch_size => 1, #max parallelisation???
     -analysis_capacity => 1000,
     -rc_name => 'normal_10gb'},

    {-logic_name => 'Run_bwa_samse_merged_chunk',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
     -batch_size => 1, #max parallelisation???
     -analysis_capacity => 1000,
     -rc_name => 'normal_10gb'},

    {-logic_name => 'Run_bwa_samse_replicate_chunk',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
     # These jobs can be run in parallell... don't put too many since it may generate many jobs...jobs!
     #-limit => 1,#what is this?
     -batch_size => 1, #max parallelisation? Although probably want to up this so we don't hit pending time
     -analysis_capacity => 1000,
     -rc_name => 'normal_10gb'},


    ### MERGE ANALYSES ###
    #These all use 2cpus as there use a double piped samtools command
    #which runs 3 samtools processes. 2cpus was the suggestion of ISG based on
    #the output of their monitoring programs
    #It may be safer and quicker to schedule if these we broken down into individual
    #processes. But this would increase footprint and require tidy up.

    #flow_mode could alternatively be flowed from the ProcessFastqs jobs
    #Leave it here for now

    #We could have had just 1 merge analysis but this would have meant selective data flow
    #for each instance. This is more flexible, and requires less hardcording.
    #We still need branch config here to handle the aligners dynamically
    #We will never use the control branch(2) in AlignmentQC,
    #This is reserved for DefineDataSet params dataflow
    #This also prevents the issue of assinging an arbitrary much higher branch number
    #to avoid clashes with any new aligner branches which might be added

    #If there is no branch_config, then we just default to flow the DataSet params
    #else we will flow the Preprocess and ReplicateFactory aligner jobs dynamically
    #using the alignment_branches config
    #This approach is implicit and may get broken if alignment_branches     #ever get batch flowed or movbed to pipeline_wide?
     #hence it would be safer to set a flow mod
     #
         #How do we fail QC without failing the job?
     #Need to update a QC report table



    {-logic_name => 'MergeControlAlignments_and_QC',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments',
     -parameters => {flow_mode => 'signal'},
     #actually, this can be inferred from the presence of result_set_groups
     #which is only set in PreprocessFastqs if run_controls is set

     #Flow to the signal alignment analyses (as we have now done the control pre-requisite).
     -flow_into =>
      {
       #2 is reserved for other Define DataSet flow (single vs multiple ResultSets).
       # Although isn't really required as the branching is handled dynamically
       #but let's keep it clean here for now.
       'A->3'  => [ 'PreprocessIDR' ],
       #alignment analyses encoded in blocks of 10
       '10'    => ['Preprocess_bwa_samse_merged'],
       '11->A' => ['Preprocess_bwa_samse_replicate'],
       },
     -batch_size => 1, #max parallelisation
     -analysis_capacity => 200,
     -rc_name => 'normal_5GB_2cpu_monitored',
    },


    {-logic_name => 'MergeAlignments_and_QC',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments',
     -parameters => {flow_mode => 'merged'},
     -flow_into => { 2 => ['DefineMergedDataSet']},
     -batch_size => 1, #max parallelisation
     -analysis_capacity => 200,
     -rc_name => 'normal_5GB_2cpu_monitored',
    },


    {-logic_name => 'MergeReplicateAlignments_and_QC',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments',
     -parameters =>
      {flow_mode => 'replicate',
       permissive_peaks         => $self->o('permissive_peaks')},
       #make permissive_peaks pipeline wide and remove from here?

     #peaks analyses encoded from 100
     -flow_into => {'100' => ['run_SWEmbl_R0005_replicate']},
     -batch_size => 1, #max parallelisation
     -analysis_capacity => 200,
     -rc_name => 'normal_5GB_2cpu_monitored',

     #This will definitely have to write some accu data
     #such that Run_IDR can pick up on any QC failures

    },






   #WARNING!
   #How do we protect again the same set details being flowed through both configs?
   #This may fail due to the presence of a duplicate input_id when flowing to the Merged analysis
   #but not before we might attempt a rollback (if rollback/recover is set)?
   #This may end up rolling back a set which is actively being processed via the other branch of the pipeline
   #Leave this for now due to safety concerns?


    ### LINK ANALYSES TO DOWNSTREAM CONFIGS

    ### WARNING
    # These would normally still run and produce output even if we haven't topped up!
    # This is normally fine for IdentifySetInput as this only creates dataflow output
    # but for analyses which create output in the DB (e.g a DataSet) then subsequent top ups
    # may create jobs which try and recreate the DataSet and hence would require some recovery/rollback
    # management.
    # To avoid this we can instruct the link analysis exit without running.
    # This is managed with a pipeline_wide_parameter run_DefineMergedOutputSet = 0
    # Then set it to 1 DefineOutputSets and RunIDR
    # This will give us successful jobs which did nothing if we add the following line
    # to the top of the fetch_input, run and write methods:
    #     return if $self->param_silent('run_DefineMergedOutputSet') == 0;
    # When trying to reseed them, the input_id will either be different
    # or we will get an error, at which point we can reset/retry them, which may be messy
    # Alternatively the jobs could fail not retry, which would be easier to pick up
    # but would be harder to differentiate between true failures jobs where
    # run_DefineMergedOutputSet == 1


     #Peaks-DefineMergedOutputSet for flowing of non-IDR merged replicate jobs i.e. CCAT
     #This is now in BaseSequenceAnalysis as it is common to all
     #either as a 'link out' analysis or as a 'link from' analysis.
     #Link from analyses, still have to be specified in the relevant config
     #as this will also host the required -flow_into spec
     #which cannot be specified for a link out analysis


     #These are in IDRPeaks.pm config,
     #This is for flowing of IDR replicate jobs  i.e. permissive SWEMBL
    #{
    # -logic_name => 'DefineMergedDataSet',
    # -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineDataSet',
    # -meadow_type => 'LOCAL',#should always be uppercase
    # -batch_size => 10,
    # -parameters => {check_analysis_can_run => 1},
    # #No flow into conf as this will be defined in the top up conf
    # },


    #Could split this out into a mixin conf (or BaseSequenceAnalysis)
    #as this is a shared analysis between ReadAlignment and IDRPeaks
     {
     -logic_name    => 'run_SWEmbl_R0005_replicate',  #SWEmbl permissive
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
     -parameters =>
      {check_analysis_can_run => 1,
       peak_analysis => $self->o('permissive_peaks'), #This will not allow batch override!
       #Would have to flow this explicity from IdentifyReplicateResultSets and MergeReplicateAlignments_and_QC
      },
     -analysis_capacity => 10,
     -rc_name => 'long_monitored_high_mem', # Better safe than sorry... size of datasets tends to increase...
     #This is to stop the beekeeper exiting when all these jobs fail
     #would it be better to have them succeed, as we will be resetting them anyway?
    },


    #This is sempahored from the ReplicateFactory
    #which flows directly DefineReplicateOutputSet and then on to the individual Peaks jobs
    {
     -logic_name => 'PreprocessIDR',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
     -batch_size => 30,
     -rc_name    => 'default',
     -parameters =>
      {check_analysis_can_run => 1,
       permissive_peaks => $self->o('permissive_peaks'), #This will not allow batch override
      },
     #No flow into conf as this will be defined in the top up conf

     #This ill have to pick up some accu data from Preprocess_BWA_replicate?
     #or MergeReplicateALignments_and_QC?
     },


   ];
}

1;

