
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


package Bio::EnsEMBL::Funcgen::Hive::Config::Peaks;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut

#sub default_options {
#  my $self = shift;
#  return {
#    %{$self->SUPER::default_options},  # inheriting database and hive tables creation
#	 };
#}

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
      #can_Link_to_Peaks       => 1,
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
 
    #This is overkill for sorting the input for the peak calling
    #But we have to use the same analysis as the Collection conf
    #otherwise we might get conflicting jobs.
 
    #Does this actually do the filtering as required?
    #i.e. is it using get_alignment_file_by_InputSets?
    #also, is it redundantly preparing the file?
    #Yes get_alignment_file_by_InputSet
    #sorts and filter MT and chrM and dups from BAM
    #Importer simply sorts bed and returns slice names
    
    #todo remove this sort duplication
    #i.e. make pre_process_file use the same code
    #and remove redudant sort?
    #Skip sort in pre_process_file, but allow it to run, so we can get the slices required.
    #This will require some monkeying aroung with the 'prepared' flags
    
    #Also, don't copy file, as this is just duplicating the footprint?
    #Will have to rename file as per expected by the InputSet importer
    #
    

    { #This is basically making sure the input file is sorted wrt genomic locations
	 -logic_name => 'PreprocessAlignments',
	 -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	 -flow_into => 
      {
       #Can we not wire all these on the same branch?
       #2 is slice jobs i.e. WriteSignalColections         
       #'1->A' => [ 'run_peaks', 'run_peaks_wide', 'run_macs' ],   
       #this would give each analysis all the input!
       #ideally we want to flow these non-reundantly
       #so we would need separate branches anyway!
       #this means CollectionWriter will need access to the default
       #analyses to know how to flow them
       #unless we create a link analysis?
     

       #Not 3->A as we don't close these funnels
	   #These branch numbers must match the corresponding default config in DefineSets
       '3'   => [ 'run_SWEmbl_R015' ],      #Normal peaks
       '4'   => [ 'run_ccat_histone' ],     #CCAT
       '5'   => [ 'run_SWEmbl_R0025' ],     #Tight peaks 
       '6'   => [ 'run_SWEmbl_R0005_IDR' ], #Permissive peaks filtered by IDR max peaks value
       '100' => [ 'run_peaks_custom' ],     #defined by feature_set_analysis
    
        
       #See notes in IdentifyInputSets config for implementation of 
       #replicate/IDR config 
         
       #This is no different to data flowing directly from the 
       #individual run_peaks analyses
       #but we can define it here once instead of in each run_peak analysis
       #'A->1'   => [ 'PeaksReport' ],
         
         #Merge this into RunPeaks now?          
      },
			 
	 -analysis_capacity => 100,
	 #Change this to hive_capacity as it may be competing with parallel peak jobs?
     -rc_name => 'default',
     #todo revise this to reserve tmp space as this is sorting the bed files
    },
    
 
  
  
    #we have an issue here as the fan changes from a slice fan, to an analysis fan
    #so branch 1 get's overloaded
    #do we need a factory here as we can't define fan on branch one and an accumulator
    #basically
    #should we just rebranch the collection confs to handle this 
    #or can we feed and wait for accumulators on the same branch?
    
    
      #Merge QC into this step? rather than flow to PeaksQC
      #Ww always want it run, but we don't really want to merge it into the PeakCaller
      #so merge it into RunPeaks?
      #We can do this for IDR as this runs on bed files i.e. before we do the write_output step
      #Until we do IDR (i.e. peak calling on separate replicates)
      #we will have to rely on peaks report
      #dowe do this for each set seprately or wait and report across all for comparison
  
      #what if only the PeaksQC failed, reload will enable us to rerun QC?
      #do we need added support for skip QC?
      #certainly need support for force_load if QC has failed
 
      #Do individual set Peak
 
    #Do not change these logic_names as they are used for dynamic dataflow/branching
    
 
	  {
	   -logic_name    => 'run_SWEmbl_R015',  #SWEmbl normal (histones/Dnase?)
	   -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
	   -analysis_capacity => 100,
	   -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...	   
	   
	  },

	  {
	   -logic_name    => 'run_ccat_histone', #CCAT
	   -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	   -parameters    =>
	    {process_file_types => ['significant_region'],
	     CCAT_parameters => {-chr_file => $self->o('data_root_dir').
	                           '/reference_files/CCAT/'.$self->o('species').'_'.$self->o('assembly').'.CCAT_chr_lengths.txt'}},
	   -analysis_capacity => 100,
	   -rc_name => 'normal_monitored_2GB', # CCAT does not need much?
	   #we were getting MEMLIMIT failures with default 200mb
	  
	  },

#fragmentSize 200 slidingWinSize 150 movingStep 20 isStrandSensitiveMode 0 minCount 10 outputNum 100000 randomSeed 123456 minScore 4.0 bootstrapPass 50



      {
	   -logic_name    => 'run_SWEmbl_R0025', #SWEmbl TF
	   -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	   -hive_capacity => 100,
	   -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...
	  # -flow_into => 
	  },
	  
	  
	  #Does this need to be shared with IDRPeaks
	  #actually making analysis names tied to branching
	  #means that we can't reuse the mechanism if we want to flow out differently
	  #Luckily this is not the case
	  
	  {
     -logic_name    => 'run_SWEmbl_R0005_IDR', #SWEmbl permissive IDR filtered
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
     -parameters    => {'filter_max_peaks' => 1},
     -hive_capacity => 100,
     -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...
    # -flow_into => 
    },
    
    
    
	  #Add others in here e.g. macs etc
	  
      { 
       -logic_name    => 'run_peaks_custom', #Peaks caller not supported by config 
       -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', # Peak caller defined by feature_set_analysis
       -hive_capacity => 100,
       -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...
       #-flow_into => 
      },
	  
  
      #{
	  # -logic_name    => 'IDR_QC',
	  # -module        => 'Bio::EnsEMBL::Funcgen::Hive::IDR_QC',
	  # -hive_capacity => 10,
	  # #Control files should be handled by setup_pipeline.
	  # -rc_name => 'long_monitored_high_mem', # Better safe than sorry... size of datasets tends to increase...

      #Do we need to pass output from here to set up MotifPipeline analysis?
      #i.e. those that fail QC should be included in MotifPipeline or infact the reg_build
      #How would we flow through to the build conf(after the MF conf)
      #Simply through IdentifySets
      #This would have to be able to do the same identify by a fetch method with constraints
      #or maybe this should simply be what is flowed
      #the constraint, rather than the set IDs?
      #constraint would simply be those marked as TO_BUILD and not FAILED_PEAKS_QC
      
	  #},



    #Merge this into RunPeaks now
    #Peaks report can only be done if a set has been loaded
    #not easy to block MotifFeatures if we don't flow directly
    #PeaksReport must branch from PreProcessAlignments

    #todo peaks report on all featuresets for comparison
    #  {
    #   -logic_name    => 'PeaksReport',
    #   -module        => 'DUMMY',#'Bio::EnsEMBL::Funcgen::Hive::PeaksReport',
    #   -hive_capacity => 10,
    #   #Control files should be handled by setup_pipeline.
    #   -rc_name => 'long_monitored_high_mem', # Better safe than sorry... size of datasets tends to increase...

    #   #Do we need to pass output from here to set up MotifPipeline analysis?
    #  #i.e. those that fail QC should be included in MotifPipeline or infact the reg_build
    #  #How would we flow through to the build conf(after the MF conf)
    #  #Simply through IdentifySets
    #  #This would have to be able to do the same identify by a fetch method with constraints
    #  #or maybe this should simply be what is flowed
    #  #the constraint, rather than the set IDs?
    #  #constraint would simply be those marked as TO_BUILD and not FAILED_PEAKS_QC
    #  
    # },

	  
	  

      #Move this to MotifFeatures config
      
	  #{
      # -logic_name => 'IdentifyFeatureSets',
      # -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',     
      # -meadow_type => 'LOCAL',#should always be uppercase
      # #general parameters to pass to all jobs, use_tracking_db?
      # -parameters => {set_type    => 'feature_set'},
      # #-input_ids #NONE! These will be data flown from PreprocessAlignments
      # #(or seeded in the next conf, if we are not topping up) 
      # #No flow into spec here, as this will be defined in the top up config(s)
      # #i.e. the motif pipeline            
      # #-flow_into => {         },
      # -analysis_capacity => 10,
      # -rc_name => 'default',
      #},
	  
	 ];
}

1;

__END__


   
=pod
   
    {
     -logic_name => 'Link_to_Peaks',
     -meadow     => 'LOCAL',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MultiConfigLinker',
     
     #Need to added in previous config as will not be updated by -analysis_topup
     #-parameters => 
     # {
     # },
     
     -flow_into =>
      { #Maintained original branches here, but could revise these
       '3'   => [ 'run_SWEmbl_R015' ],      #Normal peaks
       '4'   => [ 'run_ccat_histone' ],     #CCAT
       '5'   => [ 'run_SWEmbl_R0025' ],     #Tight peaks 
       '6'   => [ 'run_SWEmbl_R0005_IDR' ], #Permissive peaks filtered by IDR max peaks value
       '100' => [ 'run_peaks_custom' ],     #defined by feature_set_analysis
      }, 
       
     -analysis_capacity => 100,
     -rc_name => 'default',
    }, 
    
=cut   
  
    
