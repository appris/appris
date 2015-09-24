


package Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis;

use strict;
use warnings;
use base qw(Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB);

# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly

=head2 default_options

    Description : Implements default_options() interface method of 
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used
                  to initialize default options.
    Caller      : DependantOptions::process_options (init_pipeline.pl)

=cut

sub default_options {
  my $self = shift;  
  
  return 
   {
    %{$self->SUPER::default_options},
    #todo move some of these defaults move to the BaseDB or DefineOutputSets config?

     #Would these ever be needed outside the Importer context i.e. in a generic DB context?
      
     #Peaks conf will have access to this as it is topped up from DefineOutputSets which is
     #a BaseImporter, even those the Peaks pipeline doesn't use the importer
     #Actually, yes it does via PreprocessAlignments


               
     #all of the default_peaks values(logic_names) should be represented here, otherwise
     #they will be lumped into generic branch 100
     #branch number is used by CollectionWriter to flow correctly
     #from the PreprocessAlignments analysis in the Peaks conf
     #todo, can we move this to Peaks conf?
     #or will this not be updated?
     
     
    
     #These allow batch wide override of defaults
     #result_set_analysis  => undef,
     #feature_set_analysis => undef,
     #or override the following directly for fine grained controls
     
     #default feature_set_analyses
     default_peaks       => 'SWEmbl_R015', 
     default_tight_peaks => 'SWEmbl_R0025', 
     default_broad_peaks => 'ccat_histone',
     permissive_peaks    => 'SWEmbl_R0005',
     idr_peaks           => 'SWEmbl_R0005_IDR', 
     
     #default result_set_analyses
     #default_collection          => 'bwa_samse',#todo update this to reflect collection analysis
     #default_5mC                 => undef,#do we need to split this for  RRBS_merged_filtered_10 or WGBS_merged?


    #These are really DefineOutputSet specific
    #but as ReadAlignment config joins to DefineOutputSets config
    #via DefineOutputSet, we need the DefineOutputSet parameter
    #config in here(BaseDB) so it is accessible to both configs 
    
    #Now linking via IdentifyInputSets to avoid creating the output data set
    #if running ReadAlignments in isolation
    
    broad_peak_feature_types => ['H3K36me3', 'H3K27me3'],
    #Slight redundancy here wrt default_peak_analyses
    #Can we switch this around, so they are all defined like this?
    #e.g. default_peak_feature_types = ['Histone', 'Transcription Factor'],
    #This would be slightly more cumbersone wrt access
    #But woudl eradicate potential for unsync'd broud peak analyses
    #or can we map this above?   
    
    default_peak_analyses => 
     {
      #These $self->o methods interpolate above spec 
      Histone             => $self->o('default_peaks'),
      'Transcription Factor' => $self->o('default_peaks'), 
      'Polymerase'        => $self->o('default_peaks'), 
      DNase1              => $self->o('default_tight_peaks'),
      H3K36me3            => $self->o('default_broad_peaks'), 
      H3K27me3            => $self->o('default_broad_peaks'),
      #(map{ $_ => $self->o('default_broad_peaks') } @{$self->o('broad_peak_feature_types')}),
      #(map{ $_ => $self->o('default_broad_peaks') } "@{#broad_peak_feature_types#}"),
      #('#expr(map{ $_ => #default_broad_peaks# } @{#broad_peak_feature_types#})expr#'),
     }, 
   
 
    'control_feature_types' => ['Goat-IgG', 'Rabbit-IgG', 'WCE', 'rat-IgG-control', 'rabbit-IgG-control', 'mouse-IgG-control', 'GFP'], 
    # Possible add 'negative' to this?
    # was in Peaks.pm but needed for ReadAlignments 
    # Only needed for IdentifyInputSubsets, which does the control association
    # Move this to ReadAlignments as we don't need it for Peaks at all?
    
  
    'checksum_optional' => 0,   
   };
  
}

=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.
    Caller      : HiveGeneric_conf::run (init_pipeline.pl)

=cut

#$self->o here will take the cmdline option else default to what was set in default_options

#Will init_pipeline with analysis_topup, use previously stored opts so we don't have to define
#them again?


sub pipeline_wide_parameters {
  my $self = shift;
    
  #Deal with batch_params first as this may not have been subtituted yet
  #and will evaluate to a string:
  #Can't use string ("#:subst batch_params:#") as an ARRAY ref while "strict refs"
  #my $batch_params = 
  #  [ ( ref($self->o('batch_params')) ? 
  #      @{$self->o('batch_params')} : () ),
  #    'control_feature' ];
  #This doesn't work with #expr()expr# either

  return 
   {
    %{$self->SUPER::pipeline_wide_parameters}, 
    # 'alignment_root_dir' => undef, # defaults to what is set in  Base::alignment_root_dir 

    'control_feature_types'    => $self->o('control_feature_types'),    
    #'idr_analysis'            => $self->o('idr_analysis'),
    'broad_peak_feature_types' => $self->o('broad_peak_feature_types'), 
    
    
    'checksum_optional'        => $self->o('checksum_optional'),
    #DefineResultSets currently requires this in ReadAlignments, but may be able to change that
  
  
    #Todo feature_file_fromat, result_set_onluy/mode code should be moved to a new BaseSequenceAnalysis Runnable?
    #Maybe result_set stuff can stay. But certainly feature_file_format or just remove if we don't need any more?
  
    #Ad they are now in here, there is a possibility that another Runnable will inherit from Base/BaseDB but 
    #not have these as batch params
    #
  
     ### THIS NEEDS REWORKING FOR IDRPEAKS ###
           
           
    #Any of these which have defaults, also need to be specified as pipeline_wide_parameters
    #otherwise they will never be picked up
    #This is likely only required for non-null default config  
    #Currently over-ride of any default config which is not batch flown will only work 
    #when initialising the DB (is this even currently supported?)
    
  
    batch_param_names => 
      [#From Base.pm       
       'no_write', #For use with runWorker.pl -no_write, so we can write some STDOUT in run
                   #is this already available in the job, or is it just passed ot the worker?
       #'feature_file_format',
     
     
        #BaseDB.pm batch_params     
        #Generic optional params (used in Helper and elsewhere)
        'rollback',
        'full_delete',
        'slices',
        'skip_slices',
        ### More optional Helper params (currently used in DefineOutputSet)
        'result_set_only',
        'result_set_mode', #Can only be set to recover at present
        'recover',         #is this an Importer or a Helper param?
        
        
        ### Optional IdentifySetInputs parameters
        #comma separated or defined as list ref  
        #are these batch wide or just used for the seed job?
        #qw( cell_types feature_types input_sets input_set_ids
        #    experimental_groups states ),             
        #'input_analyses'     => undef,
        #'experiment'   => undef, 
        #This is actually more like study, and omit for now as 
        #input_set will handle this
        
        #DefineOutputSets.pm batch_params   
         #These allow override of defaults  
        #'result_set_analysis',#isn't this alignment analysis? Should this be in BaseImporter
        #This has now been removed in favour of more spcific analysis names as we may want >1 result_set_analysis
        #This is translated by the runnable itself, or in the case of a generic runnable
        #with have from some analysis level config
        
        'alignment_analysis',
        'peak_analysis',
        'permissive_peaks',
        #? This will require RunPeaks changes? 
        #Currently harcoded in config, but probably want 
        #to flow from IdentifyReplicateResultSets
        
        
        #Peaks.pm batch params
        'control_feature',
        
        #ReadAlignments.pm batch params
        'no_idr', #This is required
        #'bam_filtered', #Required by CollectionWriter and RunPeaks
        #Now we assume/expect the filtered files
        #this is due to potential conflicts across parallel jobs asking for the 
        #control files to be filtered.
        #
        
        
        #'alignment_analysis' #Needs to flow from IdentifyAlignInputSubsets>DefinesResultSets
        #This is just a data flow param!! But we know we always want to flow this just between
        #these two jobs, so do it explicitly    
                  
                  
       #Don't need no_idr here, as this is only used by IdentifyAlignInputSubsets
       #The rest of the IDR handling is implicit by the output and dataflow of this analysis           
        
        'indexed_ref_fasta',   
        'idr_analysis', #Not currently used as this is done outside of the DB
        'max_peaks',
        'checksum_optional'        
      ],      
    };
}


#pipeline_create_commands in Base and/or subclass config 
#pipeline_analyses in subclass config

#Move link from analyses here, so we don't risk
#having unsynced param config. Params for these analyses are likely to be in 
#here anyway.
#Is this possible? can we define dataflow and semaphores in 
#a inherited fashion?
#Likely yes, so long as we make sure we call super appropriately
#Issue here is that when using individual configs, we will get a few
#unlinked analyses



sub pipeline_analyses {
  my $self = shift;
  
  #This should only contain 'link out' analyses which are common to >1 (preferably all)
  #configs. This is to prevent redundant analysis spec whichy may be come out of sync between
  #configs. As these are link out analyses, the -flow_into spec should be omitted and specified
  #in the 'link in' analysis of the downstream config. Special care must be take to match
  #the other aspect of the config, particularly the -parameters, as they will not be updated from 
  #the 'link in' analysis spec.
  #However, multiple 'link in' analyses are permitted, and the -flow_into spec from these
  #will be updated.
  
 
  return 
   [#To use these the following must be placed in your subclass config:
    #@{$self->SUPER::pipeline_analyses}, 
  
    {
     -logic_name  => 'DefineMergedDataSet',
     -module      => 'Bio::EnsEMBL::Funcgen::Hive::DefineDataSet',     
     #-meadow_type => 'LOCAL',#should always be uppercase
     #Not a great idea as we maybe dealing with 100s, so will bottle neck 
     #on the local node/machine
     
     -parameters  => 
      {
       default_feature_set_analyses => $self->o('default_peak_analyses'),
       feature_set_analysis_type    => 'peak',
       check_analysis_can_run       => 1,
      },
       
     -analysis_capacity => 100,
     -rc_name           => 'default',
     -batch_size        => 10, 
     #None of these shoudl take > 2-3 mins unless there is some rolling back to do
     #But this will only ever be deleting annotated_feature records
     #and maybe some states or updating sets
     #Keep this fairly low, so we a getting the parallel compute running asap.

     #Not having a -failed_job_tolerance here is causing the beekeeper to 
     #exit, especially as there is no -max_retry_count set either
    
     #No flow into conf as this will be defined in the top up conf
    },
   ];
}


1;

