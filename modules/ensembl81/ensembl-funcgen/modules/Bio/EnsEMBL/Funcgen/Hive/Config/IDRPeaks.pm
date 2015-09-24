
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::IDRPeaks;

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::Hive::Config::IDRPeaks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut


sub default_options {
  my $self = shift; 
  #Define any optional params here as undef
  #If they are mandatory for one analysis, but optional for another
  #will have to catch that in the runnable
  return 
   {
    %{$self->SUPER::default_options},        
  
    ### THIS NEEDS REWORKING FOR IDRPEAKS ###
   };
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut


#Can we move some of these to Base.pm config?

sub pipeline_wide_parameters {
  my $self = shift;
               
  return 
   {
    %{$self->SUPER::pipeline_wide_parameters},
    
    #Arg! we can use this approach if we are reusing an analysis?
    #As there is no way of there is no way of the analysis knowing which param to use?
    #
    
    can_run_SWEmbl_R0005_replicate => 1,#'IDRPeaks',
    can_PreprocessIDR              => 1,#'IDRPeaks',
    #Set to the config name so the configure_pipeline script knows to reset these jobs 
    can_DefineMergedDataSet        => 0, 
   };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands
      that will create and set up the Hive database.

=cut

## WARNING!!
## Currently init_pipeline.pl doesn't run this method when a pipeline is created with the -analysis_topup option



=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.


=cut



#(2) The hive_capacity worker limiting mechanism that was in place for years is going to change slightly.
#The old meanings of the values were:
#	negative value  : checking is switched off in this particular analysis
#	zero		: no workers will be allowed to take this analysis
#	missing value   : sets the default, which is 1
#	positive value	: (including the 1 by default) this analysis is limited by the given value
#This was counter-intuitive, because by default there was already a limit, which people had to raise or switch off if they needed.
#Now, since we also have an alternative mechanism (analysis_capacity), I'd like to make both mechanisms "off" by default
#(missing value will mean "checking is switched off").
#
#So please, please, please - check whether your pipeline *RELIES* on the current default behaviour
#(which is "if you don't specify -hive_capacity, the system will set it to 1
#and will not allow any workers of other analyses to take any jobs while this analysis is running").
#
#If you know/find out that your pipeline *RELIES* on this functionality, please explicitly set -hive_capacity => 1 in the corresponding analyses.
#This will make your pipeline compatible with future releases of the Hive system.





sub pipeline_analyses {
  my $self = shift;

  return [
   @{$self->SUPER::pipeline_analyses}, #To pick up BaseSequenceAnalysis-DefineMergedOutputSet
   
   {
    -logic_name => 'IdentifyReplicateResultSets',
	  -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',	 
	  

    #Currently set to dummy until we define branch 4 output from IdentifySetInputs
    #which will group replicate ResultSet outputs based on parent merged ResultSet
    #how are we going to identify that if we don't have a parent set created?
    #Could do this simply based on a set name match?
    #by stripping off the _TR1 number
    #better way would be to re-use the code which identified them as a merged set in the first place
    #which is already in IdentifySetInputs
	   
	   
	  -meadow_type => 'LOCAL',#should always be uppercase
	  
	  #general parameters to pass to all jobs, use_tracking_db?
	  -parameters => {set_type        => 'ResultSet',
	                  only_replicates      => 1, 
	                 #This might need to take a -replicate flag
	                 #to ensure we only identify single rep InputSets
	                 #Probably need a naming convention i.e. suffix of TR_[1-9]*	                
	                 },
	             
      #This will fan into the rep peak jobs
      #and semaphore the IDR job, which will need all the input_set ids
      #including the final InputSet
      
      #This IDR jobs needs to record the analysis params and associate them with the FeatureSet
      #so we need to link on DefineOutputSets in here
      #Hence we also need all of the collection config! which I have just deleted. doh!
      
	
	#DOES THIS NEED TO BE A FACTORY Or do we need to flow into a factory?
	#We need to semaphore the IDR job based on a batch of replicates
	
	
	 #This needs to change now as we will have to fan the IDR jobs.
	 #DefineReplicateDataSet can still do the post processing tho?
	 #No wait, this is te wrong way around
	   -flow_into => 
	    {		
	     # 'A->2' => ['PreprocessSelfIDR'],
       'A->2' => ['PreprocessIDR'],  
	     '3->A' => [ 'run_SWEmbl_R0005_replicate' ], #['DefineReplicateDataSet'],
      },
			 
	   -analysis_capacity => 100, #although this runs on LOCAL
     -rc_name => 'default',
     -failed_job_tolerance => 100, 
     #We don't care about these failing, as we expect them too
    },
	
      
    {
     -logic_name    => 'run_SWEmbl_R0005_replicate',  #SWEmbl permissive
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
     -parameters => 
      {
       #peak_analysis => $self->o('permissive_peaks'), #This will not allow batch override!
       #Now batch flown, so we need to flow this explicity as peak_analysis from IdentifyReplicateResultSets and MergeReplicateAlignments_and_QC
      },
     -analysis_capacity => 300,
     -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...       
    },
  
    {
     -logic_name    => 'PreprocessSelfIDR',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
     #-meadow        => 'LOCAL', #Could potentially do this locally due to low load/runtime?
     #although it is potentially a >1 cpu job due to triple pipe command
     -analysis_capacity => 100,#Unlikely to get anywhere near this
     -rc_name    => 'default',
     -batch_size => 30, #Should really take ~1min to process each set of replicates
            
     -parameters => { permissive_peaks => $self->o('permissive_peaks') },
      #This will not allow batch override
           
     -flow_into => 
      {
       '2->A' => [ 'RunSelfIDR' ],                   # fan
       #'3->A' => [ 'RunPooledPseudoRepIDR ' ],   # fan #This doesn't need a separate analysis
       'A->3' => [ 'PostProcessSelfIDRReplicates' ], # funnel
      }, 
    },
    
    
       {
     -logic_name    => 'RunSelfIDR',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunIDR',
     -analysis_capacity => 100,
     -rc_name    => 'normal_2GB', #~10mins + (memory?)
     -batch_size => 6,#~30mins +

     #No flow into here, but this analysis should update the DB with the calculated threshold
     #This can't be in the feature_set_stats table, due to the combinations between the reps.
     #If we accu this, then we will lose the data!!! Preventing us from being able to simply drop one rep and
     #rerun the post process?
     #do we need a stand alone table for idr_thresholds, which simply has the threshold and the 2 feature_set_ids
     #should we store anything else here? 
     
    -flow_into => 
     {
      2 => [ ':////accu?self_idr_peak_counts=[accu_idx]' ],
     }

           
 
      #Or should this do all this in the same analysis
      #where are we going to cache the run_idr output for the
      #final peak calling threshold? Let's keep this in a tracking DB
      #table to prevent proliferation of analyses based on this value differing between data sets.
      
      #are we going to have problems having these two analyses together
      #if we want to rerun the analysis due to the creation of the merged rset failing?
      #would have to rerun idr too?
      
      #This is extremely unlikely to happen
      #and could do some funky job_id manipulation to set skip_idr
      
      #No we won't be able to reuse DefineResultSets here
           
    },
  
    {
     -logic_name    => 'PostProcessSelfIDRReplicates',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR',
      #-meadow        => 'LOCAL', #Could potentially do this locally due to low load/runtime?
     #although it is potentially a >1 cpu job due to triple pipe command
     -analysis_capacity => 100,
     -rc_name    => 'default', #<1mins
     -batch_size => 10,#?
     -parameters => 
      { 
       #result_set_mode              => 'none',#We never want one for an IDR DataSet
       #default_feature_set_analyses => $self->o('permissive_feature_set_analyses'), 
      },
           
     -flow_into => 
      {
       '2' => [ 'PreprocessIDR' ],
      }, 
      
      #Or should this do all this in the same analysis
      #where are we going to cache the run_idr output for the
      #final peak calling threshold? Let's keep this in a tracking DB
      #table to prevent proliferation of analyses based on this value differing between data sets.
      
      #are we going to have problems having these two analyses together
      #if we want to rerun the analysis due to the creation of the merged rset failing?
      #would have to rerun idr too?
      
      #This is extremely unlikely to happen
      #and could do some funky job_id manipulation to set skip_idr
      
      #No we won't be able to reuse DefineResultSets here
           
    },
    
    
    
  
    #SubmitIDR is a simple link job, to handle sugmitting jobs
    #as hive will not handle multi semaphores  
  
     {
     -logic_name    => 'PreprocessIDR',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
     #-meadow        => 'LOCAL', #Could potentially do this locally due to low load/runtime?
     #although it is potentially a >1 cpu job due to triple pipe command
     -analysis_capacity => 100,#Unlikely to get anywhere near this
     -rc_name    => 'default',
     -batch_size => 30, #Should really take ~1min to process each set of replicates
            
     -parameters => { permissive_peaks => $self->o('permissive_peaks') },
      #This will not allow batch override
           
     -flow_into => 
      {
       '2->A' => [ 'RunIDR' ],                   # fan
       #'3->A' => [ 'RunPooledPseudoRepIDR ' ],   # fan #This doesn't need a separate analysis
       'A->3' => [ 'PostProcessIDRReplicates' ], # funnel
      }, 
    },
  
  
    {
     -logic_name    => 'RunIDR',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunIDR',
     -analysis_capacity => 100,
     -rc_name    => 'normal_2GB', #~10mins + (memory?)
     -batch_size => 6,#~30mins +

     #No flow into here, but this analysis should update the DB with the calculated threshold
     #This can't be in the feature_set_stats table, due to the combinations between the reps.
     #If we accu this, then we will lose the data!!! Preventing us from being able to simply drop one rep and
     #rerun the post process?
     #do we need a stand alone table for idr_thresholds, which simply has the threshold and the 2 feature_set_ids
     #should we store anything else here? 
     
    -flow_into => 
     {
      2 => [ ':////accu?idr_peak_counts=[accu_idx]' ],
     }

           
 
      #Or should this do all this in the same analysis
      #where are we going to cache the run_idr output for the
      #final peak calling threshold? Let's keep this in a tracking DB
      #table to prevent proliferation of analyses based on this value differing between data sets.
      
      #are we going to have problems having these two analyses together
      #if we want to rerun the analysis due to the creation of the merged rset failing?
      #would have to rerun idr too?
      
      #This is extremely unlikely to happen
      #and could do some funky job_id manipulation to set skip_idr
      
      #No we won't be able to reuse DefineResultSets here
           
    },
  
    {
     -logic_name    => 'PostProcessIDRReplicates',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR',
      #-meadow        => 'LOCAL', #Could potentially do this locally due to low load/runtime?
     #although it is potentially a >1 cpu job due to triple pipe command
     -analysis_capacity => 100,
     -rc_name    => 'default', #<1mins
     -batch_size => 10,#?
     -parameters => 
      { 
       #result_set_mode              => 'none',#We never want one for an IDR DataSet
       #default_feature_set_analyses => $self->o('permissive_feature_set_analyses'), 
      },
           
     -flow_into => 
      {
       '2' => [ 'DefineMergedReplicateResultSet' ],
      }, 
      
      #Or should this do all this in the same analysis
      #where are we going to cache the run_idr output for the
      #final peak calling threshold? Let's keep this in a tracking DB
      #table to prevent proliferation of analyses based on this value differing between data sets.
      
      #are we going to have problems having these two analyses together
      #if we want to rerun the analysis due to the creation of the merged rset failing?
      #would have to rerun idr too?
      
      #This is extremely unlikely to happen
      #and could do some funky job_id manipulation to set skip_idr
      
      #No we won't be able to reuse DefineResultSets here
           
    },
 
 
 
    {-logic_name    => 'DefineMergedReplicateResultSet',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::DefineResultSets',
     -analysis_capacity => 100,
     -rc_name => 'default', #< 1 min probably need > 200MB for samtools merge!
     #Might this have been the source of the merged truncated file errors?
     #but hive should have caught this if it was exceeding the MEMLIMIT
     #how many cpu's are we using for merge?
     #This all depends on output format and whether rmdups is set
     #currently we do not rmdups and maintain bam format, so no piped commands
  
     -flow_into => { '2' => [ 'DefineMergedDataSet' ] }, 
    },
  
 
 
 
 
  
 
 
 
  #LINK ANALYSES TO OTHER CONFIGS ###
  
  #DefineMergedDataSet is in BaseSequenceAnalysis as it is common to all
  #either as a 'link out' analysis or as a 'link from' analysis.
  
  
  #We need this to run otherwise we will lose the association between
  #the IDR output (SWEmbl params) and the FeatureSet it is to be associated with
  #As such we need identical working analysis config here and in DefineOutputSets.pm
  #We could separate this and require/import it?
  #might be tricky with variable scoping
 
  #no longer needing to subsamble  
  #java -jar ~/tools/picard-tools-1.70/DownsampleSam.jar I=accepted_hits.bam P=0.01 R=42 O=sample.bam  


 


	 
  ];
}




1;
