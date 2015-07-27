
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

    Bio::EnsEMBL::Funcgen::Hive::Config::DefineMergedDataSets;

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::Hive::Config::DefineMergedDataSets;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly

#todo this change to Signal/Wiggle.pm config 
#with DefineSetInputs and DefineSets as shared analyses

#We really want to be able to genericise this such that we can swap out
#any given Runnable/DB. Move all specific options to parameters sections
#specific for that analysis


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut


sub default_options {
  my ($self) = @_;
 
  #Define any optional params here as undef
  #If they are mandatory for one analysis, but optional for another
  #will have to catch that in the runnable
  
  return 
   {
    %{$self->SUPER::default_options},        
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
    %{$self->SUPER::pipeline_wide_parameters}, 
    can_DefineMergedDataSet => 1,#'DefineMergedDataSets',
    #Set to the config name so the configure_pipeline script knows to reset these jobs 
    
    #These are not strictly needed for DefineMergedDataSets
    #As this is always used in conjunction with Peaks and/or Collections
    #But need these here as we are now using MutliConfigLinker 
    #which requires these params
    #can_Link_to_Peaks       => 0,
    #can_Link_to_Collections => 0,
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
   {
    -logic_name => 'IdentifyMergedResultSets',
	  -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',	  
	  -meadow_type => 'LOCAL',#should always be uppercase
	  -parameters => {set_type        => 'ResultSet'},
	  	  
  	 #These are always created, even if we have used analysis_top_up
  	 #hence at the start of the next conf, we will get input_ids which are data flown
  	 #Do we even want to redo this?
  	 # as this forces use to funnel
  	 #when what we really want is to data flow directly from batch job to batch job
  	 #Do we need a sparse shared batch job analysis?
  	 
  	 #regardless, the way to fix this is probably to remove these from here
  	 #and use seed_pipeline after init_pipeline
  	 #seed pipeline func could automatically print out the jobs created? 
	 #can seed pipeline use the existing meta table params?
	 #So init will do the validation and seed pipeline will simply
	 
	 #This is actually fine here, as this is the very first 
	 #shared IdentifySetInputs analysis
	 #Although we should probably make this consistent across the pipelines
	 #How do we catch mistyped options?
	 #Need to put these in pipeline_wide_parameters otherwise they won't be caught?
	 #This won't give us the ability to translate between input_sets and set_names
	 #This won't matter as the specificity will be implied by the logic_name
	 #although not in the meta table
	 #Not too happy about separating the input_id definition from this config
	 
	 #Putting them in pipelinewide params offers not advantage
	 #and may prevent using seed_pipeline again, as it would require 
	 #init_pipeline to be run again, probably unsafe for various reasons
	 #including over-writing of meta params
	 #let's write our own seed_pipeline script
	 #which calls hive seed_pipeline, but we catch opts better
	 #and optionally print input_id as output(all or newly seeded)
	 #this would also handle array args!
	 
	 
	
	 -flow_into => {		 
	    #2 is used for potential fan jobs from a single result set
	    #3 is used as a funnel, or for jobs with no fan   
      '3' => [ 'DefineMergedDataSet' ],
   
      
      
      
      #IDR/Replicate Implementation
      #We will need additional semaphored data_flow to run_peak_replicates here
      #to support replicate calling and IDR QC
      #This needs to wait for all replicate sets to be peak called
      #but ignore the rest(we don't want to wait for merged replicate peaks)
      #there is no way of doing this without splitting the whole analysis tree from
      #here i.e. DefineOutputSet?!!!
      #it's clear that implementing IDR optionally along side merged peak calling 
      #in the same pipeline is going to be tricky
      #This would need to flow to another link analysis(if >1 is supported)
      #Something like Run_IDR_QC
      #That link analysis would pick up the groups of input_set_ids
      #do the IDR across the linked FeatureSets and the submit another merged 
      #peak job.
      #This would mean duplicating the peak analyses as we can't loop back due to semaphore
      #IdentifyFeatureSets(or more likely PeaksReport) would also have to be semaphored from here
      #and receive the output of all peak jobs from both originally merged and IDR jobs?
      #No, we can handle this with status entries
      #Further more PreprocessALignment should not data flow to WriteCollections if it is 
      #a replicate
      
      #How can we semaphore an analysis(PeaksReport) which doesn't exist in this config?
      #would have to do this through the link analysis
      #This would wait for everything hanging of DefineOutputSet before commencing 
      #the link analysis, which would then run the IDR and merged peaks, before doing the final 
      #PeaksReport
      #This would also wait for Collections! So we would have to detach that analysis somehow
      #or split DefineOutputSet into DefineReplicateOutputSet and DefineOutputSet
      #VERY COMPLICATED!!!
     #We don't care about these failing, as we expect them too
      },
	
	   
       -failed_job_tolerance => 100, 
			 
	  -analysis_capacity => 10,
      -rc_name => 'default',
    },
	 
	 
	  #This has to mirror conf in IDRPeaks and ReadAlignments
	  #Move this to BaseSequenceAnalysis!!!
	  
	  
	  #DefineMergedDataSet is also in BaseSequenceAnalysis 
	  #as it is common to all pipelines as a link analysis
	  #but has to be here also, as we can't define the -flow_into conf
	  #from the other confs do to lack of the PreprocessAlignments analysis
	  #All spec should match apart from -flow_into, which should be blank
	  #in BaseSequenceAnalysis, but will be updated accordinly when using
	  #-analysis_top_up
	  	  
    {-logic_name => 'DefineMergedDataSet', 
	   -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineDataSet',
	   -parameters => 
	    {#Need default_feature_set_analyses/feature_set_analysis_type here, 
	     #as DefineDataSet may want some for result_set too?
	     #although this is now currently impossible since we are removing InputSet
	     #but let's keep this generic for now
	     default_feature_set_analyses => $self->o('default_peak_analyses'),
	     feature_set_analysis_type    => 'peak',
	     check_analysis_can_run       => 1,
	     #not $self->o('peak_analysis') as we want to batch flow this
	 
	     #rollback                     => $self->o('rollback'), 
	    },
	 
	   -flow_into => 
	    {
		   1 => [ 'PreprocessAlignments' ],
		  },
		 
	   -analysis_capacity => 100,
     -rc_name           => 'default',
     -batch_size        => 10, 
     #None of these shoudl take > 2-3 mins unless there is some rolling back to do
     #But this will only ever be deleting annotated_feature records
     #and maybe some states or updating sets
     #Keep this fairly low, so we a getting the parallel compute running asap.
	   
    },
    
   
   
   
   
   
   
   
    #Can't have this as this link analysis as the flow_into config is only updated
    #if it doesn't exists
    #Therefore we are just getting the flow_into spec from the first config we top up with
    #We need a ConfigLinker
    #This will simply act to separate the flow_into config, and simply pass on the data flow
    #For this purpose, it will not need any 'can_run' functionality
    #Although we will probably want to to add this in case we want to add collections 
    #in after we have run peaks
    #Currently there is no easy way to do this. Would have to blow the pipeline away and
    #init with just the config you didn't run and reseed.
    #The only reason why we need ConfigLinker, is because we are linking to multiple downstream configs
    #So this needs to be MultiConfigLinker!
   
    
    #This needs changing to take a result_set_id and not create the sets
    #shall we merge this into the next step
    #No we still need this as this batches the jobs further
    #This will simply
    
     {
	   -logic_name => 'PreprocessAlignments',#was 'SetUp',
	   #This is basically making sure the input file is sorted wrt genomic locations
	   -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	   
	   #Need to maintain this here as will not be updated by -analysis_topup
	   -parameters => 
	    {
	     feature_formats         => ['bam', 'bed'],
	     
	     #Would be nice to do something like this in each of the topup configs
	     #However, topped up analyses don't have their parameter updated!
         #-parameters => { feature_formats => '#expr([ @{#feature_formats#}, "bam"])expr#'},
	     
	     
	     #Would really like to build this array via the topped up confs
	     #i.e. each conf will add it's own formats
	     #$self->o doesn't allow inheritance style building of arrays
	     #let's try it with # variable # notation?     
	     #peak_branches           => $self->o('peak_branches'),
	     
	     
	     #This needs to use new init_branchign_by_analysis method
	     
	    },
	    
	    #-flow_into => 
	    # {'2' => [ 'Link_to_Peaks' ],
		  #  '3' => [ 'Link_to_Collections' ],
	    # },
	  	
		
	   #No -flow_into defined as this will be handled by -analysis_topup
	   # See Peaks and Collections conf.			 	 
			 
	   -analysis_capacity => 100,
     -rc_name => 'normal_high_mem_2cpu',
	   #this really need revising as this is sorting the bed files
	   #Need to change resource to reserve tmp space
	  },
   
  ];
}




1;


__END__


=pod    

    #The input id needs to be job_groups
    #Then the MultiConfigLink can simply branch_job_groups
    #according to config
    #This means that PreprocessAlignments will have
    #to flow job_groups which match config in MultiConfigLink config
    #so there is an interdependancy here
    #There is no way around this, MultiConfigLink always need to config 
    #that the previous analysis would have if it wasn't linked
    #This is the nature of the MultiConfigLinker
    #This is odd, as not only is the PreprocessAlignments 
    #data_flow spec in a seprate downstream config but it is also 
    #in a completely separate MultiConfigLink analysis
    
        
    {
     -logic_name => 'Link_to_Peaks',
     -meadow     => 'LOCAL',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MultiConfigLinker',
     
     #Need to maintain this here as will not be updated by -analysis_topup
     #-parameters => 
     # {
     # },
     #No -flow_into defined as this will be handled by -analysis_topup
     # See Peaks  conf.         
       
     -analysis_capacity => 100,
     -rc_name => 'default',
    },
    
    
     {
     -logic_name => 'Link_to_Collections',
     -meadow     => 'LOCAL',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MultiConfigLinker',
     
     #Need to maintain this here as will not be updated by -analysis_topup
     #-parameters => 
     # {
     # },
     #No -flow_into defined as this will be handled by -analysis_topup
     # See Collections conf.         
       
     -analysis_capacity => 100,
     -rc_name => 'default',
    },
    
=cut