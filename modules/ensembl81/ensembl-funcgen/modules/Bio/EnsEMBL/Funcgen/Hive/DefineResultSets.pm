
=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::DefineResultSets

=head1 DESCRIPTION

The modules creates and stores ResultSet objects given a list of input_subset_ids. Dynamic data flow is performed 
via the 'branching_by_analysis' mechanism which allows branch names to be used to data flow. This allows dataflow 
to be defined by the inputs of an analysis wrt the pre-defined analysis config.

It functions in two modes:

  1 DefineResultSets analysis
  This runs prior to any alignment analyses. For experiments which are destined for IDR analysis, 
  it creates individual replicate ResultSets, else it creates a merged replicate ResultSet. Branching by
  analysis data flows to:
    Preprocess_${alignment_analysis}_(control|merged|replicate) & PreprocessIDR (dependant on 'replciate')
    
  2 DefineMergedReplicateResultSet analysis
  This runs post IDR analysis and also merges the replicate alignments. This always dataflows to DefineMergedDataSet

=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineResultSets;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects validate_checksum);
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( merge_bams );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );


sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  $self->init_branching_by_analysis;
  
  my $isset_ids = $self->get_param_method('input_subset_ids',  'required');
  assert_ref($isset_ids, 'HASH', 'InputSubset dbIDs');
  
  $self->set_param_method('controls', delete $isset_ids->{controls});
  
  if(scalar(keys %$isset_ids) == 0){
    throw('No (signal) input_subset_ids defined. Must pass input_subset_ids hash of (optional) \'controls\''.
      ' and ResultSet name keys and input_subset_id arrayref values');
  }

  $self->get_param_method('alignment_analysis', 'required');
  #alignment_analysis is dataflowed explicitly from IdentifyAlignInputSubsets
  #to allow batch over-ride of the default/pipeline_wide value.
  #Could have put this in batch params, but it is only needed here
                  
  my $merge = $self->get_param_method('merge_idr_replicates', 'silent');

  if($merge){
    
    if( scalar(keys %$isset_ids) != 1 ){
      $self->throw_no_retry('Cannot currently specify > 1 input_subsets_ids group in merge_idr_replicate');
    }
    
    $self->get_param_method('max_peaks', 'required'); #dataflowed from PostprocessIDR
    my $ppeak_lname = $self->param_required('permissive_peaks'); #batch flown
    $self->set_param_method('idr_peak_analysis_id',
                            scalars_to_objects($self->out_db,
                                               'Analysis', 
                                               'fetch_by_logic_name',
                                               $ppeak_lname)->[0]->dbID);
  }

  #probably need to batch flow just -no_idr
  #as that will also be required in DefineMergedOutputSet? why?
  #Or can we just omit analysis to use default there?

  return;
}

#Move these to Pipeline/EFGUtils?

sub _are_controls{
  my $ctrls        = shift; 
  my $all_controls = 1;
  
  foreach my $ctrl(@$ctrls){
    
    if(! $ctrl->is_control){
      $all_controls = 0;
      last;
    }
  }
  
  return $all_controls;
}


sub _are_signals{
  my $sigs     = shift; 
  my $all_sigs = 1;
  
  foreach my $sig(@$sigs){
    
    if($sig->is_control){
      $all_sigs = 0;
      last;
    }
  }
  
  return $all_sigs;
} 


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self         = shift;
  my $helper       = $self->helper;
  my $rset_adaptor = $self->out_db->get_ResultSetAdaptor;
  my $issets       = $self->input_subset_ids;  
  my $ctrls        = [];  
  my $branch;
  
  #The branch will either be 2 for all, or may flip flop betwee 3 & for for idr and merged sets without controls
  
  my $align_lname = $self->alignment_analysis;
  my $align_anal  = &scalars_to_objects($self->out_db, 'Analysis',
                                                       'fetch_by_logic_name',
                                                       [$align_lname])->[0];
  
  #We don't use fetch_Set_input here as we are dealing with many InputSubsets
 
  my $control_branch = '';
  my $controls = $self->controls;
  
  if($controls && 
     assert_ref($controls, 'ARRAY', 'controls') &&
     @$controls){
    $ctrls = &scalars_to_objects($self->out_db, 'InputSubset',
                                                'fetch_by_dbID',
                                                $controls);
    if(! &_are_controls($ctrls)){
      throw("Found unexpected non-control InputSubsets specified as controls\n\t".
        join("\n\t", map($_->name, @$ctrls)));
    }
    
    $control_branch = 'Preprocess_'.$align_lname.'_control';
    
    
    #Test control Experiment is ALIGNED_CONTROL
    #if it has already been aligned then only submit then 
    #undef the $control_branch and things should just flow directly
    #onto the correct analyses
    
    #How are we going to handle rollback of the control alignments?
    #This would require a rollback of all dependant data!
    
    
    
    #if it is ALIGNING_CONTROL
    #Then we need to exit here
    #it would be nice to set a retry perioud for this
    #We need to output some query to enable easy polling of the 
    #DB to allow identification of when the 
    
    
    #We need to validate they are all from the same experiment, 
    #although this os done in IdentifySetInputs?
    
    my %exps;
    
    foreach my $ctrl(@$ctrls){
      $exps{$ctrl->experiment->name} = $ctrl->experiment;
    }
    
    if( scalar(keys(%exps)) != 1 ){
      throw("Failed to identify a unique control Experiment for :\n".
        join(' ', keys(%exps)));  
    }
    
    my ($exp) = values(%exps);#We only have one
    
    if ($exp->has_status('ALIGNED_CONTROL')){
        $self->helper->debug(1, 'Skipping control processing as '.$exp->name.
          ' Experiment has ALIGNED_CONTROL status');
        $control_branch = '' 
    }
    else{
      #Potential race condition here will fail on store
      
      if($exp->has_status('ALIGNING_CONTROL')){
        $self->input_job->transient_error(0); #So we don't retry  
        #Would be nice to set a retry delay of 60 mins
        throw($exp->name.' is in the ALIGNING_CONTROL state, another job may already be aligning these controls'.
          "\nPlease wait until ".$exp->name.' has the ALIGNED_CONTROL status before resubmitting this job');
      }
      
      $exp->adaptor->store_status('ALIGNING_CONTROL', $exp); 
      
      #todo check success of this in case another job has pipped us
       
    }  
  }
                                             
  
  #validate all input_subsets before we start creating anything, such that rerunning the job will be clean
  
 
  my (%rsets, %rep_bams);
       
  #setname here might not actually be set name
  #$issets is actually:
  #{input_subset_ids => {$set_name => {$replicate => [dbID1]}}
  #                       controls => have been deleted above
  
  
  #why did we need these rep keys?
  my $merge_idr_reps = $self->merge_idr_replicates;
   
  foreach my $set_name(keys %$issets){
    my $parent_set_name = $set_name.'_'.$align_anal->logic_name;
    
    my $sigs = &scalars_to_objects($self->out_db, 'InputSubset',
                                                  'fetch_by_dbID',
                                                  $issets->{$set_name});
    if(! &_are_signals($sigs)){
      throw("Found unexpected controls in signal InputSubsets\n\t".
      join("\n\t", map($_->name, @$sigs)));
    }

    my $ftype          = $sigs->[0]->feature_type;
    my $is_idr_ftype   = $self->is_idr_FeatureType($ftype);
    my $ctype          = $sigs->[0]->cell_type;
    
    #Define a single rep set with all of the InputSubsets 
    #i.e. non-IDR merged or post-IDR merged     
    my @rep_sets = ($sigs);
    my $has_reps = (scalar(@$sigs) >1) ? 1 : 0;    
    my $run_reps; 
    
    if($is_idr_ftype && $has_reps){
      
      if($merge_idr_reps){
        #Merged control file should already be present
        #Merged signal file maybe present if GeneratePseudoReplicates has been run
        $branch = 'DefineMergedDataSet';
        
        
        #Merging of fastqs is normally done in PreprocessFastqs
        #But now we need to merge bams (or other)
        #Should we move this to PreprocessAlignments?
        #This normally expects the bams to be in merged if required.
        #It seem much more natural to do the merge here
        
        #Sanity check here the control file is available?
        #Using the ResultSet below?
        
      
        $rep_bams{$parent_set_name}{rep_bams} = [];
        #my $ctrl;
        
        foreach my $rep(@$sigs){
          #Pull back the rep Rset to validate and get the alignment file for merging
          my $rep_rset_name = $parent_set_name.'_TR'.$rep->replicate;
          my $rset = $rset_adaptor->fetch_by_name($rep_rset_name);            
          #Could also fetch them with $rset_a->fetch_all_by_supporting_Sets($rep).
           
          if(! defined $rset){
            $self->throw_no_retry("Could not find ResultSet for post-IDR merged ResultSet:\t".
              $rep_rset_name); 
          }
      
          #todo validate controls are the same
          #This should already have been done in PreprocessIDR
          #but probably a good idea to do here too
          #As we may get here by means other than PreprocessIDR?
                    
          push @{$rep_bams{$parent_set_name}{rep_bams}}, 
            $self->get_alignment_files_by_ResultSet_formats($rset, ['bam'])->{bam};
        }
      }
      else{ #! $merge_idr_reps
        $run_reps = 1;
        #RunIDR semaphore handled implicitly later 
        $branch = 'Preprocess_'.$align_lname.'_replicate';       
        @rep_sets = map {[$_]} @$sigs;  #split single rep sets
      }   
    }
    else{ #single rep ID or multi-rep non-IDR
      $branch = 'Preprocess_'.$align_lname.'_merged';
      #This is a pre-alignemnt fastq merge done by PreprocessFatsqs
    }
    
    #Reset if we have controls    
    $branch = $control_branch if $control_branch;
    
    foreach my $rep_set(@rep_sets){ 
      my $rset_name = $parent_set_name;#.'_'.$align_anal->logic_name;
    
      if($is_idr_ftype && $has_reps && ! $merge_idr_reps){
        #there will be only 1 in the $rep_set 
        $rset_name .= '_TR'.$rep_set->[0]->replicate;
      }
      #else we want the parent name
    
      
      #MergeControlALignments_and_QC will also need to know which way to branch
      #This done with combination of flow_mode via the analysis config and the presence of 
      #result_set_groups in the input_id
      
      #Do not create ResultSet here as we want to validate the whole lot of 
      #input_subsets before we create any, this will help with recovery
      #Also can't flow here directly as the control and replicates branches
      #need group of ResultSets    
      #Define a cache ref to push onto
      my $cache_ref;
  
  
      #change grouping here based on $run_reps and branch?
      #This will mean we will have to break the dbIDs and set_names
      #structure, but this is fine      
      #???? What was this for? What does it mean?
  

        
      if($run_reps){  #branch can be replicate(no control) or control(with reps)  
        $rsets{$branch}->{$parent_set_name} ||= [];  
        $cache_ref                            = $rsets{$branch}->{$parent_set_name};
      }
      else{
        #branches can be
        # merged (no controls)
        # control (single rep IDR set or merged) 
        # or DefineMergedDataSet i.e. this is the IDR analysis is DefineMergedReplicateResultSet
        
        #is used of merged key here correct for DefineMergedDataSet?
        
        $rsets{$branch}{merged} ||= [];
        $cache_ref                = $rsets{$branch}{merged};
      }
      
      push @$cache_ref,
       {-RESULT_SET_NAME     => $rset_name,
        -SUPPORTING_SETS     => [@$rep_set, @$ctrls],
        -DBADAPTOR           => $self->out_db,
        -RESULT_SET_ANALYSIS => $align_anal,
        #change these to reference a specific rollback parameter
        #e.g. rollback_result_set?
        -ROLLBACK            => $self->param_silent('rollback'),
        -RECOVER             => $self->param_silent('recover'),
        -FULL_DELETE         => $self->param_silent('full_delete'),
        -CELL_TYPE           => $ctype,
        -FEATURE_TYPE        => $ftype};
    }
  }
  
  #Now do the actual ResultSet generation and cache the output_ids on the correct branch
  my %batch_params = %{$self->batch_params};
  my %branch_sets;
  my $tracking_adaptor = $self->tracking_adaptor;

  #We need to test for CONTROL_ALIGNED here
  #but this is currently only set on the ResultSets
  #and these maybe entirely new result sets
  #so this has to be set on the experiment too!
  #we can test that above and set the branch accordingly

  foreach my $branch(keys %rsets){
  
    foreach my $rset_group_name(keys %{$rsets{$branch}}){
      my $rset_group = $rsets{$branch}->{$rset_group_name};

      foreach my $rset(@$rset_group){     
        $rset = $helper->define_ResultSet(%{$rset});
        
        #How is this not failing? This is likely because there are not dependancies and 
        #everythign matches, but it should match as we should have both reps in here
        
        my %archive_files;
      
        if($merge_idr_reps){
          #Store the tracking info first!
          #This is to ensure persistance of max_peaks and permissive_peaks analysis
          #between instances of the pipeline configs.
          #what is allow_update?
          #This is to allow update of tracking info, if we have rerun
          
          $tracking_adaptor->store_tracking_info(
            $rset, 
            {allow_update => 1,
             info         => {idr_max_peaks        => $self->max_peaks,
                              idr_peak_analysis_id => $self->idr_peak_analysis_id,
             }         
            });
          
                
          #if(!exists $rep_bams{$rset->name}){#throw?}
          #Check we don't already have the file from the Generate PseudoReps step  
          #Can only do merge here, as this is the point we have access to the final ResultSet
          my $merged_file = $self->get_alignment_path_prefix_by_ResultSet($rset).'.bam'; 
      
          #This also needs to use the get_alignment_files_by_ResultSet method!
          
          #This would use and existing merged file here if Pseudo rep IDR has already created the file
          #Let's simply remove this test, and expect the file 
          #once the pseudo rep code is implemented

          # This is causing old buggy files to be used at present

          # Let's add a checsum test if it already exists
          
          
          if(! -f $merged_file || $self->param_silent('overwrite')){
            
            #In future we will treat align all replicates separately,
            #and never archive merged files
            
            
            $archive_files{to_archive} = $rep_bams{$rset->name}{rep_bams}; 
     
            
            merge_bams($merged_file, 
                       $self->sam_ref_fai($rset->cell_type->gender),
                       $rep_bams{$rset->name}{rep_bams},
                       {no_rmdups      => 1,        # Only necessary within a replicate
                        debug          => $self->debug,
                       });
                       
            $rset->adaptor->store_status('ALIGNED', $rset);
            #IMPORTED not set here, as this is used to signify whether
            #a collection file has been written.
          }
          else{
            validate_checksum($merged_file);
          }
          
          #Can't set IDR peak_analysis here, as we may lose this
          #between pipelines i.e. when using IdentifyMergedResultSets after
          #blowing away an old pipeline
        }
        
        if($branch =~ /(_merged$|^DefineMergedDataSet$)/){# (no control) job will only ever have 1 rset
          $self->branch_job_group($branch, [{%batch_params,
                                             dbID       => $rset->dbID, 
                                             set_name   => $rset->name,
                                             set_type    => 'ResultSet',
                                             %archive_files}]);
        }
        elsif($branch =~ /(_control$|_replicate$)/){
          
          $self->helper->debug(1, "Cacheing $branch branch jobs for ".$rset->name);
          
          $branch_sets{$branch}{$rset_group_name}{set_names} ||= [];
          $branch_sets{$branch}{$rset_group_name}{dbIDs}     ||= [];
          push @{$branch_sets{$branch}{$rset_group_name}{set_names}}, $rset->name;
          push @{$branch_sets{$branch}{$rset_group_name}{dbIDs}},     $rset->dbID;
                    
          if($branch =~ /_replicate$/o){
            #Just create the individual fan output_ids
            $branch_sets{$branch}{$rset_group_name}{output_ids} ||= [];
            push @{$branch_sets{$branch}{$rset_group_name}{output_ids}}, {%batch_params,
                                                                          dbID     => $rset->dbID,
                                                                          set_name => $rset->name,
                                                                          set_type => 'ResultSet'};    
          }
        }
        else{ #Sanity check
          #branch is always defined within this module
          $self->throw_no_retry("$branch is not supported by DefineResultSets");  
        }
      }
    }
  }
  
  
  #Now flow job_groups of rsets to control job/replicate & IDR 
  foreach my $branch(keys %branch_sets){            
    $self->helper->debug(1, "Processing cached branch $branch");
    #what about other branches in here
  
    if($branch =~ /_replicate$/){
      
      foreach my $rep_set(keys %{$branch_sets{$branch}}){ 
        #Add semaphore RunIDR job
        $self->branch_job_group($branch, $branch_sets{$branch}{$rep_set}{output_ids}, 'PreprocessIDR', 
                                [{%batch_params,
                                 dbIDs     => $branch_sets{$branch}{$rep_set}{dbIDs},
                                 set_names => $branch_sets{$branch}{$rep_set}{set_names},
                                 set_type  => 'ResultSet'}]);
      }           
    }       
    elsif($branch =~ /_control$/){    
      #Pick an arbitrary set for access to the controls 
      my ($any_group) = keys(%{$branch_sets{$branch}});
   
      #result_set_groups here is used by MergeControlAlignments_and_QC to flow correctly
      $self->branch_job_group($branch, 
                              [{%batch_params,
                               result_set_groups => $branch_sets{$branch},
                               set_type  => 'ResultSet',
                               dbID      => $branch_sets{$branch}{$any_group}{dbIDs}->[0],
                               set_name  => $branch_sets{$branch}{$any_group}{set_names}->[0]}]);   
    }
    #else{
 #branch not supported     
#    }
  }
  
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;