=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::CollectionWriter

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::CollectionWriter;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseImporter');

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names );
                                        
                                               #strip_param_args strip_param_flags run_system_cmd);
use Bio::EnsEMBL::Utils::Exception qw(throw);

#use Data::Dumper;

#global values for the Helper... maybe pass as parameters...
$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;


#params
#output_dir  This is a generic over ride for the default output dir

#todo
#1 default set up is that alignment will have already done the filter from bam
#  hence we woudl need to set the filter_from_format param if we ever want to change this
#2 Move a lot of the fetch_input_code to run


sub fetch_input {   # fetch parameters...
  my $self = shift;
  #Set some module defaults
  $self->param('disconnect_if_idle', 1);

  $self->SUPER::fetch_input;    
  $self->helper->debug(1, "CollectionWriter::fetch_input after SUPER::fetch_input");  
  my $rset = $self->fetch_Set_input('ResultSet');
  $self->helper->debug(1, "CollectionWriter::fetch_input got ResultSet:\t".$rset);
  
 
  $self->init_branching_by_analysis;  #Set up the branch config
  my $ftype_name = $rset->feature_type->name;
 

  
  
  $self->get_output_work_dir_methods($self->db_output_dir.'/result_feature/'.$rset->name, 1);#no work dir?
  #todo enable use of work dir in get_alignment_file_by_InputSet and Importer
  
  #This is required by sam_ref_fai which is called by get_alignment_file_by_InputSet
  #Move this to get_alignment_file_by_InputSet?
  $self->set_param_method('cell_type', $rset->cell_type, 'required'); 
  
  
  #TODO: need to dataflow 'prepared' status of bam file from alignment pipeline
  #can't data flow, but can set this in the alignment pipeline_wide params
  #then use this to perform the filtering in the first place
  
  
  #This should be in run as it can do some conversion!, also need to pass formats array through 
  #dependant on which analysis we are running bam and bed for Preprocess and just bed for WriteCollections
  #Currently no way to do this dynamically as we don't know what the analysis is (?)
  #so we have to ad as analysis_params
  #and there is no way of detecting what formats down stream analyses need
  #so again, they have to be hardcoded in analysis params
  #use all formats by default
  #Can we update -input_files in the Importer after we have created it?
  
  
  #Hack to add in states that were missed in DefineResultSets
  $rset->adaptor->store_status('ALIGNED', $rset);
  
  if($self->FeatureSet->analysis->program eq 'CCAT'){
    #throw('Need to implement control preprocessing for CCAT');  
    
    my $exp = $rset->experiment(1);#ctrl flag
    #But only in PreprocessAlignments no in WriteCollections 
    
    if(! $exp->has_status('CONTROL_CONVERTED_TO_BED')){
      #Make this status specific for now, just in case
    
      if($exp->has_status('CONVERTING_CONTROL_TO_BED')){
         $self->input_job->transient_error(0); #So we don't retry  
          #Would be nice to set a retry delay of 60 mins
          throw($exp->name.' is in the CONVERTING_CONTROL_TO_BED state, another job may already be converting these controls'.
            "\nPlease wait until ".$exp->name.' has the CONTROL_CONVERTED_TO_BED status before resubmitting this job');
      }
      else{
        #Potential race condition here will fail on store
        $exp->adaptor->store_status('CONVERTING_CONTROL_TO_BED', $exp); 
       
        $self->get_alignment_files_by_ResultSet_formats($rset,
                                                        ['bed'],
                                                        1); #control flag
      
        #todo check success of this in case another job has pipped us
        $exp->adaptor->store_status('CONTROL_CONVERTED_TO_BED',   $exp); 
        $exp->adaptor->revoke_status('CONVERTING_CONTROL_TO_BED', $exp, 1);#Validate status flag
      }   
    }
  }
  
  
  #This currently keeps the sam files! Which is caning the lfs quota                                                                    
  warn "Hardcoded PreprocessAlignments to get bed only, as sam intermediate was eating quota";
   
  my $align_files = $self->get_alignment_files_by_ResultSet_formats($rset, ['bed']);
                                                                #    $self->param_required('feature_formats'),
                                                                #    undef, #control flag
                                                                #    1);     #all formats
                                                                    
                                                                    
                                                                    
                                                                    

                                                                    
  #No need to check as we have all_formats defined? Shouldn't we just specify bed here?
  #other formats maybe required for other downstream analyses i.e. peak calls
  #but we don't know what formats yet
  #We have to do the conversion here, so we don't get parallel Collection slice jobs trying to do the conversion
  #todo review this
   
  #Now we need to convert the control file for CCAT
  #This is harcoded and needs revising, can't guarantee this will be grouped by control
  #at this point, so will have to employ status checking to avoid clashes
  #This analysis really needs to know the formats required for downstream analyses
  #These are available via the config, but not the PeakCaller modules themselves
  #So we will have to hardcode this in the config 
  
  #We can't test the filepaths, as 
  #1 We don't know wether they are finished
  #2 The code to build the file path is nested in get_alignment_files_by_ResultSet_formats 
  #  and get_files_by_formats. Probably need a no_convert mode, which just returns what's there
  #  already
  
  #This has to be on the experiment
  #Let's just do it for all rather than just non-IDR ftypes??
  
  
 
  
  
  
  $self->helper->debug(1, 'CollectionWriter::fetch_input setting new_importer_params with align_files:', 
                       $align_files);
  #Arguably, some of the should be set in default_importer_params
  #but this does the same job, and prevents having to flow two separate hashes
  


    #Default params, will not over-write new_importer_param
    #which have been dataflowed
  $self->new_Importer_params( 
   {#-prepared            => 1, #Will be if derived from BAM in get_alignment_file_by_InputSets
    #but we aren't extracting the slice names in this process yet
    -input_feature_class => 'result', #todo remove this requirement as we have output_set?
    -format              => 'SEQUENCING',#This needs changing to different types of seq
    -recover             => 1, 
    -force               => 1, #for store_window_bins_by_Slice_Parser
    -parser              => 'Bed',
    -output_set          => $rset,
    -input_files         => [$align_files->{bed}], #Can we set these after init?
    -slices              => &generate_slices_from_names
                              ($self->out_db->dnadb->get_SliceAdaptor, 
                               $self->slices, 
                               $self->skip_slices, 
                               'toplevel', 
                               0, 1),#nonref, incdups
   });
  
  #todo why do we have to set EFG_DATA?
  $ENV{EFG_DATA} = $self->output_dir;
  
  #This picks up the defaults param from the input_id
  my $imp = $self->get_Importer;
  
  if((! $self->param_silent('merge')) &&
    ($imp->prepared)){
    $self->param_required('slices');     
  }
  
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift;
  my $Imp  = $self->get_Importer;
    
  if(! $self->param_silent('merge')){ #Prepare or write slice col
   
    if( ! $Imp->prepared ){  #Preparing data...
      $Imp->read_and_import_data('prepare');
      #Dataflow only jobs for slices that the import has seen
      my %seen_slices  = %{$Imp->slice_cache};
      
      if(! %seen_slices){
        $self->throw_no_retry("Found no data after prepare step for:\n\t".
          join("\n\t", @{$Imp->input_files}));  
      }
      
      my $batch_params = $self->batch_params;
      my @slice_jobs = ();

      foreach my $slice (values %seen_slices){ 

        push @slice_jobs,
         {
          %{$batch_params},#-slices will be over-ridden by new_importer_params
          dbID           => $self->param('dbID'),
          set_name       => $self->param('set_name'),#mainly for readability
          set_type       => $self->param('set_type'),
          slices         => [$slice->seq_region_name],
          filter_from_format => undef, #Don't want to re-filter for subsequent jobs
          feature_formats    => ['bed'], #also don't need bam for Collection job
          new_importer_params => 
           {
            -input_files    => [$Imp->output_file],
            -total_features => $Imp->counts('total_features'),
            -batch_job      => 1, 
	          -prepared       => 1,  
           }
         };
      }
      
                  
      my $output_id = {%{$batch_params},#-slices will be over-ridden below
                       dbID         => $self->param('dbID'),
                       set_name     => $self->param('set_name'),#mainly for readability
                       set_type     => $self->param('set_type'),
                       filter_from_format => undef,
                       slices             => [keys %seen_slices]
                       
      }; #Don't want to re-filter for subsequent jobs
                
      #Deal with Collections and Mergecollection job group first
      #It's actually not necessary to branch like this for this analysis
      #as it is only dealing with one job group
      #so could have simply branced the funnel jobs after
      #but this makes it more explicit  
      
      #This merge job delete the input for the peak job, this could conceivable do this
      #before it get a chance to run. V unlikly tho, du eto runtime.
      #We really need to ad ad a final semaphored job, to do the clean up
      
      $self->branch_job_group(2, \@slice_jobs, 1, [{%$output_id, garbage => [$Imp->output_file]}]);
  
      #We have no way of knowing whether method was injected by fetch_Set_input
      my $fset = $self->FeatureSet;#from fetch_Set_input
      
      #todo
      #Do we need to be able to do this conditionally?
      #If we are re-running then data flowing here will create duplicate jobs
      #on any branch which has run successfully before
      #this is fine if the input_id is the same and we are using the same hive
      #as it will not create the job
      #But if we have reseeded with slightly different batch_params
      #or we are using a new hive, then this will try to rerun
      #the next analyses
      
      #Maybe we don't even won't to run the peaks, even though we have a feature set?
      #
      #see pipeline dev notes
      
      if(defined $fset){
        
        #Only flow slices to peak jobs if we have defined -slices/skip_slices   
        if(! $self->slices){ 
          #Deref $output_id so we don't inadvertantly update above reference.
          $output_id = {%$output_id};
          delete $output_id->{slices};
        }
    
        #what about peaks reports?
        #we should really flow the slices to that analysis!   
        
        $self->branch_job_group('run_'.$fset->analysis->logic_name, [{%$output_id}]);
      }
    } 
    else { #These are the fanned slice jobs  
   
      #This generates a slice .col file
      $Imp->read_and_import_data;

      #This is currently skipping tons of empty lines?!!!      
      
      #End of this branch, merge job is waiting for these
    } 
  }
  else{ #Merge
    warn "HARDCODING force_overwrite for merge_and_index_slice_collections!";
    
    $Imp->feature_adaptor->merge_and_index_slice_collections($Imp->output_set, 
                                                             $Imp->slices, 
                                                             $self->output_dir,
                                                             1);   
    #set appropriate states...
    my $rset = $self->ResultSet;     
    $rset->adaptor->set_imported_states_by_Set($rset);  
    #End of this branch & config
  }

  return;
}


sub write_output { 
  my $self = shift;    
  $self->dataflow_job_groups;
  return; 
}



1;
