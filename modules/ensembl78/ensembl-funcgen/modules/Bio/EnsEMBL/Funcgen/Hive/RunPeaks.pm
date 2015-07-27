
=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::RunPeaks

=head1 DESCRIPTION

This class runs a PeakCaller analysis and loads the results into the DB

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunPeaks;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( _init_peak_caller 
                                                    _run_peak_caller );                                               
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );



#todo rename logic_names to be generic e.g. SWEmbl_tight
#then we can have different parameter sets between species

#todo Move more to PeakCaller so it is available from SeqTools
#Or move directly to SeqTools? Max peaks filtering. write_output and store_AnnotatedFeature?
#

#input_dir and work_dir can be separate
#such that we can't point to remote fastqs
#input_dir should never be altered
#where as work_dir can be written to and cleaned up afterwards
#is this fully supported?

#params
#-output dir should never be dataflowed! ANd is really only safe when running one analysis
#as it will send all output there
#-reload



#Make this optionally take a ResultSet and a peak analysis
#to support calling without generation of Data/FeatureSet
#if there is no feature set, we should also not load the peaks
#if we define a max peaks value, then we rename to first output to unfiltered,
#then add a filter step to the final expected output

sub fetch_input {
  my $self = shift;
  #Do some thing before we SUPER::fetch_input
  $self->param('disconnect_if_idle', 1);
  $self->check_analysis_can_run;
  $self->SUPER::fetch_input;

  my $set_type = $self->param_required('set_type');  
  my $pftypes  = $self->get_param_method('process_file_types', 'silent', [undef]);
  assert_ref($pftypes, 'ARRAY', 'process_file_types');
  
  if(scalar(@$pftypes) != 1){
    #Single undef value, to enable single iteration of loops
    #below where we want default value of $peak_caller->out_file($file_type);
    #See write_output for limitations of current support
    $self->throw_no_retry("RunPeaks currently only supports a single process_file_type:\n\t".
      'process_file_type => ["'.join('", "', @$pftypes).'"]');  
  }
  
  #filter_max_peaks ensures validation of max_peaks variables for known IDR analyses
  #max_peaks also needs grabbing from result_set_tracking and dataflow from IdentifySetInputs
  
  my $max_peaks = $self->get_param_method('max_peaks', 'silent');
  
  if($self->param_silent('filter_max_peaks') &&
     (! defined $max_peaks)){
    $self->throw_no_retry('The filter_max_peaks param has been set, but no max_peaks param has been set');
  }
  
  my ($fset, $rset, $analysis);

  if($set_type eq 'ResultSet'){
    $rset = $self->fetch_Set_input('ResultSet'); 
    
    #This is likely permissive peaks for pre_IDR rep 
    
    
    #Can we not auto detect based on run_idr and is_idr_feature_type
    #No this is breaking the link between setting the analysis and the branching
    #which is based on the FeatureSet logic_name
    #This branch config is not loaded in IdentifyReplicateResultSet
    #How are we goign to do this, such that we don't risk 
    #passing the wrong analysis! i.e. It has to match the hive analysis name
    #
    
    my $peak_analysis = $self->param_required('peak_analysis');
    $analysis = scalars_to_objects($self->out_db, 'Analysis',
                                                  'fetch_by_logic_name',
                                                  [$peak_analysis])->[0];
    if(! defined $analysis){
      $self->throw_no_retry("Could not find peak_analysis in DB:\t".$peak_analysis);  
    }                            
  }
  else{ 
    $fset     = $self->fetch_Set_input('FeatureSet');
    $analysis = $fset->analysis;
    $rset     = $self->ResultSet; 
  }

  #do we need both experiment name and logic name in here?
  #logic name will be in the otufile name anyway?

  $self->get_output_work_dir_methods( $self->db_output_dir . '/peaks/' .
      $rset->experiment->name. '/' . $analysis->logic_name );

  my $align_prefix   = $self->get_alignment_path_prefix_by_ResultSet($rset, undef, 1);#validate aligned flag 
  my $control_prefix = $self->get_alignment_path_prefix_by_ResultSet($rset, 1, 1);#and control flag 
  my $sam_ref_fai = $self->sam_ref_fai($rset->cell_type->gender);  #Just in case we need to convert

  #These maybe things like extra input/reference files
  #where we don't want to store the filepath in the DB.
  my $sensitive_caller_params = $self->param_silent($analysis->program.'_parameters') || {};
  my $pfile_path = ( defined $self->bin_dir ) ?
    $self->bin_dir.'/'.$analysis->program_file : $analysis->program_file;

  my $peak_runnable = _init_peak_caller
   (-analysis => $analysis,
    -align_prefix => $align_prefix,
    -control_prefix => $control_prefix,
    -sam_ref_fai => $sam_ref_fai,
    -debug             => $self->debug,
    -peak_module_params =>
     {%$sensitive_caller_params,
      -program_file      => $pfile_path,
      -out_file_prefix   => $rset->name.'.'.$analysis->logic_name,
      -out_dir           => $self->output_dir,
      -convert_half_open => 1,
      -is_half_open      => $self->param_silent('is_half_open') || 0},    #default to closed coords
   );
   
  $self->set_param_method( 'peak_runnable', $peak_runnable );

  return;
} ## end sub fetch_input



#out_file may need a file_type spec here
#There is only partial support here for multiple output
#process_file_types array is specified in the analysis config
#all would be post processed uniformly, although this is only specified 
#for CCAT at present which has no $max_peaks specified
#write_output also lacks the ability to specify the extra FeatureSets required
#for the outher output files. Hence we only expect 1 process_file_type there

#mv max peak filtering to write_output?
#This seems more appropriate now, given that we need to do file/header manipulation
#for filtering. And sort if format specific
#Do we want to be able filter in run, without 'writing' output

sub run {
  my $self      = shift;
  my $max_peaks = $self->max_peaks;
   
  if( ! $self->param_silent('reload')){
    my $pcaller = $self->peak_runnable;

    if(! eval { _run_peak_caller(-peak_caller => $pcaller, 
                                 -max_peaks   => $self->max_peaks, 
                                 -file_types  => $self->process_file_types,
                                 -debug       => $self->debug); 1 }){
      $self->throw_no_retry('Failed to call run on '.ref($pcaller)."\n$@");                                
    }
  }
  
  return;
}


#TODO
# 1 Log counts here, in tracking DB
# 2 Dataflow to PeaksQC from here of PreprocessAlignments
# 3 Add in optional QC and PeaksReport
# 4 Handle >1 output files/formats e.g. CCAT significant.region significant.peak
#  This is already handle in run, but there is no way to handle/specify the params
#  for the process method i.e. extra feature sets?


# Move the bulk of this to SeqTools?
# This would require a pre-registered FeatureSets
# unless we import some of the define_sets code in there
# Will SeqTools ever be able to do an Feature/ResultSet import based on command line params?
# Is this overkill?

sub write_output {
  my $self = shift;
  my $fset;
  
  #Move this test to fetch_input based on -no_write status?
  #Is -no_write available to the Process? I htink not, as we have had to specify it in the input_id
  #for IdentifySetInputs
  
  
  if($self->can('FeatureSet') &&
     ($fset = $self->FeatureSet) ){
    #test assignment, as we may have the FeatureSet method from a previous
    #job in this batch     

    if ( $fset->has_status('IMPORTED') ) {
      throw( "Cannot imported feature into a \'IMPORTED\' FeatureSet:\t" .
             $fset->name );
    }
    else{ #Rollback
      #Just in case we have some duplicate records from a previously failed job
      $self->helper->rollback_FeatureSet($fset);
    }
    
    my $af_adaptor = $fset->adaptor->db->get_AnnotatedFeatureAdaptor;    
    my $params     = 
     {-file_type      => $self->process_file_types->[0],
      -processor_ref  => $self->can('store_AnnotatedFeature'),
      -processor_args => [$self,#Need to pass self as this calls a code ref
                          $fset,
                          $af_adaptor]};  

    #Need list context otherwise would get last value of list returned
    #i.e. $retvals, which we don't need here
    my ($feature_cnt) = $self->peak_runnable->process_features($params);
    
    #As CCAT is very oddly calling duplicate regions we need to skip over these in 
    #store_AnnotatedFeature
    #Hence we need to validate the feature_cnt returned matches a the DB count 
    #else throw here, rather than in store_AnnotatedFeature
    #so we can finish the load.
    #Luckily there is now flow from the CCAT analysis, so we can review these and accept them
    #throw before or after imported states?
                        
    $self->helper->debug(1, "Processed $feature_cnt features for ".$fset->name.'('.$fset->dbID.')');
    my $stored_features = $af_adaptor->generic_count('af.feature_set_id='.$fset->dbID);
    
    #Arguably this should be after the follwing test
    $fset->adaptor->set_imported_states_by_Set($fset);
    
    if($feature_cnt != $stored_features){
      #coudl change to a normal throw if we move above the satus setting
      $self->throw_no_retry('Processed feature count does not match stored '.
        "feature count:\t$feature_cnt vs $stored_features");
    }  
      
    #This is currently only setting IMPORTED_GRCh38, not IMPORTED too

    # Log counts here in tracking?
    #No data flow to PeaksQC here as this is done via semaphore from PreprocessAlignment?
  }
  else{
    warn "Failed load features as no FeatureSet is defined";
    #This is for the IDR replicate data, can we omit this error if we know that  
  }

  return;
} ## end sub write_output

#Should we attempt to skip features which overhang end of slices?
#This may be a 'feature' of the peak caller?

sub store_AnnotatedFeature {
  my ( $self, $fset, $af_adaptor, $fhash ) = @_;
  my $err;  

  if ( my $slice = $self->get_Slice( $fhash->{-seq_region} ) ) {
    delete ${$fhash}{-seq_region};
    
    if(! eval {$af_adaptor->store( Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                                    (%$fhash, -slice => $slice, -feature_set => $fset)); 1}){
      $err = 'Could not create and store ' .
        $fset->name . " AnnotatedFeature with attributes:\n\t" .
        join( "\n\t", ( map { "$_ => " . $fhash->{$_} } keys %$fhash ) ).
        "\n\t-slice => $slice\n\t-feature_set => $fset\n$@";
    }
  }
  else{
    $err = "Failed to get slice ".$fhash->{-seq_region};
  }

  warn $err if $err;
  return $err;
}

1;


