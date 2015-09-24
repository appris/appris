=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::DefineDataSet

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineDataSet;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects 
                                               get_set_prefix_from_Set );

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

#This assumes the ResultSet has been previously registered,

#params
#set_type, name, dbID
#These allow over-ride of defaults
#feature_set_analysis
#result_set_analysis
#default_feature_set_analyses (shouldn't alter)
#default_result_set_analyses (shouldn't alter)

#todo
#-slice_import_status?
#This now needs to take lists of dbIDs for InputSubsets
#These should all share the same control
#The InputSets will be created iteratively
#but the data flow will be done in a batch fashion, passing on sets of InputSet IDs
#as appropriate
#change this analysis to DefineSets? As it deals with InputSets, FeatureSets, ResultSet and DataSets
#or should we just create another analysis which is DefineInputSets?
#as it is almost entirely different?
#refactor this default_analysis code in BaseDB? As this is reused in other runnables


sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  $self->check_analysis_can_run; 
  $self->SUPER::fetch_input;  
  my $set = $self->fetch_Set_input('ResultSet');#defines result_set method
  my %default_analyses;


  #Can we move some of this into BaseSequenceAnalysis as it will need to be used
  #by Run_QC_and_ALigner to enable selective data flow into the relevant 
  #IdentifyInputSets analyses
  
  my $anal_type = $self->param_required('feature_set_analysis_type'); 
  my $set_lname = $self->get_param_method($anal_type.'_analysis', 'silent');
  
  
  ### Set post IDR FeatureSet analysis 
  if($self->is_idr_ResultSet($set)){
    #This is highly dependant on the fact we never make a DataSet for the IDR replicate peaks (only a ResultSet)
    #If this ever changes we could use a 'merged' param here or the analysis_name (DefineMergtedDataSet) as a proxy
    
    #We may not have had peak_analysis batch flown if we are starting from IdentifyMergedResultSets
    #although if defined it should match permissive_peaks here
    #This would highlights a potential conflict between no_idr and peak_analysis
    #Should peak_analysis turn on no_idr?
    #peak_analysis, is generally for non-idr data sets
    #and permissive_peaks overridesthis
    #IDR mode currently only supports SWEMBL, so redefining permissive_peaks
    #will fail in some unkown way at present?       
    
    # This is currently causing failures with only Collection config loaded
    # and will still cause failure even if we don't set peak analysis?
    # How do we re-use an existing data set?
    # The pipeline actually needs restructuring
    # as there is no need for a dataset if we are just generating bigwigs(result_set)

    # flow to unwired analysis is okay, so let's just code around this for now.

    # $set_lname = $self->get_param_method('permissive_peaks', 'required').'_IDR';  
    # $self->peak_analysis($set_lname);
    # $set_lname = $self->get_param_method('permissive_peaks', 'required').'_IDR';  

    if(my $ppeak_anal = $self->get_param_method('permissive_peaks') ){ 
      #Reset peak_analysis here for clarity as this is batch flown
      #although we woudl still use the above to create the feature_set
      #$self->peak_analysis($ppeak_anal.'_IDR');
      $set_lname = $ppeak_anal.'_IDR';
    }

    #This module has been written so that it is agnostic towards to the type of feature_set
    #it is creating, so putting permissive_peaks in here spoils that at present   
    #i.e. we pass the default_peak_analyses as the default_feature_set_analyses
  }
    
    
  if(! defined $set_lname){
    
    $default_analyses{feature_set} = $self->param_silent('default_feature_set_analyses');     
    
    if(! defined $default_analyses{feature_set}){
      throw("Please define -feature_set_analysis or add to default_feature_set_analyses".
        " in the analysis config");
    }
    
    if(exists $default_analyses{feature_set}{$set->feature_type->name}){
      $set_lname = $default_analyses{feature_set}{$set->feature_type->name}; 
    }
    elsif(exists $default_analyses{feature_set}{$set->feature_type->class}){
      $set_lname = $default_analyses{feature_set}{$set->feature_type->class};
    }
    else{
      throw("No default feature_set analysis available for ".$set->feature_type->name.
        '('.$set->feature_type->class.
        ").\n Please add FeatureType name or class to  default_feature_set_analyses ".
        "in analysis config or specify -feature_set_analysis"); 
    }
  }                                  
    
  #Catch undefs in config
  if(! defined $set_lname){
    throw('Unable to identify defined feature_set analysis in config for '.$set->feature_type->name.
    ".\nPlease define in -default_feature_set_analyses config or specify -feature_set_analysis");     
  }    
  #can't use process_params here as the param name does match the object name
  #We could over-ride this with a second arrayref of class names     
                                   
  $self->param('feature_set_analysis', 
               &scalars_to_objects($self->out_db,
                                   'Analysis', 
                                   'fetch_by_logic_name',
                                   [$set_lname])->[0]);                                              
  return;
}


#TODO
#1 Migrate this to the Importer as define_OutputSetthere is some overlap 
#  there of param validation between BaseImporter and hive 
#2 Review whether rollback handles status entries correctly
#3 Review whether we need to store some tracking info here
#4 If the dataflow ever changes from this, we will need to use branching_by_analysis
  
sub run {   # Check parameters and do appropriate database/file operations... 
  my $self   = shift;
  my $helper = $self->helper;  
  my $rset   = $self->ResultSet;
  my $fset_anal = $self->param('feature_set_analysis');
  my $set_prefix = get_set_prefix_from_Set($rset);
  my $set;    
    
  # Never set -FULL_DELETE here!
  # It is unwise to do this in a pipeline and should be handled
  # on a case by case basis using a separate rollback script    
  # Should also never really specify recover here either?
  # This bascailly ignores the fact that a ResultSet may be linked to other DataSets
  $set = $helper->define_DataSet
   (-NAME                 => $set_prefix.'_'.$fset_anal->logic_name,
    -FEATURE_CLASS        => 'annotated', #Is there overlap with rset feature_class here?
    -SUPPORTING_SETS      => [$rset],
    -DBADAPTOR            => $self->out_db,
    -FEATURE_SET_ANALYSIS => $fset_anal,
    -RESULT_SET_ANALYSIS  => $self->param('result_set_analysis'),
    -RESULT_SET_MODE      => $self->param('result_set_mode'),
    -ROLLBACK             => $self->param('rollback'),
    -RECOVER              => $self->param('recover'),
    -SLICES               => $self->slices,    
    -CELL_TYPE            => $rset->cell_type,
    -FEATURE_TYPE         => $rset->feature_type                    );  

  $self->param('output_id', 
               {%{$self->batch_params},
                dbID       => $set->dbID, 
                set_name   => $set->name,
                set_type   => 'DataSet'  });
  return ;
}


sub write_output {  # Create the relevant jobs
  my $self = shift;
  $self->helper->debug(1, 'DefineDataSet data flowing:', $self->param('output_id'));
  $self->dataflow_output_id($self->param('output_id'), 1);
  return;
}


1;
