=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::BaseImporter

=head1 DESCRIPTION

BaseImporter is a base class to support those modules which use the Importer class

=cut

package Bio::EnsEMBL::Funcgen::Hive::BaseImporter;

use warnings;
use strict;

use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');


#This defines a set of parameters based on given parameters from the pipeline:

#Removed fetch_input here as sub class normally requires access to out_db to grab
#new_importer_params before calling get_Importer


#add in param_required for mandatory config?
#or do this in subclasses

#sub fetch_input {
#  my $self = $_[0];
#  $self->SUPER::fetch_input;

#  #Add a call to get_Importer here, so we initialise it before
#  #we do anythign else, no essential 

#  $self->get_Importer;
#  return;
#}

sub get_Importer_params{
  my $self = $_[0]; 
  my %params;
  
  #All optional params are defined in default_importer_params as undef. 
  #This will also allow over-ride defaults from init_pipeline
  #based on an input_id

  #my $default_params = $self->param_required('default_importer_params');
 
 
  #This undef params need to be removed or exposed
  
  my %default_params =  
    (
     #All optional or unknowable param are set here now
     #for population by get_Importer_params later
     -db         => undef, 
     -species    => $self->species,#This should be lc already?
     -assembly   => $self->assembly, 
     -parser     => $self->get_param_method('feature_file_format'), 
     -output_dir => $self->output_dir,
     #we really want to send the intermediate files to work_dir
     #but we don't have support for this in the Importer
 
 
     #will these always be the same for a set of ids?
     #this will only be called once during init pipeline!
     #hence can't have defaults here
     #-feature_analysis => $self->param('feature_analysis'),#todo is this feature_set and/or result_set?
     #-input_feature_class => $self->o('input_feature_class'),
      
     #do we even need this registry stuff?
      
     #unfortunately yes, the registry is used to validate species
     #alias, so we are using the 'homo_sapiens' rather than say 'Human'
     #this is probably overkill in the BaseImporter
     #and most likely not required now the Funcgen DBAdaptor
     #automatically sets the species if it is present in the
     #meta table?
      
     #enable resitry param here, but we don't use them in the config
     #as the DB may not have the correct species name, hence
     #all adaptors are fetch directly from the DBAdaptor rather than using
     -registry_host    => $self->param_silent('registry_host'),
     -registry_port    => $self->param_silent('registry_port'),
     -registry_user    => $self->param_silent('registry_user'),
     -registry_version => $self->param_silent('registry_version'),
     -registry_pass    => $self->param_silent('registry_pass'),

     
     #These are generic, but may get redefined
     #todo these are not defined at present
     #-input_dir   => $self->o('root_input_dir'),
     #-output_dir  => $self->o('root_output_dir'),
     #-data_dir    => $self->o('data_dir'),#?

     -recover     => $self->param_silent('recover'),
     -verbose     => $self->param_silent('verbose'),
     -ssh         => $self->param_silent('ssh'),
     
      
     #you can only add a slice at a time... see if list of values can be passed as params...
     #-slices => undef,
     #no point in setting this to undef here 
       
        
     #will this handle []? #if so we can set this to $self->o('slices')
     #as this is defined in default_options
     #no point in setting here, unless we want to rerun only 1 slice for all of the data sets?
     #what about skip slices?
      
      
     #optional but not pipeline_wide
     #these will be populated by get_Importer_params
     #-ucsc_coords => undef,  
     #-name        => undef,
     #-format      => undef,
     #-vendor      => undef,
     #-group       => undef,
      
     #a lot of this has been replaced by -output_set in new_Importer_params 
     #todo handle the param name mapping in the data flow!
     #todo can we pass object here for some of these?
     #-feature_analysis  => #Don't need this now as we already have an output_set
     -input_set_name      => undef,# $self->param('input_set'),
     #-input_feature_class => undef,
     -result_set_name     => undef, #$self->param('input_set'), #not implemented yet
     -feature_type_name   => undef,#$self->param('feature_type'),
     #-feature_analysis    => undef,
     -cell_type_name      => undef, #$self->param('cell_type'),
  
      
     #-batch_job      => undef,
     #-prepared       => undef,
     -input_files    => undef, #[ $self->param('result_file') ],
     #-total_features => undef, 

     #Now done in caller as slice_object
     #does not inc dups? 
     #-slices         => $self->slice_objects,
     
     #-skip_slices    => $self->skip_slices,         
     #-force          => undef,
    
    
    
     #result_files here pertain to sam files, but actually we just want fastq references i.e. input_subsets
     #These are already known from the DB (not for register mode)
  
     #todo remove date from experiment
     #-exp_date not currently used 
     #-farm NEVER SET THIS AS HIVE IS NOW HANDLING BSUBING  
    );
  
    
  
  $self->helper->debug(3, "BaseImporter::get_Importer_params defaults:", \%default_params);

  #Can't we just pass them all explicitly through new_importer_params
  #This prevents having to account for all possible params
  
  #Only pass those which are defined
  #We don't actually need to set any as undef in defaults hash!
   
  foreach my $param(keys %default_params){  
    my $value = $default_params{$param};
    
    #Any over-riding will be done via new_importer_params
       
    if(defined $value){
      $self->helper->debug(1, "Importer param:\t${param} => $value"); 
      $params{${param}} = $value; 
    }   
  }
  
  #new dataflown params will over-write input_id & default params
  return {-db => $self->out_db, %params, %{$self->new_Importer_params}};
}


sub get_Importer {
  my $self = $_[0];
   $self->helper->debug(1, "BaseImporter::get_Importer");
  
  #TODO appears not to be using any dnadb params from env
  #dnadb setting reverts to ensembldb, when it shoudl revert to the registry host?
  #Is this a fix for the BaseImporter or the DBAdaptor?
      
  #todo Bring in Helper $main params here?
  #need to do some hierarchichal default_Helper_params config?   
 
  
  if(! defined $self->param_silent('importer')){ 
    my $imp_params = $self->get_Importer_params;
    $self->helper->debug(2, "BaseImporter::get_Importer params are:", $imp_params);
    my $Imp = Bio::EnsEMBL::Funcgen::Importer->new(%{$imp_params});
          
    if(! defined $Imp){ 
      throw("Could not create Importer"); 
    }
    
    $self->param('importer', $Imp);
  }
  
  return $self->param('importer');
}

#This over-writes default params passed with new_importer_params 
#which have been dataflowed, but does not undef existing values
#todo pass the defaults to get_Importer instead of here
#get_Importer would then pass them to get_Importer_params
#Then we wouldn't need this method at all?
#we would still only data flow new_importer_params
#so defaults from this analysis would not get integrated
#hence this would be more clean if we are to pass them onto another Importer 
#analysis with different defaults

sub new_Importer_params {
  my ($self, $default_params) = @_;
  
  #These will be the dataflown_params
  my $new_params = $self->param_silent('new_importer_params') || {};
  
  $self->helper->debug(3, "BaseImporter::new_Importer_params param_silent:", $new_params);
  
  if(defined $default_params){
    
    if( (! ref($default_params)) ||
        ref($default_params) ne 'HASH') {
      throw('New importer params arg must be a Hashref');     
    }
    
    $new_params = { %$default_params, %$new_params,  }; 
    $self->param('new_importer_params', $new_params);
  }
  
  $self->helper->debug(3, "BaseImporter::new_Importer_params updated:", $new_params);
   
  return $new_params;
}



1;
