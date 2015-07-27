
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB;

=head1 SYNOPSIS



=head1 DESCRIPTION

 

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB;

use strict;
use warnings;
use Data::Dumper;
use base qw(Bio::EnsEMBL::Funcgen::Hive::Config::Base);
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of 
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used
                  to initialize default options.
    Caller      : DependantOptions::process_options (init_pipeline.pl)

=cut

#These are used by init_pipeline.pl and according to the docs should the pipelinedb spec
#Anything set (e.g. test) in default_options will be accessible in other methods here 
#(e.g pipeline_wide_params, pipeline_create_commands) as $self->o(test), 
#unless specified on the cmdline in which case that will be the value of $self->o(test);
#This is basically mechanism to remove the needs for doing soemthing like this f
#or defaults in pipeline_wide_params: test => $self->o('test') || default_value
#or: test => (defined $self->o('test')) ? $self->o('test').'/ing' : '/default/ing'
#$self->o also captures undefined cmdline params, if they have not been set in defaults

#There seems to be no clear reason for not setting params in here, if we never want to set 
#a hardcoded default for them which we might want to over-write with a specific cmdlind opt
#in pipeline_wide_parameters

#It looks like it is possible to refer to a previously defined option within the same hash

sub default_options {
  my $self = $_[0];  
  
  return {
    %{$self->SUPER::default_options},
    
    #Registry params are optional
    #It's actually a bit unsafe to use the registry in pipeline code
    #as the DB may not have the correct species name, hence
    #all adaptors are fetch directly from the DBAdaptor rather than using
    #registry_host       => undef, 
	#registry_port       => undef, 
	#registry_user       => undef, 
	#registry_version    => undef,
	#registry_pass       => undef,
	
	#Have to access optional (undef) ENV params here
	#These must be set in the environment
	#Access to these from pipeline_wide_parameters is currently broken
	#passwords via env for security
	dnadb_pass          => $self->o('ENV', 'DNADB_PASS'),
	pass                => $self->o('ENV', 'DB_PASS'),
	dnadb_port          => undef,
	
	port                => undef,
	
	#Best to do this on a per analysis basis
	#disconnect_when_idle => undef, #Set on funcgen and core DBAdaptors
	ssh                 => undef, #Connect to DBs using ssh(use in Importer)
	
	
	### Optional Helper param (currently used in DefineOutputSet)
    result_set_only    => 0, #why is this 0 rather than undef?
    #todo remove this and use param_silent as this is really a batch_param?
    
 
    
   };
  
}

#Generic resource_classes are in Base



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
    my ($self) = @_;
    
    #Deal with batch_params first as this may not have been subtituted yet
    #and will evaluate to a string:
    #Can't use string ("#:subst batch_params:#") as an ARRAY ref while "strict refs"
     
   # my $batch_params = 
   #   [ ( ref($self->o('batch_params')) ? 
   #       @{$self->o('batch_params')} : () ),
   #     #Generic optional params (used in Helper and elsewhere)
    #    'rollback',
    #    'slices',
    #    'skip_slices',
    #    ### More optional Helper params (currently used in DefineOutputSet)
#    #    'result_set_only',
#        'result_set_mode', #Can only be set to recover at present
#        'recover',         #is this an Importer or a Helper param?
        
        ### Optional IdentifySetInputs parameters
        #comma separated or defined as list ref  
        #are these batch wide or just used for the seed job?
        #qw( cell_types feature_types input_sets input_set_ids
        #    experimental_groups states ),             
        #'input_analyses'     => undef,
        #'experiment'   => undef, 
        #This is actually more like study, and omit for now as 
        #input_set will handle this         
 #     ];
                            
    return {
      %{$self->SUPER::pipeline_wide_parameters}, 
               
       #todo make this optional to the extent we only specify a host
       #and let the API do the rest  
      dnadb   => 
        {
         -dnadb_host   => $self->o('dnadb_host'),
         -dnadb_pass   => $self->o('dnadb_pass'),
         -dnadb_port   => $self->o('dnadb_port'),
         -dnadb_user   => $self->o('dnadb_user'),
         -dnadb_name   => $self->o('dnadb_name'),
        },
        
      #todo CR rename to output DB we are no longer efg!
      #change the cmdline prefixes to outdb_host etc, to match the config param
      
      #do not rename this 'db', as the param wrapper method will 
      #clash with Process::db method  
      
      #These need to be made optional in Base?
      
      out_db  => 
        {
         -host   => $self->o('host'),
         -port   => $self->o('port'),
         -user   => $self->o('user'),
         -pass   => $self->o('pass'),
         -dbname => $self->o('dbname'),
        },
     
     #Can't pass optional params via the ENV as this
     #will always barf when it is absent!
     #Just define as undef in the env!
     
     
     
      #It's actually a bit unsafe to use the registry in pipeline code
      #as the DB may not have the correct species name, hence
      #all adaptors are fetch directly from the DBAdaptor rather than using
      #the Registry.
      #registry_host    => $self->o('registry_host'), 
	  #registry_port    => $self->o('registry_port'), 
	  #registry_user    => $self->o('registry_user'), 
	  #registry_version => $self->o('registry_version'),
	  #registry_user    => $self->o('registry_user'),

    
      #Now defaults in runnable
      #db_output_dir => $self->o('data_root_dir').'/output/'.$self->o('dbname'),
     
     #Optional params
     
     #These are not pipelinewide!!!! They only refer to s specific run!
     #todo move these to specific analyses requiring these
     #these should never be specified with init!
 
     #todo put in other options here
     #then validate mandatory aspects in analysis   
        
  
    #Maintain this here for now in case BaseDB is used by something other than BaseSequenceAnalysis
    #Could change this if we move the releveant code from the respective modules
  
     batch_param_names => 
      [
       #From Base.pm       
       'no_write', #For use with runWorker.pl -no_write, so we can write some STDOUT in run
                   #is this already available in the job, or is it just passed ot the worker?
 #      'feature_file_format',
   
        #BaseDB.pm batch_params     
        #Generic optional params (used in Helper and elsewhere)
        'rollback',
        'slices',
        'skip_slices',
        ### More optional Helper params (currently used in DefineOutputSet)
        'result_set_only',
        'result_set_mode', #Can only be set to recover at present
        'recover',         #is this an Importer or a Helper param?
        
 #       ### Optional IdentifySetInputs parameters
 #       #comma separated or defined as list ref  
 #       #are these batch wide or just used for the seed job?
 #       #qw( cell_types feature_types input_sets input_set_ids
 #       #    experimental_groups states ),             
 #       #'experiment'   => undef, 
#        #'input_analyses'     => undef,
 #       #This is actually more like study, and omit for now as 
 #       #input_set will handle this         
     ],
             
    };
}


#pipeline_create_commands in Base and/or subclass config 
#pipeline_analyses in subclass config


1;

