
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::BaseImporter;

=head1 SYNOPSIS

=head1 DESCRIPTION

    This is the base config file for all pipelines which have runnables which use the
    Bio::EnsEMBL::Funcgen::Importer

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut






### TO BE REMOVED!!!!!!!!







package Bio::EnsEMBL::Funcgen::Hive::Config::BaseImporter;

use strict;
use warnings;
use Data::Dumper;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB');

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
  
  #A lot of these which are only set to undef do no need structly to be here
  #Although it is good to maintain them as a reference for what can be set
  #in an input_id  
  
  return {
    %{$self->SUPER::default_options},
          
          
          
     #Todo Figure out which of these can be handled by batch_params
     #may need to change how we deal with default_importer_param
     #change these to param names, then interpolate at run time
     #instead of here, to allow over-ride of param when seeding
     #    Does this even have to be in the config?
     #can't we move it out of the config and into the code? 
     
     #We might even be able to get rid of this config then?
              
   

         
  }
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
    return {
      %{$self->SUPER::pipeline_wide_parameters},
      
      batch_params => [
      
                      ],
          
    };
}



1;

