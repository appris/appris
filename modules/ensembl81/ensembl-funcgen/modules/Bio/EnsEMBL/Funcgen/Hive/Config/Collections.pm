
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::Collections;

=head1 SYNOPSIS

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Import_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the Config file for the Import Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Import pipeline consists of several "analysis":
        * SetupImportPipeline is equivalent to the "prepare" in parse_and_import.pl
        * RunImport loads the reads per each slice...
        * WrapUpImport finalizes when all partial imports are done...

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::Hive::Config::Collections;

#Should we extend this to integrate all types of signal ResultSet generation
#or create separate configs for say methylation sets?

#This would be fine, so long as with run init_pipeline with -job_top_up, 
#using the different params required

#Or do we want to split thi sout into yet another conf?
#Currently the actual business part of the pipeline is different for meth-seq vs chip-seq
#so would need different confs, unless we have a wrapper module? 
#Too messy and unecessary coding vs sligntly redundant conf?
#If we merged confs and then specified -feature_types 5mC,DNase1 then we would
#have differing requirements for these input
#How will the pipeline know which way to flow the inputs? i.e. either a collection of a bigwig output?
#Would definately require an intelligent wrapper, 'PrepareInputAlignments' would need need to 
#flow based on ftypes or some part of the input_id specified when using job_top_up
#this could be made mandatory if the configs were merged. How can the pipeline itself identify this?
#Is this conditional flow even valid?

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly

#todo this change to Signal/Wiggle.pm config 
#with DefineSetInputs and DefineSets as shared analyses

#We really want to be able to genericise this such that we can swap out
#any given Runnable/DB. Move all specific options to parameters sections
#specific for that analysis



# WARNING
# Never set the following values in here as they should be set in a previous conf:
#   bam_filtered


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut

#Comment/Omit this method if we have none
#All of these have moved to BaseImporter as they are shared across the Peaks & Collections conf


#sub default_options {
#  my ($self) = @_;
#  return 
#    {
#     %{$self->SUPER::default_options},        
        #Define any optional params here as undef
  #If they are mandatory for one analysis, but optional for another
  #will have to catch that in the runnable
   
#  };
#}

=head2 resource_classes

    All generic resources are defined in Base.pm
  
=cut


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
  my ($self) = @_;
  return 
   {
    %{$self->SUPER::pipeline_wide_parameters}, 
    #can_Link_to_Collections       => 1,     
	 };
}


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
	   -logic_name => 'PreprocessAlignments',
	   # This is basically making sure the input file is sorted wrt genomic locations

	   -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	   #-parameters => {},
				 	 
	   -flow_into => { 2 => ['WriteBigWig'] }, 

	   -analysis_capacity => 50,
	   #Change this to hive_capacity as it may be competing with parallel peak jobs
     -rc_name => 'normal_2GB',
     #can we key this into another param for the memusage which will also be pickup by get_feature
     #file and set the same mem spec in the samtools cmd?
     #todo revise this to reserve tmp space as this is sorting the bed files
    },
    
    {
     -logic_name    => 'WriteBigWig',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools',
     -parameters    => {mode => 'RPKM'},
     #-input_ids     => [ dataflowed from PreprocessAlignments via branch 2 ]
     -analysis_capacity => 1000,
     #Change this to hive_capacity as it may be competing with parallel peak jobs
     -rc_name => 'normal_16GB_2cpu',
     # This resource usage is a product of the read depth, so it could be detected in
     # PreprocessAlignments and flowed selectively to WriteBigWig_10gb WriteBigWig_16gb
     # This would be beneficial if throughput is hampered by spec required for biggest jobs
     # ~ 2/3rds run with 10GB or less
     # ~ 1/4 run with 16GB
     # The small remainder require ~25GB
     # One clocked in over 31.1GB!
     # Actually max already reported as ~24GB (via lsf_report.pl)
     # Best current approach is to run with 16GB by default
     # Then switch the resource_class_id to the 30GB resource for the failed jobs.
     # TODO
     # Analyse the usage here wrt file size and implement the dynamic data flow as per above.
     # Likely need 2/3 cpu also due to pipe
     # Doesn't need to be monitored, as there is little to now DB work here

    },  
  ];
}

1;


__END__


=pod    
    {
     -logic_name => 'Link_to_Collectionss',
     -meadow     => 'LOCAL',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MultiConfigLinker',
     
     #Need to added in previous config as will not be updated by -analysis_topup
     #-parameters => 
     # {
     # },
     
     -flow_into =>
      { #Maintained original branches here, but could revise these
       '2->A' => [ 'WriteSignalCollection' ],
       'A->1' => [ 'MergeCollections' ],
      }, 
       
     -analysis_capacity => 100,
     -rc_name => 'default',
    },     

=cut

