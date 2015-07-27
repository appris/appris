#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::InputSet
#

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

Bio::EnsEMBL::Funcgen::Parsers::InputSet

=head1 SYNOPSIS

  use vars qw(@ISA);
  @ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::InputSet);

  
=head1 DESCRIPTION

This is a base class to support simple file format parsers. For simple imports the vendor is
set to the parser type i.e. the file format.  The generic read_and_import_simple_data assumes
a one line per feature format, other format need there own read_and_import_format_data method, 
which will need defining in the result_data config element. Features are stored either as 
ResultFeature collections or AnnotatedFeatures dependan ton the input feature class.

=cut

# To do
# Add Parsers for BAM/SAM
# Rename to InputSet
# Handle mysqlimport for large data sets e.g. reads
# Incorporate collection code
# Implement matrix storage

package Bio::EnsEMBL::Funcgen::Parsers::InputSet;

use strict;
use warnings;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar            qw( check_ref assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( species_chr_num open_file is_gzipped );
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::SegmentationFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;

use base qw( Bio::EnsEMBL::Funcgen::Parsers::BaseImporter );


#todo change input_set_name to name

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Bed class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Simple
  Exceptions : throws if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Parsers:Simple
  Status     : at risk

=cut


sub new{
  my $caller = shift;  
  my $class  = ref($caller) || $caller;
  my $self   = $class->SUPER::new(@_);

  ($self->{'input_feature_class'}, 
   $self->{total_features}, 
   $self->{force},                  #Is this generic enough to go in Importer? used by store_window_bins_by_Slice_Parser
   $self->{dbfile_data_root},       #only appropriate for result input_feature_class
   $self->{merged_replicates},
   $self->{output_set}              #Feature/ResultSet
  ) = rearrange
    (
     ['input_feature_class', 'total_features', 'force',
      'dbfile_data_root', 'merged_replicates', 'output_set'],
     @_);


  #Could potentially take fields params directly to define a custom format
  #Take direct field mappings, plus special fields which needs parsing differently
  #i.e. default is tab delimited, and GFF would define Attrs field as compound field and provide special parsing and field mapping
  

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") 
	if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
  
  $self->{'config'} =  
	{(
	  #can we omit these?
      array_data   => [],#['experiment'],
      probe_data   => [],#["probe"],
      norm_method => undef,
	  #protocols => {()},
	  'results_data' => ["and_import"], 
     )};

  #set up feature params
  $self->{'_feature_params'} = {};
  $self->{'_dbentry_params'} = [];
   
  return $self;
}


sub output_file{
  my ($self, $output_file) = @_;
  $self->{output_file} = $output_file if $output_file;  
  return $self->{output_file};
}

sub output_set { return $_[0]->{output_set}; }


sub input_file{
  my ($self, $input_file) = @_;

  $self->{'input_file'} = $input_file if $input_file;
  return $self->{'input_file'};
}



=head2 set_config

  Example    : my $self->set_config;
  Description: Sets attribute dependent config
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


#Set config assumes we are importing an InputSet, rather than providing
#an existing ResultSet!
#We shouldn't be depedant on params for config until after we have called new?

#This is not even an InputSet import really. It's a Collection/ResultSet import
#Depending on how you look at it.

#todo refine/simplify this once we separate the Imports/parsers
#This needs to move to new anyway
#input_feature_class refers to output not input!

sub set_config{
  my $self = shift;

  #Move all this to new when we fix the inheritance in Importer

  #We could set input_set_name to experiment name
  #But we would have to make warning in define_and_validate_sets mention -input_set_name
  #


  if(! defined $self->{output_set}){
    throw('Must provide an -name for a '.uc($self->vendor).' import') if ! defined $self->name();

    #Mandatory checks
    if(! defined $self->feature_analysis){
	 throw('Must define a -feature_analysis parameter for '.uc($self->vendor).' imports');
    }
  }



  #some convenience methods
  $self->{dbentry_adaptor}   = $self->db->get_DBEntryAdaptor;
  $self->{input_set_adaptor} = $self->db->get_InputSetAdaptor;
  
  my $fclass_name = $self->{input_set_adaptor}->build_feature_class_name($self->input_feature_class);
  $self->{'feature_class'} = 'Bio::EnsEMBL::Funcgen::'.$fclass_name;

  #We need to undef norm method as it has been set to the env var
  $self->{config}{norm_method} = undef;

  #dirs are not set in config to enable generic get_dir method access
  $self->{dbfile_data_root} ||= $self->get_dir('output');#Required for Collector

  my $adaptor_method         = 'get_'.$fclass_name.'Adaptor';

  #cannot validate here using 'can' as it doesn't appear to work for autoloaded methods :(
  eval { $self->{feature_adaptor} =  $self->db->$adaptor_method; };

  if($@){
    throw("$@\nFailed to $adaptor_method, ".$self->input_feature_class.
      ' is not a valid feature class');
  }


  #Validate slices
  $self->slices($self->{'slices'}) if defined $self->{'slices'};

  #Move to new when we sort out inheritance
  $self->validate_and_store_config([$self->name]);
  #Could use input_set_name here?
  #This was to support >1 input set per experiment (name)

  #This current breaks for no config imports
  #i.e. standard Bed import e.g. result_feature collections
  #segmentation imports use Bed and config
  #allow no config imports in BaseImporter?
  #or ultimately set the params as part of the user_config?

  return;
}


#Some of this can probably move to define_and_validate_sets?


#We have issues here as we need to know the rollback status of a set here
#otherwise it is unsafe, which is why we did it together in the first place
#i.e. if we rollback prior to this method, then there is no way to guarantee
#that a passed set will be rollback properly

#If we define the sets generically beforehand, we also don't have the context provided here which allows us
#to check whether the given set has the appropriate states for import
#given the rollback modes/level



sub define_sets{
  my ($self, $dset) = @_;

  my ($output_set, $eset);

  
  warn("define_sets need updating to use new set definition methods in Helper.\n".
  "i.e. needs mirror DefineOutputSet::run, then point that method here\n");

  $eset = $self->db->get_InputSetAdaptor->fetch_by_name($self->name);
    
  if(! defined $eset){

    #temporary fix until we implement input_subset_file replicate/control config
    #and merge validite_files in here
    #do not currently have non merged replicate import
    #need to define this param in all the run scripts
    my $replicate = 0 if $self->merged_replicates;

    $eset = Bio::EnsEMBL::Funcgen::InputSet->new
      (
       -name          => $self->name(),
       -experiment    => $self->experiment(),
       -feature_type  => $self->feature_type(),
       -cell_type     => $self->cell_type(),
       -analysis      => $self->feature_analysis,
       -feature_class => $self->input_feature_class,
       -replicate     => $replicate,
      );


    
    ($eset)  = @{$self->db->get_InputSetAdaptor->store($eset)};
  }

  $dset = $self->define_and_validate_sets
	 (
    -dbadaptor    => $self->db,
    -name         => $self->name,#or name?
    -feature_type => $self->feature_type,
    -cell_type    => $self->cell_type,
    -analysis     => $self->feature_analysis,
    -feature_class=> $self->input_feature_class, 
    -description  => $self->feature_set_description,
    #-append          => 1,#Omit append to ensure we only have this eset
    -recovery     => $self->recovery,
    -supporting_sets => [$eset],
    -slices        => $self->slices,
    #Can't set rollback here, as we don't know until we've validated the files
    #Can't validate the files until we have the sets.
    #So we're doing this manually in validate_files
   );

  #We are now using IMPORTED to define wheather a FeatureSet has been imported succesfully
  #However we already have IMPORTED on the InputSubSet
  #We should add it to FeatureSet to remain consistent.
  #See Helper::define_and_validate_sets for more notes on
  #potential problems with FeatureSet IMPORTED status
 
  
  $self->{'_data_set'} = $dset;



  if ($self->input_feature_class eq 'result' ||
      $self->input_feature_class eq 'dna_methylation') {
    
    $output_set = $dset->get_supporting_sets->[0];
    $eset       = $output_set->get_InputSets->[0];

    $self->result_set($output_set); #required for ResultFeature Collector and Bed Parser
  } 
  else {                      #annotated/segmentation
    $output_set = $dset->product_FeatureSet;
    $eset       = $dset->get_supporting_sets->[0];
  }

  return ($dset, $output_set, $eset);
}




#we have rollback functionality incorporated here
#why is this not integrated into define_sets?
#they are inherantly linked and for an InputSet import
#we always want to validate files.
#No way of defining replicate/control config just yet
#This should be changed to support both tracking DB and direct import.
#to correctly define the sets, we need access to the files and in the
#rep/control config
#We can do this implicitly by simply letting the order define the rep number
#-result_files rep1 rep2 rep2
#-control_files rep1 rep2
#-merged_replicates would over-ride the above and trigger a check for >1 file
#let's implement -merged_replicates for now and come back to this order based config
#rename -result_files to -input_subset_files?
#Will this support the collection import?
#Need to be mindful about the true input files and what pipeline will actually use
#as it may be some merged intermediate i.e. a sam/bed file rather than a fastq


#should we also store the dbfile_registry data here?
#this maybe done before generating the actual files

#This overlaps with new Set definition methods
#Need to use new methods which handle rollback

#todo remove InputSubset generation from here
#as this should be done by register_experiment?
#What about data types which do not have support, such as external feature sets?
#Or segmantations?


#

sub validate_files{
  my ($self, $prepare, $eset, $output_set) = @_;

  if(! (check_ref($eset, 'Bio::EnsEMBL::Funcgen::ResultSet') &&
        ($eset->table_name eq 'input_subset'))){
    throw('validate_files currently only supports processing ResultSets with InputSubset support');        
  }

  #This has now been changed to remove support for InputSet
  #and currently expects all InputSubsets to be registered before hand
  #If we are to genericse this into a generic Set importer
  #then we will need to maintain registration functionality in here


  #This no longer works, as the old code had merged InputSubset and now input_subsets 
  #represent disitnct replcicates
  
  #We simply want to change how this works to use the single input_file
  #provided.
  #We were previously setting IMPORTED on subsets, 
  #but now this must be done on ResultSets instead?
  #The match between the ResultSet name and the file name is slightly different
  #This check was to ensure cross valdiation between the file and the input_subset
  #filepath generation code is only available in the Hive at present?
    
  
  #And we will need to remove IMPORTED status setting for rsets ealier in the pipeline
  #What implications will this have?
  #The states are getting a bit messy here
  #as we can have a ResultSet which is LOADED, associated with a DataSet,
  #but itself not IMPORTED. This is what recovery mode is for.
  
  #Is there a danger of deleting a ResultSet which does not have the IMPORTED state
  #therefore invalidating the DataSet?
  #There is logic in place to prevent this
  
  #Let just simplify this massively for now to get things running
  #We just need to identify new data and rollback if we have the relevant rollback flags set

  
  #Assume we have got the right file wrt to the ResultSet
  #and we don't need to use the modified filename as a key any more, just use the full path?


  
  #Here were are tracking the import of individual files by adding them as InputSubSets
  #Recovery would never know what to delete
  #So would need to delete all, Hence no point in setting status?
  #We do not rollback IMPORTED data here.  This is done via separate scripts
  #To reduce the rick of accidentally deleting/overwriting data by leaving a stray -rollback
  #flag in the run script

  ### VALIDATE FILES ###
  #We need validate all the files first, so the import doesn't fall over half way through
  #Or if we come across a rollback halfway through
  my (%new_data);

  if ((scalar(@{$self->slices}) > 1) &&
      ! $prepare) {
    throw('Validate files does not yet support multiple Slice rollback');
  }
  
 
  #NEW CODE STARTS
  assert_ref($self->input_files, 'ARRAY', 'input_files');
  
  if(scalar(@{$self->input_files}) != 1) {
    throw("Did not find 1 input_file:\n\t".join("\n\t", @{$self->input_files}));    
  }
  
  my $filepath = $self->input_files->[0];
  (my $filename = $filepath) =~ s/.*\///;#Just used for logging now
  #$filename =~ s/(prepared\.){0,1}(.*)(?:\..*$)/$1$2/;
  
   
  $self->log('Validating '.$self->vendor." file:\t$filename");
  throw("Cannot find ".$self->vendor." file:\t$filepath") if(! -e $filepath); 
 
 
  my $rollback;
 
  #This is decoupling the STATUS from the filename(was input_subset)
  #So this is a bit risky, if we ever want to use a different file, as it may be recognised
  #as already IMPORTED, or rollback other data
  #We need to update the status here to be more specific
  #COLLECTION_IMPORTED?
  #When do we set IMPORTED for ResultSet
  #Is this even useful now
  #It was a simple status to signify whether an import of a ResultSet was complete
  #and in case of reload, either rollback incomplete imports or not
 
 
  if(! $eset->has_status('IMPORTED')){
    $self->log("Found existing ResultSet without IMPORTED status:\t$filename");
    $rollback = 1;
  }
  elsif($self->rollback){
    $self->log("Rolling back IMPORTED ResultSet:\t$filename");
    $rollback  = 1;
  }
  else{ #IMPORTED
     $new_data{$filepath} = 0;
     $self->log("Found previously IMPORTED ResultSet:\t$filename");
  }
 
  #Is recovery even valid here any more?
 
 
  if ($rollback &&         #recoverable sets i.e. exists but not IMPORTED
      ( (! $self->recovery) && (! $self->rollback) ) ) {
    throw("Found ResultSet without IMPORTED status:\n\t".$eset->name."\n".
          "You must specify -recover or -rollback to perform a full rollback");
  }
  
  
  
  if($rollback){
    $new_data{$filepath} = 1;
   
    if (! $prepare) {
      #Don't rollback when in prepare mode as it will also be done
      #redundantly in parallel slice jobs
      #and to ensure that data is definitely rolled back before import
      $self->log("Rolling back ResultSet:\t\t\t\t".$eset->name);
    
      if ($self->input_feature_class eq 'result' ||
          $self->input_feature_class eq 'dna_methylation') {
        $self->rollback_ResultSet($eset, $self->slices->[0], $self->recovery);
        #Do not rollback_InputSet here as we may have parallel Slice based imports running
      }
    }
  }

  return \%new_data;
 
 
 
  #OLD CODE STARTS
 
 
 
 
 
  #Get files from input directory
  if (! $self->input_files) {
    #todo genericise this input directory? i.e. set a defined dir, and do not concat name here
    
    #todo should change this to use get_feature_file
    #i.e. the same code which is being used in the pipeline
    
    my $list = "ls ".$self->get_dir('input').'/'.$self->name().'*.';
    my @rfiles = run_backtick_cmd($list);
    $self->input_files(\@rfiles);
  }
 
 
 
 
  
  #IMPORTED status here may prevent
  #futher slice based imports
  #so we have wait to set this until we know all the slices 
  #are loaded, unless we store slice based IMPORTED states
  #We currently get around this be never settign IMPORTED for slice based jobs
  #and always rolling back by slice before import

  #This loop supports multiple files
  my (@rollback_sets, %file_paths);
  my $auto_rollback = ($self->rollback) ? 0 : 1;

  
  #set the max seen replicate value
  my $replicate = 0;

  my (%support);  

  #foreach my $iss( @{$eset->get_InputSubsets} ){
  foreach my $iss( @{$eset->get_support} ){ 
    
    
    if(! $iss->is_control){
      $support{$iss->name} = $iss;  
    
      if($iss->replicate >  $replicate){
        $replicate = $iss->replicate;
      }
    }
  }

 
  
  


  foreach my $filepath ( @{$self->input_files} ) {
    my ($filename, $sub_set);
    chomp $filepath;
   	($filename = $filepath) =~ s/.*\///;
   	$filename =~ s/(prepared\.){0,1}(.*)(?:\..*$)/$1$2/; 
   	#Strip of any suffixes! Does this match the exact input_subset name processing?
   	
    $file_paths{$filename} = $filepath;			 
    $filename =~ s/^prepared\.// if $self->prepared; #reset filename to that originally used to create the Inputsubsets

    $self->log('Validating '.$self->vendor." file:\t$filename");
    throw("Cannot find ".$self->vendor." file:\t$filepath") if(! -e $filepath); 
    #Can deal with links
		
	
		
		
    #if ( $sub_set = $eset->get_subset_by_name($filename) ) {
    if(exists $support{$filename}){
      $sub_set = $support{$filename};
      
      
      #IMPORTED status here is just for the file
      #Any changes to analysis or coord_system should result in different InputSubset(file)
      #Will only ever be imported into one Feature|ResultSet

      #Currently conflating recover_unimported and rollback
      #as they serve the same purpose until we implement InputSubset level recovery

	
      if ( $sub_set->has_status('IMPORTED') ) {
        $new_data{$filepath} = 0;
        $self->log("Found previously IMPORTED InputSubset:\t$filename");
      } else {
        $self->log("Found existing InputSubset without IMPORTED status:\t$filename");
        push @rollback_sets, $sub_set;
      }
    } else {
      
      throw('This needs updated to support the new InputSubset data model, or removing in favour of register_experiment.pl');
      
      $replicate ++;

      $self->log("Found new InputSubset:\t${filename}");
      #throw("Should not have found new 'prepared' file:\t$filename") if $self->prepared;
      #This is no longer quite true
      #we may have preprepared the file, but simply not registered the InputSubset yet
      
      	  
      $new_data{$filepath} = 1;

      $replicate = 0 if $self->merged_replicates;
      #need to add is_control support here too.
      #tidy this up, force full creation here and remove convenience code from InputSet?
      #where else is it used?

      $sub_set = Bio::EnsEMBL::Funcgen::InputSubset->new
        (
         -name       => $filename,
         -replicate  => $replicate,
         -input_set  => $eset,
         # should also handle these here
         -is_control  => 0,#for result files, need to handle control files separately
         #to be defined
         #-display_url => $display_url,
         #-archive_id  => $archive_id,
        );

      $self->input_set_adaptor->store_InputSubsets([$sub_set]);
    }
  }


  #Does -recover allow a single extra new file to be added to an existing InputSet?
 
  if (@rollback_sets &&         #recoverable sets i.e. exists but not IMPORTED
      ( (! $self->recovery) && (! $self->rollback) ) ) {
    throw("Found partially imported InputSubsets:\n\t".join("\n\t", (map $_->name, @rollback_sets))."\n".
          "You must specify -recover or -rollback to perform a full rollback");

    if ($self->recovery) {
      #Change these to logger->warn
      $self->log("WARNING::\tCannot yet rollback for just an InputSubset, rolling back entire set? Unless slices defined");
      $self->log("WARNING::\tThis may be deleting previously imported data which you are not re-importing..list?!!!\n");
    }
  }
  
  
  if ($self->rollback) {
    #Check we have all existing InputSubsets files before we do full rollback
    #Can probably remove this if we support InputSubset(file/slice) level rollback
    $self->log('Rolling back all InputSubsets');
    @rollback_sets = @{$eset->get_InputSubsets};

    foreach my $isset (@rollback_sets) {
	  
      if (! exists $file_paths{$isset->name}) {
        throw("You are attempting a multiple InputSubset rollback without specifying an existing InputSubset:\t".$isset->name.
              "\nAborting rollback as data will be lost. Please specifying all existing InputSubset file names");
      }
    }
  }


  #This is now all handled in the define sets methods
  #todo remove all this in favour of using new define sets methods in helper
  #is there anyway we can bypass this for now if we have passed the outputset?
  #can we sub out part of this to identify those which are new and those which need rolling back?
  



  foreach my $esset (@rollback_sets) {
    #This needs to be mapped to the specified filepaths
    my $fp_key = $esset->name;
    $fp_key = 'prepared.'.$fp_key if $self->prepared;

    $new_data{$file_paths{$fp_key}} = 1;
    $self->log("Revoking states for InputSubset:\t\t\t".$esset->name);
    $eset->adaptor->revoke_states($esset);

    if (! $prepare) {
      #This was to avoid redundant rollback in prepare step
      #and to ensure that data is definitely rolled back before import
      $self->log("Rolling back InputSubset:\t\t\t\t".$esset->name);
	  
      if ($self->input_feature_class eq 'result' ||
          $self->input_feature_class eq 'dna_methylation') {
        $self->rollback_ResultSet($output_set, $self->slices->[0], $self->recovery);
        #Do not rollback_InputSet here as we may have parallel Slice based imports running
      }
      else {                  #annotated/segmentation
        $self->rollback_FeatureSet($output_set, $self->slices->[0], $self->recovery);
        #undef is full delete.
        $self->rollback_InputSet($eset);
        last;
      }			
    }
  }

  return \%new_data;
}




sub set_feature_separator{
  my ($self, $separator) = @_;

  #How do we test if something undefined was passed?
  #Rather than nothing passed at all?
  #Can't do this as this is the accessor
  #Need to split method
 
  throw('Must provide a valid feature separator') if ( (! defined $separator) || ($separator eq '') ); 

  $self->{'_feature_separator'} = $separator;

}

### SIMPLE ACCESSORS ###
# Some of these can be called for each record
# Trim the access time as much as possible

#This is to support input_set/input_subset.replicate
#Not appropriate for array experiment imports

sub merged_replicates{   return $_[0]->{merged_replicates}; }

sub input_feature_class{ return $_[0]->{input_feature_class}; }

#sub input_set_name{      return $_[0]->{input_set_name}; }  #Set in new

sub feature_adaptor{     return $_[0]->{feature_adaptor}; }

sub dbentry_adaptor{     return $_[0]->{dbentry_adaptor}; }

sub input_set_adaptor{   return $_[0]->{input_set_adaptor}; }

sub set{                 return $_[0]->{set}; } #Feature or Result, set in define_sets

sub data_set{            return $_[0]->{_data_set}; }

sub feature_separator{   return $_[0]->{_feature_separator}; }

sub feature_params{      return $_[0]->{_feature_params}; }

sub dbentry_params{      return $_[0]->{_dbentry_params}; }

sub input_gzipped{       return $_[0]->{input_gzipped}; }


sub input_file_operator{
  my ($self, $op) = @_;
  $self->{'input_file_operator'} = $op if defined $op;

  return $self->{'input_file_operator'} || '<';
}

# prepare boolean, simply stores the sets and preprocesses the file


#Still need to implement prepare in other Parsers!!

#We are currently trying to support two methods of import
# 1 Pre-registered sets
# 2 Full set import
# We would also like to support preprepared data, such that we don't have to
# make a copy of the file, and to reduce redundanct wrt sorting filtering 

#prepare (does the sorting and pre_processing) and pre_process just does the slice
#caching?


#This boils down to prepare and prepared
#Sounds counter intuitive, but prepared means it's already sorted, but prepare
#also means to identify which chrs are present



sub read_and_import_data{
  my ($self, $preprocess) = @_;

  my $action = ($preprocess) ? 'preprocessing' : 'importing';
  $action .= ' prepared' if $self->prepared;
  

  $self->log("Reading and $action ".$self->vendor()." data");
  my ($filename, $fh, $f_out, %feature_params, @lines);


  if($preprocess && 
     (! $self->prepared) &&
     (! $self->can('initialise_input_file')) ){
    throw('preprocess mode is not currently available in '.ref($self));
  }
  
 
  #Test for conflicting run modes
  if($preprocess && $self->batch_job){ 
	# ($self->batch_job || $self->prepared)){
	#prepare should be called once by the runner, not in each batch_job
	#don't prepare if already prepared
    throw('You cannot run read_and_import_data in preprocess mode with a -batch_job');
  }
  
  
  

  my ($eset);
  my $output_set = $self->output_set;
  
  #todo integrate all this into define_sets?
  
  if(! defined $output_set){
    throw('Cannot currently handle undef output set, needs updating to remove InputSet');
    (undef, $output_set, $eset) = $self->define_sets();
  }
  elsif($output_set->isa('Bio::EnsEMBL::Funcgen::ResultSet')){
     #$eset = $output_set->get_support->[0];
     $self->result_set($output_set); #required for ResultFeature Collector and Bed Parser
  }
  elsif($output_set->isa('Bio::EnsEMBL::Funcgen::FeatureSet')){
    throw('Need to update FeatureSet mode to use ResultSets instead of InputSets');
    $eset = $output_set->get_DataSet->get_supporting_sets->[0];
  }
  else{
    throw("Unsupported output set class specified:\t".$output_set); 
  }
  
  #We also need to account for bsub'd slice based import

  
  #If we can do these the other way around we can get define_sets to rollback the FeatureSet
  #Cyclical dependency for the sets :|
  
  #validate_files required that the data_set is defined
  
  #my $new_data = $self->validate_files($preprocess, $eset, $output_set);
  my $new_data = $self->validate_files($preprocess, $self->result_set, $output_set);
  my $seen_new_data = 0;
 
  
  ### READ AND IMPORT FILES ###
  foreach my $filepath(@{$self->input_files()}) {
    chomp $filepath;
    #($filename = $filepath) =~ s/.*\///;
    $filename = $filepath;

    if($new_data->{$filepath} ){	#autovivify is okay here
      $seen_new_data = 1;
      $self->input_file($filepath); #Only used by Collector::ResultFeature::reinitialise_input
    
      if( $self->can('initialise_input_file') ){ #Set input_file_operator and output_file
        $filepath = $self->initialise_input_file($filepath, $preprocess);  
      }
	  
	    $self->log_header(ucfirst($action).' '.$self->vendor);
	    $self->log("Input file:\t".$filepath);
      $self->log("Output file:\t".$self->output_file) if $self->output_file;
	        
  
  	  #We need to be able to optional open pipe to gzip | sort here
      #i.e. define open command
      $fh = open_file($filepath, $self->input_file_operator);

      #This my become way too large for some reads files
      #Currently no problems
      #This is not working as we are sorting the file!
      #$self->parse_header($fh) if $self->can('parse_header');
	  
      #For result features some times we want to run 
      #locally and just sort without dumping
      #i.e if we are not a batch job
      #as there is no need to dump if it is a single process
	  
    
      if( (($self->input_feature_class eq 'result') ||
           ($self->input_feature_class eq 'dna_methylation' )) && 
         ! $preprocess){
      
        #Use the ResultFeature Collector here
        #Omiting the 0 wsize
        #How are we going to omit 0 wsize when doing the fetch?
        #simply check table name in ResultSet?
        
        #Should we do this for multiple chrs?
        #or fail here
        # we need to pass self
        #for access to get_Features_by_Slice
        #which should be in the specific parser e.g Bed
    
    		#Will this not clash with standard ResultFeature::get_ResultFeatures_by_Slice?
    		#Could really do with separating the pure file parsers from the importer code, so these can be reused
    		#by other code. Then simply use Bed import parser for specific import functions and as wrapper to 
    		#Bed file parser
    		#So should really have
    		#Parsers::File::Bed
    		#and
    		#Parsers::Import::Bed
    		#This can probably wait until we update BioPerl and just grab the Bed parser from there?
    
        my $slices = $self->slices;
    
    		#Should this be caught in new?
        if(! @$slices){
          throw("You must define a slice to generate ResultFeature Collections from InputSet:\t".$eset->name);
        }
    		
    
        if(scalar(@$slices) > 1){
      	  throw("InputSet parser does not yet support multi-Slice import for ResultFeature collections\n"
      			."Please submit these to the farm as single slice jobs");
      	}

        #restrict to just 1 slice as we don't yet support disk seeking
        #if the slices are not in the same order as they appear in the file
        #also we want to parallelise this
        
        #Set as attr for parse_Features_by_Slice in format sepcific Parsers

	
        $self->file_handle(open_file($filepath, $self->input_file_operator));
        
        
        
        foreach my $slice(@$slices){
          $self->feature_adaptor->store_window_bins_by_Slice_Parser($slice, $self, 
        															(
        															 #Force needs reimplementing here?
        															 -force            => $self->{force},
        															 -dbfile_data_root => $self->{dbfile_data_root},
        															));
        }  

        warn "Need to update InputSubset status to IMPORTED after all slices have been loaded";
        #Do we even need to set RESULT_FEATURE_SET for input_set ResultFeatures?
        warn "Closing $filename\nDisregard the following 'Broken pipe' warning";

        #Closing the read end of a pipe before the process writing to it at the other end 
        #is done writing results in the writer receiving a SIGPIPE. If the other end can't 
        #handle that, be sure to read all the data before closing the pipe.
        #This suggests the gzip pipe has not finished reading, but it *should* be at the end of the file?
        #$SIG{PIPE} = 'IGNORE'; #Catch with subref and warn instead?
        #Or maybe an eval will catch the STDERR better?
        #sub catch_SIGPIPE{
        #  my $sig = shift @_;
        #  print " Caught SIGPIPE: $sig $1 \n";
        #  return;
        #  
        #}
        #$SIG{PIPE} = \&catch_SIGPIPE;
        #my $retval =  eval { no warnings 'all'; $fh->close };
        #if($@){
        #  warn "after eval with error $@\nretval is $retval";
        #}
        #Neither of these catch gzip: stdout: Broken pipe
        
        #IO::UnCompress::Gunzip?

        if(! eval { close($fh); 1;}){
          throw("Found errors when closing file\n$@");
        }
      }
      else{
	  
        #Revoke FeatureSet IMPORTED state here incase we fail halfway through
        #This will have already been rolled back, 
        #todo remove this when we update the define sets method
        #as stats should have been removed already?
        if ($output_set->has_status('IMPORTED') && (! $preprocess)) {
              $output_set->adaptor->revoke_status('IMPORTED', $output_set)
        }
        #What about IMPORTED_"CSVERSION"
        #This may leave us with an incomplete import which still has
        #an IMPORTED_CSVERSION state
        #We need to depend on IMPORTED for completeness of set
        #DAS currently only uses IMPORTED_CSVERSION
        #This is okayish but we also need to write HCs for any sets 
        #which do not have IMPORTED state!
        my ($line, @outlines, $out_fh);

        if($preprocess && ! $self->prepared){  
          #Assume we want gzipped output
          #filename is actually based on input, so may not have gz in file name
          $out_fh = open_file($self->output_file, "| gzip -c > %s");
        }

   	    while(defined ($line=<$fh>)){
          #Generic line processing
          #Move these to parse_line?
          $line =~ s/\r*\n//o;
          next if $line =~ /^\#/;	
          next if $line =~ /^$/;
          
          if($self->parse_line($line, $preprocess)){
            #This counts all lines as oppose to cache_slice count, which only caches
            #data lines where a slice can be cached
            $self->count('total parsed lines');

            #Cache or print to sorted file
            if($preprocess && ! $self->prepared ){

              if(scalar(@outlines) >1000){
                print $out_fh join("\n", @outlines)."\n";
                @outlines = ();
              }
              else{
                push @outlines, $line;
              }
            }
          }
        }
	 
        if(! eval { close($fh); 1;}){
          throw("Found errors when closing file\n$@");
        }
		  

        #Print last of sorted file
        if($preprocess && ! $self->prepared){
          print $out_fh join("\n", @outlines)."\n";
  
          if(! eval { close($out_fh); 1;}){
            throw("Found errors when closing file\n$@");
          }
        
          @outlines = ();
        }

        if(! $preprocess){
          #reset filename to that originally used to create the Inputsubsets
          #some files may not have this prefix
          $filename =~ s/^prepared\.// if $self->prepared;
          my $sub_set = $eset->get_subset_by_name($filename);
          $sub_set->adaptor->store_status('IMPORTED', $sub_set) if ! $self->batch_job;
        }
      }


      #Could really do with warning if seen slices don't match
      #specified slices

      foreach my $slice(@{$self->slices}) {

        if(! exists $self->{'slice_cache'}->{$slice->seq_region_name}){
          #should handle Y for female cell types here?
          warn "Found no data for specified slice:\t".$slice->seq_region_name."\n";  
        }
      }

      if($preprocess){
        $self->log("Finished preparing import:\t".$self->output_file);
      }
      else{
        #Need to tweak this for slice based import
        $self->log('Finished importing '.$self->counts('features').' '.
        $output_set->name." features from:\t$filepath");
      }
	
      foreach my $key (keys %{$self->counts}){
        $self->log('Count '.sprintf("%-25s", "$key:").$self->counts($key));
      }	  
    }
  }

  




  #Here we should set IMPORTED on the FeatureSet
  #We could also log the first dbID of each feature in a subset to facilitate subset rollback
  #in feature table
  #this would be sketchy at best
  #delete from annotated_feature where annotated_feature_id >= $first_subset_feature_id and feature_set_id=$feature_set_id
  #This may already have IMPORTED status as we don't revoke the status whilst
  #updating to protect the feature set due to lack of supportingset tracking
  #see Helper::defined_and_validate_sets for more notes.
  #Is there any point in setting it if we don't revoke it?
  #To allow consistent status handling across sets. Just need to be aware of fset status caveat.
  #Also currently happens with ResultFeatures loaded by slice jobs, as this may already be set by a parallel job

  if(! $preprocess){
    $output_set->adaptor->set_imported_states_by_Set($output_set) if $seen_new_data && ! $self->batch_job;
    #When does this actually set imported states?
    #In a post processing job, but it looks like thi swill stil try and store the collections?
   
    $self->log("No new data, skipping result parse") if ! grep /^1$/o, values %{$new_data};
    $self->log("Finished parsing and importing results"); 
  }
    
  return;
}
  

#Should be called from format parser e.g. BED, GFF, eQTL etc
#Why don't we pass feature_params and dbentry_params directly?

sub load_feature_and_xrefs{
  my $self = shift;
		  
  #warn "Loading ".($self->{_counts}{features}+1).' '.$self->feature_params->{-FEATURE_TYPE}->name."\n";
  #This now only fails once on the first run and then
  #Need to count based on feature_type?

  #new rather than new fast here as we want to validate the import
  my $feature = $self->{feature_class}->new(%{$self->feature_params});
  ($feature) = @{$self->feature_adaptor->store($feature)};
  $self->count('features');

  #Add count based on FeatureType, should be ftype name and analysis to reflect unique ftype key?



  ##This needs to be handled in caller as we are validating loci?
  #if($self->ucsc_coords){
  #	$start += 1;
  #	$end   += 1;
  #  }
  
  #This needs to be put in a separate sub and called by the caller
  #if(!  $self->cache_slice($chr)){
  #  warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
  #}
  #else{
  #grab seq if dump fasta and available		  
  #my $seq;
  #if(exists $self->feature_params->{'sequence'}){
  #	  $seq = $self->feature_params->{'sequence'};
  #	  delete $self->feature_params->{'sequence'};
  #	}
  #	else{
  #	  $self->log('No fasta sequence available for '.$self->feature_params->display_label);
  #	}
  #  }
	  
  #dump fasta here
  #if ($self->dump_fasta){
  #	$self->{'_fasta'} .= $self->generate_fasta_header($feature)."\n$seq\n";
  #  }
  
  #Store the xrefs
 
  foreach my $dbentry_hash(@{$self->{'_dbentry_params'}}){
	my $ftype = $dbentry_hash->{feature_type};
	delete $dbentry_hash->{feature_type};

	my $dbentry = Bio::EnsEMBL::DBEntry->new(%{$dbentry_hash});
	$self->dbentry_adaptor->store($dbentry, $feature->dbID, $ftype, 1);#1 is ignore release flag
	#count here? no count in caller
  }

  
  #Clean data cache
  $self->{'_feature_params'} = {};
  $self->{'_dbentry_params'} = [];

  return $feature;
}


sub total_features{
  my ($self, $total) = @_;

  $self->{'total_features'} = $total if defined $total;
  return $self->{'total_features'};
}

#Currently only required for Bed::parse_Features_by_Slice

#filehandle


sub file_handle{
  my ($self, $fh) = @_;

  $self->{'file_handle'} = $fh if defined $fh;
  return $self->{'file_handle'};
}

sub result_set{
  my ($self, $rset) = @_;

  #already tested/created by self
  
  $self->{'result_set'} = $rset if $rset;
  return $self->{'result_set'};
}


1;
