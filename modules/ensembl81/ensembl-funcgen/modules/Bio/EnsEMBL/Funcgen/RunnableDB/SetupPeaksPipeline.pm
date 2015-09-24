=pod

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::SetupPeaksPipeline

=head1 DESCRIPTION

'SetupPeakPipeline' Does all the setup before the Peaks Pipeline
Checks for existence of Cell and Feature Type
(only flags if they do not exist, does not try to create them!!)
Creates Experiment and Set Entries, Checks Analysis, Group etc...
This Runnable CANNOT be run multiple times in parallell!

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::SetupPeaksPipeline;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::SWEmbl');

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use File::Spec;

sub fetch_input {   # fetch parameters...
  my $self = shift @_;

  $self->SUPER::fetch_input();

  my $efgdba = $self->_efgdba();
  #my $dnadba = $self->param('dnadb');

  #Need to change this...

  #my $group = $self->_group_name();
  my $group = 'efg';
  $efgdba->fetch_group_details($group) || throw "Group $group does not exist in the database. Please create it";

  #Theoretically we should be getting this at the registry...
  #my $csa = $dnadba->get_adaptor('coordsystem');

  my $assembly = $self->_assembly();
  ##Check if the assembly is allowed in the current db we're working on...
  ##In certain DBs this doesn't work: e.g. dev_*_funcgen
  #if(!grep(/$assembly/, map($_->version, @{ $csa->fetch_all()}))){ throw $assembly." is not a valid assembly"; }

  #Checks if the input dir for this experiment_name exists... all input files (including controls must be in this folder)
  my $input_dir = $self->_input_dir();
  if(! -d $input_dir){ throw("Didn't find input directory $input_dir"); }

  #Sets up the output dir for this experiment_name
  my $output_dir = $self->_output_dir();
  if(! -d $output_dir){
    system("mkdir -p $output_dir") && throw("Couldn't create output directory $output_dir");
  }

  ##Better pass this as a parameter... as sometimes the db does not return the species name...
  #my $species = $dnadba->species;
  my $species = $self->_species();

  my $file_type = $self->_file_type();
  if($file_type !~ /^(bed|bam|sam)$/o) {
    throw "File type '$file_type' currently not supported";
  }
#REMOVE if above works
# if(($file_type ne 'sam') && ($file_type ne 'bed') && ($file_type ne 'bam')) { throw "File type $file_type currently not supported"; }

  #infer sam_header from species and assembly...
  if($file_type eq 'sam'){

    #Change the directory structure so it will agree with the rest, without the need to do uc()
    my $sam_header = $self->_sam_header();
    if(! -e $sam_header){
      throw " File type is sam but could not find sam header: '$sam_header'";
    }
  }

  #Check if input file exists. No pre-processing necessary as it can be done in parallel
  my $input_file = File::Spec->catfile($input_dir, $self->_input_file());
  if(! -e $input_file){
    throw " Couldn't find input file '$input_file'";
  }

  #Check if control file exists
  if(!$self->_skip_control()) {
#    my $control_file = $self->_control_file();
    my $control_file = File::Spec->catfile($input_dir, $self->_control_file());
    if(! -e $control_file ){
      #Force throw or just warn?
      #tj1: Changed to warn. Just forgetting the -s in busy times is annoying
      warning("Couldn't find control file '$control_file'. Setting '-skip_control' ");
      $self->_skip_control(1);
    }
  }
  if($self->_skip_control() && ($self->_analysis()->logic_name =~ /ccat/i) ) {
    throw "CCAT analysis requires a control. Cannot skip it";
  }

  return 1;
}

sub run {   # Check parameters and do appropriate database/file operations...
  my $self = shift @_;


  #Preprocess control file if needed. this cannot be done in parallel as several sets may require same control file
  if(!$self->_skip_control()){
    #Do not mix this "input" file with the true data input file...
    my $input_file  = File::Spec->catfile($self->_input_dir(),  $self->_control_file());
    my $output_file = File::Spec->catfile($self->_output_dir(), $self->_control_file());
    #Only do it if it hasn't already been done..
    if(! -e $output_file){
      $self->_preprocess_file(
          $input_file,
          $output_file,
          $self->_file_type(),
          ) || throw "Error processing file $input_file";
    }
  }

  # Check experiment, data set, feature set and create as appropriate...
  $self->_check_Experiment(
      $self->_analysis(),
      $self->_input_file(),
      $self->_feature_set_name());

  if(!$self->_skip_control()){
    #Add control file as a subset, if needed...
    my $input_set = $self->_efgdba()->get_InputSetAdaptor()->fetch_by_name($self->_set_name());
    if(!$input_set) { throw "Input set was not created"; }
    my $isubset_name = $self->_control_file();
    if(!$input_set->get_subset_by_name($isubset_name)){
      #Change InputSetAdaptor so one subset could be stored each time?
#$input_set->add_new_subset($isubset_name);
      #this expects a behavior where subsets already stored will just be ignored and no error is thrown
#      $input_set->adaptor->store_InputSubsets($input_set->get_InputSubsets);
      $isubset_name = Bio::EnsEMBL::Funcgen::InputSubset->new (
         -name         => $isubset_name,
         -input_set     => $input_set,
# Analysis is not being used??
#-analysis     => $self->feature_analysis,
        );
    }
  }

  return 1;
}


sub write_output {  # Create the relevant job
  my $self = shift @_;

  #These numbers need to be parameters...
  ## If feature type is H3K36me3 or H3K27me3 use broad peak caller... Maybe all histone data?
  #if($self->_feature_type()->class eq 'Histone'){
  #if(($self->_feature_type()->name eq 'H3K36me3') || ($self->_feature_type()->name eq 'H3K27me3')){
  #  $self->dataflow_output_id( $self->input_id, $self->_broad_peak_id);

  if ($self->_analysis()->logic_name =~ /ccat/i){
    $self->dataflow_output_id( $self->input_id, 4);
  }
  elsif ($self->_analysis()->logic_name =~ /macs/i) {
    $self->dataflow_output_id( $self->input_id, 5);
  }
  else {
    $self->dataflow_output_id( $self->input_id, 3);
  }

  return 1;
}

1;
