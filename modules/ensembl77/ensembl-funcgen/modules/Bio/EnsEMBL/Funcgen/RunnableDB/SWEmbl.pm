=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::SWEmbl;

=head1 DESCRIPTION

'SWEmbl' Is a base class for other classes dealing with SWEmbl
It contains virtually nothing so it may disappear and just pass to Funcgen

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::SWEmbl;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Funcgen');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;

sub fetch_input {
  my $self = shift @_;

  $self->SUPER::fetch_input();

  my $efgdba = $self->_efgdba;

  my $aa = $efgdba->get_AnalysisAdaptor;
  my $analysis = $self->param('analysis') || throw "Need to specify the analysis";
  my $analysis_obj = $aa->fetch_by_logic_name($analysis);
  #Here we do not create or modify the analysis on the fly!! - similarly as cell or feature type...
  if(!$analysis_obj){ throw "Analysis $analysis is not in the database"; }
  $self->_analysis($analysis_obj);

  $self->_feature_set_name($self->_set_name()."_".$analysis);

  my $cell_type = $self->_cell_type()->name;
  my $feature_type = $self->_feature_type()->name;
  my $file_type = $self->_file_type();
  my $experiment_name = $self->_experiment_name();

  my $work_dir = $self->_work_dir."/alignments/".$self->_species()."/".$self->_assembly()."/".$experiment_name;
  my $input_dir = $self->param('input_dir') || $work_dir;
  $self->_input_dir($input_dir);

  #TODO Change this to accept samse and sampe and others?? - add an extra parameter when needed?
  #Also not necessarily .gz ... maybe force the use of the parameter 'data_file' in the input_id
  #my $input_file =  $self->param('data_file') || $cell_type."_".$feature_type.".samse.".$file_type.".gz";
  my $input_file =  $self->param('data_file') || $self->_set_name().".samse.".$file_type;
  $input_file .=  ".gz" unless ( -e $input_dir.'/'.$input_file);
  $self->_input_file($input_file);

  my $skip_control = $self->_skip_control($self->param('skip_control'));
  if(!$skip_control){
    #May also be passed as input_id, but then input_id may need to be TEXT
    #Control must be in same dir as input file... maybe change this...
    if(!$self->param('control_file')){
      #For the moment the control file has to be the same file type as the input file, but does not need to be the case...
      #TODO Need to validate here if that's the case...
      my $control_feature = $self->param('control_feature') || throw "Need to define 'control_feature'";
      my $control_file = $cell_type."_".$control_feature."_".$experiment_name.".samse.".$file_type;
      $control_file .=  ".gz" unless ( -e $input_dir.'/'.$control_file);
      $self->_control_file($control_file);
    }
  }

  return 1;
}


# Private function only to be called by subclasses of this class
# prepares data to be used with SWEMBL...
# sorts the input, removes mythochondria and unaligned reads...
sub _preprocess_file{

  #Consider using hash to process input
  my ($self, $input, $output, $file_type) = (shift, shift, shift, shift);
  return 1  if (-e $output);
  #For the moment we always overwrite any existent file...
  #Maybe reuse previously cached files? How to check if they are corrupted?
  #Maybe create a -reuse flag?

  #running piped system commands is a potential source of untraceable errors!!
  #TODO try changing this... e.g. with Bio::DB::Sam ... ?
  my $command;
  if($file_type eq 'bed'){
    $command = "gzip -dc ${input}";

    #Remove mitochondria before SWEMBL
    warn "Excluding mytochondria reads before passing to SWEMBL. Make sure Bed file has MT for mytochondria";
    warn "Duplicates are not removed in Bed files while they are in sam files";
    #This probably won't work for many BED files... e.g chrM or sam-like "beds"...
    $command = " | grep -v '^MT' ";
    $command .= " | sort -k1,1 -k2,2n -k3,3n | gzip -c > $output";

  } elsif($file_type eq 'sam' || $file_type eq 'bam'){

    #Manual Alternative to samtools
    #$command = "gzip -dc ".$self->param($parameter);
    #$command .= ' | grep -vE \'^@\' | grep -vE "^[^[:space:]]+[[:blank:]]4[[:blank:]].*$" | sort -k3,3 -k4,4n | gzip -c > '.$work_file;

    #Sometimes the input sam file may have incorrect headings which will mess up subsequent steps...
    #$command = "gzip -dc $input | grep -v '^\@' | "; # This is not likely to happen now
    if ($file_type eq 'sam'){
      $command = "gzip -dc $input | ";
    }
    else {
      $command = "samtools view $input |";
    }

    #PREPROCESSING... remove mitochondria before SWEMBL : we need to cater for different approaches
    #TODO do this in a better, more generic way
    $command .= "grep -vE '^[^[:space:]]+[[:blank:]][^[:space:]]+[[:blank:]][^[:space:]]+\:[^[:space:]]+\:MT\:' | ";
    $command .= "grep -v '^MT' | grep -v '^chrM' | ";

    #Remove unmapped reads...
    $command .= $self->_bin_dir()."/samtools view -uSh";
    $command .= " -t ".$self->_sam_header() if $self->_sam_header();
    $command .= " -F 4 - | ";

    # This piped sort is not working!! (this problem has been reported in the samtools mailing list)
    #TODO Check why this pipe in the sort is not working...
    #$command .= "samtools sort -o - ".$self->param($parameter)."_tmp | " ;
    #$command .= "samtools view -h - | gzip -c > $work_file";

    # Create temp file for now and remove it when we figure out what's wrong with the piped sort!
    $command .= $self->_bin_dir()."/samtools sort - ${input}_tmp ; " ;

    #Add a remove duplicates step (this is not supported with BED files for the moment)
    $command .= $self->_bin_dir()."/samtools rmdup -s ${input}_tmp.bam - ";

    if ($file_type eq 'sam'){
      $command .= "| ".$self->_bin_dir()."/samtools view -h - | gzip -c > $output";
    }
    else {
      $command .= "> $output";
    }
    #Alternative with no rmdup...
    #$command .= $self->_bin_dir()."/samtools view -h ${input}_tmp.bam | gzip -c > $output";

    $command .= " ; rm -f ${input}_tmp.bam";

  }
  else{
    throw("$file_type file format not supported");
  }
  print STDERR "Preprocess cmd: $command$/";
  system($command) && throw("Failed processing $input with command $command");

  return 1;
}

# Private function only to be called by subclasses of this class
# gets the number of reads in a sam or bed file
#sub _get_number_of_reads {
#   my ($self, $file, $file_type) = (shift, shift, shift);
#   if(($file_type ne "bed") && ($file_type ne "sam")){ throw "Only bed and sam file types supported"; }
#   my $nbr_reads = 0;
#   #If needed, add an option to check if is zipped or not...
#   my $open_cmd = "gzip -dc $file |";
#   open(FILE,$open_cmd);
#   while(<FILE>){
#     if($file_type eq "sam"){
#       next if /^\@SQ/;
#     }else {
#       next if /track name=/o;
#     }
#     $nbr_reads++;
#   }
#   close FILE;
#   return $nbr_reads;
#}

# Private function only to be called by subclasses of this class
# gets the number of reads in a sam or bed file
#sub _get_slices {
#  #NOT DONE!!
#   my ($self, $file, $file_type) = (shift, shift, shift);
#   if(($file_type ne "bed") && ($file_type ne "sam")){ throw "Only bed and sam file types supported"; }
#   my $nbr_reads = 0;
#   #If needed, add an option to check if is zipped or not...
#   my $open_cmd = "gzip -dc $file |";
#   open(FILE,$open_cmd);
#   while(<FILE>){
#     if($file_type eq "sam"){
#       next if /^@SQ/;
#     }else {
#       next if /track name=/o;
#     }
#     $nbr_reads++;
#   }
#   close FILE;
#   return $nbr_reads;
#}


#Private getter / setter to the Feature Set Name
sub _feature_set_name {
  return $_[0]->_getter_setter('feature_set_name',$_[1]);
}

#Private getter / setter to the input folder
sub _input_dir {
  return $_[0]->_getter_setter('input_dir',$_[1]);
}

#Private getter / setter to the Input subset (usually a file name)
sub _input_file {
  return $_[0]->_getter_setter('input_file',$_[1]);
}

#Private getter / setter to the skip control option
sub _skip_control {
  return $_[0]->_getter_setter('skip_control',$_[1]);
}

#Private getter / setter to the control file
sub _control_file {
  return $_[0]->_getter_setter('control_file',$_[1]);
}

1;

