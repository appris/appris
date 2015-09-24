=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::RunCCAT

=head1 DESCRIPTION

'RunCCAT' Runs CCAT "broad peak" caller and stores peaks as an annotated feature set.
Assumes Files are organized and named with a specific convention 
($repository_folder)/experiment_name/cell-type_feature-type/
unless specific files are given

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::RunCCAT;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (run_system_cmd);


use base ('Bio::EnsEMBL::Funcgen::RunnableDB::SWEmbl');

sub fetch_input {   # fetch and preprocess the input file plus do some more checking
  my $self = shift @_;

  $self->SUPER::fetch_input();

  my $efgdba = $self->_efgdba();
  my $fsa = $efgdba->get_FeatureSetAdaptor();
  if(!$self->_feature_set($fsa->fetch_by_name($self->_feature_set_name()))){
    throw "Feature Set was not Created";
  }

  my $bin_dir = $self->_bin_dir();

  my $analysis = $self->_analysis();
  #Use the definitions that are on the database
  
  my $cell_type = $self->_cell_type()->name;
  my $feature_type = $self->_feature_type()->name;
  my $experiment_name = $self->_experiment_name();
  
  my $file_type = $self->_file_type();

  if(($file_type ne 'sam') && ($file_type ne 'bed') && ($file_type ne 'bam')){ 
    throw "Only sam, bam and bed currently supported for CCAT"; 
  }

  my $output_dir = $self->_output_dir();  

  my $size_file = $output_dir."/".$self->_set_name.".sizes";

  #Get the size file... similar to the sam header... 
  open(FILE, $self->_sam_header);
  #Consider having it pregenerated...
  open(SIZES,">".$size_file);
  while(<FILE>){
    chomp;
    /^(\S+)\s+(\d+)\s+/;
    my $slice = $1;
    my $slice_size = $2;
    if(!$slice || !$slice_size){ throw " Could not process sam header line $_ "; }

    #Mouse Hack!!
    next if(($self->_species eq 'mus_musculus') && !($slice =~ /chromosome/));
    print SIZES $slice."\t".$slice_size."\n";
  }
  close SIZES;
  close FILE;

  $self->_size_file($size_file);

  my $input_dir = $self->_input_dir();  
   
  my $file_name = $self->_input_file();
  my $input_file =  $input_dir."/".$file_name;
  if(-e $input_file){ 
    my $output_file = $output_dir."/".$file_name;
    if(! $self->param('reenter')){
      #TODO Validate if existent file is ok. 
      $self->_preprocess_file($input_file, $output_file, $file_type) || 
	throw "Error processing data file $input_file";
    } else {
      if(! -e $output_file){ warn "$output_file does not exist. May need to rerun from scratch."; }
    }
    $self->_input_file($output_file); 
    
    if(!$self->_results_file($self->param('results_file'))){ 
      $self->_results_file($output_file.".".$analysis->logic_name);
    }
  } else { throw "No valid data file was given: ".$input_file; }

  #Always require a control file...
  my $control_file = $output_dir."/".$self->_control_file();
  if(! -e $control_file){ throw "No valid control file was given: ".$control_file; }
  $self->_control_file($control_file); 
  
  #May need to convert it...
  if(! $self->param('reenter') && ($file_type eq 'sam' || $file_type eq 'bam')){

    my $input_file_bed = $self->_convert_file($self->_input_file);
    $self->_input_file($input_file_bed);
    
    my $control_file_bed = $self->_convert_file($self->_control_file);
    $self->_control_file($control_file_bed);
  }

  return 1;
}

sub _convert_file {
  my ($self,$file) = @_;
  my $bin_dir = $self->_bin_dir();
  my $file_type = $self->_file_type();
  
   my $bed_file = $file;
   $bed_file =~ s/\.sam/\.bed/;
   $bed_file =~ s/\.bam/\.bed/;  
    
    my @convert_cmd;
    
    if ($file_type eq 'sam') {
      push @convert_cmd, 'gunzip -c', $file;
      push @convert_cmd, '|', 'samtools view -Su - ';
    }
    else {
      push @convert_cmd, 'cat', $file;
    } 
    push @convert_cmd, '|', "${bin_dir}/bamToBed -i -";
    
    if($self->_species eq 'mus_musculus'){#Mouse hack
      push @convert_cmd, '|', 'grep \'chromosome\'';
    }
    push @convert_cmd, '>', $bed_file;
    
    my $cmd = join ' ', @convert_cmd;
    
    run_system_cmd($cmd);
    
    return $bed_file;
}

sub run {   # call SWEmbl
  my $self = shift @_;

  if($self->param('reenter')){ return 1; }

  my $analysis = $self->_analysis;
  #<CCAT path>/bin/CCAT  <ChIP library read file name>  <control library read file name> 
  # <chromosome length file name>  <config file name>  <project name>
  my $bin_dir = $self->_bin_dir();
  my $command =  $bin_dir."/".$analysis->program_file .
    " ".$self->_input_file() .  
      " ".$self->_control_file() .  
	" ". $self->_size_file() .
	  " ".$bin_dir."/ccat_config/".$analysis->parameters .
	    " " .  $self->_results_file;    
  warn "Running analysis:\t$command";
  run_system_cmd($command);
  
  return 1;
}

sub write_output {  # Store results
  my $self = shift @_;
  
  $self->_parse_result_file();

  return 1;
}

sub _parse_result_file {
  
  my $self = shift @_;
  
  ### annotated features
  my $fset = $self->_feature_set();	
  
  my $efgdba = $self->_efgdba();
  my $sa = $efgdba->get_SliceAdaptor();
  
  #Cache slices and features...
  my %slice;
  my @af;
  
  my %cache_af;

  open(FILE,$self->_results_file().".significant.region");
  while(<FILE>){
    chomp;

    #Content of CCAT output file
    my ($seqid,$summit,$start,$end,$chipreads,$ctrlreads,$fold,$fdr)= split("\t");
    #Hardcode a minimum fdr... may pass as parameter
    next if ($fdr>0.05);
    my $score = $fold;

    #This seqid may vary depending on the input given to SWEmbl... 
    # handle it manually at least for the moment... namely the sam seqid...
    #Make sure to test thoroughly to see if it works...
    #e.g. chromosome:GRCh37:15:1:102531392:1
    if($seqid =~ /^\S*:\S*:(\S+):\S+:\S+:\S/) { $seqid = $1; }
    #In case UCSC input is used... 
    $seqid =~ s/^chr//i;

    if($self->param('slice') && ($seqid ne $self->param('slice'))){
      warn "Feature being ignored as it is not in specified slice ".$self->param('slice')." : Region:".
	$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
      next;
    }

    #May have some sort of repeats(?). Since it is ordered with significance, ignore remaining hits.
    next if(defined($cache_af{$seqid."_".$start}));
    $cache_af{$seqid."_".$start} = 1;

    #next if ($seqid =~ m/^M/);

    # filtering is done as a post-processing e.g. FilterBlacklist.pm
    #$summit = int($summit);#Round up?	  
  
    unless (exists $slice{"$seqid"}) {
      $slice{"$seqid"} = $sa->fetch_by_region(undef, $seqid);
    }

    if( ($start < 1) || ($end > $slice{"$seqid"}->end)){
      warn "Feature being ignored due to coordinates out of slice: Region:".
	$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
    }
  
    #Gracefully handle errors...
    my $af;
    eval{
      $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
    	(
	 -slice         => $slice{"$seqid"},
	 -start         => $start,
	 -end           => $end,
	 -strand        => 0,
	 -score         => $score,
	 -summit        => $summit,
	 -feature_set   => $fset,
	);
    }; 
    if($@) { warn($@); next; }
    if(!$af) { warn("Could not create feature - Region:".
		    $seqid." Start:".$start." End:".$end." Score:".$score.
		    " Summit:".$summit); next; }
    
    push(@af, $af);
  } 
  close FILE;
 
  # Batch store features...
  if(scalar(@af>0)){
    $efgdba->get_AnnotatedFeatureAdaptor->store(@af);
  } else {
    warn "No significant features detected!";
  }
  #Do this on a wrapup runnable...so it will only be visible after reads are loaded...
  $fset->adaptor->set_imported_states_by_Set($fset);

  # Status should not be set at this stage

}


#Private getter / setter to the feature set
sub _feature_set {
  return $_[0]->_getter_setter('feature_set',$_[1]);
}

#Private getter / setter to the results file 
sub _results_file {
  return $_[0]->_getter_setter('results_file',$_[1]);
}

#Private getter / setter to the size file
sub _size_file {
  return $_[0]->_getter_setter('size_file',$_[1]);
}


1;

