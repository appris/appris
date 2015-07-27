=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::RunMACS

=head1 DESCRIPTION

'RunSWEmbl' Runs SWEmbl peak caller and stores peaks as an annotated feature set.
Assumes Files are organized and named with a specific convention
($repository_folder)/experiment_name/cell-type_feature-type/
unless specific files are given

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::RunMACS;

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

#use Data::Dumper;

sub fetch_input {   # fetch and preprocess the input file plus do some more checking
  my $self = shift @_;

  $self->SUPER::fetch_input();

  # Consider using: Bio::EnsEMBL::Helper to debug and/or log to a file (given parameter)
  
  my $efgdba = $self->_efgdba();
  my $fsa = $efgdba->get_FeatureSetAdaptor();
  if(!$self->_feature_set($fsa->fetch_by_name($self->_feature_set_name()))){
    throw "Feature Set was not Created";
  }

  #if(!$self->_feature_set($fsa->fetch_by_name($self->_feature_set_name()))){
  #  throw "Feature set was not created. Maybe you're skipping the Setup Step of the Pipeline";
  #}

  my $analysis = $self->_analysis();
  #Use the definitions that are on the database
  $self->_command($analysis->program_file);
  $self->_command_options($analysis->parameters);
  
  my $cell_type = $self->_cell_type()->name;
  my $feature_type = $self->_feature_type()->name;
  my $experiment_name = $self->_experiment_name();
  
  my $file_type = $self->_file_type();
  my $output_dir = $self->_output_dir();  
  my $input_dir = $self->_input_dir();  
  my $file_name = $self->_input_file();
  my $input_file =  $input_dir."/".$file_name;
  
  if(! -e $input_file){
    throw( "No valid data file was given: ".$input_file) 
  }  
     
  my $output_file = $output_dir."/".$file_name;
  #TODO Validate if existent file is ok. 
  #TODO Add flag to enable override existent file...
  if(!$self->param('reenter')){
    $self->_preprocess_file($input_file, $output_file, $file_type) || throw "Error processing data file $input_file";
  } 
  else {
    if(! -e $output_file){ warn "$output_file not found. May need to rerun from scratch"; }
  }
  $self->_input_file($output_file); 
    
  if(!$self->_results_file($self->param('results_file'))){ 
    $self->_results_file($output_file.".".$analysis->logic_name);
  }
  

  if(!$self->_skip_control()){
    my $control_file = $output_dir."/".$self->_control_file();
    if(! -e $control_file){ throw "No valid control file was given: ".$control_file; }
    $self->_control_file($control_file); 
  }

  return 1;
}




sub run {   
  my $self = shift @_;

  if($self->param('reenter')){ return 1; }

  my %formats = (
    bed   => 'BED',
    sam   => 'SAM',
    bam   => 'BAM',
  );
  
  #Ideally it will unify to sam or bam
  my $fmt = $formats{$self->_file_type()};
  throw("Could not identify valid  input format for ".$self->_file_type()) if(! $fmt);
  
  my %species_genome_size = (
    homo_sapiens => 'hs', #macs accepts initials for common species, size required for others
    mus_musculus => 'mm',
  ); 
  my $genome_size = $species_genome_size{$self->_species()};
  throw("Could not identify genome size for ".$self->_species()) if(! $genome_size);
  
  my @command = ($self->_bin_dir()."/".$self->_command());
  push @command, '-f', $fmt;
  push @command, '-t', $self->_input_file();
  push @command, '-c', $self->_control_file() if($self->_control_file() && !$self->_skip_control());
  push @command, '-g', $genome_size; 
  push @command, $self->_command_options() if ($self->_command_options());
  push @command, '-n', $self->_results_file();
  
  my $command = join(' ', @command);
  
  warn "Running analysis:\t$command";
  run_system_cmd($command);
  
  return 1;
}

sub write_output {  
  my $self = shift @_;
  
  ### annotated features
  my $fset = $self->_feature_set();	
  
  my $efgdba = $self->_efgdba();
  my $sa = $efgdba->get_SliceAdaptor();
  
  #Cache slices and features...
  my %slice;
  my @af;
  
  #Parse the output file
  my $results_file = $self->_results_file().'_peaks.xls';
  if(!open(FILE,$results_file)){ throw "Could not open results file : ".$results_file; }
  while(<FILE>){
    chomp;
    next if (! $_    # blank line
      || /^#/   # comments line
      || /^chr\tstart\tend\tlength/); # header line 

    
    my ($seqid, $start, $end, $length, $summit, $tags, $pvalue, $fold_enrichment, $fdr) = split(/\s+/);
    my $score = $fold_enrichment;   

    
    #This seqid may vary depending on the input given to SWEmbl... 
    # handle it manually at least for the moment... namely the sam seqid...
    #Make sure to test thoroughly to see if it works...
    #e.g. chromosome:GRCh37:15:1:102531392:1
    if($seqid =~ /^\S*:\S*:(\S+):\S+:\S+:\S/) { $seqid = $1; }
    #In case UCSC input is used... carefull names may not match with ensembl db!
    $seqid =~ s/^chr//i;

    if($self->param('slice') && ($seqid ne $self->param('slice'))){
      warn "Feature being ignored as it is not in specified slice ".$self->param('slice')." : Region:".
	$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
      next;
    }

    #Just in case some of the non-aligned were missed in a filtering step... though this seriously affects peak calling
    #if($seqid =~ /\*/){ next; }

    # filtering is done as a post-processing e.g. FilterBlacklist.pm
    $summit = int($summit);#Round up?	  
  
    unless (exists $slice{"$seqid"}) {
      $slice{"$seqid"} = $sa->fetch_by_region(undef, $seqid);
    }

    if( ($start < 1) || ($end > $slice{"$seqid"}->end)){
      warn "Feature being ignored due to coordinates out of slice: Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit."\n";
    }
  
    #Gracefully handle errors...
    my $af;
    eval{
      $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
        -slice       => $slice{"$seqid"},
        -start       => $start,
        -end         => $end,
        -strand      => 0,
        -score       => $score,
        -summit      => $summit,
        -feature_set => $fset,
      );

    }; 
    if($@) { warn($@); next; }
    if(!$af) { warn("Could not create feature - Region:".$seqid." Start:".$start." End:".$end." Score:".$score." Summit:".$summit); next; }

    push(@af, $af);
  } 
  close FILE;
 
  # Batch store features...
  $efgdba->get_AnnotatedFeatureAdaptor->store(@af);
  
  #Do this on a wrapup runnable...so it will only be visible after reads are loaded...
  $fset->adaptor->set_imported_states_by_Set($fset);
  
  return 1;
}

#Private getter / setter to the feature set
sub _feature_set {
  return $_[0]->_getter_setter('feature_set',$_[1]);
}

#Private getter / setter to the command to be executed
sub _command {
  return $_[0]->_getter_setter('command',$_[1]);
}

#Private getter / setter to the command options 
sub _command_options {
  return $_[0]->_getter_setter('command_options',$_[1]);
}

#Private getter / setter to the results file 
sub _results_file {
  return $_[0]->_getter_setter('results_file',$_[1]);
}


1;

