=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::CCAT

=head1 DESCRIPTION

Runs SWEmbl peak caller and optionally parses and processes output.
This is in the Hive namespace, but is not dependnant on any Hive modules and can 
be run as a stand alone job outside of the hive infrastructure.

=cut

package Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::CCAT;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd open_file);
use base qw( Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller );#Does not import


#Chromsome lengths file can be generated using create_CCAT_chr_lengths_from_sam.pl
#Note vanilla CCAT explictly support just 100 entries in this file, although gxx will probably
#allocate enough memory to handle ~130
#This has been fixed to actually use the number of entries in the file to set the 
#array size. The fixed version is available in the funcgen software area. 

sub input_formats    { return ['bed']; }  # readmet.txt suggests bed format only, sigh.
sub requires_control { return 1; }
sub out_file_types   { return [qw( significant_region 
                                   significant_peak 
                                   top100000_peak )];}
 #We can't pass through peak caller specific params, this is done via the 
  #analysis->parameters but we don't want to store internal path in the DB
  #we need to pass this from the config?
  #or at least pass a CCAT config dir, and specify just the file names
  #in the analysis params
  #Having the config file outside of the DB risk losing the config
  #let's put all the config in in the params, and get CCAT to write this file if it doesn't exist already
  #write to /tmp
  #then remove it in DESTROY
  #The chr file shoudl be specified separately
  #This can be some generic CCAT_params config
  #passed from RunPeaks
  
sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = $class->SUPER::new(@_);

  #Passed as $sensitive_caller_params from RunPeaks (CCAT_parameters in analysis config)
  my ($chr_file, $fdr_thresh) = rearrange([qw(CHR_FILE FDR_THRESHOLD)], @_);

  throw("Must provide a valid existing -chr_file parameter:\n\t$chr_file") if ! -e $chr_file;
  $self->{chr_file} = $chr_file;      
  $self->{fdr_threshold} = (defined $fdr_thresh) ? $fdr_thresh : 0.05;
  return $self;
}


sub chr_file     { return shift->{chr_file};      }
sub fdr_threshold{ return shift->{fdr_threshold}; }


sub out_file{ 
  my $self      = shift;
  my $file_type = shift;
  
  if((! defined $file_type) ||
     (! grep {/^${file_type}$/} @{$self->out_file_types})){
    throw("Must provide a valid file type argument($file_type) to build the out_file:\n\t".
      join(' ', @{$self->out_file_types}));
  }
  
  $file_type =~ s/_/./;
  return $self->out_file_prefix.'.'.$file_type;          
}


sub run {   
  my $self = shift;
  my ($signal_bed, $out_prefix, $suffix, $control_bed, 
    $gzipped_align, $gzipped_control) = @{$self->file_info};
        
  my $chr_file = $self->chr_file;
  my $conf_file;
 
  if(! $self->parameters){
    throw('No parameters have been specified, cannot write CCAT config file');  
  }
  else{
    #Dump config to local tmp file
    #Is probably as quick as validating a pre-existing file vs
    #the DB config
    #Why can't this be specified on the cmdline?
    my @config  = split(/\s+/, $self->parameters);
    $conf_file  = "/tmp/CCAT_config.$$.txt";
    my $conf_fh = open_file($conf_file, '>');
    #Slight chance we won't clean up the tmp file if we die here.
    $self->{conf_file} = $conf_file; #For DESTROY
    
    #Could do some more validation here?
    while(my($config_name, $config_val) = splice(@config, 0, 2)){
      print $conf_fh $config_name."\t".$config_val."\n";  
      print $config_name."\t".$config_val."\n" if $self->debug;
      #This is also in log file output     
    }  
    
    close($conf_fh);
  }
 
  # Usage: <library 1 tag file name> <library 2 tag file name> 
  #  <chromosome length file name> <config file name> <project name>
  #<CCAT path>/bin/CCAT  <ChIP library read file name>  <control library read file name> 
  # <chromosome length file name>  <config file name>  <project name>
  #../bin/CCAT ES_CTCF_chr19.bed ES_GFP_chr19.bed genome_length_mm8_chr19.txt config_TF.txt ES_CTCF_chr19 > ES_CTCF_chr19.log
  #This command line will generate 4 files:
  #ES_CTCF_chr19.significant.peak: The peaks that are identified to be significant;
  #ES_CTCF_chr19.significant.region: The regions that contain significant peaks. Some regions may contain multiple peaks and the most significant peak will be reported;
  #ES_CTCF_chr19.top100000: The top 100000 peaks. Sometimes you need this file to take a look at the peaks below the cut-off threshold. The number of output peaks is configurable in the config file;
  #ES_CTCF_chr19.log: the log-file containing the intermediate information.
  
  
  
  my $cmd = join(' ', $self->program_file,
                      $signal_bed,
                      $control_bed,
                      $chr_file,
                      $conf_file,
                      $out_prefix,
                      '>',
                      $out_prefix.'.log');
  #WARNING!
  #Do not remove this redirect to log file, as th CCAT binary will expect that 
  #it is directly attached to a terminal and fails with:
  #MSG: Child process died with signal 11, without coredump
  #Error:  Inappropriate ioctl for device
                     
  warn "Running:\t$cmd\n" if $self->debug;
  run_system_cmd($cmd);  #This did not cause failure when failed?
  return;
}

sub DESTROY{
  my $self = shift;
  unlink($self->{conf_file}) if $self->{conf_file};  
}

sub init_significant_region_file {
  return shift->_init_ccat_file('significant.region');  
}

sub _init_ccat_file{
  my $self   = shift;
  my $suffix = shift or throw('Must define a ccat file suffix argument');
  $self->{out_file_handle} = open_file($self->out_file_prefix.'.'.$suffix);
  return $self->out_file_handle;
}

#Filtering on FDR should be done here or in caller, so we can count
#what we are filtering out?

sub parse_significant_region_record {
  my $self = shift;
  my ($line, $fhash);
  
  if(! eval { $line = $self->out_file_handle->getline; 1}){
    throw("Failed to getline from out_file_handle, maybe you need to init_bed_file first.\n$@");
  }
 
  if(defined $line){
    chomp $line;
    #my ($seqid, $summit, $start, $end, $signal_reads, 
    #    $ctrl_reads, $fold, $fdr) = split("\t", $line);
    my ($seqid, $summit, $start, $end, undef, 
        undef, $fold_change, $fdr) = split("\t", $line);
        
    return if $fdr > $self->fdr_threshold;

    $fhash = 
     {-start      => $start,
      -end        => $end,
      -strand     => 0,
      -score      => $fold_change,
      -summit     => $summit,
      -seq_region => $seqid   };
  }
  
  return $fhash;
}



1;

