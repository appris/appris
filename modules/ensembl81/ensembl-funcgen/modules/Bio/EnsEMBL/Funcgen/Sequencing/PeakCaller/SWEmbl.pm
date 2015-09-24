=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::SWEmbl

=head1 DESCRIPTION

Runs SWEmbl peak caller and optionally parses and processes output.
This is in the Hive namespace, but is not dependnant on any Hive modules and can 
be run as a stand alone job outside of the hive infrastructure.

=cut

package Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::SWEmbl;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd 
                                               run_backtick_cmd
                                               open_file );
use base qw( Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller );#Does not import

my %fswitches = 
 (
  bed   => '-B',
  sam   => '-S',
  maq   => '-M',
  eland => '-E',
  bam   => '-F',
 );

sub requires_control { return 0;                     }
sub input_formats    { return ['bam', 'sam', 'bed']; }
sub out_file_types   { return ['txt'];               }


sub new {
  my $caller          = shift;
  my $class           = ref($caller) || $caller;
  my $self            = $class->SUPER::new(@_);  # @_ will over-ride -is_half_open
  #This needs setting here, ss we may need it for reload (i.e. without calling run)      
  $self->{out_file} ||= $self->out_file_prefix.'.'.$self->out_file_types->[0];
  return $self;
}

#todo SWEmbl handles gzipped file, so just use is_gzipped and pass -z

sub run {   
  my $self = shift;
  my ($align_file, $out_file_prefix, $suffix, $control_file, 
    $gzipped_align, $gzipped_control) = @{$self->file_info};
    
  my $format_switch = $fswitches{lc($suffix)}; # auto-vivifies
  throw("$suffix is not recognised as a valid SWEmbl input format") if(! $format_switch);
        
  # Slight work around here as is_half_open for SWEmbl is defined by input, not output format
  # This is normally set in PeakCaller::init_files via new.
  $self->{is_half_open} = 1 if lc($suffix) eq 'bed';

  #Set -z 
  my $compressed = $gzipped_align ? ' -z ' : '';
        
  if($control_file &&
     ($gzipped_control != $gzipped_align)){
    throw("-i(input_file) and -r(eference/input file cannot have mixed compression states:\n\t".
      $align_file."\n\t".$control_file);  
  }
 
  #Sometimes SWEmbl fails to open output file if it exists already
  unlink($self->out_file);
 
  my $cmd = $self->program_file." $format_switch -i $compressed ".
    $self->align_file.' '. $self->parameters.' -o '.$self->out_file;
  $cmd .= " -r ".$self->control_file if $self->control_file;
  
  warn "Running:\t$cmd\n" if $self->debug;
  #This did no cause failure when failed
  run_system_cmd($cmd);
     
  return;
}


#This is not strictly a PeakCaller thing
#But the header extraction and sort are specific enough to be in here
#as format is SWEmbl specific

#For other PeakCallers, handling the file_type and converting it to a file
#a file suffic is currently internal (see CCAT)

sub filter_max_peaks{
  my $self      = shift;
  my $max_peaks = shift;
  my $file_type = shift; 
  
  #Extract header
  my $out_file = $self->out_file;
  my $fh       = open_file($out_file);
  my $header   = $self->parse_txt_header($fh);
  close($fh);
  
  my $header_file = $out_file.'.header';
  my $header_fh   = open_file($header_file, '>');
  print $header_fh join("\n", @$header)."\n";
  close($header_fh);

  #Create headerless, sorted, filtered file
  my $cmd = "tail -n +".scalar(@$header)." $out_file | sort -k 7nr,7nr | head -n $max_peaks ".
    "| sort -k 1,2n > ".$out_file.'.filtered';
  warn $cmd if $self->debug;
  run_system_cmd($cmd);
    
  #Sanity check we have the file with the correct number of lines
  $cmd = "wc -l $out_file.filtered | awk '{print \$1'}";
  warn $cmd if $self->debug;
  my $filtered_peaks = run_backtick_cmd($cmd);
    
  if($max_peaks != $filtered_peaks){ 
    throw("Expected $max_peaks in filtered bed file, but found $filtered_peaks:\n\t".$out_file);  
  } 
  
  #Create final headered filtered file    
  $cmd = "cat $out_file.header $out_file.filtered > $out_file"; 
  warn $cmd if $self->debug;
  run_system_cmd($cmd);
      
  unlink("${out_file}.header", "${out_file}.filtered");    
  return;     
}

#TODO revise how we create and handle filehandles
# sub out get_txt_filehandle, set as attr, or exclusively pass between methods?
# out_file_handle is currently generic and does not support multiple file_types
# Also change all this to generic parse_header/record method, and pass an optional file types
# 

sub parse_txt_header{
  my $self   = shift;
  my $txt_fh = shift;
  my @header;
  #assert_ref($fh, 'FileHandle');  
  #or create file handle?

  #check we have a header
  my $prevpos = $txt_fh->getpos;
  my $line;
  
  while(($line = $txt_fh->getline) &&
         defined $line){
    chomp $line; 

    if($line =~ /^(#|Region)/){
      push @header, $line;
      $prevpos = $txt_fh->getpos;
    }
    else{
      last;  
    }
  }
  
  $txt_fh->setpos($prevpos);   
  return \@header;
}

sub init_txt_file {
  my $self   = shift;
  my $txt_fh = open_file($self->out_file);
  $self->parse_txt_header($txt_fh);  
  $self->{out_file_handle} = $txt_fh;
  return $self->out_file_handle;
}


sub parse_txt_record {
  my $self = shift;
  my ($line, $fhash);
  
  if(! eval { $line = $self->out_file_handle->getline; 1}){
    throw("Failed to getline from out_file_handle, maybe you need to init_bed_file first.\n$@");
  }
 
  if(defined $line){
    #my ($seqid, $start, $end, $cnt, $length, $uniq_pos, $score, $ref_cnt, $max_cvg, $summit) = split(/\s+/, $line);
    chomp $line;

    my ($seqid, $start, $end, undef, undef, undef, 
        $score, undef, undef, $summit) = split(/\s+/, $line);
    if( (!defined $seqid) || ($seqid eq '') ){
      throw("FILE: " . $self->out_file . "\nLINE:\n$line");
    }

    if($summit){
      $summit = int($summit + 0.5);#Rounds to nearest integer as int rounds down     
    }
     
    if($self->is_half_open && $self->convert_half_open){
      $start += 1;
    }
    
    #Handle half_open here, but slice in Caller which knows about the DB
    $fhash = 
     {-start      => $start,
      -end        => $end,
      -strand     => 0,
      -score      => $score,
      -summit     => $summit,
      -seq_region => $seqid   };
  }

  return $fhash;
}


1;

