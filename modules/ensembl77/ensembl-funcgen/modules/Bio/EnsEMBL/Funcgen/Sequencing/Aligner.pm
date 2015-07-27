=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::Aligner

=head1 DESCRIPTION

Base class for all sequence alignement analyses. Provides generic constructor method, 
along with file validation and some accessor methods.

=cut

package Bio::EnsEMBL::Funcgen::Sequencing::Aligner;

use warnings;
use strict;

use File::Basename                         qw( fileparse );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd );

#todo add force_gunzip option

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  my ($prog_file, $prog_params, $query_file, $target_file, 
      $out_dir, $format, $batch_job, $debug) =
    rearrange(['PROGRAM_FILE', 'PARAMETERS', 'QUERY_FILE', 'TARGET_FILE', 
               'OUTPUT_DIR', 'OUTPUT_FORMAT', 'BATCH_JOB', 'DEBUG'], @_);
      
  if(! ($prog_file && $query_file && $target_file)){
    throw("Some mandatory parameters are not met:\n\t".
      join("\n\t", ("-PROGRAM_FILE   => $prog_file,", 
                    "-QUERY_FILE     => $query_file,",
                    "-TARGET_FILE    => $target_file")));
  }

  #todo Move this analysis handling to BaseDB

  if(! -f $prog_file){
    my $which_file = run_backtick_cmd("which $prog_file");
    
    if(! defined $which_file){
      throw("Program file does not exist or is not a file:\n\t$prog_file"); 
    }
    else{ #redefine it to give the full path for clarity
      $prog_file = $which_file;  
      chomp $prog_file;
    } 
  }


  if($batch_job){
  
    if(! defined $ENV{LSB_JOBINDEX}){
      throw('Batch job mode is specified, '.
        'but this job does not have an LSB_JOBINDEX environment variable available');  
    }  
   
    #Now alter query file wrt LSB_INDEX
    #This is based on the default split chunk suffix length of 4
    #unlikely we will get > 9999 chunks 
    $query_file .= '_'.sprintf("%04d" ,($ENV{LSB_JOBINDEX} - 1));
  }
    
  #It looks like links are no longer files in perl 5.14?
  if(! -f $query_file){    
    my $linked_file = readlink $query_file;
        
    if($linked_file ne $query_file){
      
      if(! -f $linked_file){
        throw("Query file link target is not a file:\n".
          "\tQuery file\t$query_file\n".
          "\tLink target\t$linked_file");      
      }       
    }
    else{
      throw("Query file does not exist or is not a file:\n\t$query_file");  
    }
  }

  if($query_file =~ /\.gz/){
    throw("Aligner input is compressed:\t$query_file\n".
      "Aligner expects uncompressed input. This is a safety feature to prevent\n".
      "gunzipping of fastq files which will destroy the ability to validate the checksum.\n".
      'Please zcat your fastq.gz file or gunzip if you are not concrened about further'.
      ' checksum validation.');  
  }


  #

  if(! -f $target_file){
    throw("Target file does not exist or is not a file:\n\t$target_file");  
  }

  my $input_dir;
  ($query_file, $input_dir) = fileparse($query_file);
  warn "Input dir:\t$input_dir\n" if $debug;

  if(defined $out_dir){ 
    if(! -d $out_dir){
      throw("Output directory does not exist:\n\t".$out_dir);  
    }
    
    $self->{output_dir}      = $out_dir;
  }
  else{
    $self->{output_dir}      = $input_dir;
  }
  
  #if($format){
  #  
  #  if(! grep(/^${format}$/, @{$self->valid_output_formats}) ){
  #    throw("Output format($format) is not valid. Please specify one of:\t".
  #      join("\t", @{$self->valid_output_formats}));  
  #  }
  #  
  #  #This assumes valid_output_formats are in order of general preference
  #  #i.e. most optimsed/performant/binary first
  #  $self->{output_format} = $self->valid_output_formats->[0]; 
  #  #}
  #else{
  #  $self->{output_format} = $self->valid_input_formats
  #}
  
  $self->{input_dir}    = $input_dir;
  $self->{program_file} = $prog_file;
  $self->{query_file}   = $query_file;
  $self->{parameters}   = defined $prog_params  ? $prog_params : ''; #To avoid warnings
  $self->{target_file}  = $target_file;
  $self->{debug}        = $debug;
  
  return $self;
}

sub query_file     { return shift->{query_file};    }
sub target_file    { return shift->{target_file};      }
sub program_file   { return shift->{program_file};  }
sub parameters     { return shift->{parameters};    }
sub output_dir     { return shift->{output_dir};    }
sub input_dir      { return shift->{input_dir};     }
sub debug          { return shift->{debug};         }
#sub output_format  { return shift->{output_format}; }

#abstract methods which muct be defined in the the specific aligner sub class  
#sub run
#sub valid_output_formats

1;

