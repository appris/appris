=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller

=head1 DESCRIPTION

Base class for all peak caller analyses. Provides generic constructor method, 
along with file validation, accessor methods and a simple feature storing.

=cut

package Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped         gunzip_file 
                                               file_suffix_parse  open_file
                                               convert_strand
                                               run_backtick_cmd );
                                               
my %half_open_formats = (bed => 1);

#To be over-ridden in the sub class with hardcoded values
sub out_file_types{ return shift->{out_file_types}; }
sub input_formats{  return shift->{input_formats};  }
#sub requires_control {  }
#Currently only used in the caller to identify which files need generating
#Omitted here for safety, due to boolean return type

#output format should really only be used to define cmdline switches for
#output file if they exist

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  my ($prog_file, $prog_params, $align_file, $control_file, 
    $out_dir, $out_file_prefix, $reload, $rerun, $no_load,
    $is_half_open, $convert_half_open, $output_format, $debug) =
    rearrange(['PROGRAM_FILE', 'PARAMETERS', 'ALIGN_FILE', 'CONTROL_FILE',
      'OUT_DIR', 'OUT_FILE_PREFIX', 'RELOAD', 'RERUN', 'NO_LOAD',
      'IS_HALF_OPEN', 'CONVERT_HALF_OPEN', 'OUTPUT_FORMAT', 'DEBUG'], @_);
      
  if(! ($prog_file && (defined $prog_params) && 
        $align_file && $out_dir)){
    #this will give an undef warning
    throw("Some of the folllowing mandatory parameters are not met:\n\t".
      "-PROGRAM_FILE, -ALIGN_FILE, -PARAMETERS or -OUT_DIR\n@_");
  }

  if(! -d $out_dir){
    throw("-out_dir does not exist:\t".$out_dir);  
  }

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

  if(! defined $align_file){
    throw('-ALIGN_FILE is a mandatory parameter');  
  }

  if(! defined $control_file){
    warn "Running $prog_file without and control file\n"; #omit this warn? 
  }
 
  if((! defined $is_half_open) && 
    exists $half_open_formats{$output_format}){
    $is_half_open = $half_open_formats{$output_format};   
  }
  
  #Set default out_file_type here
  if(! defined $self->out_file_types){
    $self->{out_file_types} = [$output_format];
    #Don't test $output_format, as it is only strictly required
    #if we specify an output prefix and the subclass
    #requires the format to build the output file path
    #or we have not specified an outfile path  
  }

  $self->{control_file}      = $control_file; 
  $self->{program_file}      = $prog_file;
  $self->{align_file}        = $align_file;
  $self->{parameters}        = defined $prog_params  ? $prog_params : ''; #To avoid warnings
  $self->{out_dir}           = $out_dir;
  $self->{out_file_prefix}   = $out_file_prefix;
  $self->{output_format}     = $output_format;
  $self->{is_half_open}      = $is_half_open;
  $self->{convert_half_open} = $convert_half_open;
  $self->{debug}             = $debug;
  
  $self->init_files;
  return $self;
}

#subtlety here to the gzip status
#SWEmbl will handle both gzipd or both unzipd
#but not mixed gzip/unzipd

#todo handle $unzip, or pass back gzipped status
#what about no unzip defined unzipped files, do we then gzip?

#do we still want a reload mode, where we don't care about the input
#files, so long as we have a valid out_file_prefix

#We never really want to unzip! As this may spoil checksums!!

#Some programs require an output file param others
# an output prefix, and some have > 1 output 
#Issue here is post-processing the files, we need to know which formats we want
#and hence which methods to call
    

sub init_files {
  my ($self, $unzip) = @_;
  
  if(! defined $self->{file_info}){
  
    my ($align_file, $gzip_align) = gunzip_file($self->align_file);
    my ($file_prefix, $suffix)    = file_suffix_parse($align_file);   
    my $gzip_control              = 0;
    my $control_file              = $self->control_file;
  
    if($control_file){
      ($control_file, $gzip_control) = gunzip_file($control_file);
    }
  
    if(! defined $self->out_file_prefix){
      $self->{out_file_prefix} = $file_prefix;  
    }
    
    $self->{out_file_prefix}  = $self->out_dir.'/'.$self->out_file_prefix;#.'.'.$self->output_format;
    $self->{file_info} = [$align_file, $self->{out_file_prefix}, $suffix,
                          $control_file, $gzip_align, $gzip_control];
  }

  return $self->{file_info};
}

sub file_info         { return shift->{file_info};         }
sub control_file      { return shift->{control_file};      }
sub align_file        { return shift->{align_file};        }
sub program_file      { return shift->{program_file};      }
sub parameters        { return shift->{parameters};        }
sub out_dir           { return shift->{out_dir};           }
sub out_file_prefix   { return shift->{out_file_prefix};   }
sub output_format     { return shift->{output_format};     }
sub out_file_handle   { return shift->{out_file_handle};   }
sub out_file_header   { return shift->{out_file_header};   } 
sub is_half_open      { return shift->{is_half_open};      }
sub convert_half_open { return shift->{convert_half_open}; }   
sub debug             { return shift->{debug};             }

#This attribute should be defined directly in the subclass constructor  
#or the method should be redefined in subclass to take file_type arg if required. 
#This is requiref for reload mode. (i.e. without calling run)
sub out_file          { return shift->{out_file};          }

  
#File parsing and DB loading may be polluting this PeakCaller code
#but is not mandatory functionality, so we will allow this as is quite useful
#as it can enable DB load functionality outside of the pipeline
#This is completely agnostic toward the type of code ref
#so it could be a DB store function or something completely different!
#Move some of this to SeqTools?
#Then use SeqTools directly in RunPeaks instead.


#todo
#Change all init/parse/out_file_handle methods to be generic and take 
#an optional file_type argument?
#This would mean the arg would be superfluous for those PeakCallers
#which don't use it. None of these methods take any args other than fle_type at present
#So this is not currently an issue
#These generic methods would then despatch to the correct method?
#How much of this can be moved in here, vs PeakCaller specific functionality
#i.e. file_type to suffix conversion
#Leave thsi for now

sub process_features {
  my $self  = shift;
  my $params = shift;
  assert_ref($params, 'HASH', 'Params');
  my ($file_type, $process_cref, $process_args) = 
    rearrange([qw( file_type processor_ref processor_args )], %$params);
  
  if(! defined $file_type){
    #If there is only one, assume we want that, else die
    if(scalar @{$self->out_file_types} == 1){
      $file_type = $self->out_file_types->[0];
    }
    else{
      throw('Unable to identify a file type to process. PeakCaller does not '.
        'have a unique out file type, please specify a -file_type in the '.
        'process_features parameter hash');  
    }
  }
  
  assert_ref($process_args, 'ARRAY', 'Process feature method arguments');
  assert_ref($process_cref, 'CODE',  'Process feature method');
  my $parse_method = 'parse_'.$file_type.'_record';
  my $init_method  = 'init_'.$file_type.'_file';
  
  if(! $self->can($parse_method)){
    throw(ref($self)."::${parse_method} is not a valid method.".
    ' This class is missing a method or does not have a supported output format');      
  }
   
  if(! $self->can($init_method)){
    throw(ref($self)."::${init_method} is not a valid method.".
    ' This class is missing a method or does not have a supported output format');      
  }
  
  
  warn "Processing file:\t".$self->out_file($file_type) if $self->debug;
  $self->$init_method;
  my ($feature_hash, $feature_cnt, $retval, @retvals);
  
  while( ($feature_hash = $self->$parse_method) &&
         defined $feature_hash ){          
    $retval = $process_cref->(@$process_args, $feature_hash);     
    push @retvals, $retval if defined $retval;
    $feature_cnt++;  
  }   
   
  return ($feature_cnt, \@retvals);
}  

#For common formats there is a potential to wrap IO parsers

sub init_bed_file {
  my $self = shift;
    
  $self->{out_file_handle} = open_file($self->out_file);
  #This will validate exists
  
  #check we have a header
  my $prevpos = $self->out_file_handle->getpos;
  my $line;
  
  while(($line = $self->out_file_handle->getline) &&
         defined $line){
    
       if($line =~ /^(browser)/){
      $prevpos = $self->out_file_handle->getpos;
    }
    elsif($line =~ /^track/){
      $prevpos = $self->out_file_handle->getpos;
      $self->{out_file_header} = $line;
    }
    else{
      last;  
    }
  }
  
  $self->out_file_handle->setpos($prevpos); 
  return $self->out_file_handle;
}


#Would be nice to define a spec for what we want returning
#rather than the standard field names, i.e. what this translates to
#in our constructor of choice. Simplest way would be to pass a list of 
#of key values in the field order, to override the default ones
#This will not handle overloading fields, maybe this is something autosql can handle?

sub parse_bed_record {
  my $self = shift;
  my ($line, $fhash);
  
  if(! eval { $line = $self->out_file_handle->getline; 1}){
    throw("Failed to getline from out_file_handle, maybe you need to init_bed_file first.\n$@");
  }
 
  if(defined $line){
    my ($seqid, $start, $end, $name, $score, $strand) = split(/\s+/, $line);
    $strand = convert_strand($strand);
    
    #Need to handle seq_id conversion if the bed is UCSC format
     
    if($self->is_half_open && $self->convert_half_open){
      $start += 1;
    }
      
    #Handle half_open here, but slice in Caller which knows about the DB
    $fhash = 
     {-start      => $start,
      -end        => $end,
      -name       => $name,
      -strand     => $strand,
      -score      => $score,
      -seq_region => $seqid   };
  }

  return $fhash;
}

1;

