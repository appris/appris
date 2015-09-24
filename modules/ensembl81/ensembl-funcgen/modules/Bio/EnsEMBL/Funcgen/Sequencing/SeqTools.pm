=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::SeqTools

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;

use warnings;
use strict;

use Net::FTP;
use feature qw(say);
use DBI     qw(:sql_types);

use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::Funcgen::InputSubset;
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Utils::SqlHelper;

use File::Temp                             qw( tempfile );
use File::Basename                         qw( dirname basename );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( dump_data
                                               run_system_cmd
                                               run_backtick_cmd
                                               write_checksum
                                               validate_path
                                               check_file
                                               open_file
                                               validate_package_path
                                               which_path );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqQC; #run_QC

use base qw( Exporter );
use vars qw( @EXPORT );

@EXPORT = qw(
  _init_peak_caller
  _run_peak_caller
  convert_bam_to_sam
  convert_sam_to_bed
  create_and_populate_files_txt
  download_all_files_txt
  explode_fasta_file
  get_files_by_formats
  load_experiments_into_tracking_db
  merge_bams
  modify_files_txt_for_regulation
  post_process_IDR
  pre_process_IDR
  process_bam
  randomise_bed_file
  run_aligner
  run_IDR
  run_peak_caller
  split_fastqs
  validate_sam_header
  write_chr_length_file
  );

#This is designed to act as an object and as a standard package
#run_align_tools will call constructor which will validate that a run mode has been passed
#before cacheing all params and retuning
#then script will call run, which will call the method defined by the run mode,
#passing all the params which wil passed to the contructor
#or optionall, taking arguments, to overide those passed to the constructor?
#then will undef method, for safety


#todo add QC methods (and integrate into other methods turn off with -no_qc)
#todo add split/preprocess methods
#todo add merge method in rmdups
#todo add sort/filter methods? unmapped are normally remove in process_sam_bam
#

#process_sam_bam filters unmapped, sorts and rmdups
#merge only rmdups



#add debug and helper to all these methods


#to do, allow qc methods to be used in here too, to restrict what is run
#from the default qc_type?
#No, just make people run these separately?


#TODO
# 1 Change all paramter handling to use rearrange only? i.e. no unamed args!

sub run_aligner{
  my ($aln_pkg, $aln_params, $skip_qc, $batch_job, $debug) =
    rearrange( [qw(aligner aligner_params skip_qc batch_job debug) ], @_);
  assert_ref($aln_params, 'ARRAY', 'Aligner params');

  throw('-aligner parameter has not been specifed')  if ! defined $aln_pkg;
    #This is being run when we call require, and exiting as we have no args
    #$0 = $real_path;
    #pod2usage(-exitval => 1,
    #          -message => '-aligner parameter must be defined');
  #}

  if($aln_pkg !~ /::/){
    warn "Full $aln_pkg namespace was not specified, defaulting to:\t".
      "Bio::EnsEMBL::Funcgen::Sequencing::Aligner::$aln_pkg\n";
    $aln_pkg = "Bio::EnsEMBL::Funcgen::Sequencing::Aligner::$aln_pkg";
  }

  validate_package_path($aln_pkg);
  my $aligner = $aln_pkg->new(@$aln_params);

  #Don't really need this, although might be nice to show the args passed in the
  #error output in Aligner
  #if(! defined $aligner){
  #  throw("Failed to create $aln_pkg from params:\n",
  #    join(' ', @$aln_params));
  #}


  #Implement QC in here, wht are we going to return here?
  #return ($aligner->run, $qc_results);
  #This may cause problems if the caller does this
  #if($run_aligner($aln_pkg, $aln_params, 1)){ ... }
  #This test will return false, as this will test the last value in the array
  #in this context
  #if we return @results, where the null value is omited
  #then all is good.


  #Need to eval this, so we don't call QC on undef file if it has failed?
  #No failure should have died by now
  my @results = ($aligner->run); #This should return the output file

  if( ! $skip_qc){

  }


  return \@results;
}


#todo implement this in PreprocessFastqs
#mv md5 checking in here too? as that is integrated into check_file
#Just force use of gz files here to simplify zcat and md5 checking?
# add qc thresholds?

#is split splitting across fastq records?
#lines do no appear to be wrapped with \n?
#Or is this becuase all multi line records are 4 lines, so we just specify
#a multiple of 4!

sub split_fastqs{
  my ($files, $out_prefix, $out_dir, $work_dir,
      $check_sums, $merge, $chunk_size, $skip_qc, $debug) = rearrange(
      [qw( files out_prefix out_dir work_dir
       merge chunk_size skip_qc debug) ], @_);

  assert_ref($files, 'ARRAY', '-files');

  if(! (@$files &&
        (grep {!/fastq.gz$/} @$files) )){
    throw('-files must be an array ref of gzipped fastq files');
  }

  throw('-out_prefix is not defined') if ! defined $out_prefix;

  if(! -d $out_dir){
    throw("-out_dir $out_dir is not a valid output directory");
  }

  if(! defined $work_dir){
    $work_dir = $out_dir;
  }
  elsif(! -d $work_dir){
    throw("-work_dir $work_dir is not a valid work directory");
  }

  if($check_sums){
    assert_ref($check_sums, 'ARRAY', '-check_sums');

    if(scalar(@$check_sums) != scalar(@$files)){
      throw(scalar(@$files).' -files have been specific but only '.scalar(@$check_sums).
        " -check_sums have been specified\nTo ensure input validation these must ".
        "match, even if undef checksums have to be specified");
    }
  }

  $chunk_size ||= 16000000;#Optimised for ~ 30 mins bwa alignment bjob

  if($chunk_size % 4){
    throw("Chunk size($chunk_size) is not a multiple of 4. This is required to ensure safe splitting of 4 line fastq records.");
  }

  my (@fastqs, %params, $throw);

  foreach my $i(0..$#{$files}){
    my $found_path;
    %params = ( debug => $debug, checksum => $check_sums->[$i] );

    #Hmm, no undef checksum here means try and find one in a file

    #Look for gz files too,
    #we can't do a md5 check if we don't match the url exactly
    if(! eval { $found_path = check_file($files->[$i], 'gz', \%params); 1}){
      $throw .= "Failed to check_file:\t".$files->[$i]."\n$@";
      next;
    }
    elsif(! defined $found_path){
      $throw .= "Could not find fastq file, is either not downloaded, has been deleted or is in warehouse:\t".
        $files->[$i]."\n";
      #Could try warehouse here?
    }
    elsif($found_path !~ /\.(?:t){0,1}gz$/o){
      #use is_compressed here?
      #This will also modify the original file!
      throw("Found unzipped path, aborting as gzipping may invalidate any further md5 checking:\t$found_path");
      #run_system_cmd("gzip $found_path");
      #$found_path .= '.gz';
    }

    push @fastqs, $found_path;
  }

  throw($throw) if $throw;

  my (@results, @qc_results);



  #if qc fails here, we still split?
  #we need a way to signify QC failure easily without
  #having to test hash keys?
  #This could be an array of booleans?
  #so we would return \@new_fastqs, \@pass_fail_booleans, \@qc_hashes

 #FastQC in here


  #This currently fails as it tries to launch an X11 window!

  ### RUN FASTQC
  #18-06-10: Version 0.4 released ... Added full machine parsable output for integration into pipelines
  #use -casava option for filtering

  #We could set -t here to match the number of cpus on the node?
  #This will need reflecting in the resource spec for this job
  #How do we specify non-interactive mode???
  #I think it just does this when file args are present

  #Can fastqc take compressed files?
  #Yes, but it seems to want to use Bzip to stream the data in
  #This is currently failing with:
  #Exception in thread "main" java.lang.NoClassDefFoundError: org/itadaki/bzip2/BZip2InputStream
  #Seems like there are some odd requirements for installing fastqc
  #although this seems galaxy specific
  #http://lists.bx.psu.edu/pipermail/galaxy-dev/2011-October/007210.html

  #This seems to happen even if the file is gunzipped!
  #and when executed from /dsoftware/ensembl/funcgen
  #and when done in interative mode by loading the fastq through the File menu

  #This looks to be a problem with the fact that the wrapper script has been moved from the
  #FastQC dir to the parent bin dir. Should be able to fix this with a softlink
  #Nope, this did not fix things!

  warn "DEACTIVATED FASTQC FOR NOW:\nfastqc -f fastq -o ".$out_dir." @fastqs";
  #run_system_cmd('fastqc -o '.$self->output_dir." @fastqs");


  #todo parse output for failures
  #also fastscreen?

  warn("Need to add parsing of fastqc report here to catch module failures");

  #What about adaptor trimming? and quality score trimming?
  #FASTX? quality_trimmer, clipper (do we have access to the primers?) and trimmer?



  #For safety, clean away any that match the prefix
  run_system_cmd('rm -f '.$work_dir."/${out_prefix}.fastq_*", 1);
  #no exit flag, in case rm fails due to no old files

  my @du = run_backtick_cmd("du -ck @fastqs");
  (my $pre_du = $du[-1]) =~ s/[\s]+total//;

  my $cmd = 'zcat '.join(' ', @fastqs).' | split --verbose -d -a 4 -l '.
    $chunk_size.' - '.$work_dir.'/'.$out_prefix.'.fastq_';
  #$self->helper->debug(1, "Running chunk command:\n$cmd");
  warn "Running chunk command:\n$cmd\n" if $debug;

  my @split_stdout = run_backtick_cmd($cmd);
  (my $final_file = $split_stdout[-1]) =~ s/creating file \`(.*)\'/$1/;

  if(! defined $final_file){
    throw('Failed to parse (s/.*\`([0-9]+)\\\'/$1/) final file '.
      ' from last split output line: '.$split_stdout[-1]);
  }

  #Get files to data flow to individual alignment jobs
  my @new_fastqs = run_backtick_cmd('ls '.$work_dir."/${out_prefix}.fastq_*");
  @new_fastqs    = sort {$a cmp $b} @new_fastqs;

  #Now do some sanity checking to make sure we have all the files
  if($new_fastqs[-1] ne $final_file){
    throw("split output specified last chunk file was numbered \'$final_file\',".
      " but found:\n".$new_fastqs[-1]);
  }
  else{
    $final_file =~ s/.*_([0-9]+)$/$1/;
    $final_file  =~ s/^[0]+//;

    #$self->debug(1, "Matching final_file index $final_file vs new_fastq max index ".$#new_fastqs);
    warn "Matching final_file index $final_file vs new_fastq max index ".$#new_fastqs."\n" if $debug;


    if($final_file != $#new_fastqs){
      throw('split output specified '.($final_file+1).
        ' file(s) were created but only found '.scalar(@new_fastqs).":\n".join("\n", @new_fastqs));
    }
  }

  #and the unzipped files are at least as big as the input gzipped files
  @du = run_backtick_cmd("du -ck @new_fastqs");
  (my $post_du = $du[-1]) =~ s/[\s]+total//;

  #$self->helper->debug(1, 'Merged and split '.scalar(@fastqs).' (total '.$pre_du.'k) input fastq files into '.
  #  scalar(@new_fastqs).' tmp fastq files (total'.$post_du.')');
  warn 'Merged and split '.scalar(@fastqs).' (total '.$pre_du.'k) input fastq files into '.
    scalar(@new_fastqs).' tmp fastq files (total'.$post_du.")\n" if $debug;

  if($post_du < $pre_du){
    throw("Input fastq files totaled ${pre_du}k, but output chunks totaled only ${post_du}k");
  }

  return (\@new_fastqs);#, \%qc_results;
}




#todo
# 1 add a DESTROY method to remove any tmp sorted files which may persist after an
#   ungraceful exit. These can be added to a global $main::files_to_delete array
#   which should then also be undef'd in DESTROY so they don't persisnt to another instance

#We could use the existing Bio::SamTools package but:
#1 This will add an extra requirement
#2 This will need to be isolated in a hive/analysis only module
#3 It doesn't appear to support merge operations
#4 It wouldn't support the piping/greping we do to filter the data


sub merge_bams{
  my $outfile     = shift;
  my $sam_ref_fai = shift;
  my $bams        = shift;
  my $params      = shift || {};
  assert_ref($bams, 'ARRAY', 'bam files');

  if(! scalar(@$bams)){
    throw('Must provide an arrayref of bam files to merge');
  }

  my $out_flag = '';

  if(! defined $outfile){
    throw('Output file argument is not defined');
  }
  elsif($outfile !~ /\.(?:bam|sam)$/xo){
    #?: does not assign to $1
    $out_flag = 'b' if $1 eq 'bam';
    throw('Output file argument must have a sam or bam file suffix');
  }

  assert_ref($params, 'HASH');
  my $debug       = (exists $params->{debug})        ? $params->{debug}       : 0;
  my $no_rmdups   = (exists $params->{no_rmdups})    ? $params->{no_rmdups}   : undef;
  my $no_checksum  = (exists $params->{no_checksum}) ? $params->{no_checksum} : undef;
  warn "merge_bam_params are:\n".dump_data($params)."\n" if $debug;

  #For safety we need to validate all the bam headers are the same?
  #or at least no LN clashed for the same SN?
  #Must all be subsets of sam_header if specified, and reheader output with
  #sam_header if defined
  #else, with the merge of all the input headers?
  #This later option would permit merges of redunant headers if the SN values
  #are not identical for the same sequence
  #force sam header for safety?
  #For now, let just make sure they are identical
 
  my $view_header_opt;

  for(@$bams){
    my $tmp_opt = validate_sam_header($_, $sam_ref_fai, 1, $params);
    $view_header_opt = $tmp_opt if $tmp_opt;
  }

  #validate/convert inputs here?
  #just assume all aren bam for now
  my $cmd = '';

  #-u uncompressed BAM output for pipe (header remains in sam format)
  #-f force overwrite output
  #-h is include header in output, seems to be in sam format i.e. not binary if output is bam??
  # - To specify seding output to STDOUT for

  my $skip_merge = 0;

  if(scalar(@$bams) == 1){
    #samtools merge cannot handle a single input!
    #Instead it throws a seemingly completely unrelated error message:
    #Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users
    #  must provide the correct header with -h, or uses Picard which properly maintains
    #  the header dictionary in merging.

    #Rather than having the caller have to handle this, let's just do the expected thing here
    #and warn.
    $skip_merge = 1;
    warn 'Only 1 bam file has been specified, merge will be skipped, '.
      "otherwise file will be processed accordingly\n";
  }


  if((! $no_rmdups) || $view_header_opt){
    #-u uncompressed output for pipeing to rmdups or view
    $cmd = 'samtools merge -u - '.join(' ', @$bams).' | ' if ! $skip_merge;

    if( ! $no_rmdups ){
       #rmdup converts the header into binary format
       $cmd .= 'samtools rmdup -s ';
       $cmd .= $skip_merge ? $bams->[0].' ' : ' - ';

       if( $view_header_opt ){
         $cmd .=  ' - | '
       }
       else{
         $cmd .= $outfile;
       }
    }

    if($view_header_opt){
      #We only need to do this if the validate_sam_header
      #method identifed some of the bams without the relevant header
      #warn "Currently integrating fai header via samtools view, but it is more efficient to integrate is with sam format header in merge";
      $cmd .= "samtools view -t $sam_ref_fai -h${out_flag} - > $outfile";
    }
  }
  elsif(! $skip_merge){
    $cmd = "samtools merge $outfile ".join(' ', @$bams);
  }
  else{ #skip merge
    $cmd = 'cp '.$bams->[0].' '.$outfile;
  }


  #piping like this may cause errors downstream of the pipe to be missed
  #could we try doing an open on the piped cmd to try and catch a SIGPIPE?
  warn "Merging with:\n$cmd\n" if $debug;
  run_system_cmd($cmd);
  
  # samtools merge can create a truncated file which lacks an EOF marker
  # unfortunately this does not raise an error and so will not be caught above
  # samtools view -h will return $? == 2 here
  # however, this is non-optimal as it slows down this step for the 99% of files which merge correctly
  # samtools -H does not raise an error or a warning, 
  # samtools -H(b|u) does not raise and but does output a warning:
  # [bam_header_read] EOF marker is absent. The input is probably truncated.
  # So let's capture output here via backticks..eugh.
  # Backticks normally only capture STDOUT. So redirect STDERR to STDOUT before we discard STDOUT
  $cmd = "samtools view -Hb $outfile 2>&1 > /dev/null";
  my $uncaught_merge_error = run_backtick_cmd($cmd);

  if($uncaught_merge_error){
    throw("samtools merge appeared to create a truncated file:\n\t$uncaught_merge_error");
  }

  warn "Finished merge to $outfile" if $debug;

  if(! $no_checksum){
    write_checksum($outfile, $params);
  }

  return;
}




#TODO Implement multi-mapping filter
#ENCODE removed multimapping reads, probably by filtering based on presence of XA tag
#-n is not defined. This seems only to apply to paired reads?
#It's unclear exactly what bwa does here.
#Repetitive hits will be chosen randoml(y, and XA will be written for alternate mappings)
#This means some duplicate reads will likely be slipping through if
#they map to multiple locations
#  To filter (given bwa samse -n wasn't used)
#  -F 100 will remove non-primary mappings
#  -v XA will remove remaining primary mappings will alternative mapping present in the
#  XA field.  samtools view -F 100 -h in.bam | grep -v XA
#This only works for single end reads, and would potentially leave dangling reads if
#the other half of a pair did not have an XA tag. So you would have to grep out the QNAME (query/pair name)
#and re-filter on that.
#--> Implement and are_paired flag

#checksum in params here acts to check and write checksums
#checksum => undef tries to find a checksum file
#checksum => MD%STRING checks using string
#Probably need a new param here

# TODO Add max cpu to manage number of pipes?

sub process_sam_bam {
  my $sam_bam_path = shift;
  my $params       = shift || {};
  my $in_file;

  #undef checksum here mean try and find one to validate
  #but then we don't write one

  if(! ($in_file = check_file($sam_bam_path, undef, $params)) ){
    throw("Cannot find file:\n\t$sam_bam_path");
  }

  #Change this to use rearrange by prefixing keys with -
  assert_ref($params, 'HASH');
  my $out_file      = (exists $params->{out_file})              ? $params->{out_file}              : undef;
  my $sort          = (exists $params->{sort})                  ? $params->{sort}                  : undef;
  my $skip_rmdups   = (exists $params->{skip_rmdups})           ? $params->{skip_rmdups}           : undef;
  #Turn on checksum writing
  my $checksum      = (exists $params->{checksum})              ? 1                                : undef;
  my $fasta_fai     = (exists $params->{ref_fai})               ? $params->{ref_fai}               : undef;
  my $out_format    = (exists $params->{output_format})         ? $params->{output_format}         : undef;
  my $debug         = (exists $params->{debug})                 ? $params->{debug}                 : 0;
  my $filter_format = (exists $params->{filter_from_format})    ? $params->{filter_from_format}    : undef;
  my $force         = (exists $params->{force_process_sam_bam}) ? $params->{force_process_sam_bam} : undef;

  #sam defaults
  $out_format     ||= 'sam';
  my $in_format = 'sam';
  my $in_flag   = 'S';

  if($out_format !~ /^(?:bam|sam)$/){
    throw("$out_format is not a valid samtools output format");
  }

  if($in_file =~ /\.bam$/o){     # bam (not gzipped!)
    $in_format = 'bam';
    $in_flag   = '';
  }
  elsif($in_file !~ /\.sam(?:\.gz)*?$/o){ # sam (maybe gzipped)
    throw("Unrecognised sam/bam file:\t".$in_file);
  }

  #This is odd, we really only need a flag here
  #but we already have the filter_from_format in the params
  if(defined $filter_format &&
     ($filter_format ne $in_format) ){
    throw("Input filter_from_format($filter_format) does not match input file:\n\t$in_file");
  }


  if(! $out_file){
    ($out_file = $in_file) =~ s/\.${in_format}(?:.gz)*?$/.${out_format}/;

    if(defined $filter_format){
      $out_file =~ s/\.unfiltered//o;  #This needs doing only if is not defined
    }
  }

  #Sanity checks
  (my $unzipped_source = $in_file) =~ s/\.gz$//o;
  (my $unzipped_target = $out_file)     =~ s/\.gz$//o;

  if($unzipped_source eq $unzipped_target){
    #This won't catch .gz difference
    #so we may have an filtered file which matches the in file except for .gz in the infile
    throw("Input and output (unzipped) files are not allowed to match:\n\t$in_file");
  }

  if($filter_format){

    if($in_file !~ /unfiltered/o){
      warn("Filter flag is set but input file name does not contain 'unfiltered':\n\t$in_file");
    }

    if($out_file =~ /unfiltered/o){
      throw("Filter flag is set but output files contains 'unfiltered':\n\t$in_file");
    }
  }
  elsif(! $sort &&
        ($in_format eq $out_format) ){
    throw("Parameters would result in no change for:\n\t$in_file");
  }


  #Define and clean intermediate sorted files first
  (my $tmp_out = $in_file) =~ s/\.$in_format//;
  my $sorted_prefix = $tmp_out.'.sorted';
  $tmp_out .= ($sort) ? '.sorted' : '.tmp';

  # Simply over-write these
  #my $cmd = "rm -f $tmp_bam*";  # Is * to handle possible checksum files
  #warn $cmd."\n" if $debug;
  #run_system_cmd($cmd, 1); #no exit flag
  my $cmd;

  #Check header and define include option
  my $fasta_fai_opt = validate_sam_header($in_file, $fasta_fai, undef, $params);


  ### HANDLE SIMPLE REHEADER OR FILE COPY ###
  #my $reheader = 0;

  if((! ($filter_format || $sort)) &&
     ($out_format eq $in_format) ){
    #This could possibly be a reheader operation or simply a move

    if(! $fasta_fai_opt){

      if(! $force){
        # arguably we should just do this but it is likely the options are wrong
        # could provide a flag over-ride for this?
        throw('The options provided do not require any processing of the input file:'.
          "\n\t$in_file\nOther than copying to the output file destination:\n\t".
          $out_file."\nPlease check/revise your options or specify the force_process_sam_bam parameter");
      }
      else{
        $cmd = "cp $in_file $out_file"; 
      }
    }
    else{ #We simply want to reheader
      #in and out format are the same, so can just test in format

      if($in_format eq 'sam'){
        $cmd = "samtools view -h${in_flag} $fasta_fai_opt $in_file ";
      }
      else{ #must be bam
        throw('bam reheader is not yet supported as requires a sam format header file');
        #actually fai format is not yet being validated, so this will fail if we pass a sam header
        #as the $fasta_fai_opt will be -h (for merge) if it is in sam format
        #$cmd = 'samtools reheader $sam_fai_or_header $in_file && mv $in_file $out_file';
      }
    }
  }


  ### FILTER/SORT/CONVERT ###

  if(! $cmd){ #We want to do some filtering/sorting
    my $filter_opt = ($filter_format) ? '-F 4' : '';

    ### SIMPLE BAM TO SAM ###
    if(($out_format eq 'sam') &&
        $skip_rmdups          &&
        ! $sort){
      $tmp_out = $tmp_out.'sam';
      $cmd = "samtools view -h${in_flag} $filter_opt $fasta_fai_opt $in_file > $tmp_out";
      warn $cmd."\n" if $debug;
      run_system_cmd($cmd);
      run_system_cmd("rm -f $tmp_out");
    }
    else{ # FILTERING & SORTING 
      # Base view command to be piped to other commands
      # rmdups does not need this view unless there is a sort and filter in place
      # or a reaheader?
      # sort & rmdup only work on bam
      $cmd = "samtools view -hu${in_flag} $filter_opt $fasta_fai_opt $in_file "; 

      #dropped MT filtering, should be handled by blacklist
      #if($filter_format){
        #$cmd .= "-F 4 | ". #-F Skip alignments with bit set in flag (4 = unaligned)
        #  " grep -vE '^[^[:blank:]]+[[:blank:]][^[:blank:]]+[[:blank:]]+(MT|chrM)' "; #Filter MTs or any reference seq with an MT prefix
        #Could add blank after MT, just in case there are valid unassembed seq names with MT prefixes
        #Fairly safe to assume that all things beginning with chrM are MT or unassembled MT
      #}

      # -h include header
      #-u uncompressed bam (as we are piping)
      #-S SAM input    
      #-t  header file (could omit this if it is integrated into the sam/bam?)
      #- (dash) here is placeholder for STDIN (as samtools doesn't fully support piping).
      #This is interpreted by bash but probably better to specify /dev/stdin?
      #for clarity and as some programs can treat is as STDOUT or anything else?
      #-b output is bam
      #-m 2000000 (don't use 2G here,a s G is just ignored and 2k is used, resulting in millions of tmp files)
      #do we need an -m spec here for speed? Probably better to throttle this so jobs are predictable on the farm
      #We could also test the sorted flag before doing this?
      #But samtools sort does not set it (not in the spec)!
      #samtools view -H unsort.bam
      #@HD    VN:1.0    SO:unsorted
      #samtools view -H sort.bam
      #@HD    VN:1.0    SO:coordinate
      #We could add it here, but VN is mandatory and we don't know the version of the sam format being used?
      #bwa doesn't seem to output the HD field, not do the docs suggest which spec is used for a given version
      #mailed Heng Lee regarding this

 
      ### WRITE INTERMEDIATE BAM ###
      # This uses a tmp_bam intermediate to keep the cpu usage down
      # i.g. consider the potential pipe
      # view bam | sort | rmdups | view >sam
      # This would cause a spike in cpu usage which is likely not specified/expected
      # in the LSF resource, and so can cause failures
      # As we are not using IPC::Run we only ever get the exit status of the first command
      # so failures downstream would go uncaught
      $tmp_out = $tmp_out.'bam';
      my $rm_cmd = "rm -f $tmp_out";

      $cmd .= ($sort) ? ' | samtools sort - '.$sorted_prefix : ' > '.$tmp_out;
      warn $cmd."\n" if $debug;
      run_system_cmd($cmd);

      # TODO ENSREGULATION-200

      if($filter_format){

        if(! $skip_rmdups){
          # Removed after alignment as opposed from fastqs as we expect multiple
          # reads if they map across several loci but not necessarily at exactly
          # the same loci which indicates PCR bias
          #-s single end reads or samse (default is paired, sampe)
          $cmd = "samtools rmdup -s $tmp_out ";
        }

        if($out_format eq 'sam'){

          if($skip_rmdups){
            $cmd = "samtools view -h $tmp_out > $out_file";
          }
          else{
            $cmd .= "- | samtools view -h - > $out_file";
          }

        }
        elsif($skip_rmdups){
          $cmd = "mv $tmp_out $out_file";
          $rm_cmd = ''
        }
        else{
          $cmd .= $out_file;
        }
      }
      elsif($out_format eq 'bam'){
        # We know we have bam by now as we have done some sorting
        $cmd = "mv $tmp_out $out_file";
        $rm_cmd = '';
      }
      else{ #We need to convert to sam
        $cmd = "samtools view -h $tmp_out > $out_file";
      }

      # TODO ENSREGULATION-201

      warn $cmd."\n" if $debug;
      run_system_cmd($cmd);
      run_system_cmd($rm_cmd);
    }
  }
 
  if($checksum){  write_checksum($out_file, $params);  }
  return $out_file;
}


# Consider refactoring get_files_by_formats
# This handles g/unzipping and conversion from bam > sam > bed
# This also assumes that we only ever want to convert in this direction
# i.e. assumes bam /sam will always exist if we have bed.

# Encapsulate conversion logic in wrapper methods
# i.e. convert_bam_to_bed would call convert_bam_to_sam 
# then convert_sam_to_bed.
# This particular example is moot if we have bedtools bamToBed installed
# So check that first?
# Can still optionally keep or clean intermediates as required? No, that's 
# only possible with the more generic get_files_by_formats
# Let's just drop that functionality and let the caller handle it.

#sam params contains:
#ref_fai         => file_path
#filter_from_bam => 1
#could also support:
#include_MT   => 1,
#include_dups => 1,
#ignore header mismatch? (could do thi swith levels, which ignore supersets in fai?

#Currently hardoded for samse files name for sam and bed
#todo _validate_sam_params


#$formats should be in preference order? Although this doesn't break things, it will just return a non-optimal file format

#Slightly horrible method to manage acquisition and conversion of files between
#supported formats (not necessarily feature formats)
#all_formats is necessary such that we don't redundant process files which are on the same conversion path when we have filter_format set


#There is a possibility that the formats provided might not have the same root, and so
#filter_from_format may be invlaid for one
#In this case two method calls might be require, hence we don't want to throw here if we can't find a file

#This seems over-engineered! But we definitely need the ability to request two formats at the same time
#to prevent parallel requests for the same file

#Filtering will normally be done outside of this method, by the alignment pipeline
#however, we must support it here incase we need to refilter, or we get alignment files
#supplied outside of the pipeline


#what about if we only have the unfiltered file
#but we ask for filtered
#should we automatically filter?
#should we move handling 'unfiltered' to here from get_alignment_files_by_InputSet_formats?
#Is this too pipeline specific?
#what if some files don't use the 'unfiltered' convention?
#then we may get warnings or failures if the in_file and the out_file match
#would need to expose out_path as a parameter
#which would then need to be used as the in path for all subsequent conversion
#No this wouldn't work as it would change the in file to contain unfiltered
#which might not be the case.

#This is actually a generic method apart from the conversion paths
#which could be passed as code refs from here
#so we could move this back to EFGUtils

sub get_files_by_formats {
  my $path    = shift;
  my $formats = shift;
  my $params  = shift || {};
  assert_ref($formats, 'ARRAY');
  assert_ref($params, 'HASH');

  if(scalar(@$formats) == 0){
    throw('Must pass an Arrayref of file formats/suffixes in preference order e.g. [\'sam\', \'bed\']');
  }

  my %conversion_paths = ( bam => ['bam'],
                           sam => ['bam', 'sam'],
                           bed => ['bam', 'sam', 'bed'],
                           #we always need the target format as the last element
                           #such that we can validate the filter_format e.g. for bam
                           #if the path array only has one element, it must match the key
                           #and this constitues calling filter_bam
                           #or if filter_format not set, just grabbing the bam file

                           #This approach prevents being able to | bam sort/filters through
                           #to other cmds, so may be slower if we don't need to keep intermediate files?

                           #Could also have non-bam rooted paths in here
                           #and maybe multiple path with different roots?
                         );

  my $can_convert           = 0;
  my $clean_filtered_format = 0;
  my $filter_format         = $params->{filter_from_format};
  my $all_formats           = $params->{all_formats};
  my $done_formats          = {};

  #Add filter format if it is not in $formats
  if($filter_format){
    #Set sort for safety, but can probably remove this when we refactor this method
    #sort should always be done when doign initial file sorting/merging
    #so if standard or unfiltered bam is present, then we don't need to sort
    $params->{sort}     = 1 if ! defined $params->{sort};    #for safety

    if(!  grep { /^$filter_format$/ } @$formats){
      unshift @$formats, $filter_format;
      $clean_filtered_format = 1;
    }
  }

  #Attempt to get the first or all formats
  foreach my $format(@$formats){
    #warn "Getting $format";
    my $can_filter = 0;

    #Do this before simple file test as it is quicker
    if(grep { /^${format}$/ } keys %$done_formats){ #We have already created this format
      #warn "Skipping $format as it is already done.";
      next;
    }

    #Simple/quick file test first before we do any conversion nonsense
    #This also means we don't have to have any conversion config to get a file which

    #This is being undefd after we filter, so hence, might pick up a pre-exising file!
    if(! defined $filter_format){

       if(my $from_path = check_file($path.'.'.$format, 'gz', $params)){#we have found the required format
          #warn "Found:\t $from_path";
          $done_formats->{$format} = $from_path;
          next;
       }
    }


    ### Validate we can convert ###
    if(exists $conversion_paths{$format}){
      #warn "Found conversion path for $format";
      $can_convert = 1;

      if(defined $filter_format){

        if( ($conversion_paths{$format}->[0] ne $filter_format) &&
            ($all_formats) ){
          throw("Cannot filter $format from $filter_format for path:\n\t$path");
        }
        elsif((scalar(@{$conversion_paths{$format}}) == 1 ) ||
              (! $clean_filtered_format)){
          my $filter_method_name = 'process_'.$filter_format;
          my $filter_method;
          $can_filter = 1;

          #Sanity check we can call this
          if(! ($filter_method = Bio::EnsEMBL::Funcgen::Sequencing::SeqTools->can($filter_method_name))){
            throw("Cannot call $filter_method_name for path:\n\t$path\n".
              'Please add method or correct conversion path config hash');
          }

          #Set outfile here so we don't have to handle unfiltered in process_sam_bam
          #don't add it to $params as this will affected all convert methods
          (my $outpath = $path) =~ s/\.unfiltered$//o;

          #$format key is same as first element

          $done_formats->{$format} = $filter_method->($path.'.'.$filter_format,
                                                      {%$params,
                                                       out_file => $outpath.'.'.$filter_format} );
          #so we don't try and refilter when calling convert_${from_format}_${to_format}

          #delete $params->{filter_from_format};#Is this right?

          undef $filter_format; #Just for safety but not strictly needed
          $path = $outpath;

        }
      }
    }
    elsif($all_formats){
      throw("No conversion path defined for $format. Cannot acquire $format file for path:\n\t$path\n".
        'Please select a supported file format or add config and conversion support methods for $format');
    }
    else{
      warn "No conversion path found for $format\n";
    }

    ### Attempt conversion ###
    if($can_convert){

      if(scalar(@{$conversion_paths{$format}}) != 1){      #already handled process_${format} above
        #Go through the conversion path backwards
        #Start at last but one as we have already checked the last above i.e. the target format
        #or start at 0 if we have $filter_format defined
        my $start_i = (defined $filter_format) ? 0 : ($#{$conversion_paths{$format}} -1);

        for(my $i = $start_i; $i>=0; $i--){
          my $from_format = $conversion_paths{$format}->[$i];

          #Test for file here if we are not filtering! Else we will always go through
          #other formats and potentially redo conversion if we have tidied intermediate files
          if( (! defined $filter_format) &&
              (! grep { /^${from_format}$/ } keys %$done_formats )){
            my $from_path = $path.".${from_format}";

            if($from_path = check_file($from_path, 'gz', $params)){#we have found the required format
              $done_formats->{$from_format} = $from_path;
              #next; #next $x/$to_format as we don't want to force conversion
            }
          }


          #find the first one which has been done, or if none, assume the first is present
          if( (grep { /^${from_format}$/ } keys %$done_formats)  ||
              $i == 0){
            #then convert that to the next, and so on.

            for(my $x = $i; $x < $#{$conversion_paths{$format}}; $x++){

              my $to_format   = $conversion_paths{$format}->[$x+1];
              $from_format    = $conversion_paths{$format}->[$x];
              my $conv_method = 'convert_'.$from_format.'_to_'.$to_format;
              warn "Calling $conv_method";

              #Sanity check we can call this
              if(! ($conv_method = Bio::EnsEMBL::Funcgen::Sequencing::SeqTools->can($conv_method))){
                throw("Cannot call $conv_method for path:\n\t$path\n".
                  'Please add method or correct conversion path config hash');
              }


              $done_formats->{$to_format} = $conv_method->($path.'.'.$from_format, $params);

              #Remove '.unfiltered' from path for subsequent conversion
              if(($i==0) &&
                defined $filter_format){
                $path =~ s/\.unfiltered$//o;
              }
            }

            last; #We have finished with this $conversion_path{format}
          }
        }


        if($clean_filtered_format && ($format eq $filter_format)){
          #filter_format is not our target format, so we need to keep going
          next; #$format
        }
        elsif(! $all_formats){  #else we have found the most preferable, yay!
          last;  #$format
        }
      }
    } #end of can convert
  } #end foreach my $format


  #Now clean $done_formats

  if($clean_filtered_format){
    #actually delete filtered file here?
    delete $done_formats->{$filter_format};
  }

  foreach my $format(keys %$done_formats){
    #doesn't matter about $all_formats here

    if(! grep { /^${format}$/ } @$formats){
      delete $done_formats->{$format};
    }
  }

  return $done_formats;
}



#Is validate_checksum going to have problems as files are gunzipped
#Should validate checksum also handle .gz files i.e. check for entry without .gz, gunzip and validate?
#Maybe all checksums should be done on gunzipped files


sub process_bam{
  my $bam_file = shift;
  my $params   = shift || {};
  assert_ref($params, 'HASH');
  return process_sam_bam($bam_file, {%$params, output_format => 'bam'});
}

sub convert_bam_to_sam{
  my $bam_file = shift;
  my $params   = shift || {};
  assert_ref($params, 'HASH');
  return process_sam_bam($bam_file, {%$params, output_format => 'sam'});
}

#sub process_sam would need to check_file with gz suffix!


# Need to implement optional sort_and_filter_sam here?
# Make this use process_sam_bam by adding support for bed output format
# which would simply pipe sam through awk to reformat.

sub convert_sam_to_bed{
  my $sam_file = shift;
  my $params   = shift || {};
  my $in_file;

  if(! ($in_file = check_file($sam_file, 'gz', $params)) ){
    throw("Cannot find file:\n\t$sam_file(.gz)");
  }

  (my $bed_file = $in_file) =~ s/\.sam(\.gz)*?$/.bed/;
  run_system_cmd($ENV{EFG_SRC}."/scripts/miscellaneous/sam2bed.pl -1_based -files $in_file");

  if(exists $params->{checksum}){
    #Currently just having this exist turns on write & validation
    write_checksum($bed_file, $params);
  }

  return $bed_file;
}



#There are three use cases here:
#1 Cross validation
#2 Returning header opt in case of no samfile header
#3 Ensuring sam has header if ref header not specified

#The last two will happen by default, with 1 happening if both are present
#or the corss validate boolean is passed

#arguably the cross validate boolean could be dropped in favour of testing
#in the calling context, but convenient here


sub validate_sam_header {
  my $sam_bam_file     = shift or throw("Mandatory argument not specified:\tsam/bam file");
  my $header_or_fai    = shift;
  my $xvalidate        = shift;
  my $params           = shift || {};
  assert_ref($params, 'HASH');
  my $is_fai           = ($header_or_fai =~ /\.fai$/o) ?                1 : 0;
  my $debug            = (exists $params->{debug})     ? $params->{debug} : 0;

  if($xvalidate && ! $header_or_fai){
     throw('The cross validation boolean has been passed, but no header/fai argument has been passed');
  }

  validate_path($sam_bam_file);
  #samtools view -t
  #samtools merge -h
  my $header_opt    = ($is_fai) ? ' -t '.$header_or_fai : ' -h '.$header_or_fai;
  # Filter out non-@SQ lines
  my @infile_header = run_backtick_cmd("samtools view -H $sam_bam_file | grep '\@SQ'");

  if($!){
    #$! not $@ here which will be null string
    if(! defined $header_or_fai){
      throw('Could not find an in file header or a reference file header for:'.
        "\n$sam_bam_file\n$header_or_fai\n$!");
    }
    elsif($xvalidate){
      throw('Cross validate boolean has been passed but failed to fetch a header from the reference file:'.
        "\n\t$header_or_fai");
    }
  }
  elsif($xvalidate && ! @infile_header){
    throw('Cross validate boolean has been passed but bam/sam file has no header entries:'.
      "\n\t$sam_bam_file");
  }
  elsif($header_or_fai){
    validate_path($header_or_fai);
    my @ref_header;

    if($is_fai){
      @ref_header = run_backtick_cmd("samtools view -H $header_opt $sam_bam_file");
    }
    else{
      @ref_header = run_backtick_cmd("cat $header_or_fai");
    }


    if(! @ref_header){
       throw("Reference file has no header entries:\n\t$header_or_fai");
    }


    if(scalar(@infile_header) > scalar(@ref_header)){
      throw("Found in file header with more entries that the reference header:\n".
        scalar(@ref_header)."\t$header_or_fai\n".scalar(@infile_header)."\t$sam_bam_file");
    }

    # size difference is fine here, just so long as the file header
    # is a subset of the ref header
    my ($SN, $LN, %ref_header);
    my $hdr_cnt = 0;

    for(@ref_header){
      (undef, $SN, $LN) = split(/\s+/, $_);
      $ref_header{$SN} = $LN;
    }

    foreach my $line(@infile_header){
      (undef, $SN, $LN) = split(/\s+/, $line);
      $hdr_cnt++;

      if(! exists $ref_header{$SN} ){
        throw("$SN exist in file header but not sam header/fai\n".
          $header_or_fai."\n".$sam_bam_file);
      }
      elsif($ref_header{$SN} ne $LN){
        throw("$SN  has mismatched LN entry between file header and sam fasta index\n".
          $ref_header{$SN}."\t".$header_or_fai."\n$LN\t".$sam_bam_file);
      }
    }

    # we don't need the header file as the headers completely match
    # This may result in a feamle header bing replaced with a male header
    # due to it containing the extra @SQ SN:Y line
    # actaully there can be some gender specific top level unassembled contigs too!
    # Meaning any non-gender specific header will need to be a merge, not just the male header
    warn scalar(@ref_header)." lines in reference header:\t$header_or_fai\n".
      $hdr_cnt." lines in file header:\t$sam_bam_file\n" if $debug && ($debug >= 2);
    $header_opt = '' if $hdr_cnt == scalar(@ref_header);
  }

  return $header_opt;
}




#These are split into _init and _run methods to allow initialisation before running the
#peak caller proper. This fits well with the fetch_input/run jobs procedure of the hive.

#Can't do this easily, as we don't know exactly what we want to pass through to the
#PeakCaller. Isn't this already done with run_aligner?
#Actually, we have acces to get_files_by_format here
#so we could optionally defined the signal and control files here given the right input
#Probably to separate out get_alignment_files_by_ResultSet_formats
#If this is possible, then we can add FeatureSet support here and move all the Analysis
#processing here from RunPeaks::fetch_input

sub _init_peak_caller{
  my($peak_pkg, $params, $analysis, $align_prefix, $control_prefix,
    $sam_ref_fai,$debug) = rearrange(
     [qw(peak_module peak_module_params analysis align_prefix control_prefix
     sam_ref_fai debug)], @_);

  my (@params_array, $peak_module);

  if(defined $analysis){ #API mode

    if($peak_pkg){
      throw("Mututally exclusive paramters detected:\n\t-peak_module => $peak_pkg".
        "\n\t-analysis => $analysis(".$analysis->logic_name.')');
    }

    assert_ref($params, 'HASH', 'PeakCaller parameters');
    assert_ref($analysis, 'Bio::EnsEMBL::Analysis');

    #Identify potentially clashing params which woudl be derived from analysis
    #or FeatureSet
    #-program_file
    #-parameters
    #-align_file
    #-control_file

    #TODO identify conflict and throw
    #foreach my $conflict()

    #Let's just overwrite the params with the Analysis derived stuff for now?

    #Things we require in params but can be left to PeakCaller to validate?
    #-out_file_prefix

    #we don't expect any peak_module or analysis params!?

    #but we do require an signal prefix

    #get align file and control file

    if(! defined $align_prefix){
      throw('Must pass an -align_prefix in -analysis mode');
    }

    $peak_module = validate_package_path($analysis->module);
    my $formats = $peak_module->input_formats;
    #my $filter_format = $self->param_silent('bam_filtered') ? undef : 'bam';
    #It is currently unsafe to filter here (control clash), so expect filtered file

    #The problem here is that we are returning a hash of files keys on the format
    #This conversion may cause clashes for fan job which share the same controls
    #(e.g. peak calling jobs if they require formats other than bam)
    #Collections jobs will be pre-processed/converted individually before submitting
    #the slice job.
    #So here we really only need the first available format

    #Restrict to bam (and bed explicitly for CCAT) for now
    my $format = 'bam';

    if($formats->[0] ne 'bam'){

      if($analysis->program eq 'CCAT'){

        warn 'Hardcoding CCAT format to bed until CollectionWriter can be made PeakCaller aware wrt to formats';
        $format = 'bed';
      }
      else{
        throw("It is currently unsafe to use any non-bam format at this point.\n".
          "This is due to the possibility of filtering/format conversion clashes between parallel\n".
          "jobs which share the same control files. Please implement PreprocessAlignments to\n".
          "group jobs by controls and set/handle FILTERING_CONTROL status");
      }
    }

    $formats = [$format];

    #May need to specify checksum param separately
    my $get_files_params = {debug              => $debug,
                            ref_fai            => $sam_ref_fai,  #Just in case we need to convert
                            skip_rmdups        => 1, #This will have been done in merge_bams
                            checksum           => undef,
                            #Specifying undef here turns on file based checksum generation/validation
                           };

    my $align_file = get_files_by_formats($align_prefix,
                                         $formats,
                                         $get_files_params)->{$format};

    if(! defined $align_file){
      throw("Failed to identify file using get_files_by_formats:\n\t".$align_prefix.".${format}");
    }

    my $ctrl_file;

    if($control_prefix){
      $ctrl_file = get_files_by_formats($control_prefix,
                                       $formats,
                                       $get_files_params)->{$format};
    }

    $params->{-debug}          = $debug;
    $params->{-align_file}     = $align_file;
    $params->{-control_file}   = $ctrl_file;
    $params->{-parameters}     = $analysis->parameters;
    $params->{-program_file} ||= $analysis->program_file;
    #All the rest are passed in $params from RunPeaks
    #e.g.
    #-out_file_prefix
    #-out_dir
    #-convert_half_open
    #-is_half_open


    @params_array = %$params; #flatten hash
  }
  else{ #Script mode
    if(! defined $peak_pkg){
      throw('Must provide a -peak_module or -analysis parameter');
    }
    elsif($peak_pkg !~ /::/){
      #Might this be a path?
      warn "Full $peak_pkg namespace was not specified, defaulting to:\t".
        "Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::$peak_pkg\n";
      $peak_pkg = "Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::$peak_pkg";
    }

    assert_ref($params, 'ARRAY', 'PeakCaller params');
    @params_array = @$params;
    $peak_module = validate_package_path($peak_pkg);
  }


  return $peak_module->new(@params_array);
}

sub _run_peak_caller{
  my($pcaller, $max_peaks, $file_types)  = rearrange([qw(peak_caller max_peaks file_types)], @_);
  assert_ref($pcaller, 'Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller');

  #Do this first, so we fail early
  if( (defined $max_peaks) && (! $pcaller->can('filter_max_peaks')) ){
    throw(ref($pcaller).' cannot filter_max_peaks');
  }

  $pcaller->run; #eval in caller if required

  if(defined $max_peaks){

    foreach my $file_type(@$file_types){
      $pcaller->filter_max_peaks($max_peaks, $file_type)
    }
  }

  return;
}

sub run_peak_caller{
   my ($max_peaks, $file_types, $debug, $skip_qc) = rearrange([qw(max_peaks file_types debug skip_qc)], @_);
   my $peak_caller           = _init_peak_caller(@_);

   #Pass these explicitly, so we don't pass all the _init_peak_caller params too
   _run_peak_caller(-peak_caller => $peak_caller,
                    -max_peaks   => $max_peaks,
                    -file_types  => $file_types,
                    -skip_qc     => $skip_qc,
                    -debug       => $debug);



   #Now write/process/load as required
}

#Move peak caller post processing stuff in here too?
#This maybe overkill as it would require Set handling/defnition code


# Validate, count and define IDR threshold
#If you started with ~150 to 300K relaxed pre-IDR peaks for large genomes (human/mouse),
#then threshold of 0.01 or 0.02 generally works well.
#If you started with < 100K pre-IDR peaks for large genomes (human/mouse),
#then threshold of 0.05 is more appropriate.
#This is because the IDR sees a smaller noise component and the IDR scores get weaker.
#This is typically for use with peak callers that are unable to be adjusted to call large number of peaks (eg. PeakSeq or QuEST)
#What exactly are we counting here? Total number peaks across rep, average, or the max between reps?
#This also depends on the prevalence of the mark, it may be that a particular feature type genuinely does not have many genome wide hits

#idr_threshold is being skewed by the rep with the least peaks
#meaning more peaks will be identified from each replicate comparison
#This will however, still be low if we have had to truncate in line with
#the rep with the lowest number of pre-IDR peaks

sub pre_process_IDR{
  my $out_dir      = shift or throw('Must provide and out_dir argument');
  my $pre_idr_beds = shift;
  my $batch_name   = shift or throw('Must provide a batch_name argument');
  my $max_peaks    = shift or throw('Must provide a max_peaks argument');
  throw("out_dir does not exist:\t$out_dir") if ! -d $out_dir;
  assert_ref($pre_idr_beds, 'ARRAY');

  if(scalar(@$pre_idr_beds) < 2){
    throw("Cannot run IDR with less than 2 replicates:\n\t".$pre_idr_beds->[0]);
  }

  my ($mt_100k, $lt_100k, $log_txt);

  foreach my $bed_file(@$pre_idr_beds){

    if(! -f $bed_file){
      throw("Could not find pre-IDR bed file:\t$bed_file\n".
        'Need to dump bed from DB directly to required format!');
    }

    # Ignore comments and header
    my $cmd       = "grep -vE '(#|(^Region[[:space:]]+Start))' $bed_file | wc -l | awk '{print \$1}'";
    my $num_peaks = run_backtick_cmd($cmd);

    # $num_peaks can sometimes be 0. Throw here as it will just fail in the RunIDR step
    # We likely need to remove a replicate? 

    if($num_peaks == 0){
      throw("Found bed file with 0 peaks. This will cause IDR to fail, please remove or fix:\n\t$bed_file");
    }

    if($num_peaks > 100000){
      $mt_100k = 1;

      if($num_peaks < 150000){
        warn 'Pre-IDR peaks counts fall in threshold grey zone(100k-150k), defaulting to 0.01'.
          " but maybe consider 0.02 for:\n\t$bed_file\n";
      }
    }
    else{
      $lt_100k = 1;
    }

    if($num_peaks < $max_peaks){
      # We take the lowest number of peaks, as we need comparable numbers of peaks across all inputs
      $max_peaks = $num_peaks;
    }

    $log_txt .= $bed_file."\t".$num_peaks."\n";
  }

  # Note this does not yet support MACS yet, should prbably just ignore it as we filter to 100000
  my $idr_threshold   = ($max_peaks < 100000) ? 0.05 : 0.01;
  # Could alternatively pass all thresholds back to the caller
  my $x_thresh_adjust = 0;

  if($lt_100k && $mt_100k){
    #$self->throw_no_retry("Identified different optimal thresholds due to pre-IDR peak counts spanning the 100k limit:\n".
    #  $log_txt);
    warn 'Identified different optimal thresholds due to pre-IDR peak counts spanning the 100k limit:'.
      "\n$log_txt\nDefaulting to lower threshold:\t$idr_threshold\n";
    $x_thresh_adjust = 1;
  }

  # TODO We need some mechanism to restart this job, to force the threshold, or by dropping 1/some of the replicates.

  my $cmd = "echo -e \"Pre-IDR File\tNum Peaks\n$log_txt\nIDR Threshold = $idr_threshold\"".
    "> ${out_dir}/${batch_name}-idr-stats.txt";
  run_system_cmd($cmd);

  # TODO parallelise the filtering and reformating to speed things up, so let's semphaore than to a simple CMD job.
  # can we even do this as we already have a semaphore on the RunIDR and
  # maybe with a job factory? I think this is not possible without another analysis
  # but we could use a dummy? which then submit the RunIDR and semaphores the PostProcessIDR
  # Just do here for now
  my @np_bed_files;

  foreach my $i(0..$#{$pre_idr_beds}){
    my $bed_file     = $pre_idr_beds->[$i];
    (my $np_bed_file = $bed_file) =~ s/\.txt$/.np_idr.txt/;
    # Never re-use np_idr output file in case it is truncated due to job failure/exit. 
    # Is there a danger of an np_idr file usage clash (pseudo/self consistency IDR?

    # SWEmbl output header::
    # Input  GSE30558_GSM758361_PrEC.sorted.bam
    # Reference      GSE30558_GSM758360_LNCaP.sorted.bam
    # Sequence length        0
    # Fragment length        0
    # Background     0.000000
    # Position Background    0.036383
    # Long Background        0.181917
    # Threshold      5.000000
    # Minimum count above bg 15
    # Penalty increase       70
    # Quality cutoff 0.000000
    # Result cutoff  0.000000
    # Penalty factor 0.552834

    # and fields:
    # Region        - Part of the genome build e.g. chromosome
    # Start pos.    - Base in region where peak starts
    # End pos.      - Base in region where peak ends
    # Count         - Number of reads in experimental sample in peak
    # Length        - Length of peak (distance between start and end pos.)
    # Unique pos.   - Number of unique bases within peak at which reads begin
    # Score         - The SWEMBL score, which is basically the count of filtered thresholded reads in the peak, minus the penalties (gap distances and reference sample reads).
    # Ref. count    - Number of reads in reference sample in peak
    # Max. Coverage - Depth of reads at summit
    # Summit        - Median position of highest read coverage

    # signalValue field was being set to (Count - Ref. Count)/min
    # Min is not really required here for ranking and simply add the header requirement
  
    # TODO Get header skipping regex from PeakCaller. Currently hardcoded for SWEmbl
    # TODO Handle pipes with perl pipe or IPC::open/run

    # Strip out the header and set signal.value to score, sort on score, filter based on $max_peaks, resort based on position
    $cmd = 'awk \'BEGIN {OFS="\t"} { if($0 !~ /^(#|(Region[[:space:]]+Start))/) {print $1,$2,$3,".",$7,".",$7,-1,-1,int($9-$1)} }\' '.
      "$bed_file | sort -k 7nr,7nr | head -n $max_peaks | sort -k 1,2n > ".$np_bed_file;
    run_system_cmd($cmd);

    # Sanity check we have the file with the correct number of lines
    $cmd = "wc -l $np_bed_file | awk '{print \$1}'";
    my $filtered_peaks = run_backtick_cmd($cmd);

    if($max_peaks != $filtered_peaks){
      throw("Expected $max_peaks in filtered pre-IDR bed file, but found $filtered_peaks:\n\t".$np_bed_file);
    }

    # Need to check this is != 0?


    # TODO check the feature_set_stat or states
    # such that we know the peak calling jobs has finished and passed QC!
    # Do this for each before we submit IDR jobs, as we may have to drop some reps
    push @np_bed_files, $np_bed_file;
  }

  return (\@np_bed_files, $idr_threshold, $x_thresh_adjust);
}



sub run_IDR{
  my ($out_dir, $output_prefix, $threshold, $bed_files, $batch_name, $bam_files) =
    rearrange([ qw(out_dir output_prefix threshold bed_files batch_name bam_files) ], @_);
  #add max_npairs, this may have to be a denominator of scalar(@$idr_comparison_files)?
  defined $out_dir       or throw('Must provide an out_dir argument');
  defined $output_prefix or throw('Must provide an output_prefix argument');
  defined $threshold     or throw('Must provide an IDR threshold to count peaks');
  assert_ref($bed_files, 'ARRAY', 'bed_files');
  defined $batch_name    or throw('Must provide a batch name to log counts to idr-stats file');
  assert_ref($bam_files, 'ARRAY', 'bam_files') if defined $bam_files;

  #use a default for this? Currently defined in Preprocess_IDR

  #Check we have different files
  #This is not sensible as we will need to run self consistency IDR?

  if($bed_files->[0] eq $bed_files->[1]){
    throw("Pre-IDR peak files are identical:\t".join(' ', @$bed_files));
  }

  #Check we have 2 reps
  if(scalar (@$bed_files) != 2){
    throw("run_IDR expect 2 replicate bed files:\t".join(' ', @$bed_files));
  }

  #Check output dir exists
  if(! -d $out_dir){
    throw("Output directory does not exist:\t$out_dir");
  }

  #IDR analysis
  #TODO install idrCode in /software/ensembl/funcgen and add this an analysis?
  #my $idr_name = $self->idr_name;
  my $script_path = which_path('batch-consistency-analysis.r');
  my $cmd = "Rscript $script_path ".join(' ', @{$bed_files}).
    " -1 ${out_dir}/${output_prefix} 0 F signal.value";
  #signal.value is ranking measure here i.e. SWEmbl score
  run_system_cmd($cmd);

  #Do this here rather than in post_process so we parallelise the awk.
  $cmd = "awk '\$11 <= ".$threshold.
    " {print \$0}' ${out_dir}/${output_prefix}-overlapped-peaks.txt | wc -l";
  my $num_peaks = run_backtick_cmd($cmd);

  #Now, do we write this as an accu entry in the hive DB, or do we want it
  #in the tracking DB?
  #Probably the later, such that we can drop/add reps to an IDR set after we have dropped the hive DB.
  #There is currently no logic place to put this in the tracking DB!
  #We would have to add a result_set_idr_stats table to handle the multiplicity
  #This is probably a good place to store the other IDR stats too?
  #Just write to file for now until we know if/what we want in the table.
  #Do we need to be concerned if thresholds differ between combinations?

  my $unaltered_num_peaks = $num_peaks;

  #Temporary solution to handle differeing pre-IDR peak counts, which cross threshold boundaries
  #Do we have access to the ResultSet here? No!

  if($bam_files){
    my @align_counts;

    #- the number of peaks is strongly correlated (from the couple of examples we looked at)
    # to the number of reads.
    #- the IDR gives an estimate on the number of realistic peaks at the intersection of two files.
    #- obviously this intersection is limited by the smaller file

    #We have two files, one with N reads and the other with n < N reads.
    #From their intersection we call p peaks that look good. However, if n were greater,
    #p would also be greater. So how big should n be? At first blush, n should be each to the mean (N+n)/2.
    #p*((N+n)/2)/n

    foreach my $bam(@$bam_files){
      #Unfiltered counts are already available in the flagstat output
      #in the alignment report

      #We always want to counts of the input(filtered) bam file
      #so let's recount here for safety, even though we aren't filtering at present
      my $cmd = "samtools view $bam | wc -l";
      push @align_counts, run_backtick_cmd($cmd);
    }

    my ($smalln, $bigN) = sort { $a <=> $b } @align_counts;
    $num_peaks *= ( ($bigN + $smalln) / 2 ) / $smalln;
  }


  # Warning: Parallelised appending to file!
  # This will also cause duplicate lines if the RunIDR jobs are rerun!
  $cmd = "echo -e \"IDR Comparison\tIDR Peaks\n$output_prefix\t$num_peaks($unaltered_num_peaks)\"".
    " >> ${out_dir}/${batch_name}-idr-stats.txt";
  run_system_cmd($cmd);

  return $num_peaks;
}

# TODO
# 1 Add more IDR based QC here
# 2 Rscript will currently submit and interactive job to farm if launched on a head node
#   else will run locally if already on a farm node

sub post_process_IDR{
  my $out_dir       = shift or throw('Must provide an out_dir argument');
  my $output_prefix = shift or throw('Must provide an output_prefix argument');
  my $idr_peaks     = shift;
  my $params        = shift || {};
  assert_ref($idr_peaks, 'ARRAY', 'idr_peaks');

  #add max_npairs, this may have to be a denominator of scalar(@$idr_comparison_files)?
  my ($idr_files, $npairs, $debug) =
    rearrange([ qw(idr_files npairs debug) ], %$params);

  my $max_peaks = int((sort {$a <=> $b} @{$idr_peaks})[-1]);#Take the highest!

  my $cmd = "echo -e \"IDR Max Peaks = $max_peaks\" >> ${out_dir}/${output_prefix}-idr-stats.txt";
  run_system_cmd($cmd);

  if($output_prefix && $idr_files){
    assert_ref($idr_files, 'ARRAY');#or just handle scalar arg here?
    $npairs ||= 1;

    if($npairs > scalar(@$idr_files)){
      throw("batch-consistency-plot.r will not handle npairs($npairs)".
        ' greater than the number of idr files specified('.scalar(@$idr_files).')');
    }

    my $script_path = which_path('batch-consistency-plot.r');
    $cmd = "Rscript $script_path 1 $out_dir/$output_prefix ".join(' ', map { $out_dir.'/'.$_ } @{$idr_files});
    run_system_cmd($cmd);

    warn "IDR plots are available here:\n\t$out_dir/${output_prefix}-plot.ps\n" if $debug;
  }
  elsif($output_prefix || $idr_files){
    throw('To run batch-consistency-plot.r both the output_prefix and idr_files parameters must be set');
  }

  warn "Final IDR Max peaks threshold:\t$max_peaks\n".
    "For more details, see ${out_dir}/${output_prefix}-idr-stats.txt\n" if $debug;

  return $max_peaks;
}

# Iterates through all files.txt found in the typical goldenPath
# directory structure
# !!! Only stores fastq !!!

sub create_and_populate_files_txt {
  my ($cfg, $helper) = @_;

  # reduced to classes (potentially) present in $table
  my $class = "'Histone','RNA','Polymerase','Transcription Factor','Open Chromatin'";
  my $table = $cfg->{tables}->{registration};

  my $table_id = $table.'_id';

  my $sql_table = "
          CREATE TABLE `$table` (
            `$table_id`         INT(10) unsigned      NOT NULL  auto_increment,
            `name`              VARCHAR(100)          NOT NULL,
            `alternate_name`    VARCHAR(100)          DEFAULT NULL,
            `antibody`          VARCHAR(64)           DEFAULT NULL,
            `assembly`          enum('hg18', 'hg19')  NOT NULL,
            `cell`              VARCHAR(64)           NOT NULL,
            `compression`       VARCHAR(10)           DEFAULT NULL,
            `control`           VARCHAR(64)           DEFAULT NULL,
            `controlId`         VARCHAR(64)           DEFAULT NULL,
            `dataType`          VARCHAR(50)           DEFAULT NULL,
            `dateUnrestricted`  DATE                  DEFAULT NULL,
            `filename`          VARCHAR(100)          NOT NULL,
            `lab`               VARCHAR(100)          NOT NULL,
            `md5sum`            CHAR(32)              DEFAULT NULL,
            `objStatus`         VARCHAR(255)          DEFAULT NULL,
            `path`              VARCHAR(255)          DEFAULT NULL,
            `replicate`         INTEGER(2)            DEFAULT NULL,
            `setType`           VARCHAR(25)           DEFAULT NULL,
            `size`              VARCHAR(5)            DEFAULT NULL,
            `treatment`         VARCHAR(50)           DEFAULT NULL,
            `type`              VARCHAR(20)           NOT     NULL,
            `cell_type`         VARCHAR(120)          DEFAULT NULL,
            `class`             ENUM($class)          DEFAULT NULL,
            `ens_lab`           VARCHAR(100)          NOT NULL,
            `feature_type`      VARCHAR(40)           DEFAULT NULL,
            `logic_name`        VARCHAR(100)          DEFAULT NULL,

            PRIMARY KEY  (`$table_id`),
            UNIQUE name_assembly_idx (`name`, `assembly`)
            ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;
  ";
  # DBI->trace(2);

  $helper->execute_update(-SQL => "DROP TABLE IF EXISTS `$table`");
  $helper->execute_update(-SQL => $sql_table);


  #Make table name variable, use $cfg

  my $sql_select_by_md5sum = "
    SELECT
      name
    FROM
      $table
    WHERE
      md5sum = ?
  ";

  my $sql_insert_table = "
  INSERT INTO
    $table (
      name,
      alternate_name,
      antibody,
      assembly,
      cell,
      compression,
      control,
      controlId,
      dataType,
      dateUnrestricted,
      filename,
      lab,
      md5sum,
      objStatus,
      path,
      replicate,
      setType,
      size,
      treatment,
      type,
      cell_type,
      ens_lab,
      feature_type
    )
  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);
  ";
  my @assemblies = qw(hg18 hg19);
  for my $assembly(@assemblies) {

    my $base_dir = File::Spec->catdir('goldenPath', $assembly, 'encodeDCC');
    my $dir_dcc  = File::Spec->catdir($cfg->{directories}->{files_txt}, $base_dir);

    opendir (my $labs, $dir_dcc ) or die "Error  opening dir $dir_dcc";
    while( my $lab = readdir($labs)){
      next if($lab !~ /^wgEncode/);
      # DELETE next unless($lab eq 'wgEncodeOpenChromChip');
      my $path_table = File::Spec->catfile($dir_dcc, $lab, 'files.txt');
      if(-f $path_table){
        # !!! only fastq coming back !!!!
        my $files_txt = _read_files_txt($path_table);
        # DELETE   say dump_data($files_txt->{wgEncodeOpenChromChipHelas3CtcfRawDataRep2},1,1);

        foreach my $name (sort keys %{$files_txt}){
          my $record = $files_txt->{$name};
          # DELETE next unless ($name eq 'wgEncodeOpenChromChipHelas3CtcfRawDataRep2');
          my $path = File::Spec->catfile($base_dir,$lab);
          # DELETE   say 'Record: '. dump_data($record,1,1);
          my $md5sum = $record->{md5sum};
         # DELETE  say 'Intial value: '.$record->{md5sum};
          my $alt_name = undef;

          if(defined $md5sum){
              # DBI->trace(2);
              $alt_name =
              $helper->execute_single_result(
                -SQL      => $sql_select_by_md5sum,
                -PARAMS   => [$md5sum],
                -NO_ERROR => 1
                );
            }
         # DELETE  say 'After lookup: '.$record->{md5sum};

          # use the same values as ENCODE. Ensembl specific changes are applied later
          my $cell_type     = $record->{cell};
          my $feature_type  = $record->{antibody};
          my $ens_lab       = $record->{lab};
            # say dump_data($$table->{$name},1,1);
            # say $cell_type,
            # say $feature_type;
            # say $ens_lab;



            $helper->execute_update(
              -SQL    => $sql_insert_table,
              -PARAMS => [
              [$name,                       SQL_VARCHAR],
              [$alt_name,                   SQL_VARCHAR],
              [$record->{antibody},         SQL_VARCHAR],
              [$assembly,                   SQL_VARCHAR],
              [$record->{cell},             SQL_VARCHAR],
              [$record->{compression},      SQL_VARCHAR],
              [$record->{control},          SQL_VARCHAR],
              [$record->{controlId},        SQL_VARCHAR],
              [$record->{dataType},         SQL_VARCHAR],
              [$record->{dateUnrestricted}, SQL_VARCHAR],
              [$record->{filename},         SQL_VARCHAR],
              [$record->{lab},              SQL_VARCHAR],
              [$record->{md5sum},           SQL_VARCHAR],
              [$record->{objStatus},        SQL_VARCHAR],
              [$path,                       SQL_VARCHAR],
              [$record->{replicate},        SQL_VARCHAR],
              [$record->{setType},          SQL_VARCHAR],
              [$record->{size},             SQL_VARCHAR],
              [$record->{treatment},        SQL_VARCHAR],
              [$record->{type},             SQL_VARCHAR],
              [$cell_type,                  SQL_VARCHAR],
              [$ens_lab,                    SQL_VARCHAR],
              [$feature_type,               SQL_VARCHAR],
              ]);
  # die;
        }
      }
      else{
        say "$path_table not available";
      }
    }
    closedir($labs);
  }
}

=head2 _read_files_txt

  Argument 1  :
  Returntype  :
  Exceptions  :
  Description : Reads $table and returns a HASHREF with filename, fieldname and values, e.g.:
                $h->{wgEncodeHaibTfbsHepg2P300V0416101RawDataRep2}->{md5sum} = e4d9de1900e5cd950196d4d46f321816
                Not every files.txt in hg18 has a md5sum field, the data was provided in a
                In some cases files.txt in hg18 don't have md5sum field. In those cases the values is taken from
                the md5sum.txt file.

=cut

sub _read_files_txt {
  my ($file) = @_;
say "FILE: $file";
  my $md5sum_file = $file;
     $md5sum_file =~ s/files\.txt/md5sum\.txt/;

  my $md5sums = _read_md5sum_txt($md5sum_file);

  my $all_records = {};
  open(my $fh,'<',$file) or die "Can not open/access '$file'\n$!";
  # line: wgEncodeSydhTfbsA549Bhlhe40IggrabAlnRep1.fastq.gz  project=wgEncode; grant=Snyder; lab=Stanford; (...)
  while(my $line = <$fh>){
    chomp($line);

    # wgEncodeHudsonalphaChipSeqA549ControlPcr2xDex100nmAlnRep1.fastq.gz
    # wgEncodeHudsonalphaChipSeqA549ControlPcr2xDex100nmAlnRep1.bam

    # $1: wgEncodeBroadChipSeqAlignmentsRep1Gm12878ControlV2.tagAlign.gz
    # $2: wgEncodeBroadChipSeqAlignmentsRep1Gm12878ControlV2
    # $3: tagAlign
    # $4: gz

    $line =~ /^((\S+?)\.(\w+)\.?(\w+))?\s+/;
    # only keeping fastq's
    throw($line) if(!$3);
    next unless($3 eq 'fastq');
    my $filename = $1;
    # name is unique in input_subset
    my $name     = $2;
    #$all_records->{$name}->{format} = $2;
    $all_records->{$name}->{filename}    = $1;
    $all_records->{$name}->{compression} = $3;

    # Cut filename from line
    $line =~ s/^\S+?\s+//o;

    my @fields = split(/;/,$line);
    for my $field (@fields){
      $field =~ s/^\s//;
      $field =~ s/\s$//;
      my @values = split(/=/, $field);
      # be suspicious!
      if(exists $all_records->{$name}->{$values[0]}){
        throw "$name $values[0] duplicate";
      }
      $all_records->{$name}->{$values[0]} = $values[1];
    }
    # DELETE say "Before if: ". dump_data($all_records->{$name},1,1) if($name eq 'wgEncodeOpenChromChipHelas3CtcfRawDataRep2');

    if(! exists $all_records->{$name}->{md5sum}){

      # DELETE  say "Before if2: ";say dump_data($all_records->{$name},1,1);
      if(!exists $md5sums->{$name}){
        warn "no md5sum for $name";
      }
      else{
        $all_records->{$name}->{md5sum} = $md5sums->{$name};
      }
      # DELETE  say "After if2: ";say dump_data($all_records->{$name},1,1);


    }
  }
  close($fh);
  return $all_records;
}


=head2 modify_files_txt_for_regulation

  Arg 1  : HASH - Configuration
  Arg 2  : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Returntype  : None
  Exceptions  : Throws if non-optional arguments are missing
  Description : Downloads all files.txt and md5sum.txt from standard ENCODE directory structure.

=cut

sub modify_files_txt_for_regulation {
  my ($cfg, $helper) = @_;

  warn "\n\n++++++ These modifications should be regularly reviewed. Especially assigning class HISTONE ++++++\n";
  my $table = $cfg->{tables}->{registration};


  if($table ne 'files_txt'){
   warn "These modifications are specific to ENCODE $table";
  }

  say "\n\n+++++++++++++++++++ Modifications to $table: FeatureType +++++++++++++++++++\n";

  # Antibodies "CTCF_(SC-15914)" and "CTCF_(SC-5916)" to FeatureType "CTCF"
  # no discinct as we execute_into_hash
  # Be careful, this will also do wrong shortenings like ZNF-MIZD-CP1_(ab65767) to ZNF-MIZD-CP1
  # These are addressed individually later
  my $sql_select_antibody = "
    SELECT
      antibody,
      1
    FROM
      $table
    WHERE
      antibody LIKE '%_(%)'
  ";

  my $tmp = $helper->execute_into_hash(-SQL => $sql_select_antibody);

  my $ab_to_ft = {};
  foreach my $antibody (sort keys %{$tmp}){
    # CTCF_(SC-15914) or  Pol2(phosphoS2)
    $antibody =~ /^(.*?)_?\(/;
    $ab_to_ft->{$antibody} = $1;
  }

  my $sql_update_feature_type = "
    UPDATE
      $table
    SET
      feature_type = ?
    WHERE
      antibody = ?
  ";

  foreach my $antibody (sort keys %{$ab_to_ft}){
    my $feature_type = $ab_to_ft->{$antibody};
    say "UPDATE $table SET feature_type = $feature_type\tWHERE antibody = $antibody";
    $helper->execute_update(
      -SQL    =>  $sql_update_feature_type,
      -PARAMS => [$feature_type, $antibody],
      );
  }

  say "\n\n+++++++++++++++++++ Modifications to $table: Miscellaneous +++++++++++++++++++\n";

  my @sqls;

  # replicate
  push(@sqls, "UPDATE $table SET replicate = 1 WHERE replicate IS NULL");

  # treatment
  push(@sqls, "UPDATE $table SET treatment = 'None' WHERE treatment IS NULL");

  # setType
  push(@sqls, "UPDATE files_txt SET setType = 'input' WHERE feature_type = 'Input'   AND setType = 'exp'");
  push(@sqls, "UPDATE files_txt SET setType = 'input' WHERE antibody     = 'Control' AND setType = 'exp'");

  # antibody
  push(@sqls, "UPDATE $table SET antibody = 'DNase' WHERE name LIKE '%dnase%'");

  # cell_type
  push(@sqls, "UPDATE $table SET cell_type = 'Monocytes-CD14+' WHERE cell = 'Monocytes-CD14+_RO01746'");
  push(@sqls, "UPDATE $table SET cell_type = 'DND-41'          WHERE cell = 'Dnd41'");
  push(@sqls, "UPDATE $table SET cell_type = 'H1ESC'           WHERE cell = 'H1-hESC'");

  #ens_lab
  push(@sqls, "UPDATE $table SET ens_lab = 'UTA' WHERE lab = 'UT-A'");

  #feature_type
  push(@sqls, "UPDATE $table SET feature_type = 'H2AF'  WHERE antibody = 'H2A.Z'");
  push(@sqls, "UPDATE $table SET feature_type = 'DNase1' WHERE antibody = 'DNase'");
  push(@sqls, "UPDATE $table SET feature_type = 'CTCF'   WHERE antibody like 'CTCF_%'");
  push(@sqls, "UPDATE $table SET feature_type = 'PolII'  WHERE antibody like 'Pol2%'");
  push(@sqls, "UPDATE $table SET feature_type = 'PolIII' WHERE antibody like 'Pol3%'");

  # http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000087510;
  push(@sqls, "UPDATE $table SET feature_type = 'TFAP2C' WHERE feature_type = 'AP-2gamma' ");
  # http://moma.ki.au.dk/genome-mirror/cgi-bin/hgEncodeVocab?ra=encode%2Fcv.ra&target=%22CHD4_Mi2%22
  push(@sqls, "UPDATE $table SET feature_type = 'CHD4'   WHERE antibody = 'CHD4_Mi2' ");
  # http://www.factorbook.org/mediawiki/index.php/ERalpha_a
  push(@sqls, "UPDATE $table SET feature_type = 'ESR1'   WHERE antibody = 'ERalpha_a' ");
  push(@sqls, "UPDATE $table SET feature_type = 'EGR1'   WHERE antibody = 'Egr-1' ");
  push(@sqls, "UPDATE $table SET feature_type = 'GATA2'  WHERE antibody = 'GATA-2' ");
  # http://epigenome.cbrc.jp/cgi-bin/hgEncodeVocab?ra=encode%2Fcv.ra&target=%22HA-E2F1%22
  push(@sqls, "UPDATE $table SET feature_type = 'E2F1'   WHERE antibody = 'HA-E2F1' ");
  # http://moma.ki.au.dk/genome-mirror/cgi-bin/hgEncodeVocab?ra=encode%2Fcv.ra&target=%22NCoR%22
  push(@sqls, "UPDATE $table SET feature_type = 'NCOR1'  WHERE antibody = 'NCoR' ");
  # http://www.noncode.org/cgi-bin/hgEncodeVocab?ra=encode%2Fcv.ra&target=%22p300%22
  push(@sqls, "UPDATE $table SET feature_type = 'EP300'  WHERE antibody = 'P300_KAT3B' ");
  push(@sqls, "UPDATE $table SET feature_type = 'PAX5'   WHERE antibody like 'PAX5-%' ");
  push(@sqls, "UPDATE $table SET feature_type = 'SIN3A'  WHERE antibody = 'Sin3Ak-20' ");
  # http://epigenome.cbrc.jp/cgi-bin/hgEncodeVocab?ra=encode%2Fcv.ra&target=%22TCF7L2%22
  push(@sqls, "UPDATE $table SET feature_type = 'TCF7L2' WHERE antibody = 'TCF7L2_C9B9_(2565)' ");
  push(@sqls, "UPDATE $table SET feature_type = 'USF1'   WHERE antibody = 'USF-1' ");
  # http://www.noncode.org/cgi-bin/hgEncodeVocab?ra=encode%2Fcv.ra&target=%22ZNF-MIZD-CP1_(ab65767)%22
  push(@sqls, "UPDATE $table SET feature_type = 'ZMIZ1'  WHERE antibody = 'ZNF-MIZD-CP1_(ab65767)' ");

  # See also: http://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=control
  push(@sqls, "UPDATE $table SET feature_type = 'WCE'    WHERE setType = 'Input'");

  # logic_name
  push(@sqls, "UPDATE $table SET logic_name = 'ChIP-Seq'  WHERE dataType = 'ChipSeq'");
  push(@sqls, "UPDATE $table SET logic_name = 'DNase-Seq' WHERE dataType = 'DnaseSeq'");
  push(@sqls, "UPDATE $table SET logic_name = 'FAIRE'     WHERE dataType = 'FaireSeq'");

  # class
  # This is true for the ENCODE 2011 data freeze
  push(@sqls, "UPDATE $table SET class = 'Histone'                WHERE dataType = 'ChipSeq' AND feature_type LIKE 'H%K%'");
  push(@sqls, "UPDATE $table SET class = 'Polymerase'             WHERE dataType = 'ChipSeq' AND feature_type IN ('PolII', 'PolIII')  ");
  push(@sqls, "UPDATE $table SET class = 'Open Chromatin'         WHERE dataType = 'DnaseSeq'");
  push(@sqls, "UPDATE $table SET class = 'Transcription Factor'   WHERE dataType = 'ChipSeq' AND feature_type != 'WCE' AND class is NULL");
  push(@sqls, "UPDATE $table SET class = 'Transcription Factor'   WHERE antibody = 'H2A.Z'");


  for my $sql(@sqls){
    say $sql;
    $helper->execute_update(-SQL => $sql);
  }
}

=head2 load_experiments

  Arg 1  : HASH - Configuration
  Arg 2  : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Arg 3  : HASH - Constraints
  Arg 4  : [Optional] ARRAY of Bio::EnsEMBL::Funcgen::Experiment
  Arg 5  : [Optional] ARRAY of Bio::EnsEMBL::Funcgen::CellType
  Arg 6  : [Optional] ARRAY of Bio::EnsEMBL::Funcgen::FeatureType

  Returntype  : None
  Exceptions  : Throws if non-optional arguments are missing
  Description : Downloads all files.txt and md5sum.txt from standard ENCODE directory structure.

=cut


sub load_experiments_into_tracking_db {
  my ($cfg, $db, $constraints, $exp_data, $cell_type_data, $feature_type_data) = @_;

  assert_ref($cfg, 'HASH');
  assert_ref($db, 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  assert_ref($constraints, 'HASH');

  my $helper  = Bio::EnsEMBL::Utils::SqlHelper->new( -DB_CONNECTION => $db->dbc );

  my $exp_a = $db->get_ExperimentAdaptor;
  my $eg_a  = $db->get_ExperimentalGroupAdaptor;

  my $ct_a = $db->get_CellTypeAdaptor;
  my $ft_a = $db->get_FeatureTypeAdaptor;

  my $anal_adaptor = $db->get_AnalysisAdaptor;
  my $iss_a = $db->get_InputSubsetAdaptor;

  my $tr_a  = Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor->new(
        -user       => $db->dbc->user,
        -pass       => $db->dbc->pass,
        -host       => $db->dbc->host,
        -port       => $db->dbc->port,
        -dbname     => $db->dbc->dbname,
        -dnadb_name => $db->dnadb_name,
  );


  my $sql = "
    SELECT
      cell_type,
      dateUnrestricted,
      ens_lab,
      feature_type,
      filename,
      logic_name,
      md5sum,
      name,
      objStatus,
      path,
      replicate
    FROM
      files_txt
  ";
  # Adding constraints
  if(defined $constraints){
    $sql .= 'WHERE ';

   for my $cst(@$constraints){
     my $table = shift @$cst;
     my $values = join(', ', map { qq/"$_"/ } @$cst);
     $sql .= "$table IN ($values)";
     $sql .= " AND\n " if($cst != $constraints->[-1]);
   }
  }
  # say $sql;
  # DBI->trace(2);
  my $files = $helper->execute(
    -SQL      => $sql,
    -CALLBACK => sub {
      my @row = @{shift @_};
      return {
        cell_type         => $row[0],
        dateUnrestricted  => $row[1],
        ens_lab           => $row[2],
        feature_type      => $row[3],
        filename          => $row[4],
        logic_name        => $row[5],
        md5sum            => $row[6],
        name              => $row[7],
        objStatus         => $row[8],
        path              => $row[9],
        replicate         => $row[10],
      }
    }
  );
      # say dump_data($files,1,1);
  # die;
  foreach my $file (@$files){
    my $ft_name  = $file->{feature_type};
    my $ct_name  = $file->{cell_type};


    # type determines the style of the experiment name
    my $type = $cfg->{general}->{type};

    my $exp_name;
    if($type eq 'ENCODE'){
      $exp_name = $ct_name .'_' . $ft_name . '_' . $type . '_' . $file->{ens_lab};
    }
    else{
      throw "'$type' not implemted.";
    }

    # CellType
    my $ct  = $ct_a ->fetch_by_name($ct_name);
    if(!$ct){
      $ct =_store_cell_feature_type ($db, $ct_name, $cell_type_data);
    }

    # FeatureType
    my $ft  = $ft_a ->fetch_by_name($ft_name);
    if(!$ft){
      $ft = _store_cell_feature_type ($db, $ft_name, $feature_type_data);
    }

    # Implement store method
    my $anal = $anal_adaptor->fetch_by_logic_name($file->{logic_name});
    if(not $anal){
      warn "Analysis $file->{logic_name} not in DB. Skipping...";
      next;
    };

    # Risky me thinks
    my $control = 0;
    $control = 1 if($ft_name eq "WCE");


    my $exp = $exp_a->fetch_by_name($exp_name);
    my $iss = $iss_a->fetch_by_name($file->{name}, $exp);

    # Check if InputSubset is already linked to a different Experiment
    if(!$exp and $iss){
      if($exp_name ne $iss->experiment->name){
        throw($iss->name . ' is linked to ' . $iss->experiment->name . ' not ' . $exp_name );
      }
    }

    if(! $exp){
      my $eg = $eg_a->fetch_by_name($exp_data->{experimental_group});

      my $exp_new = Bio::EnsEMBL::Funcgen::Experiment->new
                       (
                        -cell_type           => $ct,
                        -experimental_group  => $eg,
                        -feature_type        => $ft,
                        -date                => DateTime::Format::MySQL->format_datetime(DateTime->now),
                        -description         => $exp_data->{description},
                        -name                => $exp_name,
                       );
       ($exp) = @{$exp_a->store($exp_new)};
       say "Added Experiment " . $exp->name . ' [dbID: ' . $exp->dbID .']';
    }

    if(!$iss){
      my $iss_new = Bio::EnsEMBL::Funcgen::InputSubset->new
                     (
                      -cell_type     => $ct,
                      -experiment    => $exp,
                      -feature_type  => $ft,
                      -analysis      => $anal,
                      -is_control    => $control,
                      -name          => $file->{name},
                      -replicate     => $file->{replicate},
                     );
      ($iss) = @{$iss_a->store($iss_new)};
      say "Added InputSubset " . $iss->name . ' [dbID: ' .$iss->dbID .']';

      my $web_url = File::Spec->catfile($cfg->{urls}->{base}, $file->{path}, $file->{filename});
      # catfile replaces // with /
      $web_url =~ s!:/!://!;

      if($file->{objStatus}){
        if($file->{objStatus} !~ /^[revoked|replaced]/){
          throw('Status: '. $file->{objStatus}. ' not implemted');
        }
      }

      my $tr_info->{info} = {
        availability_date => 1,
        download_url      => $web_url,
        download_date     => undef,
        local_url         => undef,
        md5sum            => $file->{md5sum},
        notes             => $file->{objStatus},
      };
      my $out = $tr_a->store_tracking_info($iss, $tr_info);
      say "Added TrackingInfo for " . $iss->name . ' [dbID: ' .$iss->dbID .']';

    }
  }
}



=head2 _download_all_files_txt

  Argument 1  : HASH - configuration
  Returntype  : None
  Exceptions  : Missing config, inaccessible remote server or local directories
  Description : Downloads all files.txt and md5sum.txt from standard ENCODE directory structure.

=cut

sub download_all_files_txt {
  my ($cfg) = @_;

  my $server = $cfg->{urls}->{base};
  my $base_data_dir = $cfg->{directories}->{data};


  my $ftp;
  $ftp = Net::FTP->new($server, Debug => 0) or throw "Cannot connect to $server: $@";
  $ftp->login("anonymous",'-anonymous@')    or throw "Cannot login ", $ftp->message;

  my @assemblies = qw(hg18 hg19);
  for my $assembly(@assemblies){
    my $dir_dcc = File::Spec->catdir('/', 'goldenPath', $assembly, 'encodeDCC');
    $ftp->cwd($dir_dcc) or throw "Cannot cd to $dir_dcc ", $ftp->message;
    my $labs = $ftp->ls;
    for my $lab(@$labs){
      next if($lab !~ /^wgEncode/);
      my $dir_lab = File::Spec->catdir($dir_dcc, $lab);
      $ftp->cwd($dir_lab) or throw "Cannot cd to $dir_lab", $ftp->message;
      my $local_dir = File::Spec->catdir($base_data_dir, $dir_lab);
      make_path($local_dir);

      my $local_files_txt  = File::Spec->catfile($local_dir, 'files.txt');
      $ftp->get('files.txt', $local_files_txt)
      or warn "No files.txt in $dir_lab\tFTP message:", $ftp->message;
      $local_files_txt =~ s/files\.txt/md5sum.txt/;
      $ftp->get('md5sum.txt', $local_files_txt)
      or warn "No md5sum.txt in $dir_lab\tFTP message:", $ftp->message;
    }
  }
}

=head2 _store_cell_feature_type

  Argument 1  : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Argument 2  : String - Name used as key in $data HASH
  Argument 3  : HASHREF - containing Cell or FeatureType objects
  Returntype  : Bio::EnsEMBL::Funcgen::CellType or Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions  : Missing arguments
  Description : PRIVATE - stores Cell or FeatureType

=cut

sub _store_cell_feature_type {
  my ($db, $name, $data) = @_;
  assert_ref($db,   'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  assert_ref($data, 'HASH', 'Cell/FeatureType cache');

  if(! ((exists $data->{$name} && defined $data->{$name})) ){
    throw("$name not defined in data.");
  }

  my $object  = $data->{$name};
  my $adaptor =  (check_ref($object, 'Bio::EsnEMBL::Funcgen::CellType')) ?
    $db->get_CellTypeAdaptor : $db->get_FeatureTypeAdaptor;

  ($object) = @{$adaptor->store($data->{$name})};
  say 'Added '.ref($object).":\t$name";
  return $object;
}

=head2 _read_md5sum_txt

  Argument 1  :
  Returntype  :
  Exceptions  : Unknown format
  Description : Only for FASTQ! Reads md5sum file, which is needed as some hg18 labs did not include
                it in files.txt

=cut

sub _read_md5sum_txt {
  my ($file) = @_;
  return undef if(!-e $file);
  my $hash;
  open(my $fh,'<',$file) or die "Can not open/access '$file'\n$!";
    while(my $line  = <$fh>){
      next unless($line =~ /fastq\.gz/);
      $line =~
      chomp($line);
      my @line = split(/\s+/,$line);
      if($line[0] =~ /^(wgEncode\w+?)\./){
        $hash->{$1} = $line[1];
      }
      elsif($line[1] =~ /^(wgEncode\w+?)\./){
        $hash->{$1} = $line[0];
      }
      else{
        throw("Unrecognised format: $line");
      }
    }
  close($fh);
  return $hash;
}

#Add no overlap mode?

#This is dependant on the following fasta header format:  >seq_region_name NA slice_name
#e.g. >18 dna:chromosome chromosome:GRCh37:18:1:78077248:1 chromosome 18

sub randomise_bed_file{
  my ($input_bed, $output_bed, $fasta_header_file, $sort_options) =
   rearrange([qw(INPUT_BED OUTPUT_BED FASTA_HEADER_FILE SORT_OPTIONS)], @_);

  $sort_options ||= '-k1,1 -k2,2n -k3,3n';
  my $ifh     = open_file($input_bed);
  my $ofh     = open_file($output_bed.'.tmp', "| sort $sort_options > \%s");
  #This sort is not working!! But the print is?

  my $headers = open_file($fasta_header_file);
  my ($line, %chrom, @orig);

  while(($line = $headers->getline) && defined $line){
    chomp $line;
    #This is new fasta header format as of release 76
    #This needs moving to SeqTools or similar, so we always use the same code for fasta header handling!

    my($sr_name, undef, $slice_name) = split(/\s+/, $line);
    $sr_name =~ s/^>//;
    (undef, undef, undef, $chrom{$sr_name}->{min}, $chrom{$sr_name}->{max}) = split(/:/, $slice_name);

    if(! $chrom{$sr_name}->{min}){
      throw("no min found for $line");
    }

    if($chrom{$sr_name}->{min} != 1){
      throw("$sr_name has high start in genome file:\t$$fasta_header_file\n".
        "Script needs updating to deal with this");
    }
  }
  $headers->close;

  my $wrote   = 0;
  my $skipped = 0;

  while(($line = $ifh->getline) && defined $line){
    chomp $line;
    #push @orig, [split("\t", $line)];
    my ($sr_name, $start, $end, $id, undef, $strand) = split("\t", $line);

    #Now assumes bed standard 0 based coords
    my $len = $end - $start; ##  +1;

    if($len >= $chrom{$sr_name}->{max}){
      warn "Length of feature($len) greater than length of genome region(".
        $chrom{$sr_name}->{max}.
        ")\nSkipping likely artefactual region:\t$sr_name\t$start\t$end\n";
      $skipped++;
      next;
    }

    my $new_start = int(rand($chrom{$sr_name}->{max} - $len +1));
    # +1 as we want to use 0 to iradicate the int rounding bias towards 0
    # if we get 0, then we set it to the under-represented max
    $new_start  ||= $chrom{$sr_name}->{max} - $len + 1;
    $new_start--; #Make 1/2 open
    my $new_end   = $new_start + $len;
    $strand = '.' if ! defined $strand; #Avoid undef warnings
    #Buffer here!
    $wrote++;
    print $ofh join("\t", ($sr_name, $new_start, $new_end, '.', '.', $strand))."\n";
  }
  $ifh->close;
  $ofh->close;
  run_system_cmd("mv ${output_bed}.tmp $output_bed");

  if(! $wrote){
    throw("Failed to write any mock peaks to:\t$output_bed");
  }
  else{
    warn "Wrote ${wrote}/".($wrote + $skipped)." mock peaks to:\t$output_bed\n";
  }

  return;
}



# get the individual sequences from the genome file and put them in files
# which have the fasta id as their name and an extension of .fa
# optionally reduce chromosome name to chr_name as in ensembl databases
# returns the list of files produced or dies

#TODO This should just delete and over-write unless recover is specified

sub explode_fasta_file{
  my $input_fasta = shift;
  my $target_dir  = shift;
  my $assembly    = shift;
  my $vlevel      = shift || 0;


  if (! (defined $input_fasta && -f $input_fasta)){
    throw("Input fasta file is not defined or does not exist:\t$input_fasta");
  }

  my @lines;

  if(! -d $target_dir){
    print "Exploding fasta file:\t$input_fasta\nTo:\t\t\t$target_dir\n" if $vlevel;
    my $tmp_dir = "${target_dir}_tmp";
    #Don't make this a tmp subdir or $target_dir, as this will obviously
    #make the $target_dir exist, leading to errors on retry
    #if this falls over

    run_system_cmd("rm -f $tmp_dir") if -d $tmp_dir;
    run_system_cmd("mkdir -p $tmp_dir");
    run_system_cmd("fastaexplode -f $input_fasta -d $tmp_dir");
    @lines = run_backtick_cmd("ls -1 ${tmp_dir}/*.fa");
    #file name expression results in paths being returned instead of basenames

    print 'Cleaning '.scalar(@lines)." exploded fasta files\n" if $vlevel;
    # now we alter the file names if short chromosome names have been requested
    if($assembly){

      foreach my $file_path (@lines){
        my $name = basename $file_path;# ls -1 always gives basename!
        # sed 's/chromosome:$assembly:/chr/'|sed 's/supercontig:$assembly://' |sed 's/supercontig:://' | sed 's/:[0-9:]*//'

        #This should all be in the first block below?
        #And make it generic to rename all slice named files to seq_region named files?
        #This is really only important if we have mismatched header formats which may result in
        #mismatched file name formats
        $name =~ s/chromosome:$assembly://o;
        $name =~ s/supercontig:$assembly://o;
        $name =~ s/supercontig:://o;
        $name =~ s/:[0-9:]*//;
        my $path = dirname $file_path;
        my $new = $path.'/'.$name;
        #This is not working as all but the Y chr are already short names
        #temp hack to get it running. This

        my $command;

        if($file_path =~ /chromosome:GRCh37:Y/){
          $command = "fastaclean -a -f $file_path > $new ; rm -f $file_path";
        }
        else{
          my $path = dirname $file_path;
          my $new = $path.'/tmp';
          $command = "fastaclean -a -f $file_path > $new ;".
          " rm -f $file_path ;mv $new $file_path";
        }

        print "Executing: $command\n" if $vlevel > 1;
        run_system_cmd($command);
      }
    }
    else{  # just do fastaclean

      foreach my $file_path (@lines){
        my $new = dirname($file_path).'/tmp';
        run_system_cmd("fastaclean -a -f $file_path > $new");
        run_system_cmd("mv $new $file_path");
      }
    }

    run_system_cmd("mkdir -p $target_dir");
    run_system_cmd("mv $tmp_dir/* $target_dir/");
    run_system_cmd("rm -rf $tmp_dir");
  }
  #else Assume this is correct if we have some files

  @lines = run_backtick_cmd("ls -1 ${target_dir}/*.fa", 1);

  if(! scalar(@lines)){
    throw("ERROR: No fasta files produced by splitting:\t".$input_fasta);
  }

  return \@lines;
}

# Environment should call ascript to write this and set path in env for pipeline config
# else will be called in tmp file mode for each analysis


sub write_chr_length_file{
  my ($slices, $out_file) = @_;
  my $fh;

  if(! defined $out_file){ 
    ($fh, $out_file) = tempfile(); 
    #DIR => '/tmp/'); #, UNLINK => 0); # Do not delete on exit? 
  }
  else{
    # Use flock here for safety?
    # over-write by default?
    $fh = open_file($out_file, '>');
  }

  foreach my $slice (@{$slices}) {
    print $fh join("\t", ($slice->seq_region_name, $slice->end - $slice->start + 1)) . "\n";
  }

  close $fh;
  return $out_file;
}  


1;

