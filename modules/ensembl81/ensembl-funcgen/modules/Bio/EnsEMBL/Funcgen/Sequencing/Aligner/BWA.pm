=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::Aligner::BWA

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Sequencing::Aligner::BWA;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );

use base qw( Bio::EnsEMBL::Funcgen::Sequencing::Aligner);#Does not import


# NOTE:
# If bwa is seg faulting with little to no useful output, it is probably running
# out of memory. Try upping the bsub -R/M spec


#todo use Rob Davies hardened bwa https://github.com/daviesrob/bwa
#The following branches are available:
#0.5.10-mt_fixes    0.5.10 with multithreaded samse/sampe
#0.6.2_fixes        0.6.2 released version
#0.6.2-mt           0.6.2 with multithreaded samse/sampe
#master_fixes       The current bwa master version
#or see if Heng has incorporated these into the main bwa github rep?

#todo look into multi-threading to speed thing up -t #nodes
#will this just reduce our submission rate as we will require an entire node
#to be free before we get a slot
#maybe try -t 2?
#What is more efficient here, multiple cpu jobs or more single cpu jobs?
#In terms of optimsed run time, curretnly 1cpu/chunk takes ~ 30mins
#Using >1cpu would likely mean upping the chunk size, or running > 1 in a batch series


sub new {
  my $caller      = shift;
  my $class       = ref($caller) || $caller;
  my $self        = $class->SUPER::new(@_);
  my ($fasta_fai) = rearrange(['SAM_REF_FAI'], @_);
  
  if(! defined $fasta_fai) {
    throw("BWA requires a -sam_ref_fai file to add a header to the output\nParams passed:\n@_");   
  }
  elsif(! -f $fasta_fai){
    throw("-sam_ref_fai parameter does not exist or is not a file:\n".$fasta_fai);
  }
  
  $self->{sam_ref_fai} = $fasta_fai; 
  
  warn "Created ".$self if $self->debug;
  return $self;
}

sub sam_ref_fai { return shift->{sam_ref_fai}; }


#todo handle bwa dir ans samtools dir better
#we could take a bin_dir and a program name/file
#todo handle samse/sampe? Default to samse
#this isn't really a parameter, and there are several commands to provide params for


#although output format is sam
#we always need to sort and convert to bam for merge

#Todo
#1 Move sort functionality to SeqTools? Careful not to have reciprocal use statements
# as SeqTools will use this will running an aligner
# Hence we may need to move the sort call back to RunAligner
# Also add @HD SO tag! Or does a more recent version
#of samtools do this? Check manifest for releases past 0.1.18, sourceforge only list 0.1.19
#which seems to be the one with the BGZD bug re-introduced
#Have the project hosting moved elsewhere?
#gut hub lists a release candidate for 0.2.0, but this is not likely

sub run {
  my $self        = shift;
  my $query_file  = $self->query_file;
  my $input_dir   = $self->input_dir;
  my $output_dir  = $self->output_dir;
  my $ref_file    = $self->target_file;
  my $bwa_bin     = $self->program_file;
  my $fasta_fai   = $self->sam_ref_fai;
  my $bin_dir     = '';
  
  if($bwa_bin =~ /\//o){
    ($bin_dir = $bwa_bin) =~ s/(.*\/)[^\/].*/$1/go;
  }
  
  #move this to Aligner as it is based on generic split output
  #only if outfile_prefix is required by other aligner
  
  
  #Handle split or whole files here
  #my $outfile_prefix;
  
  #Throw if this is gzipped, as gunzipping will invalidate any md5 checksum
  #force manual gunzipping ot safer zcating, which will maintain the original
  #file and keep the checksum valid
  
  
  
  my $outfile_prefix;
  #Strip out the fastq
  
  if($query_file =~ /\.fastq_([0-9]+)$/o){
    ($outfile_prefix = $query_file) =~ s/\.fastq_([0-9]+)$/.$1/;
  }
  else{
    ($outfile_prefix = $query_file) =~ s/\.fastq$//;   
  }
  warn "Query file:\t$query_file\nOutfile prefix:\t$outfile_prefix\n" if $self->debug;
  
  my $sai_file      = "${output_dir}/${outfile_prefix}.sai";
  my $samse_file    = "${output_dir}/${outfile_prefix}.samse.sam";
  my $unsorted_file = "${output_dir}/${outfile_prefix}.samse.bam.unsorted";
  my $bam_file      = "${output_dir}/${outfile_prefix}.bam";
  my @tmp_files = ($sai_file, $samse_file, $unsorted_file);
  
  my $rm_cmd = 'rm -f '.join(' ', @tmp_files);
  warn "Removing pre-existing intermediates:\n$rm_cmd $bam_file\n" if $self->debug;
  run_system_cmd($rm_cmd.' '.$bam_file, 1);#no exit flag in case they aren't there
  
  #assume samtools is in the same dir as bwa
  #no, just assume we are using the bin_dir param, else assume it is in the $PATH
  #so we need to pass bin_dir aswell?
  #No we can parse from the program_file as this will already have been prefixed with the bin_dir
  

  #If using -q make sure we've got the correct Sanger quality scores...
#  $bwa_cmd .= " | $bwa_bin samse $bwa_index - $input_file";
#  $bwa_cmd .= " | samtools view -uS - ";
#  $bwa_cmd .= " | samtools sort - ${input_file}.sorted";
#  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }


  #Piping all these together was causing uncaught errors

  ### FIND SUFFIX ARRAY COORDS OF SEQS
  #This seg faults if it is run locally as it runs out of memory
  my $cmd = "$bwa_bin aln $ref_file ${input_dir}/${query_file} > $sai_file";
  warn "Running:\n$cmd\n" if $self->debug;  
  run_system_cmd($cmd);
  

  ### GENERATE SINGLE END READ ALIGNMENTS
  #todo add -P option here to load index into memory which should speed execution
  #but will require ~5GB of memory
  #THIS ONLY WORKS FOR SAMPE?
    
  $cmd = "$bwa_bin samse $ref_file $sai_file ${input_dir}/${query_file} > $samse_file";
  warn "Running:\n$cmd\n" if $self->debug;  
  run_system_cmd($cmd);
  
  
  ### CONVERT TO BAM
  #-S input is sam, with header (else -t is required)
  #-u output is uncompressed bam, prefered for piping to other samtools cmds, although this has been shown to
  #be fragile, hence we have split the cmds up here
  #-h include header in output
  
  #This won't have a header output, as bwa samse doesn't output one!!!!!
  #so we need to add the fai file
  
  #we were using sam output here for collections
  #but we are deleting this at the moment anyway, 
  #so change to bam output and drop -S
  #also dropped -u as we are not piping here
  
  $cmd = "${bin_dir}/samtools view -t $fasta_fai -bh $samse_file > $unsorted_file"; 
  warn "Running:\n$cmd\n" if $self->debug;
  run_system_cmd($cmd);
  
  ### SORT BAM
  #Do this in parallel here in prep for chunk merge and sort in funnel job
  #This can give the following error output without returning an exit code?!
  #[bwa_aln_core] 4000000 sequences have been processed.
  #[samopen] SAM header is present: 194 sequences.
  #[bam_header_read] EOF marker is absent. The input is probably truncated.
  #[bam_sort_core] merging from 2 files...
  #This is actually a known bug and can be ignored
  #https://github.com/samtools/samtools/issues/18
  
  $cmd = "${bin_dir}/samtools sort $unsorted_file ${output_dir}/${outfile_prefix}";
  warn "Running:\n$cmd\n" if $self->debug;
  run_system_cmd($cmd);
  
  #This can die with the following output, but we don't exit!
  #[bam_header_read] invalid BAM binary header (this is not a BAM file).
  #[bam_sort_core] truncated file. Continue anyway.
  #Child process died with signal 11, without coredump
  #Error:
  

  
  #todo update sort flag in header, samtools should really do this
  #@HD     SO:coordinate
  
  #how are we going to handle samse here?
  #is this generic to all read aligners? nomenclature certainly isn't
  #merge job will need to know what the output format is.
  #so we need to make that available
  #although we are doing this here, as it is a pre-req for the merge step
  #just drop samse from file name for now, as we only do single end reads at present
  
  #To allow rerunning on failure, leave removing inputs to the funnel job
  warn "Removing:\n$rm_cmd\n" if $self->debug;
  run_system_cmd($rm_cmd);
  
  #return Name of output file. Will need to process output format params
  #if we implement this
  return $bam_file;
}




1;
