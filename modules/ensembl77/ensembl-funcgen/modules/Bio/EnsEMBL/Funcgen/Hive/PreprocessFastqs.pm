
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Hive::Funcgen::PrprocessFastqs

=head1 DESCRIPTION

This analysis validates the composition of the ResultSet and merges
and chunks control or signal fastq files appropriately, in preparation
for running the alignment of individual chunks.

=cut

package Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped run_system_cmd
                                               get_set_prefix_from_Set
                                               run_backtick_cmd check_file );
#use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( split_fastqs );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

#TODO... use and update the tracking database dependant on no_tracking...

#todo
# Add status tracking support for identifying which InputSubsets might already have been aligned?
# Is this actually relevant? The files may have been moved away? Options here are to force re-run
# or move files back from archive. We should enable archive access of bams?
# Normal mode would simply try and copy the files back? There maybe parallel processes doing trying to do this
# hence checksum is necessay here. Reruns of job may pick up complete copy later on.
# We need to copy the check sum file first.
# There is still a very slight risk of a race condition when files are tested and copied.
# Hence we probably need to flock here?
# dynamic file archiving should not run on controls, leave this as a manual post pipeline step
# to minimise these risks
# Add ignore_checksums flag, this should stil create them, jsut not check them?
#



sub fetch_input {  
  my $self = shift;
  #Set some module defaults
  $self->param('disconnect_if_idle', 1);
  
  $self->SUPER::fetch_input();
  my $rset = $self->fetch_Set_input('ResultSet');


  #RunAligner need not know anything about the ResultSet
  #So we can validate all that here and simply pass the files through to align
  #Need to leave module validation to RunAligner
  #but can pass the params directly, rather than in ResultSet
  

  
  #get/validate merge and run_controls
  #what about dataflow onwards, this will differ
  #dependant the next analysis 
  #e.g. DefineMergedDataSet, Run_BWA_and_QC_merged, DefineReplicateDataSet or BWA_ReplicateFactory
  #or can we do this implicitly just by check if we have set_names/ids set?
  
  #This may allow unmerged controls, if we set merge to 0 in the analysis config
  my $run_controls = $self->get_param_method('result_set_groups', 'silent') ? 1 : 0;
  $self->set_param_method('run_controls', $run_controls);
  my $merge        = $self->get_param_method('merge', 'silent', $run_controls); 
  $self->get_param_method('checksum_optional', 'silent');
  
     
  
  
  
  #we really need to define a Base Aligner class similar to the PeakCaller class
  #define a standard interface/requirements
  #location of indexes needs to be built here
  #based an analysis name
  #do we need a index_required sub?
  #Let's do all of this in BWA for now? 
  
  #Hmm this is just chunking the files! and submitting the individual jobs!
  #Hence we need more analyses:
  #1 to run the individual BWA jobs
  #2 to merge the alignments and perform QC!
  
  
  

  $self->get_param_method('fastq_chunk_size', 'silent', '16000000');#default should run in 30min-1h 
  #does this even support input sets yet?

  
  #We need to define the work dir here for the intermediate chunk/alignments files
  #output_dir here is for alignment (no need for repository)
  $self->get_output_work_dir_methods($self->alignment_dir($rset, 1, $run_controls));#default output_dir 
  return;
}

#TODO 
# 1 Extra ulatr-paranoid sanity check to test the last line of the input and the output
#   to make sure it isn't truncated

sub run {
  my $self         = shift;
  my $rset         = $self->ResultSet;
  my $run_controls = $self->run_controls;
  my $merge        = $self->merge;
 
  #Maybe we need to handle pre-aligned ResultSets here
  #check status and file
  
  #unsafe to check individual input_subsets, as these may have been processed
  #in an incomplete manner previously
  #Although we need to be mindful that a ResultSet may have had an InputSubset added
  #hence we can't re-use an old alignment file
  #this should be handled when creating/rolling back the ResultSet
  
  #This status needs to be CS specific!!
  #Actually, this should only be allowed if we are recovering
  #force should have rolled back the ResultSet ALIGNED status
  #although ideally this should be handled in the previous analysis
  #and flowed directly to DefineReplicate/MergedDataSet
  #Hence we should never reuse the merged fastq?
  #if we are recovering, we want the bam file (given the reps are the same)
  #if we are forcing or rolling back, then we should probably redo everything
  
  if($run_controls){
    my $exp = $rset->experiment(1);#control flag
    
    if($exp->has_status('ALIGNED_CONTROL')){
      throw("Need to implement force/recover_control_alignment. Found ALIGNED_CONTROL ResultSet:\t".
      $rset->name."(Control Experiment = ".$exp->name.")\n");
    }
  }
  elsif($rset->has_status('ALIGNED')){
    throw("Need to implement force/recover_alignment. Found ALIGNED ResultSet:\t".
    $rset->name."\n");
  }
  
  my @fastqs;
  my $throw = '';
  my $set_rep_suffix = '';
  my @issets = @{$rset->get_support};
  
  if((! $run_controls) && (! $merge) &&
    ($self->is_idr_FeatureType($rset->feature_type))){
  
    my @signal_issets = grep { ! $_->is_control } @issets;
  
    if(scalar(@signal_issets) != 1){
      $self->throw_no_retry('Expected 1 InputSubset(replicate) for IDR ResultSet '.
        $rset->name.' but found '.scalar(@signal_issets).' Specify merge, or restrict to 1 InputSubset');  
    }
    
    #We are not filtering for controls here!!
    $set_rep_suffix = '_TR'.$signal_issets[0]->replicate;
      
  }
  
  foreach my $isset(@issets){

    if(($isset->is_control && ! $run_controls) ||
       ($run_controls && ! $isset->is_control)){
      next;    
    }
 
    if(! $self->tracking_adaptor->fetch_tracking_info($isset)){
       $throw .= "Could not find tracking info for InputSubset:\t".
        $isset->name."\n";
      next;
    }
    
    if(! defined $isset->local_url){
      $throw .= "Found an InputSubset without a local_url, has this been downloaded?:\t".
        $isset->name."\n";
      next;  
    }

    my $found_path;
    my $params = {};#{gunzip => 1}; #NEVER DEFINE gunzip here!
    #Instead of gunzipping in the warehouse, zcat is now used to 
    #pipe directly split directly into the work area
    #This reduces tidy up and keeps footprint low, so we don't hit
    #out of space errors when running with a full warehouse
    #or a nearly full scratch space
    #This will also prevent any clashes between unzipping files in the warehouse
    #188 secs to zcat 227MB gzipped fastq
    #vs
    #10 secs to gunzip (to 1.2GB) and cat
    #This does not include rezip and tidy up time of ~90 secs 
    #(which could arguably be defered to after the pipeline run)
    #This is quite a large difference, but with and average of 2 or 3 reps 
    #this will probably make this run to ~10mins, which is negligable
    #compared to the down time from managing failed jobs due to out of space 
    #issues.
    
    #Set checksum   
    #As we already know we don't have a checksum, simply omit it here
    #which will mean validate_checksum will not be called
    #Alterntive is to set an undef checksum and checksum_optinal
    #This will cause validate_checksum to try and find a checksum from a file
    #But we know these checksums are stored in the DB
        
    if(defined $isset->md5sum || ! $self->checksum_optional ){
      $params->{checksum} = $isset->md5sum; 
    }
    
    my $local_url = $isset->local_url;
    #Look for gz files too. These would normally already be gzipped
    #if downloaded from a repository
    #But they may have been gzipped after processing if produced locally
    #add .tgz support here?
    #we can't do a md5 check if we don't match the url exactly
    eval { $found_path = check_file($local_url, 'gz', $params); };
 
    if($@){
      $throw .= "$@\n";
      next;  
    }
    elsif(! defined $found_path){
      $throw .= "Could not find fastq file, is either not downloaded, has been deleted or is in warehouse:\t".
        $local_url."\n";
      #Could try warehouse here?
    }
    elsif($found_path !~ /\.(?:t){0,1}gz$/o){
      #use is_compressed here?
      #This will also modify the original file! And potentially invalidate any checksumming
      $self->throw_no_retry("Found unzipped path, aborting as gzipping will invalidate any further md5 checking:\t$found_path");
      #run_system_cmd("gzip $found_path");
      #$found_path .= '.gz';  
    }
    
     
    push @fastqs, $found_path;  
  }
 
  throw($throw) if $throw;
  
  if((scalar(@fastqs) > 1) &&
     ! $merge){
    throw('ResultSet '.$rset->name.
      " has more than one InputSubset, but merge has not been specified:\n\t".
      join("\n\t", @fastqs));    
  }  
 
  my $set_prefix = get_set_prefix_from_Set($rset, $run_controls).
    '_'.$rset->analysis->logic_name.$set_rep_suffix; 
 
  #Will need to eval this so we can throw_no_retry 
  #split_fastqs(\@fastqs, $set_prefix, )
 
 
 
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
  
  warn "DEACTIVATED FASTQC FOR NOW:\nfastqc -f fastq -o ".$self->output_dir." @fastqs";
  #run_system_cmd('fastqc -o '.$self->output_dir." @fastqs");
  
 
  #todo parse output for failures
  #also fastscreen?

  warn("Need to add parsing of fastqc report here to catch module failures");
  
  #What about adaptor trimming? and quality score trimming?
  #FASTX? quality_trimmer, clipper (do we have access to the primers?) and trimmer?
 
  

     
        

  #For safety, clean away any that match the prefix
  run_system_cmd('rm -f '.$self->work_dir."/${set_prefix}.fastq_*", 1);
  #no exit flag, in case rm fails due to no old files
     
  my @du = run_backtick_cmd("du -ck @fastqs");   
  (my $pre_du = $du[-1]) =~ s/[\s]+total//;   
     
  my $cmd = 'zcat '.join(' ', @fastqs).' | split --verbose -d -a 4 -l '.
    $self->fastq_chunk_size.' - '.$self->work_dir.'/'.$set_prefix.'.fastq_';
  $self->helper->debug(1, "Running chunk command:\n$cmd");
  my @split_stdout = run_backtick_cmd($cmd);
  (my $final_file = $split_stdout[-1]) =~ s/creating file \`(.*)\'/$1/;
  
  if(! defined $final_file){
    $self->throw_no_retry('Failed to parse (s/.*\`([0-9]+)\\\'/$1/) final file '.
      ' from last split output line: '.$split_stdout[-1]);  
  }
  
  #Get files to data flow to individual alignment jobs
  my @new_fastqs = run_backtick_cmd('ls '.$self->work_dir."/${set_prefix}.fastq_*");
  @new_fastqs    = sort {$a cmp $b} @new_fastqs;
  
  #Now do some sanity checking to make sure we have all the files
  if($new_fastqs[-1] ne $final_file){
    $self->throw_no_retry("split output specified last chunk file was numbered \'$final_file\',".
      " but found:\n".$new_fastqs[-1]);  
  }
  else{
    $final_file =~ s/.*_([0-9]+)$/$1/;
    $final_file  =~ s/^[0]+//;
    
    $self->debug(1, "Matching final_file index $final_file vs new_fastq max index ".$#new_fastqs);
    
    if($final_file != $#new_fastqs){
      $self->throw_no_retry('split output specified '.($final_file+1).
        ' file(s) were created but only found '.scalar(@new_fastqs).":\n".join("\n", @new_fastqs));  
    }  
  }
  
  #and the unzipped files are at least as big as the input gzipped files
  @du = run_backtick_cmd("du -ck @new_fastqs");   
  (my $post_du = $du[-1]) =~ s/[\s]+total//; 
  
  $self->helper->debug(1, 'Merged and split '.scalar(@fastqs).' (total '.$pre_du.'k) input fastq files into '.
    scalar(@new_fastqs).' tmp fastq files (total'.$post_du.')');
  
  if($post_du < $pre_du){
    $self->throw_no_retry("Input fastq files totaled ${pre_du}k, but output chunks totaled only ${post_du}k");  
  }
  
  
  $self->set_param_method('fastq_files', \@new_fastqs);
  my %batch_params = %{$self->batch_params};
 
  foreach my $fq_file(@{$self->fastq_files}){
  
    #Data flow to RunAligner for each of the chunks 
    #do we need to pass result set to aligner?
    #Would need to pass gender, analysis logic_name 
    
    $self->branch_job_group(2, [{%batch_params,
                                 set_type   => 'ResultSet',
                                 set_name   => $rset->name,
                                 dbID       => $rset->dbID,
                                 output_dir => $self->work_dir, #we could regenerate this from result_set and run controls
                                 fastq_file => $fq_file}]);
  }

  my %signal_info;
  
  if($run_controls){
    %signal_info = (result_set_groups => $self->result_set_groups);
    #for flow to MergeControlAlignments_and_QC
  }

  # Data flow to the MergeQCAlignements job 
  
  #This was a config problem, we had a circular semaphore using the same branch
  
  #There is some redundancy between bam_file and garbage here
  #we could let the Merge job do the bam_file name generation replacement
  #
  
  $self->branch_job_group(3, [{%batch_params,
                             set_type   => 'ResultSet',
                             set_name   => $rset->name,
                             dbID       => $rset->dbID,
                             fastq_files => $self->fastq_files,
                             output_dir => $self->output_dir,
                             set_prefix => $set_prefix,
                             %signal_info}]);
  return;
}


sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}



1;
