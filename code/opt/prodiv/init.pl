#!/usr/bin/perl -w

#use strict;


##############################################################################################
# This script should be run after the prodiv-tmhmm program has been downloaded. Basically    #
# it sets some paths to the directory where podiv-tmhmm has been installed, which are needed #
# when running the program later.                                                            #
##############################################################################################

if($#ARGV != 0) {
  printf "Usage: init.pl <progdir>\n";
  exit;
}

my $progdir = $ARGV[0]."/";
my $utildir = "$progdir"."util/";

#create hmm_txt_files
my $hmm_file = "$progdir"."HMMS/S_TMHMM_0.92b.hmg\n";
open HMMNAMEFILE, ">"."$progdir"."HMM_FILES/S_TMHMM_0.92b.txt" or die "Could not set hmm file path"."HMM_FILES/S_TMHMM_0.92b.txt\n";
print HMMNAMEFILE "$hmm_file";
close HMMNAMEFILE;

$hmm_file = "$progdir"."HMMS/PRO_TMHMM_0.92b.hmg\n";
open HMMNAMEFILE, ">"."$progdir"."HMM_FILES/PRO_TMHMM_0.92b.txt" or die "Could not set hmm file path\n";
print HMMNAMEFILE "$hmm_file";
close HMMNAMEFILE;

$hmm_file = "$progdir"."HMMS/PRODIV_TMHMM_0.92b.hmg\n";
open HMMNAMEFILE, ">"."$progdir"."HMM_FILES/PRODIV_TMHMM_0.92b.txt" or die "Could not set hmm file path\n";
print HMMNAMEFILE "$hmm_file";
close HMMNAMEFILE;


#set path of amino_multi.pri in hmm files
my $hmm = "$progdir"."HMMS/S_TMHMM_0.92b.hmg";
my $hmm_tmp = "$progdir"."HMMS/S_TMHMM_0.92b.hmg.tmp";
open HMM_IN_FILE, "$hmm" or die "Could not open hmmfile\n";
open HMM_OUT_FILE, ">"."$hmm_tmp" or die "Could not open hmm outfile\n";
while(<HMM_IN_FILE>) {
  my $row = $_;
  if($row =~ 'amino_multi.pri') {
    $row =~ s/:.*amino_multi\.pri/: ${utildir}amino_multi\.pri/g;
  }
  print HMM_OUT_FILE "$row";
}
close HMM_IN_FILE;
close HMM_OUT_FILE;

my $res = `mv ${hmm_tmp} ${hmm}`;

$hmm = "$progdir"."HMMS/PRO_TMHMM_0.92b.hmg";
$hmm_tmp = "$progdir"."HMMS/PRO_TMHMM_0.92b.hmg.tmp";
open HMM_IN_FILE, "$hmm" or die "Could not open hmmfile\n";
open HMM_OUT_FILE, ">"."$hmm_tmp" or die "Could not open hmm outfile\n";
while(<HMM_IN_FILE>) {
  my $row = $_;
  if($row =~ 'amino_multi.pri') {
    $row =~ s/:.*amino_multi\.pri/: ${utildir}amino_multi\.pri/g;
  }
  print HMM_OUT_FILE "$row";
}
close HMM_IN_FILE;
close HMM_OUT_FILE;

$res = `mv ${hmm_tmp} ${hmm}`;

$hmm = "$progdir"."HMMS/PRODIV_TMHMM_0.92b.hmg";
$hmm_tmp = "$progdir"."HMMS/PRODIV_TMHMM_0.92b.hmg.tmp";
open HMM_IN_FILE, "$hmm" or die "Could not open hmmfile\n";
open HMM_OUT_FILE, ">"."$hmm_tmp" or die "Could not open hmm outfile\n";
while(<HMM_IN_FILE>) {
  my $row = $_;
  if($row =~ 'amino_multi.pri') {
    $row =~ s/:.*amino_multi\.pri/: ${utildir}amino_multi\.pri/g;
  }
  print HMM_OUT_FILE "$row";
}
close HMM_IN_FILE;
close HMM_OUT_FILE;

$res = `mv ${hmm_tmp} ${hmm}`;



#set program path in program run script
my $prog_script = "$progdir"."all_tmhmm_runner.pl";
my $prog_script_tmp = "$progdir"."all_tmhmm_runner.pl.tmp";
open RUNNER_IN_FILE, "$prog_script" or die "Could not open all_tmhmm_runner file\n";
open RUNNER_OUT_FILE, ">"."$prog_script_tmp" or die "Could not open all_tmhmm_runner outfile\n";

while(<RUNNER_IN_FILE>) {
  my $row = $_;
  if($row =~ 'progdir_line_is_set') {
    $row =~ s/".*"/"${progdir}"/g;
  }
  print RUNNER_OUT_FILE "$row";
}
close RUNNER_IN_FILE;
close RUNNER_OUT_FILE;

$res = `mv ${prog_script_tmp} ${prog_script}`;
