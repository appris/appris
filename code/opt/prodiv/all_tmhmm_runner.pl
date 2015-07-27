#!/usr/bin/perl -w

use strict;
use FindBin;

# inparameters are <program to run (s, pro, prodiv)> <prodiv-tmhmm directory> <modfile/seqfile directory> <ouput directory>
if ( $#ARGV != 3 ) {
  printf "Usage: all_tmhmm_runner.pl <prog> <modfiledir/seqfiledir> <outputdir> <tmpdir>\n";
  exit;
}

my ($name) = $ARGV[0];
my ($seqdir) = $ARGV[1]."/";
my ($outdir) = $ARGV[2]."/";
my ($tmpdir) = $ARGV[3]."/";
my ($progdir) = $FindBin::Bin.'/';

my ($seqnamefile) = $tmpdir."/$name.seqname.tmp";
my ($hmmnamefile) = $tmpdir."/$name.hmm.tmp";
my ($hmg_prodiv_file) = $tmpdir."/$name.prodiv";
my ($hmg_prodiv_res_file) = $outdir."/$name.prodiv.res";

# if exist prodiv result, we exist
if ( -e $hmg_prodiv_res_file and (-s $hmg_prodiv_res_file > 0) ) {
	exit 0;
}

# list all modfiles to modfile-file
eval {
	my ($modfile) = $seqdir."/$name.mod";	
	my $seqfiles = `ls $modfile`;
	open (SEQNAMEFILE, ">$seqnamefile") or die "Could not open temporary seqnamefile\n";
	print SEQNAMEFILE "$seqfiles";
	close (SEQNAMEFILE);	
};
exit 1 if($@);

# create hmmnamefile
eval {
	open (PROTNAMEFILE, ">$hmmnamefile") or die "Could not open temporary protnamefile\n";
	print PROTNAMEFILE "$hmg_prodiv_file\n";
	close (PROTNAMEFILE);	
};
exit 1 if($@);

# create prodivtmhmm from prodiv template
my ($hmg_template_file) = $progdir."/HMMS/PRODIV_TMHMM_0.92b.hmg";
$/=undef;
local(*HMG_FILE);
open(HMG_FILE, $hmg_template_file) or die "Can not open $hmg_template_file file\n";
my ($hmg_content) = <HMG_FILE>;
close(HMG_FILE);
$/='\n';

my ($amino_multi_file) = $progdir."/util/amino_multi.pri";
$hmg_content =~ s/EMISSION PRIORFILES\:([^\n]*)/EMISSION PRIORFILES\: $amino_multi_file/g;
$hmg_content =~ s/Emission prior file\:([^\n]*)/Emission prior file\: $amino_multi_file/g;

open(HMG_FILE2, ">$hmg_prodiv_file") or die "Could not open temporary protnamefile\n";
print HMG_FILE2 "$hmg_content\n";
close(HMG_FILE2);


# run prodiv
eval {
	my ($cmd) = "$progdir/modhmm0.92b/modhmms -m $hmmnamefile -s $seqnamefile -f msa -o $outdir/ -r $progdir/util/replacement_letter_multi.rpl -M GM -L -c 1 --max_d";
	print $cmd."\n";
	system ($cmd);		
};
exit 1 if($@);


#remove modile-file
#`rm $outdir/seqnamefile.tmp`;
#remove protnamefile
#`rm $outdir/protnamefile.tmp`;
