#! /usr/bin/perl

use strict;
use Getopt::Long;
use FindBin;
use Bio::SeqIO;

use APPRIS::Utils::Logger;

####################
# Input parameters #
####################
my ($db_file) = undef;
my ($name) = undef;
my ($input) = undef;
my ($out_memsat) = undef;
my ($out_chk) = undef;
my ($out_blast) = undef;
my ($out_align) = undef;
my ($tmp_dir) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'db=s'				=> \$db_file,
	'name=s'			=> \$name,
	'input=s'			=> \$input,
	'out-memsat=s'		=> \$out_memsat,
	'out-chk=s'			=> \$out_chk,
	'out-blast=s'		=> \$out_blast,
	'out-align=s'		=> \$out_align,
	'tmp-dir=s'			=> \$tmp_dir,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();
my ($logfilename) = $logger->logpath().'/'.$logger->logfile();

###################
# Internal Values #
###################
my ($local_pwd) = $FindBin::Bin;
my ($dir_program) = $local_pwd."/../Predictors/Memsat3/trunk";

##############
# MEMSAT 3.0 #
##############

if ( -e $input && (-s $input > 0) ) {
	eval {
		#my ($cmd) = "perl $dir_program/runmemsat.pl ".
		my ($cmd) = "runmemsat.pl ".
									" --db=$db_file ".
									" --name=$name ".
									" --input=$input ".
									" --out-memsat=$out_memsat ".
									" --out-chk=$out_chk ".
									" --out-blast=$out_blast ".
									" --tmp-dir=$tmp_dir ".
									" 2>&1 1>> $logfilename ";
		$logger->debug("\t-- script: $cmd\n");
		my (@out) = `$cmd`;	
	};
	$logger->debug("running runmemsat.pl\n") if($@);
	
	eval {
		my ($cmd) = "perl $local_pwd/Seq_parservX.pl --db=$db_file --in-seq=$input --in-psiblast=$out_blast --output=$out_align";
		$logger->debug("\t-- script: $cmd\n");		
		my (@out) = `$cmd`;
	};
	$logger->debug("running Seq_parservX\n") if($@);
}





