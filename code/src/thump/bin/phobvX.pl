#! /usr/bin/perl

use strict;
use Getopt::Long;
use FindBin;

use APPRIS::Utils::Logger;

####################
# Input parameters #
####################
my ($input) = undef;
my ($output) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'input=s'			=> \$input,
	'output=s'			=> \$output,
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


# Main subroutine
my ($local_pwd) = $FindBin::Bin."/../";
my ($dir_programs) = $local_pwd."/Predictors/phobius";

unless ( -e $output and (-s $output > 0) ) {
	eval {
		#my ($cmd) = "perl $dir_programs/phobius.pl $input > $output";
		my ($cmd) = "phobius.pl $input > $output";
		$logger->debug("\t-- script: $cmd\n");
		system($cmd);	
	};
	exit 1 if($@);	
}
