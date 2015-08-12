#!/usr/bin/perl -w

use strict;
use warnings;
use threads;

use Getopt::Long;
use FindBin;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_ENSEMBL_DB_FILE
	$CONFIG_INI_APPRIS_DB_FILE
	$LOGLEVEL
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_ENSEMBL_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/ensembldb.ini';
$LOGLEVEL					= 'INFO';

# Input parameters
my ($species) = undef;
my ($conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'species|s=s'		=> \$species,	
	'conf|c=s'			=> \$conf_file,	
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
if ( !defined $species && !defined $conf_file ) {
	print `perldoc $0`;
	exit 1;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();

#################
# Method bodies #
#################


# Main subroutine
sub main()
{
	# get the list of input genes with translation
	my (@data_genes);
	eval {
		my ($cmd) = "appris_check_ls_report -c $species -g";
		$logger->info("\n** script: $cmd\n");
		@data_genes = `$cmd`; map { s/\s+$// } @data_genes;
	};
	$logger->error("deleting log files of appris") if($@);
	
	# get the list of input transcripts with translation
	my (@data_transc);
	eval {
		my ($cmd) = "appris_check_ls_report -c $species -t";
		$logger->info("\n** script: $cmd\n");
		@data_transc = `$cmd`; map { s/\s+$// } @data_transc;
	};
	$logger->error("deleting log files of appris") if($@);	

	# get the list of annotation files
	my (@annot_appris);
	eval {
		my ($cmd) = "appris_check_ls_annots -c $species";
		$logger->info("\n** script: $cmd\n");
		@data_transc = `$cmd`; map { s/\s+$// } @data_transc;
	};
	$logger->error("deleting log files of appris") if($@);	

#$logger->info("DATA_GENES:\n".Dumper(@data_genes));
#$logger->info("DATA_TRANS:\n".Dumper(@data_transc));

	$logger->finish_log();
	
	exit 0;	
}

main();


1;

__END__

=head1 NAME

check_files

=head1 DESCRIPTION

Compare number the APPRIS result files and GENCODE/Ensembl data geneset  

=head1 SYNOPSIS

check_files

=head2 Required arguments (data input):

  -s, --species {string} <Name of species>
	* Hsap  - Homo sapiens -\n
	* Mmus  - Mus musculus -\n
	* Rnor  - Rattus norvegicus -\n
	* Drer  - Danio rerio -\n
	* Sscr  - Sus scrofa -\n
	* Dmel  - Drosophila melanogaster -\n
	* Cele  - Caenorhabditis elegans -\n
	* Lpar  - Lynx pardinus -\n  
  
 Or
 
  -c, --conf {file} <Config file name>

=head2 Optional arguments (log arguments):
	
  --loglevel=LEVEL <define log level (default: NONE)>
  
  --logfile=FILE <Log to FILE (default: *STDOUT)>
  
  --logpath=PATH <Write logfile to PATH (default: .)>
  
  --logappend= <Append to logfile (default: truncate)>

=head1 EXAMPLE

	check_files.pl -s Hsap

=head1 EXAMPLE2

	check_files.pl -c conf/script/apprisrc.XXX
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
