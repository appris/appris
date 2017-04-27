#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile prepare_workspace );

use lib "$FindBin::Bin";
use matador3d;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$WSPACE_TMP
	$WSPACE_CACHE
	$RUN_PROGRAM
	$PROG_DB
	$PROG_DB_DIR
);
       
# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($input_file) = undef;
my ($output_file) = undef;
my ($appris) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf|c=s'			=> \$config_file,
	'input|i=s'			=> \$input_file,
	'output|o=s'		=> \$output_file,
	'loglevel|ll=s'		=> \$loglevel,
	'logfile|lf=s'		=> \$logfile,
	'logpath|lp=s'		=> \$logpath,
	'logappend|la'		=> \$logappend,
);

# Required arguments
unless ( defined $config_file and defined $input_file and defined $output_file )
{
    print `perldoc $0`;
    exit 1;
}

# Get conf vars
my ($cfg) 			= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
our $CACHE_FLAG		= $cfg->val('MATADOR3D2_VARS', 'cache');
$RUN_PROGRAM		= $ENV{APPRIS_PRG_OPT_HMMER_31b2_BIN_DIR}.'/'.$cfg->val('MATADOR3D2_VARS', 'program');
our $PROG_DB		= $cfg->val('MATADOR3D2_VARS', 'db');
$PROG_DB_DIR		= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$PROG_DB;


# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);


#####################
# Method prototypes #
#####################


#################
# Method bodies #
#################
# Main subroutine
sub main()
{

	# Load input file, $sequences is a hash, key1 is gene_id, key 2 is transcript_id, the value is the sequence of the correspondng transcript
	$logger->info("-- read fasta file\n");
	my $sequences = matador3d::read_fasta_file($input_file);

	# Run each hmmscan for each sequence and store results in temporary directory
	$logger->info("-- run hmmscan\n");
	matador3d::run_hmmscan($RUN_PROGRAM, $sequences, $WSPACE_CACHE, $PROG_DB_DIR);
	
	# Score isoforms
	$logger->info("-- score isoforms\n");
	my $scoring = matador3d::run_scoring($sequences, $WSPACE_CACHE);
		
	# Print records by transcript ---------------
	my ($print_out) = printStringIntoFile($scoring, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}
	
	$logger->finish_log();
	
	exit 0;	
}


main();



__END__

=head1 NAME

matador3d

=head1 DESCRIPTION

Run Matador3D (version2) program

=head1 SYNOPSIS

matador3d

=head2 Required arguments:

    --conf <Config file>

    --input <Fasta sequence file>

    --output <Annotation output file>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
    

=head1 EXAMPLE

perl matador3d.pl

	--conf=../conf/pipeline.ini

	--input=examples/ENSG00000016864.faa
	
	--output=examples/ENSG00000016864.output


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
