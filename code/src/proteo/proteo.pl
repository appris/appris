#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile );

###################
# Global variable #
###################
use vars qw(
	$PROG_DB
	$PROG_CUTOFF
);

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($data_file) = undef;
my ($output_file) = undef;
my ($appris) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf=s'			=> \$config_file,
	'data=s'			=> \$data_file,	
	'output=s'			=> \$output_file,
	'appris'			=> \$appris,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $config_file and defined $data_file and defined $output_file )
{
    print `perldoc $0`;
    exit 1;
}

# Get conf vars
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$PROG_DB				= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('PROTEO_VARS', 'db');
$PROG_CUTOFF			= $cfg->val('PROTEO_VARS', 'cutoff');

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
	# Declare vars
	my ($output_content) = '';
	
	# obtain the gene id
	my ($gene_id);
	eval
	{
		$logger->info("-- obtain gene_id\n");
		my ($cmd) = "awk '{if (\$3 == \"gene\") {print \$10}}' $data_file | sed 's/[\"|;]//g' | sed 's/\\.[0-9]*\$//'";
		$logger->debug("$cmd\n");						
		my (@gene_id_out) = `$cmd`;
		if ( scalar(@gene_id_out) == 0 ) {
			$logger->error("we can not retrieve gene_id\n")	
		}
		else {
			$gene_id = $gene_id_out[0];
			$gene_id =~ s/\s*//g;
		}		
	};
	$logger->error("obtaining gene_id: $!\n") if($@);
	
	# retrieve peptides from text file
	if ( defined $gene_id and ($gene_id ne '') ) {
		$logger->info("-- retrieve peptides from a gene\n");
		my ($cmd) = "awk -F ',' '\$2 ~ /$gene_id/ {print \$0}' $PROG_DB"; # version with commas
		#my ($cmd) = "awk '\$2 ~ /$gene_id/ {print \$0}' $PROG_DB"; # version with tabs
		$logger->debug("$cmd\n");						
		my (@proteo_txt_out) = `$cmd`;
		if ( scalar(@proteo_txt_out) == 0 ) {
			$logger->warning("there are not peptides\n")	
		}
		else {
			$output_content = join("",@proteo_txt_out);
		}		
	}
		
	# print output
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}
	
	$logger->finish_log();
		
	exit 0;	
}

main();

__END__

=head1 NAME

proteo

=head1 DESCRIPTION

Run the main script of proteo pipeline

=head1 SYNOPSIS

proteo

=head2 Required arguments:

    --conf <Config file>
    
    --data=  <Gencode data file>

    --output <Annotation output file>

=head2 Optional arguments:

	--appris <Flag that enables the output for APPRIS (default: NONE)>

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
    

=head1 EXAMPLE

perl proteo.pl
	
	--conf=../conf/pipeline.ini
	
	--data=examples/ENSG00000198563.annot.gtf
	
	--output=examples/ENSG00000198563.output


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
