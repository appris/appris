#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Utils::File qw( printStringIntoFile getTotalStringFromFile );
use APPRIS::Utils::Logger;

use lib "$FindBin::Bin/lib";
use common qw( get_main_report get_label_report );

###################
# Global variable #
###################

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($input_old_main_file) = undef;
my ($input_old_seq_file) = undef;
my ($input_main_file) = undef;
my ($input_label_file) = undef;
my ($input_seq_file) = undef;
my ($input_data_file) = undef;
my ($output_file) = undef;
my ($output_ens_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'input-main=s' 			=> \$input_main_file,
	'input-label=s'			=> \$input_label_file,
	'input-seq=s'    		=> \$input_seq_file,
	'input-data=s'  	  	=> \$input_data_file,
	'input-old-main=s' 		=> \$input_old_main_file,
	'input-old-seq=s'    	=> \$input_old_seq_file,
	'output=s'				=> \$output_file,
	'output-ens=s'			=> \$output_ens_file,
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);
unless ( defined $input_main_file and defined $input_label_file and defined $input_seq_file and defined $output_file and defined $output_ens_file )
{
	print `perldoc $0`;
	exit 1;
}

# Optional arguments

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
	# Get data backup from file
	my ($old_main_report, $old_seq_report) = (undef,undef);
	if ( defined $input_old_main_file and defined $input_old_seq_file ) {
		$logger->debug("-- get old main data that contains CCDS -------\n");	
		($old_main_report, $old_seq_report) = common::get_main_report($input_old_main_file, $input_old_seq_file, undef);			
	}
	#$logger->debug("OLD_MAIN_REPORT:\n".Dumper($old_main_report)."\n");
	#$logger->debug("OLD_SEQ_REPORT:\n".Dumper($old_seq_report)."\n");
	
	# Get data from file
	$logger->debug("-- get main data from files -------\n");
	my ($main_report, $seq_report) = common::get_main_report($input_main_file, $input_seq_file, $input_data_file);	
	$logger->debug("MAIN_REPORT:\n".Dumper($main_report)."\n");
	$logger->debug("SEQ_REPORT:\n".Dumper($seq_report)."\n");
	
	# Get label from file
	$logger->debug("-- get label data from files -------\n");
	my ($label_report) = common::get_label_report($input_label_file, $input_seq_file, $input_data_file);	
	$logger->debug("LABEL_REPORT:\n".Dumper($label_report)."\n");
			
	# Get annots
	$logger->debug("-- get anntos -------\n");
	eval {
		my ($cmd) = "awk -F \"\\t\" 'NF >= 20 {if (\$20 ~ /^PRINCIPAL:[0-9,M]+\$/ || \$20 ~ /^ALTERNATIVE:[0-9,M]+\$/ ) {print \$2\"\t\"\$1\"\t\"\$3\"\t\"\$8\"\t\"\$20}}' $input_main_file > $output_file";
		$logger->debug("** script: $cmd\n");
		system($cmd);		
	};
	die ("ERROR: getting principal annots") if ($@);
		
	$logger->debug("-- get anntos (for Ensembl) -------\n");
	eval {
		my ($cmd) = "awk -F \"\\t\" 'NF >= 20 {if (\$20 ~ /^PRINCIPAL:[0-9,M]+\$/ || \$20 ~ /^ALTERNATIVE:[0-9,M]+\$/ ) {print \$1\"\t\"\$3\"\t\"\$20}}' $input_main_file > $output_ens_file";
		$logger->debug("** script: $cmd\n");
		system($cmd);		
	};
	die ("ERROR: getting principal annots (for Ensembl)") if ($@);

	$logger->finish_log();
	
	exit 0;	
	
}

main();


1;

__END__

=head1 NAME

retrieve_main_data

=head1 DESCRIPTION

Get the main list of transcripts that have been labeled as main isoform

=head1 SYNOPSIS

retrieve_main_data

=head2 Required arguments:

	--input-main <Result file of APPRIS's scores>

	--input-label <Result file of APPRIS's labels>
	
	--input-seq <Input sequence comment file>

	--output <Output file that has the main isoforms>
	
	--output-ens <Output file that has the main isoforms FOR ESNSEMBL>	
	
=head2 Optional arguments (log arguments):

	--input-data=  <Gene annotation file>
	
	--input-old-main <Result file of APPRIS's scores>

	--input-old-seq <Input sequence comment file>

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_main_data.pl
	
		--input-main=data/appris_data.appris.txt
		
		--input-label=data/appris_data.appris_label.txt

		--input-seq=features/gencode.v7.pc_translations.fa

		--input-data=features/gencode.v7.annotation.gtf
		
		--input-old-main=data_gen20/appris_data.appris.txt
		
		--input-old-seq=features_gen20/appris_data.transl.fa		

		--output=data/appris_data.principal.txt
		
		--output-ens=data/appris_data.principal.ForENSEMBL.txt

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
