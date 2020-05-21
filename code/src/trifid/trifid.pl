#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Config::IniFiles;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile );


###################
# Global variable #
###################
use vars qw(
	$PROG_DB
);

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($data_file) = undef;
my ($output_file) = undef;
my ($loglevel) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;


&GetOptions(
	'conf=s'		=> \$config_file,
	'data=s'		=> \$data_file,
	'output=s'		=> \$output_file,
	'loglevel=s'	=> \$loglevel,
	'logfile=s'		=> \$logfile,
	'logpath=s'		=> \$logpath,
	'logappend'		=> \$logappend,
);

# Required arguments
unless ( defined $config_file and defined $data_file and defined $output_file )
{
	print `perldoc $0`;
	exit 1;
}

# Get conf vars
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$PROG_DB = $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('TRIFID_VARS', 'db');

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE	=> $logfile,
	-LOGPATH	=> $logpath,
	-LOGAPPEND	=> $logappend,
	-LOGLEVEL	=> $loglevel,
);
$logger->init_log($str_params);

# Main subroutine
sub main()
{
	# Declare vars
	my ($output_content) = "";

	# obtain the gene id
	my ($gene_id);
	eval
	{
		$logger->info("-- obtain gene_id\n");
		open(my $fh, $data_file) or $logger->error("failed to open file: '$data_file'\n");
		LINE: while (<$fh>) {
			chomp;
			my ($feature, $attr_field) = (split(/\s*\t\s*/))[2, 8];
			if ( $feature eq "gene" ) {
				my @attrs = split(/\s*;\s*/, $attr_field);
				foreach my $attr (@attrs) {
					if ( $attr =~ /^(?<tag>\S+)\s+"(?<value>.+?)"$/ ) {
						if ( $+{"tag"} eq "gene_id" ) {
							$gene_id = $+{"value"};
							last LINE;
						}
					}
				}
			}
		}
		close($fh) or $logger->error("failed to close file: '$data_file'\n");

		if ( ! defined($gene_id) ) {
			$logger->error("can not retrieve gene_id\n");
		}
	};
	$logger->error("obtaining gene_id: $!\n") if($@);

	# retrieve Trifid predictions from tabix-indexed tab file
	if ( defined($gene_id) ) {
		$logger->info("-- retrieve Trifid predictions for gene $gene_id\n");
		my ($cmd) = "tabix -h $PROG_DB $gene_id";
		$logger->debug("$cmd\n");

		my (@trifid_lines) = `$cmd`;
		if ( scalar(@trifid_lines) > 1 ) {  # first line is header
			$output_content = join("", @trifid_lines);
		} else {
			$logger->warning("no Trifid predictions found\n");
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

trifid

=head1 DESCRIPTION

Obtain Trifid results for the given gene

=head2 Required arguments:

	--conf=FILE <Config file>

	--data=FILE <Gencode data file>

	--output=FILE <Annotation output file>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>

	--logfile=FILE <Log to FILE (default: *STDOUT)>

	--logpath=PATH <Write logfile to PATH (default: .)>

	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl trifid.pl

	--conf=../conf/pipeline.ini

	--data=examples/ENSG00000254647/annot.gtf

	--output=examples/ENSG00000254647/trifid

=head1 AUTHORS

Original APPRIS method script by Jose Manuel Rodriguez Carrasco,
adapted for TRIFID by Thomas Walsh.

For contact details see the L<APPRIS website|http://appris-tools.org>.

=cut
