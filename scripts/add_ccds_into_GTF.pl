#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$GENCODE_TSL_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;

# Input parameters
my ($data_file) = undef;
my ($old_data_file) = undef;
my ($out_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'data=s'			=> \$data_file,
	'old-data=s'		=> \$old_data_file,
	'outfile=s'			=> \$out_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $data_file and defined $old_data_file and defined $out_file )
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
$logger->init_log();

#####################
# Method prototypes #
#####################

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	my ($output) = ""; 

	# Get Gene Report from older GTF file
	$logger->debug("-- get gene report from older GTF file -------\n");
	my ($old_genedata) = APPRIS::Parser::_parse_indata($old_data_file);
#	$logger->debug(Dumper($old_genedata)."\n");
	
	
	# Add CCDS if not already exits
	$logger->debug("-- add CCDS into current GTF file if not already exits -------\n");	
	open (IN_FILE, $data_file) or throw('Can not open file');
	while ( my $line = <IN_FILE> ) {
		#ignore header		
		if ( $line =~ /^#/ ) { $output .= $line; next; }
		
		my ($fields) = APPRIS::Parser::_parse_dataline($line);
		unless ( defined $fields ) { $output .= $line; next; }
		if ( ($fields->{'type'} eq 'transcript') ) {
			unless ( $line =~ /ccds_id/ or $line =~ /ccdsid/ ) {
				my ($gene_id) = $fields->{'attrs'}->{'gene_id'};
				my ($transc_id) = $fields->{'attrs'}->{'transcript_id'};
				if ( exists $old_genedata->{$gene_id} and exists $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id} and exists $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'} ) {
					my ($ccds_id) = $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'};
					$line =~ s/\n*$//g;
					$output .= $line . " ccds_id \"$ccds_id\"\n";
				}
			}
			else { $output .= $line }
		}
		else { $output .= $line }
	}
	
	# Print output
	if ($output ne '') {
		$logger->debug("-- print output -------\n");
		my ($printing_file_log) = printStringIntoFile($output, $out_file);
		$logger->error("printing") unless(defined $printing_file_log);		
	}
	
	$logger->finish_log();
	
	exit 0;	
}

main();


1;

__END__

=head1 NAME

add_ccds_into_GTF

=head1 DESCRIPTION

Script that add the CCDS ids of older version when Ensembl give not us the updated ones.
 
=head1 SYNOPSIS

retrieve_exon_data

=head2 Required arguments:

	--data=  <Current Gene annotation file>
	
	--old-data=  <Older Gene annotation file>
	
	--outfile <Output Gene annotation file>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl add_ccds_into_GTF.pl
	
	--data={ABSOLUTE_PATH}/Homo_sapiens.GRCh38.83.gtf
	
	--old-data={ABSOLUTE_PATH}/Homo_sapiens.GRCh38.82.gtf

	--outfile={ABSOLUTE_PATH}/Homo_sapiens.GRCh38.83-fixed.gtf


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
