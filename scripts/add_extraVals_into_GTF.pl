#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( printStringIntoFile getTotalStringFromFile );
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
my ($tsl_data_file) = undef;
my ($out_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'data=s'			=> \$data_file,
	'old-data=s'		=> \$old_data_file,
	'tsl=s'				=> \$tsl_data_file,
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

	$logger->info("-- get gene report from older GTF file -------\n");
	my ($old_genedata) = APPRIS::Parser::_parse_indata($old_data_file);
	$logger->debug(Dumper($old_genedata)."\n");
	
	$logger->info("-- extract TSL list -------\n");
	my ($tsl_list) = get_tsl_annots($tsl_data_file);
	$logger->debug(Dumper($tsl_list)."\n");

	$logger->info("-- add extra values into GTF -------\n");	
	open (IN_FILE, $data_file) or throw('Can not open file');
	while ( my $line = <IN_FILE> ) {
		my ($new_values) = "";
		#ignore header		
		unless ( $line =~ /^#/ ) {
			my ($fields) = APPRIS::Parser::_parse_dataline($line);
			if ( defined $fields ) {
				if ( ($fields->{'type'} eq 'transcript') ) {
					my ($gene_id) = $fields->{'attrs'}->{'gene_id'};
					my ($transc_id) = $fields->{'attrs'}->{'transcript_id'};
					$gene_id =~ s/\.[0-9]*$//g;
					$transc_id =~ s/\.[0-9]*$//g;
					
					# add CCDS
					unless ( $line =~ /ccds_id/ or $line =~ /ccdsid/ ) {
						if ( exists $old_genedata->{$gene_id} and exists $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id} and exists $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'} ) {
							my ($ccds_id) = $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'};
							$new_values .= "; ccds_id \"$ccds_id\"";
						}
					}					
					# add RT
					unless ( $line =~ /"readthrough_transcript"/ ) {
						if ( exists $old_genedata->{$gene_id} and exists $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id} and exists $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'tag'} ) {
							my ($old_tags) = $old_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'tag'};
							foreach my $old_tag ( @{$old_tags} ) {
								if ( $old_tag eq 'readthrough_transcript' ) {
									$new_values .= '; tag "readthrough_transcript"';									
								}
								
							}
						}						
					}
					# add TSL
					unless ( $line =~ /"tsl[0-9]"/ ) {
						if ( exists $tsl_list->{$transc_id} ) {
							my ($tsl) = $tsl_list->{$transc_id};
							$new_values .= "'; tag \"$tsl\"";
						}
					}
					
				}				
			}
		}		
		$line =~ s/\;*\n*$//g;
		$output .= $line . $new_values . "\n";
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

sub get_tsl_annots($)
{
	my ($file) = @_;
	my ($report);
	my ($flines) = getTotalStringFromFile($file);
	foreach my $line (@{$flines}) {
			my (@txt_cols) = split('\t', $line);
			my ($g_id) = $txt_cols[0];
			my ($t_id) = $txt_cols[1];
			my ($tsl_annot) = $txt_cols[2];
			my ($appris_annot) = $txt_cols[3];
			if ( $tsl_annot =~ /(tsl[0-9])/ ) {
				$report->{$t_id} = $1;
			}			
	}
	return $report;
}

main();


1;

__END__

=head1 NAME

add_extraVals_into_GTF

=head1 DESCRIPTION

Script that add the extra values into GTF dataset. Values as readthrought tags, CCDS ids, TSL id.
 
=head1 SYNOPSIS

retrieve_exon_data

=head2 Required arguments:

	--data=  <Current Gene annotation file>
	
	--old-data=  <Older Gene annotation file with CCDS, RT values>
	
	--tsl=  <file with TSL annotations>
	
	--outfile <Output Gene annotation file>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl add_extraVals_into_GTF.pl
	
	--data={ABSOLUTE_PATH}/Homo_sapiens.GRCh38.84.gtf
	
	--old-data={ABSOLUTE_PATH}/gencode.v24.annotation.gtf
	
	--tsl=appris/db/TSL.annots.e80_g22_gM5.txt

	--outfile={ABSOLUTE_PATH}/Homo_sapiens.GRCh38.83.TSL-CCDS-RT.gtf


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
