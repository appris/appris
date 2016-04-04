#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( getTotalStringFromFile printStringIntoFile );
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
my ($gen_datafile) = undef;
my ($rs_datafile) = undef;
my ($xref_file) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'gdata=s'			=> \$gen_datafile,
	'rdata=s'			=> \$rs_datafile,
	'xref=s'			=> \$xref_file,
	'outfile=s'			=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $gen_datafile and defined $rs_datafile and defined $xref_file and defined $outfile )
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

	$logger->debug("-- get gene report from GENCODE GTF file -------\n");
	my ($gen_data) = APPRIS::Parser::_parse_indata($gen_datafile);
	$logger->debug(Dumper($gen_data)."\n");
	
	$logger->debug("-- extract readthoutgh list -------\n");
	my ($rt_list) = create_rt($gen_data);
	$logger->debug(Dumper($rt_list)."\n");

	$logger->debug("-- create xref report between RefSeq and Ensembl -------\n");
	my ($xref_data) = create_xref_ensembl($xref_file);
	$logger->debug(Dumper($xref_data)."\n");
	
	# Add CCDS if not already exits
	$logger->debug("-- add CCDS into current GTF file if not already exits -------\n");
	my $_extract_dbXref = sub {
		my ($data,$patt) = @_;
		my ($match);
		if ( $data =~ /$patt\:([^\,]*)/ ) {
			$match = $1;
		}
		return $match;
	};
	open (IN_FILE, $rs_datafile) or throw('Can not open file');
	while ( my $line = <IN_FILE> ) {
		my ($add_txt) = "";
		#ignore header		
		unless ( $line =~ /^#/ ) {
			my ($fields) = APPRIS::Parser::_parse_dataline($line);			
			my ($source) = $fields->{'source'};
			my ($type) = $fields->{'type'};
			my ($attribs) = $fields->{'attrs'};
			next unless ( ($source =~ /BestRefSeq/) or ($source =~ /Curated Genomic/) or ($source =~ /Gnomon/) );			
			if (exists $attribs->{'ID'} ) {
				my ($gene_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'GeneID');
				my ($accesion_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'Genbank');				
				if (defined $gene_id and defined $accesion_id ) # Any type of molecule: mRNA, transcript, ncRNA, etc.
				{
					$accesion_id =~ s/\.[0-9]*$//g;
					if ( exists $xref_data->{$accesion_id} and exists $rt_list->{ $xref_data->{$accesion_id} }) { # exits cross-refenence with Ensembl and its RT
						if ( $type eq 'transcript' or $type eq 'mRNA' or $type eq 'ncRNA' or $type eq 'primary_transcript' or $type eq 'rRNA' or $type eq 'tRNA' ) {
							$add_txt = ";tag=readthrough_transcript";
						}
					}
				}
								
			}
		}		
		$line =~ s/\n*$//g;
		$output .= $line . $add_txt . "\n";
	}
	
	# Print output
	if ($output ne '') {
		$logger->debug("-- print output -------\n");
		my ($printing_file_log) = printStringIntoFile($output, $outfile);
		$logger->error("printing") unless(defined $printing_file_log);		
	}
	
	$logger->finish_log();
	
	exit 0;	
}

sub create_rt($)
{
	my ($data) = @_;
	my ($report);
	while ( my ($gene_id, $gene_features) = each(%{$data}) )
	{
		while (my ($transcript_id, $transcript_features) = each(%{$gene_features->{'transcripts'}}) )
		{
			if ( exists $transcript_features->{'tag'} ) {
				if (ref($transcript_features->{'tag'}) eq 'ARRAY') {
					foreach my $tag (@{$transcript_features->{'tag'}}) {
						if ( $tag eq 'readthrough_transcript' ) {
							$report->{$transcript_id} = 1;
							last;						
						}
					}
				}
				else {
					if ( $transcript_features->{'tag'} eq 'readthrough_transcript' ) {
						$report->{$transcript_id} = 1;
					}
				}
			}
		}
	}
	return $report;
}

sub create_xref_ensembl($)
{
	my ($file) = @_;
	my ($report);
	
	#9606    5893    ENSG00000002016 NM_001297420.1  ENST00000545564 NP_001284349.1  ENSP00000440268
	#9606    5893    ENSG00000002016 NM_001297422.1  ENST00000536177 NP_001284351.1  ENSP00000440486
	#9606    5893    ENSG00000002016 NM_134424.3     ENST00000358495 NP_602296.2     ENSP00000351284
	my ($flines) = getTotalStringFromFile($file);
	for (my $i=0; $i < scalar(@{$flines}); $i++ ) {
		my (@cols) = split('\t', $flines->[$i]);
		my ($tax_id) = $cols[0];
		my ($entrez_id) = $cols[1];
		my ($gene_id) = $cols[2];
		my ($ref_transc_id) = $cols[3];
		my ($transc_id) = $cols[4];
		$entrez_id =~ s/\s*//g;
		$ref_transc_id =~ s/\s*//g;
		
		if ( defined $gene_id and defined $transc_id and defined $entrez_id and defined $ref_transc_id and $gene_id =~ /^ENS/ and $transc_id =~ /^ENS/ and $entrez_id ne '' and $ref_transc_id ne '' ) {
			$ref_transc_id =~ s/\.[0-9]*$//g;
			#$report->{$gene_id}->{$transc_id}->{$ref_transc_id} = 1;
			$report->{$ref_transc_id} = $transc_id;
		}
	}
	return $report;
}

main();


1;

__END__

=head1 NAME

add_RT_into_GFF3

=head1 DESCRIPTION

Script that add the ReadThrought into RefSeq GFF3 from Xref in Ensembl
 
=head1 SYNOPSIS

add_RT_into_GFF3

=head2 Required arguments:

	--data=  <GENCODE Gene annotation file>
	
	--xref=  <Xref file between RefSeq and Ensembl>
	
	--outfile <Output Gene annotation file>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl add_RT_into_GFF3.pl
	
	--gdata=e81_g23/gencode.v23.annotation.gtf
	
	--rdata=rs107/ref_GRCh38.p2_top_level.gff3
	
	--xref=rs107/xref.rs107_ensembl.txt

	--outfile=rs107/ref_GRCh38.p2_top_level.RT.gff3


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
