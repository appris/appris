#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( printStringIntoFile getTotalStringFromFile );
use APPRIS::Utils::Logger;
use lib "/local/jmrodriguez/appris/scripts/lib";
use common;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$GENCODE_TSL_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;

# Input parameters
my ($is_refseq) = undef;
my ($xref_file) = undef;
my ($data_file) = undef;
my ($seq_file) = undef;
my ($extra_data_file) = undef;
my ($out_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'ref|r'				=> \$is_refseq,
	'xref=s'			=> \$xref_file,
	'data=s'			=> \$data_file,
	'seq-data=s'		=> \$seq_file,
	'extra-data=s'		=> \$extra_data_file,
	'outfile=s'			=> \$out_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $data_file and defined $extra_data_file and defined $out_file )
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

	$logger->info("-- get gene report from extra GTF file -------\n");
	my ($extra_genedata) = APPRIS::Parser::_parse_indata($extra_data_file);
	$logger->debug(Dumper($extra_genedata)."\n");
	
	my ($xref_data);
	my $_extract_dbXref;
	if ( defined $is_refseq ) {
		$logger->info("-- create xref report between RefSeq and Ensembl -------\n");
		$xref_data = create_xref_ensembl($xref_file);
		$logger->debug(Dumper($xref_data)."\n");
		
		$_extract_dbXref = sub {
			my ($data,$patt) = @_;
			my ($match);
			if ( $data =~ /$patt\:([^\,]*)/ ) {
				$match = $1;
			}
			return $match;
		};
		
		#$logger->info("-- export CCDS for RefSeq -------\n");
		#my ($seq_report) = common::get_seq_report($seq_file);
		
	}	

	$logger->info("-- add extra values into GTF -------\n");	
	open (IN_FILE, $data_file) or throw('Can not open file');
	while ( my $line = <IN_FILE> ) {
		my ($new_values) = "";
		#ignore header		
		unless ( $line =~ /^#/ ) {
			my ($fields) = APPRIS::Parser::_parse_dataline($line);
			if ( defined $fields ) {
				
				# for RefSeq dataset
				if ( defined $is_refseq ) {
					my ($source) = $fields->{'source'};
					my ($type) = $fields->{'type'};
					my ($attribs) = $fields->{'attrs'};
					next unless ( ($source =~ /BestRefSeq/) or ($source =~ /Curated Genomic/) or ($source =~ /Gnomon/) );			
					if (exists $attribs->{'ID'} ) {
						my ($gene_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'GeneID');
						my ($accesion_id) = $_extract_dbXref->($attribs->{'Dbxref'}, 'Genbank');			
						# Any type of molecule: mRNA, transcript, ncRNA, etc.
						if (defined $gene_id and defined $accesion_id and ( $type eq 'transcript' or $type eq 'mRNA' or $type eq 'ncRNA' or $type eq 'primary_transcript' or $type eq 'rRNA' or $type eq 'tRNA') ) {
							$accesion_id =~ s/\.[0-9]*$//g;
							# exits cross-refenence with Ensembl
							if ( exists $xref_data->{$accesion_id} ) {
								my ($ens_gene_id) = $xref_data->{$accesion_id}->{'gene_id'};
								my ($ens_transc_id) = $xref_data->{$accesion_id}->{'transc_id'};
								
								# add CCDS
								unless ( $line =~ /ccds_id/ or $line =~ /ccdsid/ ) {
									if ( exists $extra_genedata->{$ens_gene_id} and exists $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id} and exists $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id}->{'tsl'} ) {
										my ($ccds_id) = $extra_genedata->{$gene_id}->{'transcripts'}->{$ens_transc_id}->{'ccdsid'};
										$new_values .= "; ccds_id \"$ccds_id\"";
									}
								}								
								# add RT
								if ( exists $extra_genedata->{$ens_gene_id} and exists $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id} and exists $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id}->{'tag'} ) {
									my ($extra_tags) = $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id}->{'tag'};
									if ( $extra_tags =~ /readthrough_transcript/ ) {
										$new_values .= ";tag=readthrough_transcript";									
									}							
								}
								# add TSL
								if ( exists $extra_genedata->{$ens_gene_id} and exists $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id} and exists $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id}->{'tsl'} ) {
									my ($extra_val) = $extra_genedata->{$ens_gene_id}->{'transcripts'}->{$ens_transc_id}->{'tsl'};
									if ( $extra_val ne '' ) {
										$new_values .= ";tsl=$extra_val";
									}
								}
							}
						}
					}
				}
				# for Ensembl/GENCODE dataset
				else {
					if ( ($fields->{'type'} eq 'transcript') ) {
						my ($gene_id) = $fields->{'attrs'}->{'gene_id'};
						my ($transc_id) = $fields->{'attrs'}->{'transcript_id'};
						$gene_id =~ s/\.[0-9]*$//g;
						$transc_id =~ s/\.[0-9]*$//g;
						
						# add CCDS
						unless ( $line =~ /ccds_id/ or $line =~ /ccdsid/ ) {
							if ( exists $extra_genedata->{$gene_id} and exists $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id} and exists $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'} ) {
								my ($ccds_id) = $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'ccdsid'};
								$new_values .= "; ccds_id \"$ccds_id\"";
							}
						}					
						# add RT
						unless ( $line =~ /"readthrough_transcript"/ ) {
							if ( exists $extra_genedata->{$gene_id} and exists $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id} and exists $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'tag'} ) {
								my ($extra_tags) = $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'tag'};
								if ( $extra_tags =~ /readthrough_transcript/ ) {
									$new_values .= '; tag "readthrough_transcript"';									
								}							
							}				
						}					
						# add TSL
						unless ( $line =~ /transcript_support_level/ or $line =~ /\;\s*tsl/ ) {
							if ( exists $extra_genedata->{$gene_id} and exists $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id} and exists $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'tsl'} ) {
								my ($extra_val) = $extra_genedata->{$gene_id}->{'transcripts'}->{$transc_id}->{'tsl'};
								if ( $extra_val ne '' ) {
									$new_values .= "; tsl \"$extra_val\"";									
								}
							}				
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
			$report->{$ref_transc_id} = {
				'gene_id'	=> $gene_id,
				'transc_id'	=> $transc_id
			};
			#$report->{$ref_transc_id} = $transc_id;
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

	--ref  <Flag that indicated RefSeq dataset>
	
	--xref=  <Xref file between RefSeq and Ensembl>
	
	--data=  <Current Gene annotation file>
	
	--extra-data=  <Extra Gene annotation file with CCDS, RT values>
	
	--outfile <Output Gene annotation file>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

=head2 GENCODE/Ensembl example

perl add_extraVals_into_GTF.pl
	
	--data=e84_g24/Homo_sapiens.GRCh38.84.gtf
	
	--extra-data=e83_g24/gencode.v23.annotation.gtf
	
	--outfile=e84_g24/Homo_sapiens.GRCh38.84.extra.gtf

=head2 RefSeq example

perl add_extraVals_into_GTF.pl
	
	-r
	
	--xref=rs107/xref.rs107_ensembl.txt
	
	--data=rs107/ref_GRCh38.p2_top_level.gff3
	
	--seq-data=rs107/protein.fa

	--extra-data=e84_g24/gencode.v24.annotation.gtf
	
	--outfile=rs107/ref_GRCh38.p2_top_level.extra.gff3
	
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
