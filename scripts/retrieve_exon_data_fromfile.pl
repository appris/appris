#!/usr/bin/perl

use strict;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;

use lib $ENV{'APPRIS_HOME'}.'/scripts/lib';
use common qw( get_main_report get_label_report );
use appris qw( create_indata );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$GENCODE_TSL_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$GENCODE_TSL_FILE			= '/Users/jmrodriguez/projects/APPRIS/data/homo_sapiens/ens76.v7.9Feb2015/gencode.tsl.e78.txt';

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($data_file) = undef;
my ($transcripts_file) = undef;
my ($translations_file) = undef;
my ($input_main_file) = undef;
my ($input_label_file) = undef;
my ($input_seq_file) = undef;
my ($input_data_file) = undef;
my ($outpath) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'data=s'			=> \$data_file,
	'transc=s'			=> \$transcripts_file,
	'transl=s'			=> \$translations_file,
	'input-label=s'		=> \$input_label_file,
	'outpath=s'			=> \$outpath,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $data_file and defined $transcripts_file and defined $translations_file and defined $input_label_file and defined $outpath )
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
sub extract_gencode_tsl($);
sub get_exon_data($$);

# extract GENCODE/TSL info
my ($GEN_TSL_REPORT) = extract_gencode_tsl($GENCODE_TSL_FILE);
#$logger->debug("GENCODE_TSL_REPORT:\n".Dumper($GEN_TSL_REPORT)."\n");

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->info("-- get gene dataset from files -------\n");
	my ($data_report) = appris::create_indata($data_file, $transcripts_file, $translations_file);
	$logger->debug("DATA_REPORT:\n".Dumper($data_report)."\n");

	# Get data from file
	$logger->info("-- get label data from files -------\n");
	my ($label_report) = common::get_label_report($input_label_file, $translations_file, $data_file);	
	$logger->debug("MAIN_REPORT:\n".Dumper($label_report)."\n");
	
	# Get data by region
	$logger->info("-- get exon data -------\n");
	my ($outfile) = $outpath.'/appris_data.exons.gff';
	my ($output_content) = get_exon_data($data_report, $label_report);	
	if ($output_content ne '') {
		my ($printing_file_log) = printStringIntoFile($output_content, $outfile);
		$logger->error("printing") unless(defined $printing_file_log);		
	}
		
	# get prin transc
	$logger->info("-- get principal exons -------\n");
	my ($outfile_prin) = $outpath.'/appris_data.exons.prin.gff';
	eval {
		my ($cmd) = "grep 'appris_annot \"PRINCIPAL' $outfile > $outfile_prin";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	throw("getting principal exons") if($@);
	
	# get alternative transc
	$logger->info("-- get alternative exons -------\n");
	my ($outfile_alt) = $outpath.'/appris_data.exons.alt.gff';
	eval {
		my ($cmd) = "grep 'appris_annot \"ALTERNATIVE' $outfile > $outfile_alt";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	throw("getting alternative exons") if($@);
	
	# get the specific principal exons (that do not overlap with alternative exons)
	$logger->info("-- get the specific principal exons -------\n");
	my ($outfile_PI) = $outpath.'/appris_data.exons.PI.gff';
	eval {
		my ($cmd) = "subtractBed -a $outfile_prin -b $outfile_alt > $outfile_PI";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	throw("getting prin exons that do not overlap") if($@);
		
	# get the specific alternative exons (that do not overlap with principal exons)
	$logger->info("-- get the specific alternative exons -------\n");
	my ($outfile_NPI) = $outpath.'/appris_data.exons.NPI.gff';
	eval {
		my ($cmd) = "subtractBed -a $outfile_alt -b $outfile_prin > $outfile_NPI";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	throw("getting alt exons that do not overlap") if($@);

	# get the overlap regions PI+NPI
	$logger->info("-- get the overlapping regions PI+NPI -------\n");
	my ($outfile_OVER) = $outpath.'/appris_data.exons.OVER.gff';
	eval {
		my ($cmd) = "intersectBed -a $outfile_prin -b $outfile_alt -wb > $outfile_OVER";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	throw("getting over regions") if($@);
	
	
#	# Add start/stop codons for alt exons
#	$logger->info("-- add start/stop codons for alt exons -------\n");
#	my ($outcont_NPI_wCodons) = add_codons($data_report, $outfile_NPI);
#	my ($outfile_NPI_wCodons) = $outpath.'/appris_data.exons.NPIwithCodons.gff';
#	if ($outcont_NPI_wCodons ne '') {
#		my ($printing_file_log) = printStringIntoFile($outcont_NPI_wCodons, $outfile_NPI_wCodons);
#		$logger->error("printing") unless(defined $printing_file_log);		
#	}
		
	
#	# get the specific principal exons (that do not overlap with alternative exons)
#	$logger->info("-- get the specific principal exons -------\n");
#	my ($outfile_PI) = $outpath.'/appris_data.exons.PI.gff';
#	eval {
#		my ($cmd) = "intersectBed -v -a $outfile_prin -b $outfile_alt > $outfile_PI";
#		$logger->debug("\n** script: $cmd\n");
#		system ($cmd);
#	};
#	throw("getting prin exons that do not overlap") if($@);
#		
#	# get the specific alternative exons (that do not overlap with principal exons)
#	$logger->info("-- get the specific alternative exons -------\n");
#	my ($outfile_NPI) = $outpath.'/appris_data.exons.NPI.gff';
#	eval {
#		my ($cmd) = "intersectBed -v -a $outfile_alt -b $outfile_prin > $outfile_NPI";
#		$logger->debug("\n** script: $cmd\n");
#		system ($cmd);
#	};
#	throw("getting alt exons that do not overlap") if($@);


#	$logger->info("-- merge over regions -------\n");
#	my ($outfile_OVER) = $outpath.'/appris_data.exons.OVER.gff';
#	eval {
#		#my ($cmd) = "bedtools merge -i $outfile_aOVER -c 4,5,6,7,8,9 -o collapse -delim '|' > $outfile_OVER ";
#		my ($cmd) = "bedtools merge -i $outfile_aOVER -c 9 -o collapse -delim '|' > $outfile_OVER ";
#		$logger->debug("\n** script: $cmd\n");
#		system ($cmd);
#	};
#	throw("mergging over regions") if($@);
	
	
#	# get prin transc
#	$logger->info("-- get count bp of principal exons -------\n");
#	eval {
#		my ($cmd) = "bedtools coverage -a $outfile_prin -b $outfile_prin | cut -f 12 | awk '{s+=$1} END {printf \"%.0f\", s}'";
#		$logger->debug("\n** script: $cmd\n");
#		system ($cmd);
#	};
#	throw("getting principal exons") if($@);
#		
#	# get prin transc
#	$logger->info("-- get count bp of alternative exons -------\n");
#	eval {
#		my ($cmd) = "bedtools coverage -a $outfile_alt -b $outfile_alt | cut -f 12 | awk '{s+=$1} END {printf \"%.0f\", s}'";
#		$logger->debug("\n** script: $cmd\n");
#		system ($cmd);
#	};
#	throw("getting principal exons") if($@);

	
	$logger->finish_log();
	
	exit 0;	
}

sub count_bp() {
	my ($file) = @_;
	
	
}
sub get_exon_data($$)
{
	my ($dataset, $labelset) = @_;
	#my ($c_report);
	my ($output) = '';
	
	foreach my $gene (@{$dataset}) {
		my ($gene_id) = $gene->stable_id;
		my ($gene_name) = $gene->external_name;
		my ($chr) = $gene->chromosome;

		#my ($g_report);		
		my ($e_report);
		my ($num_transc) = 0;
		my ($g_appris_annot) = 'REJECTED';
		my ($g_appris_transc_list) = '';
		
		foreach my $transcript (@{$gene->transcripts}) {
			my ($t_report);
			my ($transcript_id) = $transcript->stable_id;
			my ($transcript_name) = $transcript->external_name;

			if ($transcript->translate and $transcript->translate->cds) {
				my ($translate) = $transcript->translate;
				#my ($cds) = $translate->cds;			
				my ($ccds_id) = '-';				
				if ( $transcript->xref_identify ) {
					foreach my $xref_identify (@{$transcript->xref_identify}) {								
						if ($xref_identify->dbname eq 'CCDS') {
							$ccds_id = $xref_identify->id;
						}
					}					
				}
				my ($t_biotype) = $transcript->biotype;
				# get appris annotation
				#my ($a_analysis) = $registry->fetch_analysis_by_stable_id($transcript_id,'appris');
				my ($a_analysis) = $labelset->{$gene_id}->{'transcripts'}->{$transcript_id}->{'annotations'};
				if ( defined $a_analysis and (scalar(@{$a_analysis}) >= 10) ) {
					#my ($appris_annot) = $a_analysis->[9];
					#my ($a_annot);
					#if ( $appris_annot eq 'YES' ) {
					#	$a_annot = 'PRINCIPAL';
					#	$g_appris_annot = 'PRINCIPAL';
					#	$g_appris_transc_list .= $transcript_id.',';
					#}
					#elsif ( $appris_annot eq 'UNKNOWN' ) {
					#	$a_annot = 'ALTERNATIVE';
					#	$g_appris_annot = 'ALTERNATIVE';
					#	$g_appris_transc_list .= $transcript_id.',';	
					#}
					#elsif ( $appris_annot eq 'NO' ) {
					#	$a_annot = 'REJECTED';					
					#}					
					my ($appris_annot) = $a_analysis->[10];
					my ($a_annot);
					if ( $appris_annot =~ /PRINCIPAL/ ) {
						$a_annot = $appris_annot;
						$g_appris_annot = $appris_annot;
						$g_appris_transc_list .= $transcript_id.',';
					}
					elsif ( $appris_annot =~ /ALTERNATIVE/ ) {
						$a_annot = $appris_annot;
						#$g_appris_annot = $appris_annot;
						#$g_appris_transc_list .= $transcript_id.',';
					}
					elsif ( $appris_annot eq 'MINOR' ) {
						$a_annot = $appris_annot;					
					}					
					if ( defined $a_annot ) {
						
						# get isoform annotation (appris)
						$t_report->{'name'}				= $transcript_name;
						$t_report->{'ccds_id'}			= $ccds_id;
						$t_report->{'biotype'}			= $t_biotype;
						$t_report->{'appris_annot'}		= $a_annot;
						#$g_report->{'cds'}				= $cds;
						
						# get specie conservation (corsair)
						#my ($c_analysis) = $registry->fetch_analysis_by_stable_id($transcript_id,'corsair');						
						#if ( $c_analysis and $c_analysis->corsair and defined $c_analysis->corsair->score ) {
						#	$t_report->{'corsair_score'} = $c_analysis->corsair->score;								
						#}
						my ($corsair_annot) = $a_analysis->[2];;
						my ($c_annot);
						if ( ($corsair_annot eq 'YES') or ($corsair_annot eq 'UNKNOWN') ) {
							$c_annot = 'CONSERVE';						
						}
						elsif ( $corsair_annot eq 'NO' ) {
							$c_annot = 'NO_CONSERVE';
						}
						if ( defined $c_annot ) {
							$t_report->{'corsair_annot'} = $c_annot;
						}						
						
						# cds annotation
						for (my $icds = 0; $icds < scalar(@{$translate->cds}); $icds++) {
							my ($cds) = $translate->cds->[$icds];	
							my ($exon_id) = $cds->stable_id;
							my ($cds_start) = $cds->start;
							my ($cds_end) = $cds->end;
							my ($cds_strand) = $cds->strand;
							my ($cds_phase) = $cds->phase;
							my ($cds_index) = $cds_start.'-'.$cds_end.':'.$cds_strand.':'.$cds_phase;
							
							if ( defined $exon_id ) {
								$e_report->{$cds_index}->{'exon_id'} = $exon_id;
							}
							if ( exists $t_report->{'appris_annot'} and defined $t_report->{'appris_annot'} ) {
								$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'appris_annot'} = $t_report->{'appris_annot'};
							}
							#if ( exists $t_report->{'corsair_score'} and defined $t_report->{'corsair_score'} ) {
							#	$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'corsair_score'} = $t_report->{'corsair_score'};
							#}							
							if ( exists $t_report->{'corsair_annot'} and defined $t_report->{'corsair_annot'} ) {
								$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'corsair_annot'} = $t_report->{'corsair_annot'};
							}
							if ( exists $t_report->{'biotype'} and defined $t_report->{'biotype'} ) {
								$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'biotype'} = $t_report->{'biotype'};
							}
						}							
																		
						# save trans report
						#$g_report->{'chr'}								= $chr;					
						#$g_report->{'transcripts'}->{$transcript_id}	= $t_report if ( defined $t_report );
						$num_transc++;
					}											
				}
			}
		}
		#$c_report->{$gene_id} = $g_report if ( defined $g_report );
		
		# print sorted exons per gene
		foreach my $cds_index (sort { $a cmp $b } keys %{$e_report} )
		{
			my ($cds_report) = $e_report->{$cds_index};
			my ($cds_start,$cds_end,$cds_strand,$cds_phase) = (undef,undef,undef,undef);
			if ( $cds_index =~ /^([^\-]*)\-([^\:]*)\:([^\:]*)\:([^\$]*)$/) {
				($cds_start,$cds_end,$cds_strand,$cds_phase) = ($1,$2,$3,$4);
			}
			my ($exon_id) = '-';
			if ( exists $cds_report->{'exon_id'} and defined $cds_report->{'exon_id'} ) {
				$exon_id = $cds_report->{'exon_id'};
			}
			my ($t_appris_annot) = '-';
			my ($c_annot) = '-';
			my ($transc_list) = '';
			my ($biotype_rep);
			my ($biotype_list) = '';
			my ($g_annot) = '-';
			my ($t_annot) = '-';
			while ( my ($transc_id, $trans_report) = each(%{$cds_report->{'trans'}}) ) {
				if ( $trans_report->{'appris_annot'} =~ /^PRINCIPAL/ ) {
					$t_appris_annot = $trans_report->{'appris_annot'};					
				}
				$transc_list .= $transc_id.',';
				#if ( exists $trans_report->{'corsair_score'} and defined $trans_report->{'corsair_score'} and ($trans_report->{'corsair_score'} > 0) ) {
				#	$c_annot = $trans_report->{'corsair_score'};
				#}				
				if ( exists $trans_report->{'corsair_annot'} and defined $trans_report->{'corsair_annot'} ) {
						$c_annot = $trans_report->{'corsair_annot'};
				}				
				if ( exists $trans_report->{'biotype'} and defined $trans_report->{'biotype'} ) {
					my ($t) = $trans_report->{'biotype'};
					unless ( exists $biotype_rep->{$t} ) {
						$biotype_rep->{$t} = 1;
					}
				}				
				# extract GENCODE/TSL info
				if ( exists $GEN_TSL_REPORT->{$gene_id}->{$transc_id}->{'gencode'} ) {
					if ( $g_annot eq '-' ) {
						$g_annot = $GEN_TSL_REPORT->{$gene_id}->{$transc_id}->{'gencode'};
					}
				}
				if ( exists $GEN_TSL_REPORT->{$gene_id}->{$transc_id}->{'tsl'} ) {
					# get the best TSL
					if ( $t_annot eq '-' ) {
						$t_annot = $GEN_TSL_REPORT->{$gene_id}->{$transc_id}->{'tsl'};
					}
					else {
						my ($new) = $GEN_TSL_REPORT->{$gene_id}->{$transc_id}->{'tsl'};
						if ( ( ($new ne 'NA') and ($t_annot eq 'NA') ) or ( ($new ne 'NA') and ($t_annot ne 'NA') and ($new < $t_annot) ) ) {
							$t_annot = $GEN_TSL_REPORT->{$gene_id}->{$transc_id}->{'tsl'};
						}
					}
				}				
			}
			$transc_list =~ s/\,$//mg;
			$biotype_list = join(',', keys(%{$biotype_rep}) );

			$t_appris_annot = 'ALTERNATIVE' if ( $t_appris_annot eq '-');
			
			if ( $transc_list ne '' ) {
 				$output .=
 						$chr."\t".
 						'APPRIS'."\t".
 						'CDS'."\t".
						$cds_start."\t".
						$cds_end."\t".
						'.'."\t".
						$cds_strand."\t".
						$cds_phase."\t".
						'exon_id "'.$exon_id.'"'."; ".
						'gene_id "'.$gene_id.'"'."; ".
						'gene_name "'.$gene_name.'"'."; ".
						'transc_list "'.$transc_list.'"'."; ".
						'biotype_list "'.$biotype_list.'"'."; ".
						'appris_gene_annot "'.$g_appris_annot.'"'."; ".
						'appris_annot "'.$t_appris_annot.'"'."; ".
						'gencode "'.$g_annot.'"'."; ".
						'tsl "'.$t_annot.'"'."\n";
			}			
		}
	}
	return $output;
}

sub extract_gencode_tsl($) {
	my ($file) = @_;
	my ($report);
	
	local(*FILE);
	open(FILE,$file) or return undef;
	my(@string)=<FILE>;
	close(FILE);
	
	#Ensembl Gene ID Ensembl Transcript ID   GENCODE basic annotation        Transcript Support Level (TSL)
	#ENSG00000276385 ENST00000618935 GENCODE basic   tslNA
	#ENSG00000197468 ENST00000508957 GENCODE basic   tslNA
	#ENSG00000275151 ENST00000614589 GENCODE basic   tslNA
	#ENSG00000231049 ENST00000435337 GENCODE basic   tslNA
	#ENSG00000280296 ENST00000624531 GENCODE basic	
	for (my $i = 1; $i <= scalar(@string); $i++) {
		my ($line) = $string[$i];
		if ( defined $line and ($line ne '') ) {
			my (@cols) = split("\t", $line);
			my ($g_id) = $cols[0];
			my ($t_id) = $cols[1];
			my ($g_label) = $cols[2];
			my ($t_label) = $cols[3];
			$g_id =~ s/\s*//mg; $g_id =~ s/\.\d*$//;
			$t_id =~ s/\s*//mg; $t_id =~ s/\.\d*$//;			
			if ( ($g_label ne '') and ($g_label =~ /GENCODE basic/) ) {
				$report->{$g_id}->{$t_id}->{'gencode'} = 'GENCODE';				
			}
			if ( ($t_label ne '') and ($t_label =~ /tsl([^\s]*)/) ) {
				$report->{$g_id}->{$t_id}->{'tsl'} = $1;				
			}
		}
	}
	return $report;
}
#
#sub _add_codon($$) {
#	my ($gendata,$i_transc_id) = @_;
#	my ($output) = '';
#	
#	foreach my $gene (@{$gendata}) {
#		my ($gene_id) = $gene->stable_id;
#		my ($chr) = $gene->chromosome;		
#		if ( exists $gene->{'_index_transcripts'}->{$i_transc_id} ) {
#			my ($index) = $gene->{'_index_transcripts'}->{$i_transc_id};
#			my ($transcript) = $gene->transcripts->[$index];
#			if ( $transcript->translate and $transcript->translate->codons ) {
#				my ($translation) = $transcript->translate;
#				foreach my $codon (@{$translation->codons}) {
#					my ($type) = $codon->type.'_codon';
#	 				$output .=
#	 						$chr."\t".
#	 						'APPRIS'."\t".
#	 						$type."\t".
#							$codon->start."\t".
#							$codon->end."\t".
#							'.'."\t".
#							$codon->strand."\t".
#							$codon->phase."\t".
#							'gene_id "'.$gene_id.'"'."; ".
#							'transc_id"'.$i_transc_id.'"'."\n";
#				}
#			}
#			last;
#		}
#	}
#	return $output;
#}

sub get_codon_data($)
{
	my ($dataset) = @_;
	my ($report);	
	foreach my $gene (@{$dataset}) {
		foreach my $transcript (@{$gene->transcripts}) {
			my ($transcript_id) = $transcript->stable_id;
			if ($transcript->translate and $transcript->translate->codons) {
				my ($translate) = $transcript->translate;
				$report->{$transcript_id} = $translate->codons;
			}
		}
	}
	return $report;
}

sub add_codons($$) {
	my ($gendata,$file) = @_;
	my ($output) = '';
	my ($trans_codons_list);

	my ($codondata) = get_codon_data($gendata);
	
	local(*FILE);
	open(FILE,$file) or return undef;
	my(@string)=<FILE>;
	close(FILE);
	
	for (my $i = 0; $i <= scalar(@string); $i++) {
		my ($line) = $string[$i];
		my ($fields) = APPRIS::Parser::_parse_dataline($line);
		if ( defined $fields and defined $fields->{'attrs'} ) {
			if ( exists $fields->{'chr'} and defined $fields->{'chr'} and
				 exists $fields->{'attrs'}->{'gene_id'} and defined $fields->{'attrs'}->{'gene_id'} and
				 exists $fields->{'attrs'}->{'transc_list'} and defined $fields->{'attrs'}->{'transc_list'}
			) {
				my ($chr) = $fields->{'chr'};
				my ($gene_id) = $fields->{'attrs'}->{'gene_id'};
				# add codons if it is not exist yet
				foreach my $transc_id ( split(',', $fields->{'attrs'}->{'transc_list'}) ) {
					unless ( $trans_codons_list->{$transc_id} ) {
						if ( exists $codondata->{$transc_id} and defined $codondata->{$transc_id} ) {
							my ($codons) = $codondata->{$transc_id};
							foreach my $codon (@{$codons}) {
								my ($type) = $codon->type.'_codon';
				 				$output .=
				 						$chr."\t".
				 						'APPRIS'."\t".
				 						$type."\t".
										$codon->start."\t".
										$codon->end."\t".
										'.'."\t".
										$codon->strand."\t".
										$codon->phase."\t".
										'gene_id "'.$gene_id.'"'."; ".
										'transc_id"'.$transc_id.'"'."\n";
							}
							
						}
						$trans_codons_list->{$transc_id} = 1;
					}
				}
			}
			$output .= $line;
		}
	}
	return $output;
}

main();


1;

__END__

=head1 NAME

retrieve_exon_data_fromfile

=head1 DESCRIPTION

Get the sorted list of exons per gene.
Each exon is labeled which the following column:
	PRINCIPAL/POTENTIAL -> the exon belongs to principal isoform or the variant is possible principal isoform.
	CONSERVE/NO_CONSERVE -> the exon is conserve agains vertebrates.
	OVERLAP/NO_OVERLAP -> the exon is not in the whole variants labaled as principal isoform. 

=head1 SYNOPSIS

retrieve_exon_data

=head2 Required arguments:

	--data=  <Gene annotation file>
	
	--transc= <Transcript sequences file>
	
	--transl= <Translation sequences file>
	
	--input-label <Result file of APPRIS's labels>
	
	--outpath <Output path to save files>	
	
=head2 Optional arguments:

	--position= <Genome position> (chr21 or chr1:109102711-109187522)

	--apprisdb-conf <Config file of APPRIS database (default: 'conf/apprisdb.ini' file)>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl retrieve_exon_data_fromfile.pl
	
	--data={ABSOLUTE_PATH}/gencode.v20.annotation.gtf
	
	--transc={ABSOLUTE_PATH}/gencode.v20.pc_transcripts.fa
	
	--transl={ABSOLUTE_PATH}/gencode.v20.pc_translations.fa
	
	--input-label=data/appris_data.appris_label.txt

	--outpath=../data/
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
