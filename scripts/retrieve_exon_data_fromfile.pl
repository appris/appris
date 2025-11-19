#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use File::Path qw( make_path );
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;

# the order of package it's important. Conflict of names
use lib "$FindBin::Bin/lib";
use appris;
use common;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($data_file) = undef;
my ($translations_file) = undef;
my ($input_main_file) = undef;
my ($outpath) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'data=s'			=> \$data_file,
	'transl=s'			=> \$translations_file,
	'main=s'			=> \$input_main_file,
	'outpath=s'			=> \$outpath,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $data_file and defined $translations_file and defined $input_main_file and defined $outpath )
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
sub get_exon_data($$);
sub add_codons($$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->info("-- get gene dataset from files -------\n");
	my ($data_report) = appris::create_indata($data_file, undef, $translations_file);
	$logger->debug("DATA_REPORT:\n".Dumper($data_report)."\n");

	# Get data from file
	$logger->info("-- get main data from files -------\n");
	my ($main_report) = common::get_main_report($input_main_file, $translations_file, $data_file);	
	$logger->debug("MAIN_REPORT:\n".Dumper($main_report)."\n");
	
	# Ensure output directory exists.
	$logger->info("-- ensure output directory exists -------\n");
	make_path($outpath);

	# Get data by region
	$logger->info("-- get exon data -------\n");
	my ($outfile) = $outpath.'/appris_data.exons.gff';
	my ($output_content) = get_exon_data($data_report, $main_report);	
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
	my ($outfile_OVER_txt) = $outpath.'/appris_data.exons.OVER.txt';
	eval {
		my ($cmd) = "intersectBed -a $outfile_prin -b $outfile_alt -wb > $outfile_OVER_txt";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	throw("getting over regions") if($@);
	my ($outfile_OVER_gff) = $outpath.'/appris_data.exons.OVER.gff';
	eval {
		my ($cmd) = "cut -f 1-9 $outfile_OVER_txt > $outfile_OVER_gff";
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
	
	$logger->finish_log();
	
	exit 0;	
}

sub get_exon_data($$)
{
	my ($dataset, $mainset) = @_;
	my ($output) = '';
	
	foreach my $gene (@{$dataset}) {
		my ($gene_id) = $gene->stable_id;
		my ($gene_name) = $gene->external_name;
		my ($chr) = $gene->chromosome;

		my ($e_report);
		my ($num_transc) = $mainset->{$gene_id}->{'num_trans'};
		my ($num_isof) = $mainset->{$gene_id}->{'num_isof'};
		my ($g_appris_annot) = 'REJECTED';
		my ($g_appris_transc_list) = '';
		
		if ( $gene->transcripts ) {
		foreach my $transcript (@{$gene->transcripts}) {
			my ($t_report);
			my ($transcript_id) = $transcript->stable_id;
			my ($transcript_name) = $transcript->external_name;

			if ($transcript->translate and $transcript->translate->cds) {
				my ($translate) = $transcript->translate;
				my ($ccds_id) = '-';				
				if ( $transcript->xref_identify ) {
					foreach my $xref_identify (@{$transcript->xref_identify}) {								
						if ($xref_identify->dbname eq 'CCDS') {
							$ccds_id = $xref_identify->id;
						}
					}					
				}
				my ($t_biotype) = $transcript->biotype;
				my ($t_tsl) = $transcript->tsl;
				# get appris annotation
				my ($a_analysis) = $mainset->{$gene_id}->{'transcripts'}->{$transcript_id}->{'annotations'};
				if ( defined $a_analysis and (scalar(@{$a_analysis}) >= 10) ) {
					my ($appris_annot) = $a_analysis->[10];
					my ($a_annot);
					if ( $appris_annot =~ /PRINCIPAL/ ) {
						$a_annot = $appris_annot;
						$g_appris_annot = $appris_annot;
						$g_appris_transc_list .= $transcript_id.',';
					}
					elsif ( $appris_annot =~ /ALTERNATIVE/ ) {
						$a_annot = $appris_annot;
					}
					elsif ( $appris_annot =~ /MINOR/ ) {
						$a_annot = $appris_annot;					
					}					
					if ( defined $a_annot ) {						
						# get isoform annotation (appris)
						$t_report->{'name'}				= $transcript_name;
						$t_report->{'ccds_id'}			= $ccds_id;
						$t_report->{'biotype'}			= $t_biotype;
						$t_report->{'tsl'}				= $t_tsl;
						$t_report->{'appris_annot'}		= $a_annot;
						
						# get specie conservation (corsair)
						my ($corsair_annot) = $a_analysis->[2];
						if ( defined $corsair_annot ) {
							$t_report->{'corsair_score'} = $corsair_annot;
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
							if ( exists $t_report->{'corsair_score'} and defined $t_report->{'corsair_score'} ) {
								$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'corsair_score'} = $t_report->{'corsair_score'};
							}
							if ( exists $t_report->{'biotype'} and defined $t_report->{'biotype'} ) {
								$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'biotype'} = $t_report->{'biotype'};
							}
							if ( exists $t_report->{'tsl'} and defined $t_report->{'tsl'} ) {
								$e_report->{$cds_index}->{'trans'}->{$transcript_id}->{'tsl'} = $t_report->{'tsl'};
							}
						}							
					}											
				}
			}
		}
		}
		
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
			my ($appris_type) = '-';
			my ($t_appris_annot) = '-';
			my ($c_annot) = '-';
			my ($transc_list) = '';
			my ($biotype_rep);
			my ($biotype_list) = '';
			my ($tsl_rep);
			my ($tsl_list) = '';
			while ( my ($transc_id, $trans_report) = each(%{$cds_report->{'trans'}}) ) {
				$transc_list .= $transc_id.',';
				if ( $trans_report->{'appris_annot'} =~ /^PRINCIPAL/ ) {
					$t_appris_annot = $trans_report->{'appris_annot'};
					if    ( $appris_type eq '-' ) { $appris_type = 'uniqPI' }
					elsif ( $appris_type ne 'uniqPI' ) { $appris_type = 'overPI-NPI' }
				}
				else {
					if    ( $appris_type eq '-' ) { $appris_type = 'uniqNPI' }
					elsif ( $appris_type ne 'uniqNPI' ) { $appris_type = 'overPI-NPI' }
				}
				if ( exists $trans_report->{'corsair_score'} and defined $trans_report->{'corsair_score'} and ($trans_report->{'corsair_score'} > 0) ) {
					$c_annot = $trans_report->{'corsair_score'};
				}
				if ( exists $trans_report->{'biotype'} and defined $trans_report->{'biotype'} ) {
					my ($t) = $trans_report->{'biotype'};
					unless ( exists $biotype_rep->{$t} ) {
						$biotype_rep->{$t} = 1;
					}
				}				
				if ( exists $trans_report->{'tsl'} and defined $trans_report->{'tsl'} ) {
					my ($t) = $trans_report->{'tsl'};
					unless ( exists $tsl_rep->{$t} ) {
						$tsl_rep->{$t} = 1;
					}
				}				
			}
			$transc_list =~ s/\,$//mg;
			$biotype_list = join(',', keys(%{$biotype_rep}) );
			$tsl_list = join(',', keys(%{$tsl_rep}) );

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
						'num_transc "'.$num_transc.'"'."; ".
						'num_isof "'.$num_isof.'"'."; ".
						'transc_list "'.$transc_list.'"'."; ".
						'biotype_list "'.$biotype_list.'"'."; ".
						'appris_gene_annot "'.$g_appris_annot.'"'."; ".
						'appris_annot "'.$t_appris_annot.'"'."; ".
						'appris_type "'.$appris_type.'"'."; ".
						'corsair_sc "'.$c_annot.'"'."; ".
						'tsl "'.$tsl_list.'"'."\n";
			}			
		}
	}
	return $output;
}

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

	- Get the annotations for all exons '/appris_data.exons.gff'

	- Get prin transc '/appris_data.exons.prin.gff'
	
	- Get alternative transc '/appris_data.exons.alt.gff';
	
	- Get the specific principal exons (that do not overlap with alternative exons) '/appris_data.exons.PI.gff';
	
	- Get the specific alternative exons (that do not overlap with principal exons) '/appris_data.exons.NPI.gff';

	- get the overlap regions PI+NPI '/appris_data.exons.OVER.gff';

=head1 SYNOPSIS

retrieve_exon_data

=head2 Required arguments:

	--data=  <Gene annotation file>
	
	--transl= <Translation sequences file>
	
	--main <Result file of APPRIS's scores>
	
	--outpath <Output path to save files>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl retrieve_exon_data_fromfile.pl
	
	--data=data/appris_data.annot.gtf
	
	--transl=data/appris_data.transl.fa
	
	--main=data/appris_data.appris.txt

	--outpath=data/
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
