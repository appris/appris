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
$logger->init_log();


#####################
# Method prototypes #
#####################
sub get_appris_annot($$$;$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Get data backup from file
	$logger->debug("-- get old main data that contains CCDS -------\n");
	my ($old_main_report, $old_seq_report) = common::get_main_report($input_old_main_file, $input_old_seq_file, undef);	
	#$logger->debug("OLD_MAIN_REPORT:\n".Dumper($old_main_report)."\n");
	#$logger->debug("OLD_SEQ_REPORT:\n".Dumper($old_seq_report)."\n");
	
	# Get data from file
	$logger->debug("-- get main data from files -------\n");
	my ($main_report, $seq_report) = common::get_main_report($input_main_file, $input_seq_file, $input_data_file);	
	#$logger->debug("MAIN_REPORT:\n".Dumper($main_report)."\n");
	#$logger->debug("SEQ_REPORT:\n".Dumper($seq_report)."\n");
	
	# Get label from file
	$logger->debug("-- get label data from files -------\n");
	my ($label_report) = common::get_label_report($input_label_file, $input_seq_file, $input_data_file);	
	#$logger->debug("LABEL_REPORT:\n".Dumper($label_report)."\n");
			
	# Get annots
	$logger->debug("-- get anntos -------\n");
	my ($output_content, $output_ens_content) = get_appris_annot($main_report, $label_report, $seq_report, $old_main_report);
	if ( ($output_content ne '') and ($output_ens_content ne '') ) {
		my ($printing_file_log) = printStringIntoFile($output_content,$output_file);
		die ("ERROR: printing ") unless(defined $printing_file_log);		
		my ($printing_file_log2) = printStringIntoFile($output_ens_content,$output_ens_file);
		die ("ERROR: printing ") unless(defined $printing_file_log2);		
	}	

	$logger->finish_log();
	
	exit 0;	
	
}
sub get_appris_annot($$$;$)
{
	my ($main_report, $label_report, $seq_report, $old_main_report) = @_;
	my ($output,$output_ens) = ('','');
	
	# get appris annot
	while (my ($gene_id, $g_report) = each(%{$label_report}) )
	{
		$logger->debug("-- $gene_id: ");
		
		my ($gene_name) = '-';
		if (exists $g_report->{'gene_name'}) {
			$gene_name = $g_report->{'gene_name'};	
		}
		
		my ($g_output);
		my ($len_report, $len_transcs);
		#my ($sc_report, $sc_transcs);
		#my ($ccds_seq_report,  $ccds_longest_report, $ccds_id_report, $ccds_transcs);
		my ($ccds_seq_report,  $ccds_id_report, $ccds_transcs);
		while (my ($transcript_id, $t_report) = each(%{$g_report->{'transcripts'}}) ) {
			$logger->debug("$transcript_id ");
			
			# get sequence
			my ($translation_seq) = '';
			my ($translation_length) = 0;
			if (exists $seq_report->{$gene_id}->{'transcripts'}->{$transcript_id}) {
				$translation_seq = $seq_report->{$gene_id}->{'transcripts'}->{$transcript_id};
				$translation_length = length($translation_seq);
			}
						
			# get CCDS ids (from the current version or from an older version)
			my ($ccds_id) = '-';
			if (exists $t_report->{'ccds_id'}) {
				$ccds_id = $t_report->{'ccds_id'};	
			}
			else {
				if ( defined $old_main_report and exists $old_main_report->{$gene_id}->{'transcripts'}->{$transcript_id}) {
					my ($old_g_report) = $old_main_report->{$gene_id};
					if (exists $old_g_report->{'transcripts'}->{$transcript_id}) {
						my ($old_t_report) = $old_g_report->{'transcripts'}->{$transcript_id};
						if (exists $old_t_report->{'ccds_id'}) {
							my ($old_ccds_id) = $old_t_report->{'ccds_id'};
							if ( exists $old_g_report->{'ccds_id'}->{$old_ccds_id} ) {
								my ($old_translation_seq) = $old_g_report->{'ccds_id'}->{$old_ccds_id};
								if ( $old_translation_seq eq $translation_seq ) {
									$ccds_id = $old_ccds_id;
								}
							}
						}
					}
				}				
			}
			
			# get 'appris' & 'corsair' score
			my ($corsair_score) = 0;
			my ($appris_score) = 0;
			if (exists $main_report->{$gene_id}->{'transcripts'}->{$transcript_id}) {
				my ($t_m_report) = $main_report->{$gene_id}->{'transcripts'}->{$transcript_id};
				if ( exists $t_m_report->{'annotations'} and (scalar(@{$t_m_report->{'annotations'}}) > 0) ) {
					my ($analysis) = $t_m_report->{'annotations'};
					$appris_score = $analysis->[9];
					$corsair_score = $analysis->[2];
				}
			}

			# get annotations
			if ( exists $t_report->{'annotations'} and (scalar(@{$t_report->{'annotations'}}) > 0) ) {
				my ($analysis) = $t_report->{'annotations'};				
				my ($appris_annot) = $analysis->[9];
				if ( $appris_annot eq 'YES' ) {
					my ($t_rep) = {
								'gene_id'				=> $gene_id,
								'gene_name'				=> $gene_name,
								'transcript_id'			=> $transcript_id,
								'ccds_id'				=> $ccds_id,
								'translation_length'	=> $translation_length,
								'appris_score'			=> $appris_score,
								'appris_annot'			=> $appris_annot,				
								'corsair_score'			=> $corsair_score,				
					};
					$g_output->{$transcript_id}		= $t_rep;
				}
				elsif ( $appris_annot eq 'UNKNOWN' ) {
					my ($t_rep) = {
								'gene_id'				=> $gene_id,
								'gene_name'				=> $gene_name,
								'transcript_id'			=> $transcript_id,
								'ccds_id'				=> $ccds_id,
								'translation_length'	=> $translation_length,
								'appris_score'			=> $appris_score,
								'appris_annot'			=> $appris_annot,
								'corsair_score'			=> $corsair_score,		
					};
					$g_output->{$transcript_id}		= $t_rep;
					push(@{$len_report->{$translation_length}}, $transcript_id);
					#push(@{$sc_report->{$appris_score}}, $transcript_id);
					if ( ($ccds_id ne '-') and ($translation_seq ne '') ) {
						push(@{$ccds_seq_report->{$translation_seq}}, $transcript_id);
						#push(@{$ccds_longest_report->{$translation_length}}, $transcript_id);
						unless ( exists $ccds_id_report->{$translation_seq} ) {
							$ccds_id_report->{$translation_seq} = $ccds_id;
						}
						else {
							if ( $ccds_id_report->{$translation_seq} ne $ccds_id ) {
								my ($old_ccds_id) = $ccds_id_report->{$translation_seq};
								my ($new_ccds_id) = $ccds_id;
								$old_ccds_id =~ s/CCDS//; $new_ccds_id =~ s/CCDS//;
								if ( $new_ccds_id < $old_ccds_id ) {
									$ccds_id_report->{$translation_seq} = $ccds_id;
									print STDERR "\nWARNING: $gene_id has duplicate CCDS ids for the same sequence. We change to the eldest id: CCDS$old_ccds_id -> CCDS$new_ccds_id\n";
								}
								else {
									print STDERR "\nWARNING: $gene_id has duplicate CCDS ids for the same sequence. We keep the eldest id: CCDS$new_ccds_id -> CCDS$old_ccds_id \n";
								}
							}
						}
						#push(@{$ccds_id_report->{$ccds_id}}, $translation_seq);
					}
				}
				#elsif ( $appris_annot eq 'NO' ) {
				#	my ($t_rep) = {
				#				'gene_id'				=> $gene_id,
				#				'transcript_id'			=> $transcript_id,
				#				'ccds_id'				=> $ccds_id,
				#				'translation_length'	=> $translation_length,
				#				'appris_annot'			=> 'alternative',				
				#	};
				#	$g_output->{$transcript_id}		= $t_rep;
				#}
			}
		}
		
$logger->debug("\n");
$logger->debug("CCDS_SEQ_REPORT:\n".Dumper($ccds_seq_report)."\n");
#$logger->debug("CCDS_LONGER_REPORT:\n".Dumper($ccds_longest_report)."\n");
$logger->debug("CCDS_ID_REPORT:\n".Dumper($ccds_id_report)."\n");
#$logger->debug("SC_REPORT:\n".Dumper($sc_report)."\n");
$logger->debug("LEN_REPORT:\n".Dumper($len_report)."\n");

		# step 1: save the candidates with ccds, and select the eldest one.
		my ($num_ccds) = 0;		
		my ($num_age_ccds) = 0;
		if ( defined $ccds_id_report and defined $ccds_seq_report ) {
			my (@sort_ccds_transc) = sort {
				my ($anum,$bnum);
				$ccds_id_report->{$a} =~ /CCDS(\d*)/;
				$anum = $1;
				$ccds_id_report->{$b} =~ /CCDS(\d*)/;
				$bnum = $1;
				$anum <=> $bnum
			} keys %{$ccds_id_report};
			#my (@sort_ccds_transc) = sort { $ccds_id_report->{$a} <=> $ccds_id_report->{$b} } keys %{$ccds_id_report};
			$num_ccds = scalar(@sort_ccds_transc);
			if ( $num_ccds >= 2 ) {
				my ($eldest_seq) = $sort_ccds_transc[0];
				my ($eldest_seq2) = $sort_ccds_transc[1];
				my ($eldest_seq_ccds) = $ccds_id_report->{$eldest_seq};
				my ($eldest_seq_ccds2) = $ccds_id_report->{$eldest_seq2};
				$eldest_seq_ccds =~ s/CCDS//; $eldest_seq_ccds2 =~ s/CCDS//;
$logger->debug("CCDS_AGE_TRANSC:\n".Dumper(@sort_ccds_transc)."\n");
$logger->debug("ELDEST_ID: ($eldest_seq_ccds2 - $eldest_seq_ccds): $eldest_seq\n");
$logger->debug("ELDEST_SEQ: $eldest_seq\n");
				if ( abs($eldest_seq_ccds2 - $eldest_seq_ccds) >= 10 ) {
					foreach my $transc_id (@{$ccds_seq_report->{$eldest_seq}}) {
						$ccds_transcs->{$transc_id} = $eldest_seq;
					}				
					$num_age_ccds = scalar(keys(%{$ccds_id_report}));						
				}				
			}				
		}

$logger->debug("CCDS_TRANSC:\n".Dumper($ccds_transcs)."\n");

		# step 2: save the candidates with ccds, and select the longest one.
		my ($num_length_ccds) = 0;
		#if ( defined $ccds_longest_report and defined $ccds_seq_report and !(defined $ccds_transcs) ) {
		if ( defined $ccds_seq_report and !(defined $ccds_transcs) ) {
			my (@sort_ccds_transc) = sort {length($b) <=> length($a)} keys %{$ccds_seq_report};
			my ($longest_seq) = $sort_ccds_transc[0];			
$logger->debug("CCDS_LONG_TRANSC:\n".Dumper(@sort_ccds_transc)."\n");
$logger->debug("LONG_SEQ:\n".Dumper($longest_seq)."\n");
			foreach my $transc_id (@{$ccds_seq_report->{$longest_seq}}) {
				$ccds_transcs->{$transc_id} = $longest_seq;
			}				
			#$num_ccds = scalar(keys(%{$ccds_seq_report}));
			#$num_length_ccds = scalar(keys(%{$ccds_longest_report}));
			$num_length_ccds = scalar(@sort_ccds_transc);
		}
						
		# step 3: from the list of step before, get longest transc
		if ( defined $len_report ) {
			my (@long_transc) = sort {$b <=> $a} keys %{$len_report};
			my ($long_seq) = $long_transc[0];
			foreach my $transc_id (@{$len_report->{$long_seq}}) {
				$len_transcs->{$transc_id} = $long_seq;
			}
		}
		
		# step 4: from the list of step before, get highest score transc
		#my ($num_a_scores) = 0;
		#my ($num_t_high_scores) = 0;
		#if ( defined $sc_report ) {
		#	my (@high_transc) = sort {$b <=> $a} keys %{$sc_report};
		#	if ( scalar(keys(%{$sc_report})) > 0 ) {
		#		$num_a_scores = scalar(keys(%{$sc_report}));
		#		my ($high_sc) = $high_transc[0];
		#		foreach my $transc_id (@{$sc_report->{$high_sc}}) {
		#			$sc_transcs->{$transc_id} = $high_sc;
		#		}				
		#	}
		#	$num_t_high_scores = scalar(keys(%{$sc_transcs}));
		#}
		
		
#$logger->debug("CCDS_TRANS:\n".Dumper($ccds_transcs)."\n");
##$logger->debug("SC_TRANS:\n".Dumper($sc_transcs)."\n");
#$logger->debug("LON_TRANS:\n".Dumper($len_transcs)."\n");

		# print annots
		if ( defined $g_output ) {
			while (my ($transc_id, $t_report) = each(%{$g_output}) ) {
				my ($appris_annot) = $t_report->{'appris_annot'};
				my ($annot) = '-';
				#my ($score) = $t_report->{'appris_score'};
				if ( $appris_annot eq 'YES' ) {
					#$annot = 'PRINCIPAL:1,APPRIS Principal Isoform';
					$annot = 'PRINCIPAL:1';
				}
				elsif ( $appris_annot eq 'UNKNOWN' ) {
					# retrieves eldest CCDS when there are different CCDS
					if ( $num_ccds == 1 and defined $ccds_transcs and exists $ccds_transcs->{$transc_id} ) {
							#$annot = 'PRINCIPAL:2,APPRIS candidate isoform with unique CCDS';
							$annot = 'PRINCIPAL:2';
					}
					# retrieves eldest CCDS when there are different CCDS
					elsif ( ($num_ccds > 1) and ($num_age_ccds > 1) and ($num_ccds == $num_age_ccds) and defined $ccds_transcs and exists $ccds_transcs->{$transc_id} ) {
						#$annot = 'PRINCIPAL:3,APPRIS candidate isoform with lower CCDS';
						$annot = 'PRINCIPAL:3';
					}
					# retrieves longest CCDS when there are different CCDS
					elsif ( ($num_ccds > 1) and ($num_length_ccds > 1) and ($num_ccds == $num_length_ccds) and defined $ccds_transcs and exists $ccds_transcs->{$transc_id} ) {
						#$annot = 'PRINCIPAL:4,APPRIS candidate isoform with longest CCDS';
						$annot = 'PRINCIPAL:4';
					}
					##elsif ( $num_a_scores == 1 ) { # when there is unique appris score
					#elsif ( $num_t_high_scores == 1 ) { # when there is unique transcript with hightest score of appris
					#	if ( defined $sc_transcs and exists $sc_transcs->{$transc_id} ) {
					#		$annot = 'appris_candidate_highest_score';
					#	}
					#}
					# retrieves longest prot. sequence when there are not CCDS
					elsif ( defined $len_transcs and exists $len_transcs->{$transc_id} and !(defined $ccds_transcs) ) {
						#$annot = 'PRINCIPAL:5,APPRIS candidate isoform with longest coding sequence';
						$annot = 'PRINCIPAL:5';
					}
					else {
						if ( exists $t_report->{'corsair_score'} and ($t_report->{'corsair_score'} >= 4.5) ) {
							#$annot = 'ALTERNATIVE:1,APPRIS candidate isoform that is conserved in at least three tested non-primate species';
							$annot = 'ALTERNATIVE:1';
						}
						else {
							#$annot = 'ALTERNATIVE:2,APPRIS candidate isoform that appears to be conserved in fewer than three tested non-primate species';
							$annot = 'ALTERNATIVE:2';
						}
					}
				}
				$output .=	$t_report->{'gene_name'}."\t".
							$t_report->{'gene_id'}."\t".
							$transc_id."\t".
							$t_report->{'ccds_id'}."\t".
							$annot."\n";
				$output_ens .=	$t_report->{'gene_id'}."\t".
								$transc_id."\t".
								$annot."\n";
			}
		}		
		$logger->debug("\n");
	}
	return ($output,$output_ens);
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
