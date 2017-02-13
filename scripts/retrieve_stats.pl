#!/usr/bin/perl -w

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use Data::Dumper;

use APPRIS::Utils::File qw( printStringIntoFile getTotalStringFromFile );
use APPRIS::Utils::Logger;

use lib "$FindBin::Bin/lib";
use common qw(
	get_seq_report
	get_longest_seq
	get_main_report
	get_label_report
	get_prin_report
);

###################
# Global variable #
###################
use vars qw(
	@METHODS
	@PRIN_LABELS
);

@METHODS = (
	'firestar',
	'matador3d',
	'corsair',
	'spade',
	'inertia',
	'cexonic',
	'thump',
	'crash_sp',
	'crash_tp',
	'appris'
);

@PRIN_LABELS = (
	'principal',
	'single',
	'uncertain_unique_ccds',
	'uncertain_earliest_ccds',
	'uncertain_longest_ccds',		
	'uncertain_longest_seq',
);

# Input parameters
my ($input_main_file) = undef;
my ($input_label_file) = undef;
my ($input_prin_file) = undef;
my ($input_seq_file) = undef;
my ($output_file) = undef;
my ($output_gdisagree_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'input-main=s' 		=> \$input_main_file,
	'input-label=s'		=> \$input_label_file,
	'input-prin=s' 		=> \$input_prin_file,	
	'input-seq=s'    	=> \$input_seq_file,	
	'output=s'			=> \$output_file,
	'out-gdisa=s'		=> \$output_gdisagree_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless(
	defined $input_main_file and
	defined $input_label_file and
	defined $input_prin_file and
	defined $input_seq_file and
	defined $output_file
) {
	print `perldoc $0`;
	exit 1;
}

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
sub _get_coverage($$);
sub _get_prin($);
sub _get_ccds_disagreement($);
sub _print_coverage($);
sub _print_ccds($);
sub _print_gene_disagree($$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Declare vars
	my ($output_content) = '';
	
	# Get data from file
	$logger->debug("-- get seq data from files -------\n");
	my ($seq_report) = common::get_seq_report($input_seq_file);	
	my ($lon_seq_report) = common::get_longest_seq($seq_report);
	#$logger->debug("SEQ_REPORT:\n".Dumper($seq_report)."\n");
	#$logger->debug("LONSEQ_REPORT:\n".Dumper($lon_seq_report)."\n");

	# Get data from file
	$logger->debug("-- get main data from files -------\n");
	my ($main_report) = common::get_main_report($input_main_file, $input_seq_file);	
	$logger->debug("MAIN_REPORT:\n".Dumper($main_report)."\n");

	# Get data from file
	$logger->debug("-- get detailed data from files -------\n");
	my ($label_report) = common::get_label_report($input_label_file, $input_seq_file);	
	#$logger->debug("LABEL_REPORT:\n".Dumper($label_report)."\n");

	# Get principal from file
	$logger->debug("-- get principal list from files -------\n");
	my ($prin_report) = common::get_prin_report($input_prin_file);	
	$logger->debug("PRIN_REPORT:\n".Dumper($prin_report)."\n");

	# Get data from file
	$logger->debug("-- get coverages of methods -------\n");
	my ($cov_report) = _get_coverage($main_report,$label_report);	
	#$logger->debug("COVERAGE_REPORT:\n".Dumper($cov_report)."\n");
	
	# Get data from file
	$logger->debug("-- get principal stats -------\n");
	my ($prin_stats_report) = _get_prin_stats($main_report,$prin_report);	
	#$logger->debug("PRIN_STATS_REPORT:\n".Dumper($prin_stats_report)."\n");

	# Get data from file
	$logger->debug("-- get agreement with genes with unique CCDS -------\n");
	my ($ccds_report, $genes_disagree) = _get_ccds_disagreement($label_report);	
	$logger->debug("CCDS_REPORT:\n".Dumper($ccds_report)."\n");
	$logger->debug("GENES_DISAGREE_CCDS:\n".Dumper($genes_disagree)."\n");
	
	# print stats
	$logger->debug("-- print coverages of methods -------\n");
	$output_content .= _print_coverage($cov_report);
	$logger->debug("-- print principal stats -------\n");
	$output_content .= _print_prin_stats($prin_stats_report);
	$logger->debug("-- print ccds agreements -------\n");
	$output_content .= _print_ccds($ccds_report);
	if ( $output_content ne '' ) {
		my ($printing_file_log) = printStringIntoFile($output_content, $output_file);
		$logger->error("Printing _get_ccds_stats_content\n") unless(defined $printing_file_log);			
	}
	
	# print disagreements
	$logger->debug("-- print ccds disagreements (genes/methods) -------\n");
	my ($output_content2) = _print_gene_disagree($genes_disagree,$main_report,$prin_report);
	if ( $output_content2 ne '' and defined $output_gdisagree_file ) {
		my ($printing_file_log) = printStringIntoFile($output_content2, $output_gdisagree_file);
		$logger->error("Printing _print_gene_disagree_content\n") unless(defined $printing_file_log);			
	}
		
	$logger->finish_log();
	
	exit 0;
		
}

sub _get_coverage($$)
{
	my ($main_report, $label_report) = @_;
	my ($report) = {
		'genes'		=> 0,
		'firestar'	=> 0,
		'matador3d'	=> 0,
		'corsair'	=> 0,
		'spade'		=> 0,
		'thump'		=> 0,
		'crash_sp'	=> 0,
		'crash_tp'	=> 0,
		'inertia'	=> 0,
		'proteo'	=> 0,
		'appris'	=> 0,
	};
	
	# scan genes
	foreach my $gene_id (keys %{$main_report}) {
		my ($gene_report) = $main_report->{$gene_id};
		my ($gene_label_report) = $label_report->{$gene_id};
		my ($gene_flags) = {
			'firestar'	=> 0,
			'matador3d'	=> 0,
			'corsair'	=> 0,
			'spade'		=> 0,
			'thump'		=> 0,
			'crash_sp'	=> 0,
			'crash_tp'	=> 0,
			'inertia'	=> 0,
			'proteo'	=> 0,
			'appris'	=> 0,
		};		
		foreach my $transcript_id (keys %{$gene_report->{'transcripts'}}) {
			my ($trans_report) = $gene_report->{'transcripts'}->{$transcript_id};
			my ($trans_label_report) = $gene_label_report->{'transcripts'}->{$transcript_id};
			if ( exists $trans_report->{'annotations'} ) {
				my ($annotation_list) = $trans_report->{'annotations'};
				my ($annotation_label_list) = $trans_label_report->{'annotations'};
				for ( my $i_met=0; $i_met < scalar(@{$annotation_list}); $i_met++) {
					if ( $i_met == 0 ) {
						my ($annot) = $annotation_list->[$i_met];
						if ( ($annot ne '-') and ($annot > 0) ) {
							my ($met) = 'firestar';
							$gene_flags->{$met} = 1;
						}
					}
					elsif ( $i_met == 1 ) {
						my ($annot) = $annotation_list->[$i_met];
						if ( ($annot ne '-') and ($annot > 0) ) {
							my ($met) = 'matador3d';
							$gene_flags->{$met} = 1;
						}
					}					
					elsif ( $i_met == 2 ) {
						my ($annot) = $annotation_list->[$i_met];
						if ( ($annot ne '-') and ($annot >= 1.5) ) {
							my ($met) = 'corsair';
							$gene_flags->{$met} = 1;
						}
					}					
					elsif ( $i_met == 3 ) {
						my ($annot) = $annotation_list->[$i_met];
						$annot =~ s/\-[^\$]*$//g;
						if ( ($annot ne '-') and ($annot > 0) ) {
							my ($met) = 'spade';
							$gene_flags->{$met} = 1;
						}
					}					
					elsif ( $i_met == 4 ) {
						my ($annot) = $annotation_list->[$i_met];
						$annot =~ s/\-[^\$]*$//g;
						if ( ($annot ne '-') and ($annot > 0) ) {
							my ($met) = 'thump';
							$gene_flags->{$met} = 1;
						}
					}					
					elsif ( $i_met == 5 ) {
						my ($annot) = $annotation_label_list->[$i_met];
						if ( ($annot ne '-') and ($annot eq 'YES') ) {
							my ($met) = 'crash_sp';
							$gene_flags->{$met} = 1;
						}
					}
					elsif ( $i_met == 6 ) {
						my ($annot) = $annotation_label_list->[$i_met];
						if ( ($annot ne '-') and ($annot eq 'YES') ) {
							my ($met) = 'crash_tp';
							$gene_flags->{$met} = 1;
						}
					}
					elsif ( $i_met == 7 ) {
						my ($annot) = $annotation_list->[$i_met];
						if ( ($annot ne '-') and ($annot > 0) ) {
							my ($met) = 'inertia';
							$gene_flags->{$met} = 1;
						}
					}					
					elsif ( $i_met == 8 ) {
						my ($annot) = $annotation_list->[$i_met];
						if ( ($annot ne '-') and ($annot > 0) ) {
							my ($met) = 'proteo';
							$gene_flags->{$met} = 1;
						}
					}					
					elsif ( $i_met == 9 ) {
						my ($annot) = $annotation_list->[$i_met];
						if ( ($annot ne '-') ) {
							my ($met) = 'appris';
							$gene_flags->{$met} = 1;
						}
					}					
				}
			}
		}		
		$report->{'genes'}++;
		if ( $gene_flags->{'firestar'} == 1 ) {
			$report->{'firestar'}++;
		}
		if ( $gene_flags->{'matador3d'} == 1 ) {
			$report->{'matador3d'}++;
		}
		if ( $gene_flags->{'corsair'} == 1 ) {
			$report->{'corsair'}++;
		}
		if ( $gene_flags->{'spade'} == 1 ) {
			$report->{'spade'}++;
		}
		if ( $gene_flags->{'thump'} == 1 ) {
			$report->{'thump'}++;
		}
		if ( $gene_flags->{'crash_sp'} == 1 ) {
			$report->{'crash_sp'}++;
		}
		if ( $gene_flags->{'crash_tp'} == 1 ) {
			$report->{'crash_tp'}++;
		}
		if ( $gene_flags->{'inertia'} == 1 ) {
			$report->{'inertia'}++;
		}
		if ( $gene_flags->{'proteo'} == 1 ) {
			$report->{'proteo'}++;
		}
		if ( $gene_flags->{'appris'} == 1 ) {
			$report->{'appris'}++;
		}
	}
	return $report;
	
} # end _get_coverage

sub _get_prin_stats($$)
{

	my ($main_report, $prin_report) = @_;
	my ($report) = {
		'genes'						=> 0,
		'principal'					=> 0,
		'single'					=> 0,
		'uncertain_unique_ccds'		=> 0,
		'uncertain_earliest_ccds'	=> 0,
		'uncertain_longest_ccds'	=> 0,		
		'uncertain_longest_seq'		=> 0,
	};
	
	# scan genes
	foreach my $gene_id (keys %{$main_report}) {
		my ($gene_report) = $main_report->{$gene_id};
		my ($gene_prin_report) = $prin_report->{$gene_id}->{'appris_num'};
		my ($num_trans) = scalar(keys %{$gene_report->{'varsplic'}} );
		# the order is important!!!
		if ( $num_trans == 1 ) {
			$report->{'single'}++;
		}
		elsif ( exists $gene_prin_report->{'PRINCIPAL:1'} and ($gene_prin_report->{'PRINCIPAL:1'} >= 1) ) {
			$report->{'principal'}++;
		}
		elsif ( (exists $gene_prin_report->{'PRINCIPAL:2'} and ($gene_prin_report->{'PRINCIPAL:2'} >= 1) ) ) {
			$report->{'uncertain_unique_ccds'}++;
		}
		elsif ( (exists $gene_prin_report->{'PRINCIPAL:3'} and ($gene_prin_report->{'PRINCIPAL:3'} >= 1) ) ) {
			$report->{'uncertain_earliest_ccds'}++;
		}
		elsif ( (exists $gene_prin_report->{'PRINCIPAL:4'} and ($gene_prin_report->{'PRINCIPAL:4'} >= 1) ) ) {
			$report->{'uncertain_longest_ccds'}++;
		}
		elsif ( (exists $gene_prin_report->{'PRINCIPAL:5'} and ($gene_prin_report->{'PRINCIPAL:5'} >= 1) ) ) {
			$report->{'uncertain_longest_seq'}++;
		}
		$report->{'genes'}++;		
	}
	return $report;
	
} # end _get_prin_stats

sub _get_ccds_disagreement($)
{
	my ($main_report) = @_;
	my ($report) = {
		'genes'		=> 0,
		'firestar'	=> 0,
		'matador3d'	=> 0,
		'corsair'	=> 0,
		'spade'		=> 0,
		'thump'		=> 0,
		'crash'		=> 0,
		'inertia'	=> 0,
		'proteo'	=> 0,
		'appris'	=> 0,
	};
	my ($g_disagre) = {
		'firestar'	=> {},
		'matador3d'	=> {},
		'corsair'	=> {},
		'spade'		=> {},
		'thump'		=> {},
		'crash'		=> {},
		'inertia'	=> {},
		'proteo'	=> {},
		'appris'	=> {},
	};
	
	# scan genes
	foreach my $gene_id (keys %{$main_report}) {
		my ($gene_report) = $main_report->{$gene_id};
		my ($gene_flags) = {
			'firestar'	=> 0,
			'matador3d'	=> 0,
			'corsair'	=> 0,
			'spade'		=> 0,
			'thump'		=> 0,
			'crash'		=> 0,
			'inertia'	=> 0,
			'proteo'	=> 0,
			'appris'	=> 0,
		};		
		my ($num_ccds);
		
		if ( exists $gene_report->{'ccds_id'} and defined $gene_report->{'ccds_id'} ) {
			my ($ccds_report) = $gene_report->{'ccds_id'};
			$num_ccds = scalar(keys(%{$ccds_report}));
		}
		unless ( defined $num_ccds and ($num_ccds == 1) ) {
			print STDERR "MORE THAN ONE: $gene_id\n";
			next;
		}
		
		foreach my $transcript_id (keys %{$gene_report->{'transcripts'}}) {
			my ($trans_report) = $gene_report->{'transcripts'}->{$transcript_id};
			my ($ccds_annot);
			if ( exists $trans_report->{'ccds_id'} and defined $trans_report->{'ccds_id'} and $trans_report->{'ccds_id'} ne '-' ) {
				$ccds_annot = $trans_report->{'ccds_id'};
				print STDERR "UNIQUE_CCDS\t$gene_id\t$transcript_id\t$ccds_annot\n";
			}
			if ( defined $ccds_annot and exists $trans_report->{'annotations'} ) {
				my ($annotation_list) = $trans_report->{'annotations'};
				for ( my $i_met=0; $i_met < scalar(@{$annotation_list}); $i_met++) {
					if ( $i_met == 0 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'firestar';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
#print STDERR "UNIQUE_CCDS_REJECT_".uc($met).":\t$gene_id\t$transcript_id\t$ccds_annot\n";
						}
					}
					elsif ( $i_met == 1 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'matador3d';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
#print STDERR "UNIQUE_CCDS_REJECT_".uc($met).":\t$gene_id\t$transcript_id\t$ccds_annot\n";							
						}
					}					
					elsif ( $i_met == 2 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'corsair';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
#print STDERR "UNIQUE_CCDS_REJECT_".uc($met).":\t$gene_id\t$transcript_id\t$ccds_annot\n";							
						}
					}					
					elsif ( $i_met == 3 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'spade';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
#print STDERR "UNIQUE_CCDS_REJECT_".uc($met).":\t$gene_id\t$transcript_id\t$ccds_annot\n";							
						}
					}	
					elsif ( $i_met == 4 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'thump';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
#print STDERR "UNIQUE_CCDS_REJECT_".uc($met).":\t$gene_id\t$transcript_id\t$ccds_annot\n";							
						}
					}					
					elsif ( $i_met == 7 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'inertia';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
						}
					}					
					elsif ( $i_met == 8 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'proteo';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
						}
					}					
					elsif ( $i_met == 9 ) {
						my ($annot) = $annotation_list->[$i_met];
						my ($met) = 'appris';
						if ( $annot eq 'NO' ) {
							$gene_flags->{$met} = 1;
							#$g_disagre->{$met}->{$gene_id} = 1 unless (exists $g_disagre->{$met}->{$gene_id} );
							$g_disagre->{$met}->{$gene_id}->{$transcript_id} = $ccds_annot;
#print STDERR "UNIQUE_CCDS_REJECT_".uc($met).":\t$gene_id\t$transcript_id\t$ccds_annot\n";
						}
					}					
				}
			}
			
		}		
		$report->{'genes'}++;
		if ( $gene_flags->{'firestar'} == 1 ) {
			$report->{'firestar'}++;
		}
		if ( $gene_flags->{'matador3d'} == 1 ) {
			$report->{'matador3d'}++;
		}
		if ( $gene_flags->{'corsair'} == 1 ) {
			$report->{'corsair'}++;
		}
		if ( $gene_flags->{'spade'} == 1 ) {
			$report->{'spade'}++;
		}
		if ( $gene_flags->{'thump'} == 1 ) {
			$report->{'thump'}++;
		}
		if ( $gene_flags->{'inertia'} == 1 ) {
			$report->{'inertia'}++;
		}
		if ( $gene_flags->{'proteo'} == 1 ) {
			$report->{'proteo'}++;
		}
		if ( $gene_flags->{'appris'} == 1 ) {
			$report->{'appris'}++;
		}
	}
	return ($report,$g_disagre);
	
} # end _get_ccds_disagreement

sub _print_coverage($)
{
	my ($report) = @_;
	my ($content) = '';
	
	if ( exists $report->{'genes'} and ($report->{'genes'} != 0) ) {
		foreach my $method ( @METHODS ) {
			my ($value) = '';
			$value = "0[0%]";
			if ( exists $report->{$method} and ($report->{$method} != 0) ) {
				my ($per) = sprintf("%.2f", ($report->{$method}/$report->{'genes'})*100);
				$value = $report->{$method}."[$per%]";
			}				
			$content .= $value."\t";
		}		
		if ( $content ne '' ) {
			$content = "# Coverage of methods/protein-coding genes\n".
						"genes"."\t".join("\t",@METHODS)."\n".
						$report->{'genes'}."\t".$content."\n";		
		}
	}
	return $content;
		
} # end _print_coverage

sub _print_prin_stats($)
{
	my ($report) = @_;
	my ($content) = '';
	
	if ( exists $report->{'genes'} and ($report->{'genes'} != 0) ) {
		foreach my $method ( @PRIN_LABELS ) {
			my ($value) = '';
			$value = "0[0%]";
			if ( exists $report->{$method} and ($report->{$method} != 0) ) {
				my ($per) = sprintf("%.2f", ($report->{$method}/$report->{'genes'})*100);
				$value = $report->{$method}."[$per%]";
			}				
			$content .= $value."\t";
		}		
		if ( $content ne '' ) {
			$content = "# Principal stats genes\n".
						"genes"."\t".join("\t",@PRIN_LABELS)."\n".
						$report->{'genes'}."\t".$content."\n";		
		}
	}
	return $content;
		
} # end _print_prin_stats

sub _print_ccds($)
{
	my ($report) = @_;
	my ($content) = '';
	
	if ( exists $report->{'genes'} and ($report->{'genes'} != 0) ) {
		foreach my $method ( @METHODS ) {
			my ($value) = '';
			$value = "0[0%]";
			if ( exists $report->{$method} and ($report->{$method} != 0) ) {
				my ($per) = sprintf("%.2f", ($report->{$method}/$report->{'genes'})*100);
				$value = $report->{$method}."[$per%]";
			}				
			$content .= $value."\t";
		}		
		if ( $content ne '' ) {
			$content = "# Percentage of rejected genes with unique CCDS in one or more isoform\n".
						"genes"."\t".join("\t",@METHODS)."\n".
						$report->{'genes'}."\t".$content."\n";	
		}
	}
	return $content;
	
} # end _print_ccds

sub _print_gene_disagree($$$)
{
	my ($report, $main_report, $prin_report) = @_;
	my ($content) = 'gene_id'."\t".
					'gene_name'."\t".
					'transc_id'."\t".
					'ccds_id'."\t".
					'firestar'."\t".
					'matador3d'."\t".
					'corsair'."\t".
					'spade'."\t".
					'thump'."\t".
					#'inertia'."\t".
					#'proteo'."\t".
					'appris'."\n";
	
	#foreach my $gene_id (keys %{$main_report}) {
	while (my ($gene_id,$gene_report) = each(%{$main_report}) ) {
		if (
			exists $report->{'firestar'}->{$gene_id} or
			exists $report->{'matador3d'}->{$gene_id} or
			exists $report->{'corsair'}->{$gene_id} or
			exists $report->{'spade'}->{$gene_id} or
			exists $report->{'thump'}->{$gene_id} or
			#exists $report->{'inertia'}->{$gene_id} or
			#exists $report->{'proteo'}->{$gene_id} or
			exists $report->{'appris'}->{$gene_id}
		) {
			my ($cont) = undef;
			if ( exists $report->{'firestar'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			if ( exists $report->{'matador3d'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			if ( exists $report->{'corsair'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			if ( exists $report->{'spade'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			if ( exists $report->{'thump'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			#if ( exists $report->{'inertia'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			#if ( exists $report->{'proteo'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			if ( exists $report->{'appris'}->{$gene_id} ) { push(@{$cont}, 1) } else { push(@{$cont}, 0) }
			#$content .= $gene_id."\t".join("\t", @{$cont})."\n";
			my ($ccds_disagree_cont) = join("\t", @{$cont});
			
			#foreach my $transc_id (keys %{$gene_report->{'transcripts'}}) {
			while (my ($transc_id,$transc_report) = each(%{$gene_report->{'transcripts'}}) ) {
print STDERR "TRANS_ID: $transc_id\n".Dumper($report->{'firestar'}->{$gene_id})."\n";
				if (
					exists $report->{'firestar'}->{$gene_id}->{$transc_id} or
					exists $report->{'matador3d'}->{$gene_id}->{$transc_id} or
					exists $report->{'corsair'}->{$gene_id}->{$transc_id} or
					exists $report->{'spade'}->{$gene_id}->{$transc_id} or
					exists $report->{'thump'}->{$gene_id}->{$transc_id} or
					#exists $report->{'inertia'}->{$gene_id}->{$transc_id} or
					#exists $report->{'proteo'}->{$gene_id}->{$transc_id} or
					exists $report->{'appris'}->{$gene_id}->{$transc_id}
				) {
					my ($ccds_id) = ( exists $transc_report->{'ccds_id'} ) ? $transc_report->{'ccds_id'} : '-';
					my ($gene_name) = ( exists $prin_report->{$gene_id} and exists $prin_report->{$gene_id}->{'gene_name'}) ? $prin_report->{$gene_id}->{'gene_name'} : '-';
					$content .= $gene_id."\t".
									$gene_name."\t".
									$transc_id."\t".
									$ccds_id."\t".
									$ccds_disagree_cont."\n";					
				}
			}			
		}
	}
	return $content;
		
} # end _print_gene_disagree


main();

1;


__END__

=head1 NAME

retrive_stats

=head1 DESCRIPTION

Count the coverage of methods and the percentage of disagrement with CCDS ids.

=head1 VERSION

0.1

=head2 Options

	--input-main <Result file of APPRIS's scores>

	--input-label <Result file of APPRIS's labels>
	
	--input-prin <Result file of principal isoforms>
	
	--input-seq <Input sequence comment file>

	--output <Output file that has the rejections of CCDS>
	
=head2 Examples

	perl stats_appris.pl

		--input-main=data/appris_data.appris.txt
		
		--input-label=data/appris_data.appris_label.txt
		
		--input-prin=data/appris_data.principal.txt

		--input-seq=data/appris_data.transl.fa
		
		--output=data/appris_stats.txt
		
		--out-gdisa=data/appris_stats.ccds_g.txt
		
		
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut

