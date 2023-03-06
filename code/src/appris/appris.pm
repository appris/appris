package appris;

use strict;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use Storable qw(dclone);

use APPRIS::Utils::Exception qw( info throw warning deprecate );

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	get_method_metric_names
	get_method_scores
	get_method_annots
	get_normalized_method_scores
	get_final_annotations
	get_score_output
	get_nscore_output
	get_label_output
);

###################
# Global variable #
###################
use vars qw(
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL

	$METHOD_METRICS
	$UNPHASED_METRICS
	$METRIC_LABELS
	$METRIC_PHASES
	$METRIC_EFF_RANGES
	$METRIC_WEIGHTED
);

$OK_LABEL				= 'YES';
$UNKNOWN_LABEL			= 'UNKNOWN';
$NO_LABEL				= 'NO';

$METHOD_METRICS = {
	'firestar'	=>  ['firestar'],
	'matador3d'	=> ['matador3d'],
	'matador3d2' => ['matador3d2'],
	'corsair'	=> ['corsair'],
	'spade' => [
		'spade_integrity',
		'spade'  # Spade bitscore
	],
	'thump' => ['thump'],
	'crash' => ['crash'],
	'inertia' => ['inertia'],
	'proteo' => ['proteo'],
	'appris' => ['appris']
};
$UNPHASED_METRICS = [
	'firestar',
	'matador3d',
	'matador3d2',
	'corsair',
	'spade',  # Spade bitscore
	'thump',
	'crash',
	'inertia',
	'proteo',
	'appris'
];
$METRIC_LABELS = {
	'firestar'	=> [ 'functional_residue',			'num_functional_residues'							],
	'matador3d'	=> [ 'conservation_structure',		'score_homologous_structure'						],
	'matador3d2'=> [ 'conservation_structure',		'score_homologous_structure'						],
	'corsair'	=> [ 'vertebrate_signal',			'score_vertebrate_signal'							],
	'spade_integrity' => [ 'domain_integrity_signal',	'score_domain_integrity_signal'],
	'spade' => [ 'domain_signal',	'score_domain_signal'],  # Spade bitscore
	'thump'		=> [ 'transmembrane_signal',		'score_transmembrane_signal'						],
	'crash'		=> [ 'peptide_signal',				'score_peptide_signal',
					 'mitochondrial_signal',		'score_mitochondrial_signal'						],
	'inertia'	=> [ 'unusual_evolution',			'num_unusual_exons'									],
	'proteo'	=> [ 'peptide_evidence',			'num_peptides'										],
	'appris'	=> [ 'principal_isoform',			'score_principal_isoform',			'reliability'	]
};
$METRIC_PHASES = {
	'phased_bi' => {
		'firestar'	=>  [1, 2],
		'matador3d'	=> [1, 2],
		'matador3d2' => [1, 2],
		'corsair'	=> [1, 2],
		'spade' => [1],  # Spade bitscore
		'spade_integrity' => [2],
		'thump' => [1, 2],
		'crash' => [1, 2],
		'inertia' => [1, 2],
		'proteo' => [1, 2],
		'appris' => [1, 2]
	},
	'phased_ib' => {
		'firestar'	=>  [1, 2],
		'matador3d'	=> [1, 2],
		'matador3d2' => [1, 2],
		'corsair'	=> [1, 2],
		'spade_integrity' => [1],
		'spade' => [2],  # Spade bitscore
		'thump' => [1, 2],
		'crash' => [1, 2],
		'inertia' => [1, 2],
		'proteo' => [1, 2],
		'appris' => [1, 2]
	}
};
$METRIC_EFF_RANGES = {
	'firestar'	=> 9,
	'matador3d'	=> 3,
	'spade_integrity' => 2.5,
	'spade' => 100  # Spade bitscore
};
$METRIC_WEIGHTED = {
	'firestar'	=> [{
					  'max'    => 2,
					  'weight' => 2
					},
					{
					  'max'    => 3,
				  	  'weight' => 3	
					},
					{
					  'max'    => 4,
				  	  'weight' => 4	
					},
					{
					  'max'    => 5,
				  	  'weight' => 5	
					},
					{
					  'max'    => 6,
					  'weight' => 6	
					}],
	'matador3d'	=> 6,
	'matador3d2'=> 6,
	'spade_integrity' => 2,
	'spade' => 6,  # Spade bitscore
	'thump'		=> 0,
	'crash'		=> 0,
	'inertia'	=> 0,
	'proteo'	=> 0,
};


sub appris_decider_classic($$$$$$$$)
{
	my ($gene, $tag, $princ_list, $isof_report, $scores, $s_scores, $nscores, $annots) = @_;

	# 2. from preserved transcripts, we keep those with a single unique CCDS
	$princ_list = step_ccds($princ_list, $isof_report, $nscores);
#	warning("CLASSIC_PRINC_LIST_2: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 2;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 3-1. from preserved transcripts, we keep those with eldest CCDS
	$princ_list = step_ccds_eldest($princ_list, $isof_report, $nscores);
#	warning("CLASSIC_PRINC_LIST_3: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 3;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 3-2. from preserved transcripts, we keep those with TSL1
	$princ_list = step_tsl($princ_list, $isof_report, $nscores);
#	warning("CLASSIC_PRINC_LIST_3: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 3;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 4. from preserved transcript, we keep those with longest seq with CCDS
	$princ_list = step_ccds_longest($princ_list, $isof_report, $nscores);
#	warning("CLASSIC_PRINC_LIST_4: \n".Dumper($princ_list)."\n");
#	warning("CLASSIC_ISOF_REPORT_4: \n".Dumper($isof_report)."\n");
#	warning("CLASSIC_NSCORES_4: \n".Dumper($nscores)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 4;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 5-1. from preserved transcripts, we keep those that have been validated manually
	$princ_list = step_validated($princ_list, $isof_report, $nscores);
#	warning("CLASSIC_PRINC_LIST_5_val: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 5;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 5-2. from preserved transcripts, we keep those with longest sequence
	$princ_list = step_longest($princ_list, $isof_report, $nscores);
#	warning("CLASSIC_PRINC_LIST_5_len: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 5;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}
	else {  # 5-3. we take the transcript with the first-sorting ID
		$princ_list = step_smaller_id($princ_list, $isof_report, $nscores);
#		warning("CLASSIC_PRINC_LIST_5_id: \n".Dumper($princ_list)."\n");
		step_tags(5, $scores, $princ_list, $isof_report, \$annots);
	}

	return $tag;

} # End appris_decider_classic

sub appris_decider_trifid($$$$$$$$$)
{
	my ($gene, $tag, $princ_list, $isof_report, $scores, $s_scores, $nscores, $annots,
	    $trifid_report) = @_;

	# 2. from preserved transcripts, we keep the dominant TRIFID transcript, if available
	my $min_lead = $main::EXP_CFG->val('appris', 'trifid_min_lead', $main::TRIFID_MIN_LEAD);
	$princ_list = step_trifid($princ_list, $scores, $gene, $trifid_report, $min_lead);
#	warning("TRIFID_PRINC_LIST_2 ".$main::TRIFID_MIN_LEAD.": \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 2;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 3. from preserved transcripts, we keep those with proteomics support
	if ( exists $s_scores->{'proteo'} ) {
		$princ_list = step_proteo($princ_list, $s_scores->{'proteo'});
#		warning("TRIFID_PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 3;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}
	}

	# 4. from preserved transcripts, we keep the best-scoring TRIFID transcript, if available
	$princ_list = step_trifid($princ_list, $scores, $gene, $trifid_report, 0);
#	warning("TRIFID_PRINC_LIST_4: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 4;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}

	# 5. from preserved transcripts, we keep those with longest sequence
	$princ_list = step_longest($princ_list, $isof_report, $nscores);
#	warning("TRIFID_PRINC_LIST_5_len: \n".Dumper($princ_list)."\n");
	if ( is_unique($princ_list, $isof_report) ) {
		$tag = 5;
		step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
		return $tag;
	}
	else {  # 5. we take the transcript with the first-sorting ID
		$princ_list = step_smaller_id($princ_list, $isof_report, $nscores);
#		warning("TRIFID_PRINC_LIST_5_id: \n".Dumper($princ_list)."\n");
		step_tags(5, $scores, $princ_list, $isof_report, \$annots);
	}

	return $tag;

} # End appris_decider_trifid

sub get_method_metric_names($) {
	my ($methods_str) = @_;
	my @methods = split(',', $methods_str);
	my @metrics = map { @{$METHOD_METRICS->{$_}} } @methods;
	return join(',', @metrics);
}

sub unphased_metric_filter($) {
	my ($metrics_str) = @_;
	my @metrics = split(',', $metrics_str);
	my @unphased_metrics;
	foreach my $metric (@metrics) {
		if ( grep {$_ eq $metric} @{$UNPHASED_METRICS} ) {
			push(@unphased_metrics, $metric);
		}
	}
	return join(',', @unphased_metrics);
}

sub phase_metric_filter($$$) {
	my ($metrics_str, $phasing_type, $phase) = @_;
	my @metrics = split(',', $metrics_str);
	my @phase_metrics;
	foreach my $metric (@metrics) {
		if ( grep {$_ == $phase} @{$METRIC_PHASES->{$phasing_type}->{$metric}} ) {
			push(@phase_metrics, $metric);
		}
	}
	return join(',', @phase_metrics);
}


# get the main functional isoform from methods of appris
sub get_method_scores($$$)
{
	my ($gene, $reports, $involved_methods) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($scores, $s_scores);
		
	# get annotations for each method (or group of method) -------------------
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($index) = $reports->{'_index_transcripts'}->{$transcript_id};
			my ($result) = $reports->transcripts->[$index];
			my ($aa_length) = length($transcript->translate->sequence);
			if ( !exists $s_scores->{'aa_length'} ) {
				$s_scores->{'aa_length'}->{'max'} = $aa_length;
				$s_scores->{'aa_length'}->{'min'} = $aa_length;
			}
			elsif ( $aa_length > $s_scores->{'aa_length'}->{'max'} ) {
				$s_scores->{'aa_length'}->{'max'} = $aa_length;
			}
			elsif ( $aa_length < $s_scores->{'aa_length'}->{'min'} ) {
				$s_scores->{'aa_length'}->{'min'} = $aa_length;
			}
			foreach my $m ( split(',', $involved_methods) ) {				
				if ( $m eq 'firestar' and $result and $result->analysis and $result->analysis->firestar ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->firestar;
					if ( defined $analysis->num_residues ) {
						my ($sc) = $analysis->num_residues;
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;
					}
				}
				elsif ( $m eq 'matador3d' and $result and $result->analysis and $result->analysis->matador3d ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->matador3d;
					if ( defined $analysis->score ) {
						my ($sc) = $analysis->score;
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;					
					}				
				}
				elsif ( $m eq 'matador3d2' and $result and $result->analysis and $result->analysis->matador3d2 ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->matador3d2;
					if ( defined $analysis->score ) {
						my ($sc) = $analysis->score;
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;					
					}				
				}
				elsif ( $m eq 'corsair' and $result and $result->analysis and $result->analysis->corsair ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->corsair;
					if ( defined $analysis->score ) {
						my ($sc) = $analysis->score;
						# if variant has 'start/stop codon not found', the score is 0
						#if ( $transcript->translate->codons ) {
						#	my ($codons) = '';
						#	foreach my $codon (@{$transcript->translate->codons}) {
						#		if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) { $codons .= $codon->type.',' }
						#	}
						#	unless ( $codons =~ /start/ and $codons =~ /stop/ ) { $sc = 0 }
						#}
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;
					}			
				}
				elsif ( $m eq 'spade' and $result and $result->analysis and $result->analysis->spade ) {
					my ($analysis) = $result->analysis->spade;
					my %metric_map = (
						'spade_integrity' => $analysis->domain_integrity,
						'spade' => $analysis->bitscore  # Spade bitscore
					);
					while (my ($metric_key, $metric_val) = each(%metric_map) ) {
						my ($annot_sc) = $METRIC_LABELS->{$metric_key}->[1];
						if ( defined $metric_val ) {
							my ($sc) = $metric_val;
							if ( !exists $s_scores->{$metric_key} ) {
								$s_scores->{$metric_key}->{'max'} = $sc;
								$s_scores->{$metric_key}->{'min'} = $sc;
							}
							elsif ( $sc > $s_scores->{$metric_key}->{'max'} ) {
								$s_scores->{$metric_key}->{'max'} = $sc;
							}
							elsif ( $sc < $s_scores->{$metric_key}->{'min'} ) {
								$s_scores->{$metric_key}->{'min'} = $sc;
							}
							push(@{$s_scores->{$metric_key}->{'scores'}->{$sc}}, $transcript_id);
							$scores->{$transcript_id}->{$annot_sc} = $sc;
						}
					}
				}
				elsif ( $m eq 'thump' and $result and $result->analysis and $result->analysis->thump ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->thump;
					if ( defined $analysis->num_tmh and 
						 defined $analysis->num_damaged_tmh ) {
						my ($sc) = $analysis->num_tmh;
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;					
					}				
				}
				elsif ( $m eq 'crash' and $result and $result->analysis and $result->analysis->crash ) {
					my ($annot_sc_sp) = $METRIC_LABELS->{$m}->[1];
					my ($annot_sc_tp) = $METRIC_LABELS->{$m}->[3];
					my ($analysis) = $result->analysis->crash;
					if ( defined $analysis->sp_score and defined $analysis->tp_score ) {
						my ($sc) = $analysis->sp_score.','.$analysis->tp_score;
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc_sp} = $analysis->sp_score;
						$scores->{$transcript_id}->{$annot_sc_tp} = $analysis->tp_score;
					}
				}
				elsif ( $m eq 'inertia' and $result and $result->analysis and $result->analysis->inertia ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->inertia;
					if ( defined $analysis->regions ) {
						 my ($sc) = get_num_unusual_exons($analysis);
						 push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						 $scores->{$transcript_id}->{$annot_sc} = $sc;
					}
				}
				elsif ( $m eq 'proteo' and $result and $result->analysis and $result->analysis->proteo ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->proteo;
					if ( defined $analysis->num_peptides ) {
						my ($sc) = $analysis->num_peptides;
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;					
					}
				}				
			}
		}
	}
	
	return ($scores, $s_scores);
	
} # End get_method_scores

# get the main functional isoform from methods of appris
sub get_method_scores_from_appris_rst($$$)
{
	my ($gene, $reports, $involved_methods) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($scores, $s_scores);
	
	# get annotations for each method (or group of method) -------------------
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($index) = $reports->{'_index_transcripts'}->{$transcript_id};
			my ($result) = $reports->transcripts->[$index];
			my ($aa_length) = length($transcript->translate->sequence);
			if ( !exists $s_scores->{'aa_length'} ) {
				$s_scores->{'aa_length'}->{'max'} = $aa_length;
				$s_scores->{'aa_length'}->{'min'} = $aa_length;
			}
			elsif ( $aa_length > $s_scores->{'aa_length'}->{'max'} ) {
				$s_scores->{'aa_length'}->{'max'} = $aa_length;
			}
			elsif ( $aa_length < $s_scores->{'aa_length'}->{'min'} ) {
				$s_scores->{'aa_length'}->{'min'} = $aa_length;
			}
			foreach my $m ( split(',', $involved_methods) ) {
				if ( $m eq 'firestar' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->functional_residues_score ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->functional_residues_score;
					if ( !exists $s_scores->{$m} ) {
						$s_scores->{$m}->{'max'} = $sc;
						$s_scores->{$m}->{'min'} = $sc;
					}
					elsif ( $sc > $s_scores->{$m}->{'max'} ) {
						$s_scores->{$m}->{'max'} = $sc;
					}
					elsif ( $sc < $s_scores->{$m}->{'min'} ) {
						$s_scores->{$m}->{'min'} = $sc;
					}
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc} = $sc;
				}
				elsif ( $m eq 'matador3d' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->homologous_structure_score ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->homologous_structure_score;
					if ( !exists $s_scores->{$m} ) {
						$s_scores->{$m}->{'max'} = $sc;
						$s_scores->{$m}->{'min'} = $sc;
					}
					elsif ( $sc > $s_scores->{$m}->{'max'} ) {
						$s_scores->{$m}->{'max'} = $sc;
					}
					elsif ( $sc < $s_scores->{$m}->{'min'} ) {
						$s_scores->{$m}->{'min'} = $sc;
					}
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc} = $sc;					
				}
				elsif ( $m eq 'matador3d2' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->homologous_structure_score ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->homologous_structure_score;
					if ( !exists $s_scores->{$m} ) {
						$s_scores->{$m}->{'max'} = $sc;
						$s_scores->{$m}->{'min'} = $sc;
					}
					elsif ( $sc > $s_scores->{$m}->{'max'} ) {
						$s_scores->{$m}->{'max'} = $sc;
					}
					elsif ( $sc < $s_scores->{$m}->{'min'} ) {
						$s_scores->{$m}->{'min'} = $sc;
					}
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc} = $sc;					
				}
				elsif ( $m eq 'corsair' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->vertebrate_conservation_score ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->vertebrate_conservation_score;
					if ( !exists $s_scores->{$m} ) {
						$s_scores->{$m}->{'max'} = $sc;
						$s_scores->{$m}->{'min'} = $sc;
					}
					elsif ( $sc > $s_scores->{$m}->{'max'} ) {
						$s_scores->{$m}->{'max'} = $sc;
					}
					elsif ( $sc < $s_scores->{$m}->{'min'} ) {
						$s_scores->{$m}->{'min'} = $sc;
					}
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc} = $sc;
				}
				elsif ( $m eq 'spade' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->domain_score ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->domain_score;
					if ( !exists $s_scores->{$m} ) {
						$s_scores->{$m}->{'max'} = $sc;
						$s_scores->{$m}->{'min'} = $sc;
					}
					elsif ( $sc > $s_scores->{$m}->{'max'} ) {
						$s_scores->{$m}->{'max'} = $sc;
					}
					elsif ( $sc < $s_scores->{$m}->{'min'} ) {
						$s_scores->{$m}->{'min'} = $sc;
					}
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc} = $sc;						
				}
				elsif ( $m eq 'thump' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->transmembrane_helices_score ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->transmembrane_helices_score;
					if ( !exists $s_scores->{$m} ) {
						$s_scores->{$m}->{'max'} = $sc;
						$s_scores->{$m}->{'min'} = $sc;
					}
					elsif ( $sc > $s_scores->{$m}->{'max'} ) {
						$s_scores->{$m}->{'max'} = $sc;
					}
					elsif ( $sc < $s_scores->{$m}->{'min'} ) {
						$s_scores->{$m}->{'min'} = $sc;
					}
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc} = $sc;					
				}
				elsif ( $m eq 'crash' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->peptide_score and defined $result->analysis->appris->mitochondrial_score ) {
					my ($annot_sc_sp) = $METRIC_LABELS->{$m}->[1];
					my ($annot_sc_tp) = $METRIC_LABELS->{$m}->[3];
					my ($sc_sp) = $result->analysis->appris->peptide_score;
					my ($sc_tp) = $result->analysis->appris->mitochondrial_score;
					if ( $sc_sp ne '-' and $sc_tp ne '-' ) {
						my ($sc) = $sc_sp.','.$sc_tp;					
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc_sp} = $sc_sp;
						$scores->{$transcript_id}->{$annot_sc_tp} = $sc_tp;
					}
				}
				elsif ( $m eq 'inertia' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->peptide_score and defined $result->analysis->appris->mitochondrial_score ) {
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($analysis) = $result->analysis->inertia;
					if ( defined $analysis->regions ) {
						 my ($sc) = get_num_unusual_exons($analysis);
						 push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						 $scores->{$transcript_id}->{$annot_sc} = $sc;
					}
				}
				elsif ( $m eq 'proteo' and $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->peptide_evidence_score ) {				
					my ($annot_sc) = $METRIC_LABELS->{$m}->[1];
					my ($sc) = $result->analysis->appris->peptide_evidence_score;
					if ( $sc ne '-' ) {
						if ( !exists $s_scores->{$m} ) {
							$s_scores->{$m}->{'max'} = $sc;
							$s_scores->{$m}->{'min'} = $sc;
						}
						elsif ( $sc > $s_scores->{$m}->{'max'} ) {
							$s_scores->{$m}->{'max'} = $sc;
						}
						elsif ( $sc < $s_scores->{$m}->{'min'} ) {
							$s_scores->{$m}->{'min'} = $sc;
						}
						push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
						$scores->{$transcript_id}->{$annot_sc} = $sc;					
					}
				}				
			}
		}
	}
	
	return ($scores, $s_scores);
	
} # End get_method_scores_from_appris_rst

# get the main functional isoform from methods of appris
sub get_method_annots($$$)
{
	my ($gene, $scores, $involved_metrics) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($annots);
		
	# get annotations for each method (or group of method) -------------------
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {			
			foreach my $m ( split(',', $involved_metrics) ) {
				if ( $m eq 'firestar' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_firestar_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'matador3d' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_matador3d_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'matador3d2' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_matador3d2_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'spade_integrity' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_spade_integrity_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'spade' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {
					# Spade bitscore
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_spade_bitscore_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'corsair' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_corsair_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'thump' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
					my ($annot_label) = $METRIC_LABELS->{$m}->[0];
					get_thump_annots($gene, $scores->{$m}, $annot_label, \$annots);
				}
				elsif ( $m eq 'crash' and defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
					my ($annot_label_sp) = $METRIC_LABELS->{$m}->[0];
					my ($annot_label_tp) = $METRIC_LABELS->{$m}->[2];
					get_crash_annots($gene, $scores->{$m}, $annot_label_sp, $annot_label_tp, \$annots);
				}				
			}
		}
	}
	
	return $annots;
	
} # End get_method_annots

# get normalized scores of methods
sub get_normalized_method_scores($$$\$\$)
{
	my ($gene, $annots, $involved_metrics, $ref_scores, $ref_s_scores) = @_;
	my ($nscores);

	# obtain normalized scores for each method
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			foreach my $metric ( split(',', $involved_metrics) ) {
				my ($n_sc) = 0;
				my ($sc) = 0;
				my ($max) = (exists $$ref_s_scores->{$metric} and exists $$ref_s_scores->{$metric}->{'max'})? $$ref_s_scores->{$metric}->{'max'} : 0;
				my ($min) = (exists $$ref_s_scores->{$metric} and exists $$ref_s_scores->{$metric}->{'min'})? $$ref_s_scores->{$metric}->{'min'} : 0;
				my ($label) = $METRIC_LABELS->{$metric}->[1];
				my ($alabel) = $METRIC_LABELS->{$metric}->[0];
				my ($appris_label) = $annots->{$transcript_id}->{$alabel};

				# create normalize scores
				if ( exists $$ref_scores->{$transcript_id}->{$label} ) {
					$sc = $$ref_scores->{$transcript_id}->{$label};

					if ( grep { $_ eq $metric } keys %{$METRIC_EFF_RANGES} ) {

						if ( $max > 0.0 ) {
							my $eff_range = $main::EXP_CFG->val( $metric, 'effective_range', $METRIC_EFF_RANGES->{$metric} );
							if ( looks_like_number($eff_range) && $eff_range > 0.0 ) {
								$eff_range = $eff_range <= $max ? $eff_range : $max ;
							} else {
								die("invalid effective range: $eff_range");
							}

							$n_sc = $max - $sc < $eff_range ? ($eff_range - ($max - $sc)) / $eff_range : 0.0;
						} else {
							$n_sc = 0;
						}
					} elsif ( $metric eq 'corsair' ) {
						my ($delta) = 1.0;
						my ($d_sc) = $sc > $delta ? $sc - $delta : 0.0 ;  # decrement Corsair score
						my ($d_max) = $max > $delta ? $max - $delta : 0.0 ;  # decrement maximum Corsair score
						$n_sc = $d_max > 0.0 ? $d_sc / $d_max : 0.0 ;
					} else {
						if ( defined $appris_label && $appris_label ne $NO_LABEL ) { $sc = $max } # give the max value if it pass the method filters (method annotations)
						if ( $max != 0 and ($max - $min != 0) ) { $n_sc = $sc/$max } # normalize when there are differences between the max and min
						else { $n_sc = 0 }
					}
				}
				if ( $n_sc < 0 ) { $n_sc = 0 }
				$n_sc = sprintf("%.3f",$n_sc);
				$nscores->{$transcript_id}->{$metric} = $n_sc;
			}
		}
	}

	return ($nscores);

} # End get_normalized_method_scores

# get gene scores of methods for appris
sub get_appris_scores($$\$\$\$)
{
	my ($gene, $involved_metrics, $ref_scores, $ref_s_scores, $ref_n_scores) = @_;

	# clear any existing APPRIS scores, so that we can get APPRIS score multiple times
	foreach my $transcript_id ( keys(%{$$ref_scores}) ) {
		if ( exists $$ref_scores->{$transcript_id}->{'score_principal_isoform'} ) {
			delete $$ref_scores->{$transcript_id}->{'score_principal_isoform'};
		}
	}
	if ( exists $$ref_s_scores->{'appris'} ) {
		delete $$ref_s_scores->{'appris'};
	}
	foreach my $transcript_id ( keys(%{$$ref_n_scores}) ) {
		if ( exists $$ref_n_scores->{$transcript_id}->{'appris'} ) {
			delete $$ref_n_scores->{$transcript_id}->{'appris'};
		}
	}

	# obtain normalize scores (weighted normalize score) foreach method
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my $is_x = index($transcript->translate->sequence, "X");
			my ($appris_score) = 0;
			foreach my $metric ( split(',', $involved_metrics) ) {
				my ($n_sc) = $$ref_n_scores->{$transcript_id}->{$metric};
				my ($max) = (exists $$ref_s_scores->{$metric} and exists $$ref_s_scores->{$metric}->{'max'})? $$ref_s_scores->{$metric}->{'max'} : 0;
				my ($min) = (exists $$ref_s_scores->{$metric} and exists $$ref_s_scores->{$metric}->{'min'})? $$ref_s_scores->{$metric}->{'min'} : 0;

				# apply weights to normalize scores
				my ($weight) = 0;
				if ( $metric eq 'firestar' ) {
					if    ( $max >= $METRIC_WEIGHTED->{$metric}->[4]->{'max'} ) { $weight = $METRIC_WEIGHTED->{$metric}->[4]->{'weight'} }
					elsif ( $max >= $METRIC_WEIGHTED->{$metric}->[3]->{'max'} ) { $weight = $METRIC_WEIGHTED->{$metric}->[3]->{'weight'} }
					elsif ( $max >= $METRIC_WEIGHTED->{$metric}->[2]->{'max'} ) { $weight = $METRIC_WEIGHTED->{$metric}->[2]->{'weight'} }
					elsif ( $max >= $METRIC_WEIGHTED->{$metric}->[1]->{'max'} ) { $weight = $METRIC_WEIGHTED->{$metric}->[1]->{'weight'} }
					elsif ( $max >= $METRIC_WEIGHTED->{$metric}->[0]->{'max'} ) { $weight = $METRIC_WEIGHTED->{$metric}->[0]->{'weight'} }
				}
				elsif ( $metric eq 'corsair' ) {
					if    ( $max >= 13  )  { $weight = 4; }
					elsif ( $max >= 8   )  { $weight = 3; }
					elsif ( $max >= 4   )  { $weight = 2; }
					elsif ( $max >= 1   )  { $weight = 1; }
					else                   { $weight = 0; }
				}
				elsif ( $metric eq 'spade_integrity' ) {
					$weight = $main::EXP_CFG->val( 'spade_integrity', 'score_weight', $METRIC_WEIGHTED->{'spade_integrity'} );
				}
				else {
					$weight = $METRIC_WEIGHTED->{$metric};
				}
				$appris_score += $weight*$n_sc;
			}
			
			# Penalize X
			if ( $is_x != -1 ) {
				$appris_score = -1;
			}			

			# filter by biotype
			if ( $transcript->biotype and ( ($transcript->biotype eq 'nonsense_mediated_decay') or ($transcript->biotype eq 'polymorphic_pseudogene') ) ) {
				$appris_score = -1;
			}
			# filter by readthrough_transcript
			if ( $transcript->tag and $transcript->tag =~ /readthrough_transcript/ ) {
				$appris_score = -1;
			}
			# filter by start/stop codon not found
			if ( $transcript->translate->codons ) {
				my ($codons) = '';
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) { $codons .= $codon->type.',' }
				}
				unless ( $codons =~ /start/ and $codons =~ /stop/ ) { $appris_score = -1 }
			} else {
				$appris_score = -1;
			}
			# save appris score and normalize score
			my ($m) = 'appris';
			my ($label) = $METRIC_LABELS->{$m}->[1];
			if ( !exists $$ref_s_scores->{$m} ) {
				$$ref_s_scores->{$m}->{'max'} = $appris_score;
				$$ref_s_scores->{$m}->{'min'} = $appris_score;
			}
			elsif ( $appris_score > $$ref_s_scores->{$m}->{'max'} ) {
				$$ref_s_scores->{$m}->{'max'} = $appris_score;
			}
			elsif ( $appris_score < $$ref_s_scores->{$m}->{'min'} ) {
				$$ref_s_scores->{$m}->{'min'} = $appris_score;
			}
			push(@{$$ref_s_scores->{$m}->{'scores'}->{$appris_score}}, $transcript_id);
			$$ref_scores->{$transcript_id}->{$label} = $appris_score;
		}
	}

	# get normalize scores for appris
	my ($m) = 'appris';
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($n_sc) = 0;
			my ($max) = $$ref_s_scores->{$m}->{'max'};
			my ($min) = $$ref_s_scores->{$m}->{'min'};
			if ( $max != 0 and ($max - $min != 0) ) {
				my ($label) = $METRIC_LABELS->{$m}->[1];
				my ($sc) = $$ref_scores->{$transcript_id}->{$label};
				$n_sc = ($sc - $min)/($max - $min);
				$n_sc = sprintf("%.3f",$n_sc);
			}
			$$ref_n_scores->{$transcript_id}->{$m} = $n_sc;
		}
	}

} # End get_appris_scores

# get the final annotation
sub get_final_annotations($$$$$$;$)
{
	my ($gene, $scores, $s_scores, $nscores, $annots, $involved_metrics, $trifid_report) = @_;
	my ($method) = 'appris';
	my ($tag) = 0;

	my $appris_score_mode = $main::EXP_CFG->val( 'appris', 'score_mode', 'phased_bi' );
	if ( $appris_score_mode =~ /^phased_/ ) {
		my $phase_1_metrics = phase_metric_filter($involved_metrics, $appris_score_mode, 1);
		get_appris_scores($gene, $phase_1_metrics, $scores, $s_scores, $nscores);
	} elsif ( $appris_score_mode eq 'all' ) {
		# include all metrics independently, even those from the same method
		get_appris_scores($gene, $involved_metrics, $scores, $s_scores, $nscores);
	} elsif ( $appris_score_mode eq 'unphased' ) {
		my $unphased_metrics = unphased_metric_filter($involved_metrics);
		get_appris_scores($gene, $unphased_metrics, $scores, $s_scores, $nscores);
	} else {
		die("unknown APPRIS score mode: ${appris_score_mode}");
	}

	# scan transcripts sorted by appris score
	if ( defined $s_scores and exists $s_scores->{$method} and exists $s_scores->{$method}->{'scores'} and scalar(keys(%{$s_scores->{$method}->{'scores'}})) > 0 )
	{
		# 0. create isoform list and report for gene.
		my ($princ_list, $isof_report) = init_isoform_report($gene);
#		warning("GENE_0: \n".Dumper($gene)."\n");
#		warning("PRINC_LIST_0: \n".Dumper($princ_list)."\n");
#		warning("PRINC_REP_0: \n".Dumper($isof_report)."\n");

		# 1_1 acquire the dominant transcripts from first-phase appris score.
		# They have to pass the cutoff to be added into "principal" list
		$princ_list = step_appris($princ_list, $isof_report, $s_scores->{$method});
#		warning("PRINC_LIST_1_1: \n".Dumper($princ_list)."\n");
#		warning("PRINC_REP_1_1: \n".Dumper($isof_report)."\n");

		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 1;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return ($scores, $nscores, $annots);
		}

		if ( $appris_score_mode =~ /^phased_/ ) {
			# 1_2 acquire the dominant transcripts from second-phase appris score,
			# if there is overlap with the dominant transcripts from the first phase;
			# otherwise revert to using the dominant transcripts from the first phase.
			my $phase2_metrics = phase_metric_filter($involved_metrics, $appris_score_mode, 2);

			my $phase2_scores = dclone($scores);
			my $phase2_s_scores = dclone($s_scores);
			my $phase2_nscores = dclone($nscores);
			get_appris_scores($gene, $phase2_metrics, $phase2_scores, $phase2_s_scores, $phase2_nscores);

			my $phase2_isof_report = dclone($isof_report);
			my $phase2_princ_list = step_appris($princ_list, $phase2_isof_report, $phase2_s_scores->{$method});

			if ( scalar keys %{$phase2_princ_list} > 0 ) {
				$princ_list = $phase2_princ_list;
				$isof_report = $phase2_isof_report;
				$scores = $phase2_scores;
				$nscores = $phase2_nscores;
				$s_scores = $phase2_s_scores;
				# warning("PRINC_LIST_1_2: \n".Dumper($princ_list)."\n");
				# warning("PRINC_REP_1_2: \n".Dumper($isof_report)."\n");
			}

			if ( is_unique($princ_list, $isof_report) ) {
				$tag = 1;
				step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
				return ($scores, $nscores, $annots);
			}
		}

		if ( defined($trifid_report) ) {
			$tag = appris_decider_trifid($gene, $tag, $princ_list, $isof_report,
				$scores, $s_scores, $nscores, $annots, $trifid_report);
		} else {
			$tag = appris_decider_classic($gene, $tag, $princ_list, $isof_report,
				$scores, $s_scores, $nscores, $annots);
		}
	}

	return ($scores, $nscores, $annots);

} # End get_final_annotations

sub is_unique($$)
{
	my ($princ_list, $isof_report) = @_;
	my ($unique) = 0;
	
	my ($same_seq) = 1;
	my ($num_prin) = 0;
	my ($seq) = '';
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $princ_list->{$transc_id} ) {
			$num_prin++;
			my ($transl_seq) = $princ->{'seq'};
			if ( $seq eq '' ) {
				$seq = $transl_seq
			}
			else {
				# skip when we find one sequences is different
				unless ( $seq eq $transl_seq ) {
					$same_seq = 0; 
					last;
				}					
			}
		}
	}
	if ( ($num_prin >= 1) and ($same_seq == 1) ) { $unique = 1 }		
	
	return $unique;
	
} # end is_unique


sub init_isoform_report($)
{
	my ($gene) = @_;
	my ($princ_list, $isof_report);

	# generate initial report for each transcript of gene
	while( my($transc_id, $transc_idx) = each %{$gene->{'_index_transcripts'}} ) {
		my ($transcript) = $gene->transcripts->[$transc_idx];

		if ( $transcript->translate and $transcript->translate->sequence ) {
			$princ_list->{$transc_id} = 1;

			my ($transc_name) = $transcript->external_name;
			my ($transl_seq) = $transcript->translate->sequence;
			my ($transc_rep) = {
					'id'		=> $transc_id,
					'name'		=> $transc_name,
					'seq'		=> $transl_seq,
					'length'	=> length($transl_seq)
			};
			# save CCDS id
			if ( $transcript->xref_identify ) {
				foreach my $xref_identify (@{$transcript->xref_identify}) {
					if ($xref_identify->dbname eq 'CCDS') {
						my ($ccds_id) = $xref_identify->id;
						if ( defined $ccds_id ) {
							$ccds_id =~ s/CCDS//; $ccds_id =~ s/\.\d+$//;
							$transc_rep->{'ccds'} = $ccds_id;
						}
					}
				}
			}
			# save TSL(1) annot
			if ( $transcript->tsl and $transcript->tsl eq '1' ) {
				$transc_rep->{'tsl'} = 1;
			}

			push(@{$isof_report}, $transc_rep);
		}
	}
	return ($princ_list, $isof_report);

} # end init_isoform_report

sub step_appris($$$)
{
	my ($i_princ_list, $isof_report, $s_scores) = @_;
	my ($princ_list) = {};

	# map transcript IDs to report indices, clear any 'principal'
	# labels that may have been added in a previous APPRIS step
	my (%index_transcripts);
	foreach my $transc_idx ( 0 .. $#{ $isof_report } ) {
		my $transc_id = $isof_report->[$transc_idx]{'id'};
		$index_transcripts{$transc_id} = $transc_idx;
		if ( exists $isof_report->[$transc_idx]{'principal'} ) {
			delete $isof_report->[$transc_idx]{'principal'};
		}
	}

	# scan the transcripts from the sorted APPRIS scores
	my ($highest_score) = $s_scores->{'max'};
	my (@sorted_ap_scores) = sort { $b <=> $a } keys (%{$s_scores->{'scores'}});
	my ($num_distinct_scores) = scalar(@sorted_ap_scores);
	for ( my $i = 0; $i < scalar(@sorted_ap_scores); $i++ ) {
		my ($ap_score) = $sorted_ap_scores[$i];

		# filter by input principal transcript list, to take account of previous APPRIS steps
		my @score_transc_ids = grep { exists $i_princ_list->{$_} } @{$s_scores->{'scores'}->{$ap_score}};

		foreach my $transc_id (@score_transc_ids) {
			# APPRIS says could be a principal when:
			# 1.	the Core of pipeline has not rejected the transcript (IT DOES NOT APPLY YET)
			# APPRIS rejected a transcript when:
			# 2.	it is a NMD: app_score is -1. Exception: we accept the last terms when the transcript is unique
			#
			my ($core_flag) = 1;
			my ($unique_transc) = ( $num_distinct_scores == 1 and scalar(@score_transc_ids) >= 1 ) ? 1 : 0 ;
			if ( $core_flag == 1 ) {
				if ( ( ($highest_score - $ap_score) <=  $main::APPRIS_CUTOFF) and ( ($ap_score >= 0) or ($unique_transc == 1) ) ) {
					$isof_report->[$index_transcripts{$transc_id}]{'principal'} = 1;
					$princ_list->{$transc_id} = 1;
				}
			}
		}
	}
	return ($princ_list);
	
} # end step_appris

sub step_trifid($$$$$)
{
	my ($i_princ_list, $scores, $gene, $trifid_report, $min_lead) = @_;
	my ($report);

	my (@scoring_principals) = grep {
		$scores->{$_}{'score_principal_isoform'} >= 0  # i.e. not -1
	} keys(%{$i_princ_list});

	if (@scoring_principals) {

		if ( defined($trifid_report) ) {

			my (%transc_to_seq);
			my ($num_missing_seqs) = 0;
			foreach my $transc_id (@scoring_principals) {
				my ($transcript) = $gene->transcript($transc_id);
				if ( $transcript->translate and $transcript->translate->sequence ) {
					$transc_to_seq{$transc_id} = $transcript->translate->sequence;
				} else {
					$num_missing_seqs++;
				}
			}

			if ( $num_missing_seqs == 0 ) {

				my (%seq_to_score);
				foreach my $transc_report (@{$trifid_report->transcripts}) {

					my ($transc_id) = $transc_report->stable_id;
					next unless( grep { $transc_id eq $_ } @scoring_principals );
					next unless( defined($transc_report->analysis) &&
								 defined($transc_report->analysis->trifid) );

					my ($score) = $transc_report->analysis->trifid->norm_trifid_score;
					my ($seq) = $transc_to_seq{$transc_id};
					if ( ! exists($seq_to_score{$seq}) ||
							$score > $seq_to_score{$seq} ) {
						$seq_to_score{$seq} = $score;
					}
				}

				my (%score_to_seqs);
				while ( my ($seq, $score) = each(%seq_to_score) ) {
					push(@{$score_to_seqs{$score}}, $seq);
				}
				my @dec_scores = sort { $b <=> $a } keys(%score_to_seqs);
				# get the best transcripts comparing with the best score
				my $num_dec_scores = scalar(@dec_scores);
				my $best_score = $dec_scores[0];
				if ( $num_dec_scores == 1 || ($num_dec_scores >= 2 &&
						($best_score - $dec_scores[1]) >= $min_lead) ) {
					# save the transc with the best score
					my (@best_seqs) = @{$score_to_seqs{$best_score}};
					foreach my $transc_id (@scoring_principals) {
						if ( grep( /^$transc_to_seq{$transc_id}$/, @best_seqs) ) {
							$report->{$transc_id} = 1;
						}
					}
					# save the rest of transcripts that pass the threshold comparing with the best
					for ( my $i=1; $i < scalar(@dec_scores); $i++ ) {
						if ($best_score - $dec_scores[$i] <= $min_lead) {
							my (@best_seqs) = @{$score_to_seqs{$dec_scores[$i]}};
							foreach my $transc_id (@scoring_principals) {
								if ( grep( /^$transc_to_seq{$transc_id}$/, @best_seqs) ) {
									$report->{$transc_id} = 1;
								}
							}
						}
					}
				}
			}
		}
	}

	if ( ! defined($report) ) {
		$report = $i_princ_list;
	}

	return $report;

} # end step_trifid

sub step_proteo($$)
{
	my ($i_princ_list, $s_scores) = @_;
	my ($report);

	# the minumum of scores is 2
	if ( $s_scores->{'max'} >= 2) {

		# get the list of transcrips and peptides (min peptides 2)
		my ($t_scores);
		while ( my ($score, $trans) = each(%{$s_scores->{'scores'}}) ) {
			foreach my $transc_id (@$trans) {
				if ( exists $i_princ_list->{$transc_id} ) {
					push(@{$t_scores->{$score}}, $transc_id);
				}	
			}
		}
		if ( defined($t_scores) ) { 
			# get the biggest score (num. peptides)
			my @dec_scores = sort { $b <=> $a } keys(%$t_scores);
			my $num_dec_scores = scalar(@dec_scores);
			# get the principal with the biggest num. peptides
			if ( $num_dec_scores >= 1) {
				# save the best proteins
				my $best_score = $dec_scores[0];
				foreach my $transc_id (@{$t_scores->{$best_score}}) {
					$report->{$transc_id} = 1;
				}
				# save the rest of transc that pass the threshold
				if ( $num_dec_scores >= 2 && ($best_score - $dec_scores[1] < $main::PROTEO_CUTOFF) ) {
					foreach my $transc_id (@{$t_scores->{$dec_scores[1]}}) {
						$report->{$transc_id} = 1;
					}
				}
			}
		}
		else {
			$report = $i_princ_list;
		}
	}
	else {
		$report = $i_princ_list;
	}

	return $report;

} # end step_proteo

sub step_ccds($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# preliminar report
	# discarding the transcripts are not protein coding (NMD): app_score is -1
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			if ( exists $nscores->{$transc_id} and defined $nscores->{$transc_id} and $nscores->{$transc_id}->{'appris'} != '-1' ) {
				if ( exists $princ->{'ccds'} ) { push(@{$princ_isof}, $princ) }				
			}
		}
	}	

	# print princ isoforms with CCDS ids
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
		foreach my $princ ( @{$princ_isof} ) {
			my ($transc_id) = $princ->{'id'};
			$report->{$transc_id} = 1;
		}				
	}
	else { $report = $i_princ_list }

	return $report;
		
} # end step_ccds

sub step_tsl($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# preliminar report
	# discarding the transcripts are not protein coding (NMD): app_score is -1
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			if ( exists $nscores->{$transc_id} and defined $nscores->{$transc_id} and $nscores->{$transc_id}->{'appris'} != '-1' ) {
				if ( exists $princ->{'tsl'} ) { push(@{$princ_isof}, $princ) }
			}
		}
	}
		
	# print princ isoforms with TSL(1)
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
		foreach my $princ ( @{$princ_isof} ) {
			my ($transc_id) = $princ->{'id'};
			if ( $princ->{'tsl'} == 1 ) { $report->{$transc_id} = 1 }
		}				
	}
	else { $report = $i_princ_list }

	return $report;
		
} # end step_tsl

sub step_ccds_eldest($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# preliminar report
	# discarding the transcripts are not protein coding (NMD): app_score is -1
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {			
			if ( exists $nscores->{$transc_id} and defined $nscores->{$transc_id} and $nscores->{$transc_id}->{'appris'} != '-1' ) {
				if ( exists $princ->{'ccds'} ) { push(@{$princ_isof}, $princ) }
			}
		}
	}	

	# print princ isoforms with eldests CCDS
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
		# the first elem is the eldests (a <=> b - ascending)
		my (@sort_transc) = sort { $a->{'ccds'} <=> $b->{'ccds'} } @{$princ_isof};
		my ($eldest_ccds) = $sort_transc[0]->{'ccds'};
		foreach my $princ ( @{$princ_isof} ) {
			my ($transc_id) = $princ->{'id'};
			my ($ccds) = $princ->{'ccds'};
			if ( abs($ccds - $eldest_ccds) < 10 ) {
				$report->{$transc_id} = 1;
			}			
		}				
	}
	else { $report = $i_princ_list }

	return $report;
		
} # end step_ccds_eldest

sub step_ccds_longest($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			
			if ( exists $princ->{'ccds'} ) { push(@{$princ_isof}, $princ) }
			
		}
	}	

	# print princ isoforms with longest CCDS
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {		
		# the first elem is the longest (b <=> a - descending)
		my (@sort_transc) = sort { $b->{'length'} <=> $a->{'length'} } @{$princ_isof};
		my ($longest_ccds) = $sort_transc[0]->{'length'};
		foreach my $princ ( @{$princ_isof} ) {
			my ($transc_id) = $princ->{'id'};
			my ($len) = $princ->{'length'};
			if ( abs($len - $longest_ccds) == 0 ) {
				$report->{$transc_id} = 1;
			}			
		}				
	}
	else { $report = $i_princ_list }

	return $report;
		
} # end step_ccds_longest

sub step_validated($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			if ( exists $princ->{'name'} ) {
				my ($name) = $princ->{'name'};
				if ( $name =~ /NM\_/ or $name =~ /\-0[0-9]*$/ ) { # Not automatic transcripts
					push(@{$princ_isof}, $princ);
				}
			}			
		}
	}
	
	# print princ isoforms for validated manually
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
		foreach my $princ ( @{$princ_isof} ) {
			my ($transc_id) = $princ->{'id'};
			$report->{$transc_id} = 1;
		}				
	}
	else { $report = $i_princ_list }

	return $report;
		
} # end step_validated

sub step_longest($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	my ($princ_isof_seqs);
	my ($princ_isof_scores);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			if ( exists $princ->{'length'} and exists $princ->{'seq'} ) {
				my ($len) = $princ->{'length'};
				my ($seq) = $princ->{'seq'};
				push(@{$princ_isof->{$len}}, $transc_id);
				$princ_isof_seqs->{$transc_id} = $seq;
			}			
			if ( exists $nscores->{$transc_id} and exists $nscores->{$transc_id}->{'appris'} ) {
				my ($sc) = $nscores->{$transc_id}->{'appris'};
				push(@{$princ_isof_scores->{$sc}}, $transc_id);
			}
		}
	}
	
	# print princ isoforms with longest seq
	if ( defined $princ_isof ) {
		# the first elem is the longest (b <=> a - descending)
		my (@sort_lengths) = sort { $b <=> $a } keys (%{$princ_isof});
		my ($longest) = $sort_lengths[0];
		if ( scalar(@{$princ_isof->{$longest}}) == 1 ) {
			my ($transc_id) = ($princ_isof->{$longest})->[0];
			$report->{$transc_id} = 1;
		}
		else {
			# check if they are the same seqs
			my ($equal) = 1;
			my ($transc_id_0) = ($princ_isof->{$longest})->[0];
			my ($seq_0) = $princ_isof_seqs->{$transc_id_0};
			my ($higher_sc_longest) = $nscores->{$transc_id_0}->{'appris'};
			for (my $i=1; $i < scalar(@{$princ_isof->{$longest}}); $i++) {
				my ($transc_id_i) = ($princ_isof->{$longest})->[$i];
				my ($seq_i) = $princ_isof_seqs->{$transc_id_i};
				my ($sc_i) = $nscores->{$transc_id_i}->{'appris'};
				if ( $seq_i ne $seq_0 ) { $equal = 0 }
				if ( $sc_i > $higher_sc_longest ) { $higher_sc_longest = $sc_i }
			}
			if ( $equal == 1 ) {
				foreach my $transc_id ( @{$princ_isof->{$longest}} ) {
					$report->{$transc_id} = 1;
				}
			}
			else {
				# if the seqs are different, get the one with higger appris-score. otherwise, all of them
				foreach my $transc_id ( @{$princ_isof->{$longest}} ) {
					my ($sc) = $nscores->{$transc_id}->{'appris'};
					if ( $sc == $higher_sc_longest ) {
						$report->{$transc_id} = 1;
					}
				}
			}
		}
	}
	else { $report = $i_princ_list }

	return $report;
		
} # end step_longest

sub step_smaller_id($$$)
{
	my ($i_princ_list, $isof_report, $nscores) = @_;
	my ($report);
	
	# sort the transc ids and selects the first one
	my (@princ_isof) = sort { $a cmp $b } keys(%{$i_princ_list});
	my ($transc_id) = $princ_isof[0];
	$report->{$transc_id} = 1;

	return $report;
		
} # end step_smaller_id

sub step_tags($$$$\$)
{
	my ($step, $scores, $princ_list, $isof_report, $ref_annots) = @_;
	
	# add tags (reliability) and annotations for
	my ($m) = 'appris';
	my ($label) = $METRIC_LABELS->{$m}->[0];
	my ($label_flag) = $METRIC_LABELS->{$m}->[2];
	foreach my $isof_rep ( @{$isof_report} ) {
		my ($transc_id) = $isof_rep->{'id'};
		my ($status,$annot);		
		if ( exists $princ_list->{$transc_id} and exists $isof_rep->{'principal'} ) { # principal isoforms
			if ( $step == 1 ) {
				$status = $OK_LABEL;
				$annot = 'PRINCIPAL:'.$step;
			}
			else {
				$status = $UNKNOWN_LABEL;
				$annot = 'PRINCIPAL:'.$step;
			}			
		}
		elsif ( !(exists $princ_list->{$transc_id}) and exists $isof_rep->{'principal'} ) { # alternative isoforms
			if ( defined $scores and exists $scores->{$transc_id} and exists $scores->{$transc_id}->{'score_vertebrate_signal'} and ($scores->{$transc_id}->{'score_vertebrate_signal'} >= 4.5) ) {
				$status = $UNKNOWN_LABEL;
				$annot = 'ALTERNATIVE:1';
			}
			else {
				$status = $UNKNOWN_LABEL;
				$annot = 'ALTERNATIVE:2';
			}			
		}
		elsif ( !(exists $princ_list->{$transc_id}) and !(exists $isof_rep->{'principal'}) ) { # minor isoforms
			$status = $NO_LABEL;
			$annot = 'MINOR';		
		}
		$$ref_annots->{$transc_id}->{$label} = $status;
		$$ref_annots->{$transc_id}->{$label_flag} = $annot;		
	}
		
} # step_tags

sub get_score_output($$$)
{
	my ($gene, $scores, $annots) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($gene_name) = $gene->external_name;
	my ($content) = '';
	
	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		my ($tranlation_id) = '-';
		my ($biotype) = '-';
		my ($flags) = '-';
		my ($flag_transl) = 'TRANSLATION';
		my ($ccds_id) = '-';
		my ($tsl) = '-';
		my ($transl_len) = 0;		
		my ($firestar_annot) = '-';
		my ($matador3d_annot) = '-';
		my ($corsair_annot) = '-';
		my ($spade_annot) = '-';
		my ($thump_annot) = '-';
		my ($crash_sp_annot) = '-';
		my ($crash_tp_annot) = '-';
		my ($proteo_annot) = '-';
		my ($inertia_annot) = '-';
		my ($appris_annot) = '-';
		my ($appris_relia) = '-';
		$biotype = $transcript->biotype if ( defined $transcript->biotype);
		$flags = $transcript->biotype if ( defined $transcript->biotype);
		if ( $transcript->translate and $transcript->translate->sequence ) {
			$tranlation_id = $transcript->translate->stable_id;
			
			$firestar_annot = $scores->{$transcript_id}->{'num_functional_residues'} if ( exists $scores->{$transcript_id}->{'num_functional_residues'} );
			$matador3d_annot = $scores->{$transcript_id}->{'score_homologous_structure'} if ( exists $scores->{$transcript_id}->{'score_homologous_structure'} );
			$corsair_annot = $scores->{$transcript_id}->{'score_vertebrate_signal'} if ( exists $scores->{$transcript_id}->{'score_vertebrate_signal'} );		
			$spade_annot = $scores->{$transcript_id}->{'score_domain_signal'} if ( exists $scores->{$transcript_id}->{'score_domain_signal'} );
			$thump_annot = $scores->{$transcript_id}->{'score_transmembrane_signal'} if ( exists $scores->{$transcript_id}->{'score_transmembrane_signal'} );
			$crash_sp_annot = $scores->{$transcript_id}->{'score_peptide_signal'} if ( exists $scores->{$transcript_id}->{'score_peptide_signal'} );
			$crash_tp_annot = $scores->{$transcript_id}->{'score_mitochondrial_signal'} if ( exists $scores->{$transcript_id}->{'score_mitochondrial_signal'} );
			$inertia_annot = $scores->{$transcript_id}->{'num_unusual_exons'} if ( exists $scores->{$transcript_id}->{'num_unusual_exons'} );
			$proteo_annot = $scores->{$transcript_id}->{'num_peptides'} if ( exists $scores->{$transcript_id}->{'num_peptides'} );
			$appris_annot = $scores->{$transcript_id}->{'score_principal_isoform'} if ( exists $scores->{$transcript_id}->{'score_principal_isoform'} );
						
			$appris_relia = $annots->{$transcript_id}->{'reliability'} if ( exists $annots->{$transcript_id}->{'reliability'} );
			
			$transl_len = length($transcript->translate->sequence);
			
			if ( $transcript->xref_identify ) {
				foreach my $xref (@{$transcript->xref_identify}) {
					if ( $xref->dbname eq 'CCDS') {
						$ccds_id = $xref->id;
					}
				}
			}		
			my ($no_codons) = 'start/stop';
			if ( $transcript->translate->codons ) {
				my ($start_found) = 0;
				my ($stop_found) = 0;
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( $codon->type eq 'start' ) {
						$start_found = 1;
					} elsif ( $codon->type eq 'stop' ) {
						$stop_found = 1;
					}
				}
				if ( $start_found && $stop_found ) {
					$no_codons = '-';
				} elsif ( ! $start_found ) {
					$no_codons = 'start';
				} elsif ( ! $stop_found ) {
					$no_codons = 'stop';
				}
			}
			$tsl = $transcript->tsl if ( defined $transcript->tsl);
			if ( defined $transcript->tag and $transcript->tag =~ /readthrough_transcript/ ) { $flags .= ',RT' }		
						
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$tranlation_id."\t".
						$flag_transl."\t".						
						$flags."\t".
						$no_codons."\t".
						$ccds_id."\t".
						$tsl."\t".
						$transl_len."\t".
						$firestar_annot."\t".
						$matador3d_annot."\t".
						$corsair_annot."\t".
						$spade_annot."\t".
						$thump_annot."\t".
						$crash_sp_annot.",".$crash_tp_annot."\t".
						$inertia_annot."\t".
						$proteo_annot."\t".
						$appris_annot."\t".
						$appris_relia."\n";
		}
		else {
			$flag_transl = 'NO_TRANSLATION';
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$tranlation_id."\t".
						$flag_transl."\t".
						$flags."\n";
		}
	}
	
	return $content;
	
} # End get_score_output

sub get_nscore_output($$)
{
	my ($gene, $nscores) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($gene_name) = $gene->external_name;
	my ($content) = '';

	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		my ($tranlation_id) = '-';
		my ($flags) = '-';
		my ($ccds_id) = '-';
		my ($tsl) = '-';
		my ($transl_len) = 0;		
		my ($firestar_annot) = '-';
		my ($matador3d_annot) = '-';
		my ($corsair_annot) = '-';
		my ($spade_annot) = '-';
		my ($thump_annot) = '-';
		my ($crash_annot) = '-';
		my ($proteo_annot) = '-';
		my ($inertia_annot) = '-';
		my ($appris_annot) = '-';
		$flags = $transcript->biotype if ( defined $transcript->biotype);		
		if ( $transcript->translate and $transcript->translate->sequence ) {
			$tranlation_id = $transcript->translate->stable_id;
						
			$firestar_annot = $nscores->{$transcript_id}->{'firestar'} if ( exists $nscores->{$transcript_id}->{'firestar'} );
			$matador3d_annot = $nscores->{$transcript_id}->{'matador3d'} if ( exists $nscores->{$transcript_id}->{'matador3d'} );
			$corsair_annot = $nscores->{$transcript_id}->{'corsair'} if ( exists $nscores->{$transcript_id}->{'corsair'} );		
			$spade_annot = $nscores->{$transcript_id}->{'spade'} if ( exists $nscores->{$transcript_id}->{'spade'} );
			$thump_annot = $nscores->{$transcript_id}->{'thump'} if ( exists $nscores->{$transcript_id}->{'thump'} );
			$crash_annot = $nscores->{$transcript_id}->{'crash'} if ( exists $nscores->{$transcript_id}->{'crash'} );
			$inertia_annot = $nscores->{$transcript_id}->{'inertia'} if ( exists $nscores->{$transcript_id}->{'inertia'} );
			$proteo_annot = $nscores->{$transcript_id}->{'proteo'} if ( exists $nscores->{$transcript_id}->{'proteo'} );
			$appris_annot = $nscores->{$transcript_id}->{'appris'} if ( exists $nscores->{$transcript_id}->{'appris'} );
						
			$transl_len = length($transcript->translate->sequence);
			
			if ( $transcript->xref_identify ) {
				foreach my $xref (@{$transcript->xref_identify}) {
					if ( $xref->dbname eq 'CCDS') {
						$ccds_id = $xref->id;
					}
				}
			}		
			my ($no_codons) = 'start/stop';
			if ( $transcript->translate->codons ) {
				my ($start_found) = 0;
				my ($stop_found) = 0;
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( $codon->type eq 'start' ) {
						$start_found = 1;
					} elsif ( $codon->type eq 'stop' ) {
						$stop_found = 1;
					}
				}
				if ( $start_found && $stop_found ) {
					$no_codons = '-';
				} elsif ( ! $start_found ) {
					$no_codons = 'start';
				} elsif ( ! $stop_found ) {
					$no_codons = 'stop';
				}
			}
			$tsl = $transcript->tsl if ( defined $transcript->tsl);
			if ( defined $transcript->tag and $transcript->tag =~ /readthrough_transcript/ ) { $flags .= ',RT' }
						
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$tranlation_id."\t".
						$flags."\t".
						$no_codons."\t".
						$ccds_id."\t".
						$tsl."\t".
						$transl_len."\t".
						$firestar_annot."\t".
						$matador3d_annot."\t".
						$corsair_annot."\t".
						$spade_annot."\t".
						$thump_annot."\t".
						$crash_annot."\t".
						$inertia_annot."\t".
						$proteo_annot."\t".
						$appris_annot."\n";
		}
	}
	
	return $content;
	
} # End get_nscore_output

sub get_label_output($$)
{
	my ($gene, $annots) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($gene_name) = $gene->external_name;
	my ($content) = '';

	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		my ($tranlation_id) = '-';
		my ($flags) = '-';
		my ($flag_transl) = 'TRANSLATION';
		my ($ccds_id) = '-';
		my ($tsl) = '-';
		my ($transl_len) = 0;		
		my ($firestar_annot) = '-';
		my ($matador3d_annot) = '-';
		my ($corsair_annot) = '-';
		my ($spade_annot) = '-';
		my ($thump_annot) = '-';
		my ($crash_sp_annot) = '-';
		my ($crash_tp_annot) = '-';
		my ($inertia_annot) = '-';
		my ($proteo_annot) = '-';
		my ($appris_annot) = '-';
		my ($appris_relia) = '-';
		$flags = $transcript->biotype if ( defined $transcript->biotype);
		if ( $transcript->translate and $transcript->translate->sequence ) {
			$tranlation_id = $transcript->translate->stable_id;
			
			$firestar_annot = $annots->{$transcript_id}->{'functional_residue'} if ( exists $annots->{$transcript_id}->{'functional_residue'} );
			$matador3d_annot = $annots->{$transcript_id}->{'conservation_structure'}if ( exists $annots->{$transcript_id}->{'conservation_structure'} );
			$corsair_annot = $annots->{$transcript_id}->{'vertebrate_signal'} if ( exists $annots->{$transcript_id}->{'vertebrate_signal'} );
			$spade_annot = $annots->{$transcript_id}->{'domain_signal'} if ( exists $annots->{$transcript_id}->{'domain_signal'} );
			$thump_annot = $annots->{$transcript_id}->{'transmembrane_signal'} if ( exists $annots->{$transcript_id}->{'transmembrane_signal'} );
			$crash_sp_annot = $annots->{$transcript_id}->{'peptide_signal'} if ( exists $annots->{$transcript_id}->{'peptide_signal'} );
			$crash_tp_annot = $annots->{$transcript_id}->{'mitochondrial_signal'} if ( exists $annots->{$transcript_id}->{'mitochondrial_signal'} );
			$inertia_annot = $annots->{$transcript_id}->{'unusual_evolution'} if ( exists $annots->{$transcript_id}->{'unusual_evolution'} );
			$proteo_annot = $annots->{$transcript_id}->{'peptide_evidence'} if ( exists $annots->{$transcript_id}->{'peptide_evidence'} );			
			
			$appris_annot = $annots->{$transcript_id}->{'principal_isoform'} if ( exists $annots->{$transcript_id}->{'principal_isoform'} );
			$appris_relia = $annots->{$transcript_id}->{'reliability'} if ( exists $annots->{$transcript_id}->{'reliability'} );
			
			$transl_len = length($transcript->translate->sequence);
			
			if ( $transcript->xref_identify ) {
				foreach my $xref (@{$transcript->xref_identify}) {
					if ( $xref->dbname eq 'CCDS') {
						$ccds_id = $xref->id;
					}
				}
			}		
			my ($no_codons) = 'start/stop';
			if ( $transcript->translate->codons ) {
				my ($start_found) = 0;
				my ($stop_found) = 0;
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( $codon->type eq 'start' ) {
						$start_found = 1;
					} elsif ( $codon->type eq 'stop' ) {
						$stop_found = 1;
					}
				}
				if ( $start_found && $stop_found ) {
					$no_codons = '-';
				} elsif ( ! $start_found ) {
					$no_codons = 'start';
				} elsif ( ! $stop_found ) {
					$no_codons = 'stop';
				}
			}
			$tsl = $transcript->tsl if ( defined $transcript->tsl);
			if ( defined $transcript->tag and $transcript->tag =~ /readthrough_transcript/ ) { $flags .= ',RT' }
						
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$tranlation_id."\t".
						$flag_transl."\t".						
						$flags."\t".
						$no_codons."\t".
						$ccds_id."\t".
						$tsl."\t".
						$transl_len."\t".				
						$firestar_annot."\t".
						$matador3d_annot."\t".
						$corsair_annot."\t".
						$spade_annot."\t".
						$thump_annot."\t".
						$crash_sp_annot.",".$crash_tp_annot."\t".
						$inertia_annot."\t".
						$proteo_annot."\t".
						$appris_annot."\t".
						$appris_relia."\n";
		}
		else {
			$flag_transl = 'NO_TRANSLATION';
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$tranlation_id."\t".
						$flag_transl."\t".
						$flags."\n";
		}
	}
	
	return $content;
	
} # End get_label_output

sub get_firestar_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);
	
	# the maximim number of residues has to pass the threshold. Otherwise, we have not a decision
	my ($max) = $scores->{'max'};
	if ( $max <= $main::FIRESTAR_MINRES ) {
		foreach my $transcript (@{$gene->transcripts}) {
			if ( $transcript->translate and $transcript->translate->sequence ) {	
				my ($transc_id) = $transcript->stable_id;
				$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
			}
		}
	}
	else {
		# sort by descending order
		my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}}); 
		for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
			if ( $i == 0 ) { # biggest score
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@unknows,$transc_id);
				}
			}
			else {
				if ( $sorted_scores[0] - $sorted_scores[$i] <= $main::FIRESTAR_CUTOFF ) {
					my ($sc) = $sorted_scores[$i];
					foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
						push(@unknows,$transc_id);
					}
				}
				else {
					my ($sc) = $sorted_scores[$i];
					foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
						push(@nos,$transc_id);
					}
				}
			}
		}
		if ( scalar(@unknows) == 1 ) {
			foreach my $transc_id (@unknows) {
				$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
			}
		}
		else {
			foreach my $transc_id (@unknows) {
				$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
			}
		}
		foreach my $transc_id (@nos) {
			$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
		}		
	}
		
} # end get_firestar_annots

sub get_matador3d_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);
	
	# sort by descending order
	my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}}); 
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		if ( $i == 0 ) { # biggest score
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@unknows,$transc_id);
			}
		}
		else {
			if ( $sorted_scores[0] - $sorted_scores[$i] < $main::MATADOR3D_CUTOFF ) {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@unknows,$transc_id);
				}
			}
			else {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@nos,$transc_id);
				}
			}
		}
	}
	if ( scalar(@unknows) == 1 ) {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
		}
	}
	else {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
		}
	}
	foreach my $transc_id (@nos) {
		$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
	}
		
} # end get_matador3d_annots

sub get_matador3d2_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);
	
	# sort by descending order
	my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}}); 
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		if ( $i == 0 ) { # biggest score
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@unknows,$transc_id);
			}
		}
		else {
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $main::MATADOR3D2_CUTOFF ) {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@unknows,$transc_id);
				}
			}
			else {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@nos,$transc_id);
				}
			}
		}
	}
	if ( scalar(@unknows) == 1 ) {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
		}
	}
	else {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
		}
	}
	foreach my $transc_id (@nos) {
		$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
	}
		
} # end get_matador3d2_annots

sub get_spade_integrity_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);

	# sort by descending order
	my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}});
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		if ( $i == 0 ) { # biggest score
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@unknows,$transc_id);
			}
		}
		else {
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $main::SPADE_INTEGRITY_CUTOFF ) {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@unknows,$transc_id);
				}
			}
			else {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@nos,$transc_id);
				}
			}
		}
	}
	if ( scalar(@unknows) == 1 ) {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
		}
	}
	else {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
		}
	}
	foreach my $transc_id (@nos) {
		$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
	}

} # end get_spade_integrity_annots

sub get_spade_bitscore_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);

	# sort by descending order
	my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}});
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		if ( $i == 0 ) { # biggest score
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@unknows,$transc_id);
			}
		}
		else {
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $main::SPADE_CUTOFF ) {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@unknows,$transc_id);
				}
			}
			else {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@nos,$transc_id);
				}
			}
		}
	}
	if ( scalar(@unknows) == 1 ) {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
		}
	}
	else {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
		}
	}
	foreach my $transc_id (@nos) {
		$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
	}

} # end get_spade_bitscore_annots

sub get_corsair_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);
	
	# sort by descending order
	my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}}); 
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		if ( $i == 0 ) { # biggest score
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@unknows,$transc_id);
			}
		}
		else {
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $main::CORSAIR_CUTOFF ) {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@unknows,$transc_id);
				}
			}
			else {
				my ($sc) = $sorted_scores[$i];
				foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
					push(@nos,$transc_id);
				}
			}
		}
	}
	if ( scalar(@unknows) == 1 ) {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
		}
	}
	else {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
		}
	}
	foreach my $transc_id (@nos) {
		$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
	}
		
} # end get_corsair_annots

sub get_thump_annots($$$\$)
{
	my ($gene, $scores, $annot_label, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);
	
	# sort by descending order
	my (@sorted_scores) = sort { $b <=> $a } keys (%{$scores->{'scores'}}); 
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		if ( $i == 0 ) { # biggest score
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@unknows,$transc_id);
			}
		}
		else {
			my ($sc) = $sorted_scores[$i];
			foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
				push(@nos,$transc_id);
			}
		}
	}
	if ( scalar(@unknows) == 1 ) {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $OK_LABEL;
		}
	}
	else {
		foreach my $transc_id (@unknows) {
			$$ref_annots->{$transc_id}->{$annot_label} = $UNKNOWN_LABEL;
		}
	}
	foreach my $transc_id (@nos) {
		$$ref_annots->{$transc_id}->{$annot_label} = $NO_LABEL;
	}
		
} # end get_thump_annots

sub get_crash_annots($$$$\$)
{
	my ($gene, $scores, $annot_label_sp, $annot_label_tp, $ref_annots) = @_;
	my ($annots);
	my (@unknows);
	my (@nos);
	
	# sort by descending order
	my (@sorted_scores) = sort { $b cmp $a } keys (%{$scores->{'scores'}}); 
	for ( my $i=0; $i < scalar(@sorted_scores); $i++ ) {
		my ($sc) = $sorted_scores[$i];
		my ($sc_sp,$sc_tp) = split(',', $sc);
		my ($label_sp, $label_tp);
		if ($sc_sp >= 2) {
			$label_sp = $OK_LABEL;
		} elsif (($sc_sp == 0) or ($sc_sp == 1)) {
			$label_sp = $UNKNOWN_LABEL;
		} elsif ($sc_sp <= -1) {
			$label_sp = $NO_LABEL;
		}
		if ($sc_tp >= 2) {
			$label_tp = $OK_LABEL;
		} elsif (($sc_tp == 0) or ($sc_tp == 1)) {
			$label_tp = $UNKNOWN_LABEL;
		} elsif ($sc_tp <= -1) {
			$label_tp = $NO_LABEL;
		}
		foreach my $transc_id (@{$scores->{'scores'}->{$sc}}) {
			$$ref_annots->{$transc_id}->{$annot_label_sp} = $label_sp;
			$$ref_annots->{$transc_id}->{$annot_label_tp} = $label_tp;
		}
	}
		
} # end get_crash_annots

# INERTIA: get the number of unusual exons
sub get_num_unusual_exons($)
{
	my ($analysis) = @_;
	
	my (@types) = ('inertia','maf', 'prank', 'kalign', 'compara');	
	my ($num_un_exons);
	
	foreach my $type (@types) {	
					
		$num_un_exons->{$type} = 0;
		my ($regions);

		if ( ($type eq 'inertia') and $analysis->regions ) { # consensus
			$regions = $analysis->regions;
		}
		elsif ( ($type eq 'maf') and $analysis->mafft_alignment and $analysis->mafft_alignment->regions ) {
			$regions = $analysis->mafft_alignment->regions;
		}
		elsif ( ($type eq 'prank') and $analysis->prank_alignment and $analysis->prank_alignment->regions ) {
			$regions = $analysis->prank_alignment->regions;
		}
		elsif ( ($type eq 'kalign') and $analysis->kalign_alignment and $analysis->kalign_alignment->regions ) {
			$regions = $analysis->kalign_alignment->regions;
		}
		elsif ( ($type eq 'compara') and $analysis->compara_alignment and $analysis->compara_alignment->regions ) {
			$regions = $analysis->compara_alignment->regions;
		}

		foreach my $residue (@{$regions}) {
			if ( $residue->unusual_evolution and defined $residue->unusual_evolution and ($residue->unusual_evolution eq $NO_LABEL) ) {
				$num_un_exons->{$type} += 1
			}		
		}
	}
	
	return $num_un_exons;
	
} # End get_num_unusual_exons

1;
