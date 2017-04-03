package appris;

use strict;
use Data::Dumper;

use APPRIS::Utils::Exception qw( info throw warning deprecate );

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	get_method_scores
	get_method_annots
	get_final_scores
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

	$METHOD_LABELS
	$METHOD_WEIGHTED
);

$OK_LABEL				= 'YES';
$UNKNOWN_LABEL			= 'UNKNOWN';
$NO_LABEL				= 'NO';

$METHOD_LABELS = {
	'firestar'	=> [ 'functional_residue',			'num_functional_residues'							],
	'matador3d'	=> [ 'conservation_structure',		'score_homologous_structure'						],
	'corsair'	=> [ 'vertebrate_signal',			'score_vertebrate_signal'							],
	'spade'		=> [ 'domain_signal',				'score_domain_signal'								],
	'thump'		=> [ 'transmembrane_signal',		'score_transmembrane_signal'						],
	'crash'		=> [ 'peptide_signal',				'score_peptide_signal',
					 'mitochondrial_signal',		'score_mitochondrial_signal'						],
	'inertia'	=> [ 'unusual_evolution',			'num_unusual_exons'									],
	'proteo'	=> [ 'peptide_evidence',			'num_peptides'										],
	'appris'	=> [ 'principal_isoform',			'score_principal_isoform',			'reliability'	]
};
$METHOD_WEIGHTED = {
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
	'spade'		=> 6,
	'corsair'	=> [{
					  'max'    => 3,
				  	  'weight' => 1.5	
					},
					{
					  'max'    => 4,
				  	  'weight' => 2	
					},
					{
					  'max'    => 5,
					  'weight' => 3	
					}],
	'thump'		=> 0,
	'crash'		=> 0,
	'inertia'	=> 0,
	'proteo'	=> 0,
};


# get the main functional isoform from methods of appris
sub get_method_scores($$)
{
	my ($gene, $reports) = @_;
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
			
			my ($m) = 'firestar';
			if ( $result and $result->analysis and $result->analysis->firestar ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			$m = 'matador3d';
			if ( $result and $result->analysis and $result->analysis->matador3d ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			# get corsair
			$m = 'corsair';
			if ( $result and $result->analysis and $result->analysis->corsair ) {
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->corsair;
				if ( defined $analysis->score ) {
					my ($sc) = $analysis->score;
					# if variant has 'start/stop codon not found', the score is 0
					if ( $transcript->translate->codons ) {
						my ($codons) = '';
						foreach my $codon (@{$transcript->translate->codons}) {
							if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) { $codons .= $codon->type.',' }
						}
						unless ( $codons =~ /start/ and $codons =~ /stop/ ) { $sc = 0 } 
					}
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
			$m = 'spade';
			if ( $result and $result->analysis and $result->analysis->spade ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->spade;
				if ( defined $analysis->bitscore ) {
					my ($sc) = $analysis->bitscore;
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
			$m = 'thump';
			if ( $result and $result->analysis and $result->analysis->thump ) {
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			# crash: signalp + targetp
			$m = 'crash';
			if ( $result and $result->analysis and $result->analysis->crash ) {
				my ($annot_sc_sp) = $METHOD_LABELS->{$m}->[1];
				my ($annot_sc_tp) = $METHOD_LABELS->{$m}->[3];	
				my ($analysis) = $result->analysis->crash;
				if ( defined $analysis->sp_score and defined $analysis->tp_score ) {
					my ($sc) = $analysis->sp_score.','.$analysis->tp_score;
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc_sp} = $analysis->sp_score;
					$scores->{$transcript_id}->{$annot_sc_tp} = $analysis->tp_score;
				}
			}
			$m = 'inertia';
			if ( $result and $result->analysis and $result->analysis->inertia ) {
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];				
				my ($analysis) = $result->analysis->inertia;
				if ( defined $analysis->regions ) {
					 my ($sc) = get_num_unusual_exons($analysis);
					 push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					 $scores->{$transcript_id}->{$annot_sc} = $sc;
				}
			}
			# (WE DON'T USE PROTEO FOR APPRIS DECISION)
			$m = 'proteo';
			if ( $result and $result->analysis and $result->analysis->proteo ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
	
	return ($scores, $s_scores);
	
} # End get_method_scores

# get the main functional isoform from methods of appris
sub get_method_scores_from_appris_rst($$)
{
	my ($gene, $reports) = @_;
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
			
			my ($m) = 'firestar';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->functional_residues_score ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			$m = 'matador3d';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->homologous_structure_score ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			$m = 'corsair';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->vertebrate_conservation_score ) {
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			$m = 'spade';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->domain_score ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			$m = 'thump';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->transmembrane_helices_score ) {
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
			# crash: signalp + targetp
			$m = 'crash';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->peptide_score and defined $result->analysis->appris->mitochondrial_score ) {
				my ($annot_sc_sp) = $METHOD_LABELS->{$m}->[1];
				my ($annot_sc_tp) = $METHOD_LABELS->{$m}->[3];
				my ($sc_sp) = $result->analysis->appris->peptide_score;
				my ($sc_tp) = $result->analysis->appris->mitochondrial_score;
				if ( $sc_sp ne '-' and $sc_tp ne '-' ) {
					my ($sc) = $sc_sp.','.$sc_tp;					
					push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
					$scores->{$transcript_id}->{$annot_sc_sp} = $sc_sp;
					$scores->{$transcript_id}->{$annot_sc_tp} = $sc_tp;
				}
			}
			#$m = 'inertia';
			#if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->peptide_score and defined $result->analysis->appris->mitochondrial_score ) {
			#	my ($annot_sc) = $METHOD_LABELS->{$m}->[1];				
			#	my ($analysis) = $result->analysis->inertia;
			#	if ( defined $analysis->regions ) {
			#		 my ($sc) = get_num_unusual_exons($analysis);
			#		 push(@{$s_scores->{$m}->{'scores'}->{$sc}}, $transcript_id);
			#		 $scores->{$transcript_id}->{$annot_sc} = $sc;
			#	}
			#}
			# (WE DON'T USE PROTEO FOR APPRIS DECISION)
			$m = 'proteo';
			if ( $result and $result->analysis and $result->analysis->appris and defined $result->analysis->appris->peptide_evidence_score ) {				
				my ($annot_sc) = $METHOD_LABELS->{$m}->[1];
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
	
	return ($scores, $s_scores);
	
} # End get_method_scores_from_appris_rst

# get the main functional isoform from methods of appris
sub get_method_annots($$)
{
	my ($gene, $scores) = @_;	
	my ($stable_id) = $gene->stable_id;
	my ($annots);
		
	# get annotations for each method (or group of method) -------------------
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			
			my ($m) = 'firestar';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_firestar_annots($gene, $scores->{$m}, $annot_label, \$annots);
			}
			$m = 'matador3d';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_matador3d_annots($gene, $scores->{$m}, $annot_label, \$annots);
			}
			$m = 'spade';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_spade_annots($gene, $scores->{$m}, $annot_label, \$annots);
			}
			$m = 'corsair';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_corsair_annots($gene, $scores->{$m}, $annot_label, \$annots);
			}
			$m = 'thump';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_thump_annots($gene, $scores->{$m}, $annot_label, \$annots);
			}
			$m = 'crash';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label_sp) = $METHOD_LABELS->{$m}->[0];
				my ($annot_label_tp) = $METHOD_LABELS->{$m}->[2];
				get_crash_annots($gene, $scores->{$m}, $annot_label_sp, $annot_label_tp, \$annots);
			}
		}
	}
	
	return $annots;
	
} # End get_method_annots

# get gene scores of methods for appris
sub get_final_scores($$\$\$)
{
	my ($gene, $annots, $ref_scores, $ref_s_scores) = @_;
	my ($nscores);
	
	# obtain normalize scores (weighted normalize score) foreach method
	my ($max_appris_score) = 0;
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($appris_score) = 0;
			foreach my $method ( keys(%{$METHOD_WEIGHTED}) ) {
				my ($n_sc) = 0;
				my ($sc) = 0;
				my ($max) = (exists $$ref_s_scores->{$method} and exists $$ref_s_scores->{$method}->{'max'})? $$ref_s_scores->{$method}->{'max'} : 0;
				my ($min) = (exists $$ref_s_scores->{$method} and exists $$ref_s_scores->{$method}->{'min'})? $$ref_s_scores->{$method}->{'min'} : 0;
				my ($label) = $METHOD_LABELS->{$method}->[1];
				
				# create normalize scores
				if ( exists $$ref_scores->{$transcript_id}->{$label} ) {
					$sc = $$ref_scores->{$transcript_id}->{$label};
					if ( $max != 0 and ($max - $min != 0) ) { $n_sc = $sc/$max } #$n_sc = ($sc - $min)/($max - $min)
					else { $n_sc = 0 }
				}
				if ( $n_sc < 0 ) { $n_sc = 0 }
				$n_sc = sprintf("%.3f",$n_sc);
				$nscores->{$transcript_id}->{$method} = $n_sc;
				
				# apply weights to normalize scores
				if ( $method eq 'firestar' ) {
					if ( $max >= $METHOD_WEIGHTED->{$method}->[4]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[4]->{'weight'}*$n_sc;
					}
					elsif ( $max >= $METHOD_WEIGHTED->{$method}->[3]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[3]->{'weight'}*$n_sc;
					}
					elsif ( $max >= $METHOD_WEIGHTED->{$method}->[2]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[2]->{'weight'}*$n_sc;
					}
					elsif ( $max >= $METHOD_WEIGHTED->{$method}->[1]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[1]->{'weight'}*$n_sc;
					}
					elsif ( $max >= $METHOD_WEIGHTED->{$method}->[0]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[0]->{'weight'}*$n_sc;
					}
				}
				elsif ( $method eq 'corsair' ) {
					if ( $max >= $METHOD_WEIGHTED->{$method}->[2]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[2]->{'weight'}*$n_sc;
					}
					elsif ( $max >= $METHOD_WEIGHTED->{$method}->[1]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[1]->{'weight'}*$n_sc;
					}
					elsif ( $max >= $METHOD_WEIGHTED->{$method}->[0]->{'max'} ) {
						$appris_score += $METHOD_WEIGHTED->{$method}->[0]->{'weight'}*$n_sc;
					}
				}
				else {
					$appris_score += $METHOD_WEIGHTED->{$method}*$n_sc;
				}
			}
			
			# filter by biotype
			if ( $transcript->biotype and ($transcript->biotype eq 'nonsense_mediated_decay') ) {
				$appris_score = -1;
			}
			# filter by readthrough_transcript
			if ( $transcript->tag and $transcript->tag =~ /readthrough_transcript/ ) {
				$appris_score = -1;
			}
						
			# save appris score and normalize score
			my ($m) = 'appris';
			my ($label) = $METHOD_LABELS->{$m}->[1];			
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
				my ($label) = $METHOD_LABELS->{$m}->[1];
				my ($sc) = $$ref_scores->{$transcript_id}->{$label};
				$n_sc = ($sc - $min)/($max - $min);
				$n_sc = sprintf("%.3f",$n_sc);
			}
			$nscores->{$transcript_id}->{$m} = $n_sc;			
		}
	}
		
	return ($nscores);
	
} # End get_final_scores

# get the final annotation
sub get_final_annotations($$$$$)
{
	my ($gene, $scores, $s_scores, $nscores, $annots) = @_;
	my ($method) = 'appris';
	my ($annotations);
	my ($tag) = 0;
		
	# scan transcripts sorted by appris score
	if ( defined $s_scores and exists $s_scores->{$method} and exists $s_scores->{$method}->{'scores'} and scalar(keys(%{$s_scores->{$method}->{'scores'}})) > 0 )
	{
		# 1. acquire the dominant transcripts from appris score.
		# They have to pass the cutoff to be added into "principal" list
		my ($princ_list, $isof_report) = step_appris($gene, $s_scores->{$method}, $annots);
		warning("GENE_2: \n".Dumper($gene)."\n");
		warning("PRINC_LIST_1: \n".Dumper($princ_list)."\n");
		warning("PRINC_REP_1: \n".Dumper($isof_report)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 1;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}

		# 2. from preserved transcript, we keep transcripts that they have got CCDS
		$princ_list = step_ccds($princ_list, $isof_report, $nscores);
		warning("PRINC_LIST_2: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 2;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}

		# 3_1. from preserved transcript, we keep transcripts that they have got eldest CCDS
		$princ_list = step_ccds_eldest($princ_list, $isof_report, $nscores);
		warning("PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 3;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}

		# 3_2. from preserved transcript, we keep transcripts that they have got TSL1
		$princ_list = step_tsl($princ_list, $isof_report, $nscores);
		warning("PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 3;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}

		# 4. from preserved transcript, we keep transcripts that they have got longest seq with CCDS
		$princ_list = step_ccds_longest($princ_list, $isof_report, $nscores);
		warning("PRINC_LIST_4: \n".Dumper($princ_list)."\n");
		warning("ISOF_REPORT_4: \n".Dumper($isof_report)."\n");
		warning("NSCORES_4: \n".Dumper($nscores)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 4;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}
		
		# 5. from preserved transcript, we keep transcripts that they have been validated manually
		$princ_list = step_validated($princ_list, $isof_report, $nscores);
		warning("PRINC_LIST_5_val: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 5;
			step_tags($tag, $scores, $princ_list, $isof_report, \$annots);
			return $tag;
		}
		
		# 5. from preserved transcript, we keep transcripts that they have got longest seq
		$princ_list = step_longest($princ_list, $isof_report, $nscores);
		warning("PRINC_LIST_5: \n".Dumper($princ_list)."\n");
		$annotations = step_tags(5, $scores, $princ_list, $isof_report, \$annots);		
	}
	
	return $tag;
	
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

sub step_appris($$$)
{
	my ($gene, $s_scores, $annots) = @_;
	my ($gene_id) = $gene->stable_id;
	my ($princ_list, $isof_report);
	my ($ccds_ids);
	
	# scan the transcripts from the sorted APPRIS scores
	my ($highest_score) = $s_scores->{'max'};
	my (@sorted_ap_scores) = sort { $b <=> $a } keys (%{$s_scores->{'scores'}});			
	for ( my $i = 0; $i < scalar(@sorted_ap_scores); $i++ ) {
		
		my ($ap_score) = $sorted_ap_scores[$i];		
		foreach my $transc_id (@{$s_scores->{'scores'}->{$ap_score}}) {
			my ($index) = $gene->{'_index_transcripts'}->{$transc_id};
			my ($transcript) = $gene->transcripts->[$index];
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
							# Only print warning message if CCDS is duplicated.
							unless ( exists $ccds_ids->{$transl_seq} ) {
								$ccds_ids->{$transl_seq} = $ccds_id;
							}
							else {
								if ( $ccds_ids->{$transl_seq} ne $ccds_id ) {
									my ($old_ccds_id) = $ccds_ids->{$transl_seq};
									my ($new_ccds_id) = $ccds_id;
									if ( $new_ccds_id < $old_ccds_id ) {
#										warning("$gene_id has duplicate CCDS ids for the same sequence. We change to the eldest id: CCDS$old_ccds_id -> CCDS$new_ccds_id\n");
									}
									else {
#										warning("$gene_id has duplicate CCDS ids for the same sequence. We keep the eldest id: CCDS$new_ccds_id -> CCDS$old_ccds_id \n");
									}
								}
							}						
						}
					}
				}
			}			
			# save TSL(1) annot
			if ( $transcript->tsl and $transcript->tsl eq '1' ) {
				$transc_rep->{'tsl'} = 1;
			}
					
			# APPRIS says could be a principal when:
			# 1.	the Core of pipeline has not rejected the transcript (IT DOES NOT APPLY YET)
			# APPRIS rejected a transcript when:
			# 2.	it is a NMD: app_score is -1. Exception: we accept the last terms when the transcript is unique
			#
			my ($core_flag) = 1;
			my ($unique_transc) = ( scalar(keys(%{$s_scores->{'scores'}})) == 1 and scalar(@{$s_scores->{'scores'}->{$sorted_ap_scores[0]}}) >= 1 ) ? 1 : 0;
			if ( $core_flag == 1 ) {
				if ( ( ($highest_score - $ap_score) <=  $main::APPRIS_CUTOFF) and ( ($ap_score >= 0) or ($unique_transc == 1) ) ) {
					$transc_rep->{'principal'} = 1;
					$princ_list->{$transc_id} = 1;
				}				
			}
			
			push(@{$isof_report}, $transc_rep);

		}			
	}
	return ($princ_list, $isof_report);
	
} # end step_appris

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

sub step_tags($$$$\$)
{
	my ($step, $scores, $princ_list, $isof_report, $ref_annots) = @_;
	
	# add tags (reliability) and annotations for
	my ($m) = 'appris';
	my ($label) = $METHOD_LABELS->{$m}->[0];
	my ($label_flag) = $METHOD_LABELS->{$m}->[2];	
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
		my ($no_codons) = '-';
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
			if ( $transcript->translate->codons ) {
				my ($aux_codons) = '';
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) {
						$aux_codons .= $codon->type.',';							
					}
				}
				$no_codons = 'start/' unless ( $aux_codons =~ /start/ );
				$no_codons .= 'stop' unless ( $aux_codons =~ /stop/ );
				$no_codons =~ s/^\-// if ($no_codons ne '-');
				$no_codons =~ s/\/$// if ($no_codons ne '-');
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
		my ($no_codons) = '-';
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
			if ( $transcript->translate->codons ) {
				my ($aux_codons) = '';
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) {
						$aux_codons .= $codon->type.',';							
					}
				}
				$no_codons = 'start/' unless ( $aux_codons =~ /start/ );
				$no_codons .= 'stop' unless ( $aux_codons =~ /stop/ );
				$no_codons =~ s/^\-// if ($no_codons ne '-');
				$no_codons =~ s/\/$// if ($no_codons ne '-');
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
		my ($no_codons) = '-';
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
			if ( $transcript->translate->codons ) {
				my ($aux_codons) = '';
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) {
						$aux_codons .= $codon->type.',';							
					}
				}
				$no_codons = 'start/' unless ( $aux_codons =~ /start/ );
				$no_codons .= 'stop' unless ( $aux_codons =~ /stop/ );
				$no_codons =~ s/^\-// if ($no_codons ne '-');
				$no_codons =~ s/\/$// if ($no_codons ne '-');
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
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $main::MATADOR3D_CUTOFF ) {
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

sub get_spade_annots($$$\$)
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
	
} # end get_spade_annots

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
