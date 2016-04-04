#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Parser qw(
	parse_infiles
	parse_transl_data
	parse_appris_methods
);
use APPRIS::Utils::Logger;
use APPRIS::Utils::File;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$APPRIS_CUTOFF
	$FIRESTAR_MINRES
	$FIRESTAR_CUTOFF
	$MATADOR3D_CUTOFF
	$SPADE_CUTOFF
	$CORSAIR_AA_LEN_CUTOFF
	$CORSAIR_CUTOFF
	$THUMP_CUTOFF
	$TSL_DB
	$TSL_CUTOFF

	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL

	$METHOD_WEIGHTED
	$METHOD_LABELS
);

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($data_file) = undef;
my ($transcripts_file) = undef;
my ($translations_file) = undef;
my ($firestar_file) = undef;
my ($matador3d_file) = undef;
my ($corsair_file) = undef;
my ($spade_file) = undef;
my ($thump_file) = undef;
my ($crash_file) = undef;
my ($inertia_file) = undef;
my ($proteo_file) = undef;
my ($output_main_file) = undef;
my ($output_nscore_file) = undef;
my ($output_label_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf=s'			=> \$config_file,
	'data=s'			=> \$data_file,
	'transcripts=s'		=> \$transcripts_file,
	'translations=s'	=> \$translations_file,
	'firestar=s'		=> \$firestar_file,
	'matador3d=s'		=> \$matador3d_file,
	'corsair=s'			=> \$corsair_file,
	'spade=s'			=> \$spade_file,
	'thump=s'			=> \$thump_file,
	'crash=s'			=> \$crash_file,
	'inertia=s'			=> \$inertia_file,
	'proteo=s'			=> \$proteo_file,
	'output=s'			=> \$output_main_file,
	'output_nscore=s'	=> \$output_nscore_file,	
	'output_label=s'	=> \$output_label_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Get conf vars
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD				= $FindBin::Bin;
$APPRIS_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'cutoff');
$FIRESTAR_MINRES		= $cfg->val( 'APPRIS_VARS', 'firestar_minres');
$FIRESTAR_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'firestar_cutoff');
$MATADOR3D_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'matador3d_cutoff');
$SPADE_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'spade_cutoff');
$CORSAIR_AA_LEN_CUTOFF	= $cfg->val( 'APPRIS_VARS', 'corsair_aa_cutoff');
$CORSAIR_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'corsair_cutoff');
$THUMP_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'thump_cutoff');

$TSL_DB					= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('TSL_VARS', 'db');
$TSL_CUTOFF				= $cfg->val( 'TSL_VARS', 'cutoff');

$OK_LABEL				= 'YES';
$UNKNOWN_LABEL			= 'UNKNOWN';
$NO_LABEL				= 'NO';
$METHOD_WEIGHTED = {
	'firestar'	=> 6,
	'matador3d'	=> 6,
	'spade'		=> 6,
	'corsair'	=> [{
					  'max'    => 3,
				  	  'weight' => 1	
					},
					{
					  'max'    => 4,
				  	  'weight' => 1.5	
					},
					{
					  'max'    => 6,
				  	  'weight' => 2	
					},
					{
					  'max'    => 7,
					  'weight' => 3	
					}],
	'thump'		=> 0,
	'crash'		=> 0,
	'inertia'	=> 0,
	'proteo'	=> 0,
};
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
 

# Required arguments
unless ( defined $config_file and 
		#defined $data_file and 
		#defined $transcripts_file and 
		defined $translations_file and 
		defined $firestar_file and defined $matador3d_file and defined $corsair_file and defined $spade_file and
		#defined $inertia_file and defined $thump_file and defined $crash_file and defined $proteo_file and
		defined $output_main_file and defined $output_label_file )
{
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
$logger->init_log($str_params);


#####################
# Method prototypes #
#####################
sub get_method_scores($$);
sub get_method_annots($$);
sub get_final_scores($$\$\$);
sub get_final_annotations($$$$$);

sub is_unique($$);
sub step_appris($$$);
sub step_ccds($$$);
sub step_tsl($$$);
sub step_ccds_eldest($$$);
sub step_ccds_longest($$$);
sub step_longest($$$);
sub step_tags($$$$\$);

sub get_score_output($$$);
sub get_nscore_output($$);
sub get_label_output($$);

sub get_firestar_annots($$$\$);
sub get_matador3d_annots($$$\$);
sub get_spade_annots($$$\$);
sub get_corsair_annots($$$\$);
sub get_thump_annots($$$\$);
sub get_crash_annots($$$$\$);
sub get_num_unusual_exons($);
sub get_tsl_annots($);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Local variables
	my ($firestar_result);
	my ($matador3d_result);
	my ($corsair_result);
	my ($spade_result);
	my ($thump_result);
	my ($crash_result);
	my ($inertia_result);
	my ($proteo_result);

	# get sequence data
	$logger->info("-- get data files\n");
	my ($gencode_data);
	my ($gene);
	if ( defined $data_file and defined $transcripts_file and defined $translations_file ) {
		($gencode_data) = parse_infiles($data_file, $transcripts_file, $translations_file);
	}
	elsif ( defined $translations_file ) {
		($gencode_data) = parse_transl_data($translations_file);
	}
	if ( defined $gencode_data and UNIVERSAL::isa($gencode_data, 'ARRAY') and (scalar(@{$gencode_data}) > 0) ) {
			$gene = $gencode_data->[0];
	}
	else {
		$logger->error("can not get gene result: $!\n");
	}
#	$logger->debug("GENE:\n".Dumper($gene)."\n");	
	
	# get reports
	$logger->info("get firestar result\n");
	if ( -e $firestar_file and (-s $firestar_file > 0) ) {
		$firestar_result = getStringFromFile($firestar_file);
		unless ( defined $firestar_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}
	$logger->info("get matador3d result\n");
	if ( -e $matador3d_file and (-s $matador3d_file > 0) ) {
		$matador3d_result = getStringFromFile($matador3d_file);
		unless ( defined $matador3d_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}
	$logger->info("get corsair result\n");
	if ( -e $corsair_file and (-s $corsair_file > 0) ) {
		$corsair_result = getStringFromFile($corsair_file);
		unless ( defined $corsair_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}	
	$logger->info("get spade result\n");
	if ( -e $spade_file and (-s $spade_file > 0) ) {
		$spade_result = getStringFromFile($spade_file);
		unless ( defined $spade_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}
	$logger->info("get thump result\n");
	if ( -e $thump_file and (-s $thump_file > 0) ) {
		$thump_result = getStringFromFile($thump_file);
		unless ( defined $thump_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}	
		
	$logger->info("get crash result\n");
	if ( -e $crash_file and (-s $crash_file > 0) ) {
		$crash_result = getStringFromFile($crash_file);
		unless ( defined $crash_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}
	$logger->info("get inertia result\n");
	if ( -e $inertia_file and (-s $inertia_file > 0) ) {
		$inertia_result = getStringFromFile($inertia_file);
		unless ( defined $inertia_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}		
	$logger->info("get proteo result\n");
	if ( -e $proteo_file and (-s $proteo_file > 0) ) {
		$proteo_result = getStringFromFile($proteo_file);
		unless ( defined $proteo_result ) {
			$logger->error("can not open method result: $!\n");
		}
	}
	else {
		$logger->info("file does not exit\n");
	}	

	# get object of reports
	$logger->info("-- create reports\n");
	my ($reports) = parse_appris_methods($gene, $firestar_result, $matador3d_result, $spade_result, $corsair_result, $crash_result, $thump_result, $inertia_result, $proteo_result);
	$logger->debug("REPORTS:\n".Dumper($reports)."\n");

	# get scores of methods for each transcript
	$logger->info("-- get scores of methods for each variant\n");
	my ($scores,$s_scores) = get_method_scores($gene, $reports);
	$logger->debug("PRE_SCORES:\n".Dumper($scores)."\n");
	$logger->debug("PRE_S_SCORES:\n".Dumper($s_scores)."\n");
	
	# get annots of methods for each transcript
	$logger->info("-- get annots of methods for each variant\n");
	my ($annots) = get_method_annots($gene, $s_scores);
	$logger->debug("PRE_ANNOTS:\n".Dumper($annots)."\n");

	# get scores/annots of appris for each transcript
	$logger->info("--  get scores/annots of appris for each transcript\n");
	my ($nscores) = get_final_scores($gene, $annots, $scores, $s_scores);
	$logger->debug("PRE2_SCORES:\n".Dumper($scores)."\n");
	$logger->debug("PRE2_S_SCORES:\n".Dumper($s_scores)."\n");
	$logger->debug("PRE2_NSCORES:\n".Dumper($nscores)."\n");

	# get annotations indexing each transcript
	$logger->info("-- get final annotations\n");
	get_final_annotations($gene, $scores, $s_scores, $nscores, $annots);
	$logger->debug("ANNOTS:\n".Dumper($annots)."\n");
	
	# print outputs
	my ($score_content) = get_score_output($gene, $scores, $annots);
	my ($p_score_out) = APPRIS::Utils::File::printStringIntoFile($score_content, $output_main_file);
	unless( defined $p_score_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($nscore_content) = get_nscore_output($gene, $nscores);
	my ($p_nscore_out) = APPRIS::Utils::File::printStringIntoFile($nscore_content, $output_nscore_file);
	unless( defined $p_nscore_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($label_content) = get_label_output($gene, $annots);
	my ($p_label_out) = APPRIS::Utils::File::printStringIntoFile($label_content, $output_label_file);
	unless( defined $p_label_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	
	$logger->finish_log();
	
	exit 0;	
	
}

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
				get_firestar_annots($gene, $scores->{$m}, $annot_label, $annots);
			}
			$m = 'matador3d';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_matador3d_annots($gene, $scores->{$m}, $annot_label, $annots);
			}
			$m = 'spade';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_spade_annots($gene, $scores->{$m}, $annot_label, $annots);
			}
			$m = 'corsair';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_corsair_annots($gene, $scores->{$m}, $annot_label, $annots);
			}
			$m = 'thump';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label) = $METHOD_LABELS->{$m}->[0];
				get_thump_annots($gene, $scores->{$m}, $annot_label, $annots);
			}
			$m = 'crash';
			if ( defined $scores and exists $scores->{$m} and exists $scores->{$m}->{'scores'} ) {				
				my ($annot_label_sp) = $METHOD_LABELS->{$m}->[0];
				my ($annot_label_tp) = $METHOD_LABELS->{$m}->[2];
				get_crash_annots($gene, $scores->{$m}, $annot_label_sp, $annot_label_tp, $annots);
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
				if ( $method eq 'corsair' ) {
					if ( $max >= $METHOD_WEIGHTED->{$method}->[3]->{'max'} ) {
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
				else {
					$appris_score += $METHOD_WEIGHTED->{$method}*$n_sc;
				}
			}
			
			# filter by biotype
			if ( $transcript->biotype and ($transcript->biotype eq 'nonsense_mediated_decay') ) {
				$appris_score = -1;
			}
			# filter by readthrough_transcript
			if ( $transcript->tag and exists $transcript->tag->{'readthrough_transcript'} ) {
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
		$logger->debug("GENE_2: \n".Dumper($gene)."\n");
		$logger->debug("PRINC_LIST_1: \n".Dumper($princ_list)."\n");
		$logger->debug("PRINC_REP_1: \n".Dumper($isof_report)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 1;
			step_tags($tag, $scores, $princ_list, $isof_report, $annots);
			return $tag;
		}

		# 2. from preserved transcript, we keep transcripts that they have got CCDS
		$princ_list = step_ccds($princ_list, $isof_report, $nscores);
		$logger->debug("PRINC_LIST_2: \n".Dumper($princ_list)."\n");			
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 2;
			step_tags($tag, $scores, $princ_list, $isof_report, $annots);
			return $tag;
		}

		# 3. from preserved transcript, we keep transcripts that they have got TSL1
		$princ_list = step_tsl($princ_list, $isof_report, $nscores);
		$logger->debug("PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 3;
			step_tags($tag, $scores, $princ_list, $isof_report, $annots);
			return $tag;
		}

		# 3_2. from preserved transcript, we keep transcripts that they have got eldest CCDS
		$princ_list = step_ccds_eldest($princ_list, $isof_report, $nscores);
		$logger->debug("PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 3;
			step_tags($tag, $scores, $princ_list, $isof_report, $annots);
			return $tag;
		}

		# 4. from preserved transcript, we keep transcripts that they have got longest seq with CCDS
		$princ_list = step_ccds_longest($princ_list, $isof_report, $nscores);
		$logger->debug("PRINC_LIST_4: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$tag = 4;
			step_tags($tag, $scores, $princ_list, $isof_report, $annots);
			return $tag;
		}

		# 5. from preserved transcript, we keep transcripts that they have got longest seq
		$princ_list = step_longest($princ_list, $isof_report, $nscores);
		$logger->debug("PRINC_LIST_5: \n".Dumper($princ_list)."\n");
		$annotations = step_tags(5, $scores, $princ_list, $isof_report, $annots);
		
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
	
	# get TSL(1) annot
	my ($tsl_transcs) = get_tsl_annots($gene_id);

	# scan the transcripts from the sorted APPRIS scores
	my ($highest_score) = $s_scores->{'max'};
	my (@sorted_ap_scores) = sort { $b <=> $a } keys (%{$s_scores->{'scores'}});			
	for ( my $i = 0; $i < scalar(@sorted_ap_scores); $i++ ) {
		
		my ($ap_score) = $sorted_ap_scores[$i];		
		foreach my $transc_id (@{$s_scores->{'scores'}->{$ap_score}}) {
			my ($index) = $gene->{'_index_transcripts'}->{$transc_id};
			my ($transcript) = $gene->transcripts->[$index];
			my ($transl_seq) = $transcript->translate->sequence;
			my ($transc_rep) = {
					'id'		=> $transc_id,
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
										$logger->warning("$gene_id has duplicate CCDS ids for the same sequence. We change to the eldest id: CCDS$old_ccds_id -> CCDS$new_ccds_id\n");
									}
									else {
										$logger->warning("$gene_id has duplicate CCDS ids for the same sequence. We keep the eldest id: CCDS$new_ccds_id -> CCDS$old_ccds_id \n");
									}
								}
							}						
						}
					}
				}
			}			
			# save TSL(1) annot
			if ( exists $tsl_transcs->{$transc_id} and ($tsl_transcs->{$transc_id} eq 'tsl1') ) {
				$transc_rep->{'tsl'} = 1;
				my ($xref_identities) = $transcript->xref_identify();
				push(@{$xref_identities},
						APPRIS::XrefEntry->new
						(
							-id				=> '1',
							-dbname			=> 'TSL'
						)
				);
				$transcript->xref_identify($xref_identities) if (defined $xref_identities);
			}
					
			# APPRIS says could be a principal when:
			# 1.	the Core of pipeline has not rejected the transcript (IT DOES NOT APPLY YET)
			# APPRIS rejected a transcript when:
			# 2.	it is a NMD: app_score is -1. Exception: we accept the last terms when the transcript is unique
			#
			my ($core_flag) = 1;
			#while (my ($method, $weight) = each(%{$METHOD_WEIGHTED}) ) {
			#	if ( $weight > 1 ) { # the core of pipeline has to be a weight... (Note: We have discarded THUMP with weight1)
			#		my ($label) = $METHOD_LABELS->{$method}->[0];
			#		if ( $annots->{$transc_id}->{$label} eq $NO_LABEL ) {
			#			$core_flag = 0;
			#			last;					
			#		}
			#	}
			#}
			my ($unique_transc) = ( scalar(keys(%{$s_scores->{'scores'}})) == 1 and scalar(@{$s_scores->{'scores'}->{$sorted_ap_scores[0]}}) >= 1 ) ? 1 : 0;
			if ( $core_flag == 1 ) {
				if ( ( ($highest_score - $ap_score) <=  $APPRIS_CUTOFF) and ( ($ap_score >= 0) or ($unique_transc == 1) ) ) {
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
		}
		if ( exists $nscores->{$transc_id} ) {
			if ( exists $nscores->{$transc_id}->{'appris'} ) {
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
			for (my $i=1; $i < scalar(@{$princ_isof->{$longest}}); $i++) {
				my ($transc_id_i) = ($princ_isof->{$longest})->[$i];
				my ($seq_i) = $princ_isof_seqs->{$transc_id_i};
				if ( $seq_i ne $seq_0 ) { $equal = 0 }
			}
			if ( $equal == 1 ) {
				foreach my $transc_id ( @{$princ_isof->{$longest}} ) {
					$report->{$transc_id} = 1;
				}
			}
			else { # the decision base on biggest appris-score
				my (@sort_sc) = sort { $b <=> $a } keys (%{$princ_isof_scores});
				my ($longest_sc) = $sort_sc[0];
				if ( scalar(@{$princ_isof_scores->{$longest_sc}}) == 1 ) {
					my ($transc_id) = $princ_isof_scores->{$longest_sc}[0];
					$report->{$transc_id} = 1;
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
		my ($biotype) = '-';
		my ($translation) = 'TRANSLATION';
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
		if ( $transcript->translate and $transcript->translate->sequence ) {
			
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
					elsif ( $xref->dbname eq 'TSL') {
						$tsl = $xref->id;
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
						
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$translation."\t".						
						$biotype."\t".
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
			$translation = 'NO_TRANSLATION';
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$translation."\n";
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
		my ($biotype) = '-';
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
		$biotype = $transcript->biotype if ( defined $transcript->biotype);		
		if ( $transcript->translate and $transcript->translate->sequence ) {
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
					elsif ( $xref->dbname eq 'TSL') {
						$tsl = $xref->id;
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
						
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$biotype."\t".
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
		my ($biotype) = '-';
		my ($translation) = 'TRANSLATION';
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
		$biotype = $transcript->biotype if ( defined $transcript->biotype);
		if ( $transcript->translate and $transcript->translate->sequence ) {
			
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
					elsif ( $xref->dbname eq 'TSL') {
						$tsl = $xref->id;
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
						
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$translation."\t".						
						$biotype."\t".
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
			$translation = 'NO_TRANSLATION';
			$content .= $stable_id."\t".
						$gene_name."\t".
						$transcript_id."\t".
						$translation."\n";					
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
	if ( $max <= $FIRESTAR_MINRES ) {
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
				if ( $sorted_scores[0] - $sorted_scores[$i] <= $FIRESTAR_CUTOFF ) {
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
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $MATADOR3D_CUTOFF ) {
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
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $SPADE_CUTOFF ) {
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
			if ( $sorted_scores[0] - $sorted_scores[$i] <= $CORSAIR_CUTOFF ) {
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

sub get_tsl_annots($)
{
	my ($gene_id) = @_;
	my ($report);
	
	# retrieve tsl annots from text file
	if ( defined $gene_id and ($gene_id ne '') ) {
		$logger->info("-- retrieve tsl annots from text file\n");
		my ($cmd) = "awk '{ if ( (\$1 ~ /$gene_id/) && (\$3 ~ /tsl1/) ) { print \$0 } }' $TSL_DB";
		my (@txt_out) = `$cmd`;
		foreach my $txt_line (@txt_out) {
			my (@txt_cols) = split('\t', $txt_line);
			my ($g_id) = $txt_cols[0];
			my ($t_id) = $txt_cols[1];
			my ($tsl_annot) = $txt_cols[2];
			my ($appris_annot) = $txt_cols[3];
			if ( $tsl_annot =~ /^tsl1/ ) {
				$report->{$t_id} = 'tsl1';
			}			
		}
	}	
	return $report;
	
} # End get_tsl_annots


main();



__END__

=head1 NAME

appris

=head1 DESCRIPTION

Run APPRIS program

=head1 SYNOPSIS

appris

=head2 Required arguments:

	--conf <Config file>
	
	--data=  <Gencode data file>
	
	--transcripts=  <Gencode transcript file>
	
	--translations=  <Gencode translations file>
	
	--firestar <firestar results for a gene>
	
	--matador3d <Matador3D results for a gene>

	--corsair <CORSAIR results for a gene>

	--spade <SPADE results for a gene>
	
	--thump <THUMP results for a gene>

	--crash <CRASH results for a gene>
	
	--inertia <INERTIA results for a gene>

	--proteo <PROTEO results for a gene>

	--output <Annotation output file>
	
	--output_nscore <Output file of normalized scores>
    
	--output_label <Output file of labels>
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
    

=head1 EXAMPLE

perl appris.pl

	--conf=../conf/pipeline.ini
	
	--data=examples/ENSG00000142185/ENSG00000142185.gencode.v7.annotation.gtf
	
	--transcripts=examples/ENSG00000142185/ENSG00000142185.gencode.v7.pc_transcripts.fa
	
	--translations=examples/ENSG00000142185/ENSG00000142185.gencode.v7.pc_translations.fa

	--firestar=examples/ENSG00000142185/ENSG00000142185.firestar
	
	--matador3d=examples/ENSG00000142185/ENSG00000142185.matador3d
	
	--corsair=examples/ENSG00000142185/ENSG00000142185.corsair

	--spade=examples/ENSG00000142185/ENSG00000142185.spade
	
	--thump=examples/ENSG00000142185/ENSG00000142185.thump

	--crash=examples/ENSG00000142185/ENSG00000142185.crash
	
	--inertia=examples/ENSG00000142185/ENSG00000142185.inertia
	
	--proteo=examples/ENSG00000142185/ENSG00000142185.proteo	
	
	--output=examples/ENSG00000142185/ENSG00000016864.appris
	
	--output_nscore=examples/ENSG00000142185/ENSG00000016864.appris.nscore
	
	--output_label=examples/ENSG00000142185/ENSG00000016864.appris.label
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut