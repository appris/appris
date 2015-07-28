#!/usr/bin/perl -W

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
	$FIRESTAR_CUTOFF
	$MATADOR3D_CUTOFF
	$CORSAIR_CUTOFF
	$CORSAIR_AA_LEN_CUTOFF
	$SPADE_CUTOFF
	$THUMP_CUTOFF
	$TSL_DB
	$TSL_CUTOFF
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
	$FIRESTAR_ACCEPT_LABEL
	$FIRESTAR_REJECT_LABEL
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
$FIRESTAR_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'firestar_cutoff');
$MATADOR3D_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'matador3d_cutoff');
$CORSAIR_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'corsair_cutoff');
$CORSAIR_AA_LEN_CUTOFF	= $cfg->val( 'APPRIS_VARS', 'corsair_aa_cutoff');
$SPADE_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'spade_cutoff');
$THUMP_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'thump_cutoff');

$TSL_DB					= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('TSL_VARS', 'db');
$TSL_CUTOFF				= $cfg->val( 'TSL_VARS', 'cutoff');

$OK_LABEL				= 'YES';
$UNKNOWN_LABEL			= 'UNKNOWN';
$NO_LABEL				= 'NO';
$FIRESTAR_ACCEPT_LABEL	= 'ACCEPT';
$FIRESTAR_REJECT_LABEL	= 'REJECT';
$METHOD_WEIGHTED = {
	'firestar'	=> 6,
	'matador3d'	=> 6,
	'spade'		=> 4,
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
	'thump'		=> 1,
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
sub get_scores($$);
sub get_final_scores($$$);
sub filter_by_class_codons($$\$\$);
sub get_annotations($$$);
#sub is_discarted_by_nmd_codons($);
sub is_unique($$);
sub step_appris($$);
sub step_ccds($$);
sub step_tsl($$);
sub step_ccds_eldest($$);
sub step_ccds_longest($$);
sub step_longest($$);
#sub step_num_peptides($$$$$);
sub step_tags($$$$);
sub get_score_output($$$);
sub get_nscore_output($$$$);
sub get_label_output($$$);
sub get_strange_output($$);
sub get_crash_annots($$$);
sub get_crash_annot($$);
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
	$logger->debug("GENE:\n".Dumper($gene)."\n");	
	
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

	# get scores of each transcript
	$logger->info("-- get scores for each variant\n");
	my ($pre_scores, $max_scores) = get_scores($gene, $reports);
	$logger->debug("PRE_SCORES:\n".Dumper($pre_scores)."\n");
	$logger->debug("PRE_MAX_SCORES:\n".Dumper($max_scores)."\n");
	
#	# get scores of each transcript
#	$logger->info("-- get scores for each variant\n");
#	my ($scores, $s_scores, $nscores) = get_final_scores($gene, $pre_scores, $max_scores);
#	$logger->debug("SCORES:\n".Dumper($scores)."\n");
#	$logger->debug("N_SCORES:\n".Dumper($nscores)."\n");
#	$logger->debug("G_SCORES:\n".Dumper($s_scores)."\n");

	# get scores of each transcript
	$logger->info("-- get scores for each variant\n");
	my ($scores, $pre_s_scores, $nscores) = get_final_scores($gene, $pre_scores, $max_scores);
	$logger->debug("PRE2_SCORES:\n".Dumper($scores)."\n");
	$logger->debug("PRE2_S_SCORES:\n".Dumper($pre_s_scores)."\n");
	$logger->debug("PRE2_NSCORES:\n".Dumper($nscores)."\n");

	# filter list of transcripts/scores by CLASS (and CODONS is disabled)
	my ($s_scores) = filter_by_class_codons($gene, $pre_s_scores, $scores, $nscores);
	$logger->debug("SCORES:\n".Dumper($scores)."\n");
	$logger->debug("N_SCORES:\n".Dumper($nscores)."\n");
	$logger->debug("G_SCORES:\n".Dumper($s_scores)."\n");

	# get annotations indexing each transcript
	$logger->info("-- get final annotations\n");
	my ($annots) = get_annotations($gene, $scores, $s_scores);
	$logger->debug("ANNOTS:\n".Dumper($annots)."\n");
	
	# print outputs
	my ($score_content) = get_score_output($gene, $scores, $annots);
	my ($p_score_out) = APPRIS::Utils::File::printStringIntoFile($score_content, $output_main_file);
	unless( defined $p_score_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($nscore_content) = get_nscore_output($gene, $scores, $nscores, $annots);
	my ($p_nscore_out) = APPRIS::Utils::File::printStringIntoFile($nscore_content, $output_nscore_file);
	unless( defined $p_nscore_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($label_content) = get_label_output($gene, $scores, $annots);
	my ($p_label_out) = APPRIS::Utils::File::printStringIntoFile($label_content, $output_label_file);
	unless( defined $p_label_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($strange_content) = get_strange_output($gene, $annots);
	if ( $strange_content ne '' ) {
		my ($output_strange_file) = $output_main_file.'.strange';
		my ($p_strange_out) = APPRIS::Utils::File::printStringIntoFile($strange_content, $output_strange_file);
		unless( defined $p_strange_out ) {
			$logger->error("Can not create output strange file: $!\n");
		}		
	}
	
	$logger->finish_log();
	
	exit 0;	
	
}

# get the main functional isoform from methods of appris
sub get_scores($$)
{
	my ($gene, $reports) = @_;

	my ($stable_id) = $gene->stable_id;
	my ($scores, $max_scores);
	
	# init max scores
	foreach my $method ( keys(%{$METHOD_LABELS}) ) {
		$max_scores->{$method} = 0;
	}
	$max_scores->{'aa_length'} = 0;
	
	# get annotations for each method (or group of method) -------------------
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($index) = $reports->{'_index_transcripts'}->{$transcript_id};
			my ($result) = $reports->transcripts->[$index];
			my ($aa_length) = length($transcript->translate->sequence);
			if ( $aa_length > $max_scores->{'aa_length'} ) { $max_scores->{'aa_length'} = $aa_length }
			
			# get firestar
			my ($m) = 'firestar';
			if ( $result and $result->analysis and $result->analysis->firestar ) {				
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->firestar;
				if ( $analysis->functional_residue and ($analysis->functional_residue eq $FIRESTAR_REJECT_LABEL) ) {	
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->functional_residue and ($analysis->functional_residue eq $FIRESTAR_ACCEPT_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				if ( defined $analysis->num_residues ) {
					my ($sc) = $analysis->num_residues;
					$scores->{$transcript_id}->{$k_annot2} = $sc;
					if ( !(exists $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc; 
					}
					elsif ( exists $max_scores->{$m} and ($sc > $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc;
					}
				}
			}
			# get matador3d
			$m = 'matador3d';
			if ( $result and $result->analysis and $result->analysis->matador3d ) {				
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];				
				my ($analysis) = $result->analysis->matador3d;
				if ( $analysis->conservation_structure and ($analysis->conservation_structure eq $NO_LABEL) ) {	
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->conservation_structure and ($analysis->conservation_structure eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->conservation_structure and ($analysis->conservation_structure eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $OK_LABEL;
				}
				if ( defined $analysis->score ) {
					my ($sc) = $analysis->score;
					$scores->{$transcript_id}->{$k_annot2} = $sc;					
					if ( !(exists $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc; 
					}
					elsif ( exists $max_scores->{$m} and ($sc > $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc;
					}
				}				
			}
			# get corsair
			$m = 'corsair';
			if ( $result and $result->analysis and $result->analysis->corsair ) {
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->corsair;
				if ( $analysis->vertebrate_signal and ($analysis->vertebrate_signal eq $NO_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->vertebrate_signal and ($analysis->vertebrate_signal eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->vertebrate_signal and ($analysis->vertebrate_signal eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $OK_LABEL;
				}
				if ( defined $analysis->score ) {
					my ($sc) = $analysis->score;
					$scores->{$transcript_id}->{$k_annot2} = $sc;					
					if ( !(exists $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc; 
					}
					elsif ( exists $max_scores->{$m} and ($sc > $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc;
					}
				}			
			}
			# get spade
			$m = 'spade';
			if ( $result and $result->analysis and $result->analysis->spade ) {				
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->spade;
				if ( $analysis->domain_signal and ($analysis->domain_signal eq $NO_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->domain_signal and ($analysis->domain_signal eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->domain_signal and ($analysis->domain_signal eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $OK_LABEL;
				}
				if ( defined $analysis->num_domains and 
					 defined $analysis->num_possibly_damaged_domains and
					 defined $analysis->num_damaged_domains and
					 defined $analysis->num_wrong_domains ) {
					 	my ($sum_domains) = $analysis->num_domains + $analysis->num_possibly_damaged_domains;
					 	my ($sum_damaged_domains) = $analysis->num_damaged_domains + $analysis->num_wrong_domains;
						my ($sc) = $sum_domains;
						#$scores->{$transcript_id}->{$k_annot2} = $sum_domains . '-' . $sum_damaged_domains;
						$scores->{$transcript_id}->{$k_annot2} = $sum_domains;						
						if ( !(exists $max_scores->{$m}) ) {
							$max_scores->{$m} = $sc; 
						}
						elsif ( exists $max_scores->{$m} and ($sc > $max_scores->{$m}) ) {
							$max_scores->{$m} = $sc;
						}
				}				
			}
			# get thump
			$m = 'thump';
			if ( $result and $result->analysis and $result->analysis->thump ) {
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->thump;
				if ( $analysis->transmembrane_signal and ($analysis->transmembrane_signal eq $NO_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->transmembrane_signal and ($analysis->transmembrane_signal eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->transmembrane_signal and ($analysis->transmembrane_signal eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $OK_LABEL;
				}
				if ( defined $analysis->num_tmh and 
					 defined $analysis->num_damaged_tmh ) {
					my ($sc) = $analysis->num_tmh;
					#$scores->{$transcript_id}->{$k_annot2} = $analysis->num_tmh.'-'.$analysis->num_damaged_tmh;
					$scores->{$transcript_id}->{$k_annot2} = $analysis->num_tmh;
					if ( !(exists $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc; 
					}
					elsif ( exists $max_scores->{$m} and ($sc > $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc;
					}
				}				
			}
			# get crash: signalp + targetp
			$m = 'crash';
			if ( $result and $result->analysis and $result->analysis->crash ) {
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot1) = $METHOD_LABELS->{$m}->[1];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[2];
				my ($k_annot3) = $METHOD_LABELS->{$m}->[3];		
				my ($analysis) = $result->analysis->crash;
				# signalp
				if ( $analysis->peptide_signal and ($analysis->peptide_signal eq $NO_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->peptide_signal and ($analysis->peptide_signal eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->peptide_signal and ($analysis->peptide_signal eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $OK_LABEL;
				}
				my ($crash_sp_annot, $crash_sp_signal) = get_crash_annots('crash_sp', $transcript_id, $reports);
				if ( defined $crash_sp_annot ) {
					if ( defined $analysis->sp_score ) {
						$scores->{$transcript_id}->{$k_annot1} = $analysis->sp_score;							
					}
				}
				# targetp
				if ( $analysis->mitochondrial_signal and ($analysis->mitochondrial_signal eq $NO_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot2} = $NO_LABEL;
				}
				elsif ( $analysis->mitochondrial_signal and ($analysis->mitochondrial_signal eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot2} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->mitochondrial_signal and ($analysis->mitochondrial_signal eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot2} = $OK_LABEL;
				}
				my ($crash_tp_annot, $crash_tp_signal) = get_crash_annots('crash_tp', $transcript_id, $reports);				
				if ( defined $crash_tp_annot ) {
					if ( defined $analysis->tp_score ) {
						$scores->{$transcript_id}->{$k_annot3} = $analysis->tp_score;							
					}					
				}
			}
			# get inertia
			$m = 'inertia';
			if ( $result and $result->analysis and $result->analysis->inertia ) {				
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];				
				my ($analysis) = $result->analysis->inertia;
				if ( $analysis->unusual_evolution and ($analysis->unusual_evolution eq $NO_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $NO_LABEL;
				}
				elsif ( $analysis->unusual_evolution and ($analysis->unusual_evolution eq $UNKNOWN_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $UNKNOWN_LABEL;
				}
				elsif ( $analysis->unusual_evolution and ($analysis->unusual_evolution eq $OK_LABEL) ) {
					$scores->{$transcript_id}->{$k_annot} = $OK_LABEL;
				}
				if ( defined $analysis->regions ) {
					 my ($num_un_exons) = get_num_unusual_exons($analysis);
					 $scores->{$transcript_id}->{$k_annot2} = $num_un_exons->{'inertia'};
				}
			}
			# get proteo (WE DON'T USE PROTEO FOR APPRIS DECISION)
			$m = 'proteo';
			if ( $result and $result->analysis and $result->analysis->proteo ) {				
				my ($k_annot) = $METHOD_LABELS->{$m}->[0];
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
				my ($analysis) = $result->analysis->proteo;
				if ( defined $analysis->num_peptides ) {
					my ($sc) = $analysis->num_peptides;
					$scores->{$transcript_id}->{$k_annot2} = $sc;
					if ( !(exists $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc; 
					}
					elsif ( exists $max_scores->{$m} and ($sc > $max_scores->{$m}) ) {
						$max_scores->{$m} = $sc;
					}
				}
			}
		}
	}
	
	return ($scores, $max_scores);
	
} # End get_scores

# get gene scores of methods for appris
sub get_final_scores($$$)
{
	my ($gene, $scores, $max_scores) = @_;
	my ($s_scores, $nscores);
	
	# obtain normalize scores (weighted normalize score) foreach method
	my ($max_appris_score) = 0;
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($appris_score) = 0;
			foreach my $method ( keys(%{$METHOD_WEIGHTED}) ) {
				my ($n_sc) = 0;
				my ($sc) = 0;
				my ($max) = $max_scores->{$method};
				my ($k_annot2) = $METHOD_LABELS->{$method}->[1];
				
				# create normalize scores
				if ( exists $scores->{$transcript_id}->{$k_annot2} ) {
					$sc = $scores->{$transcript_id}->{$k_annot2};
					if ( $method eq 'firestar' ) {
						if ( $max != 0 ) { $n_sc = 1 + ( ($sc - $max)/15 ) }
						else { $n_sc = 0 }
					}
					elsif ( $method eq 'matador3d' ) {
						if ( $max != 0 ) { $n_sc = 1 + ( ($sc - $max)*(2/9) ) }
						else { $n_sc = 0 }
					}
					elsif ( $method eq 'spade' ) {
						if ( $max != 0 ) { $n_sc = 1 + ( ($sc - $max)/3 ) }
						else { $n_sc = 0 }
					}
					elsif ( $method eq 'thump' ) {
						if ( $max != 0 ) { $n_sc = 1 + ( ($sc - $max)/2 ) }
						else { $n_sc = 0 }
					}					
					elsif ( $method eq 'corsair' ) {
						# cutoff of AA length. 0 score when is bigger or equal than AA_LEN_CUTOFF
						my ($max_aa_legnth) = $max_scores->{'aa_length'};
						if ( $max_aa_legnth < $CORSAIR_AA_LEN_CUTOFF ) { $sc = $scores->{$transcript_id}->{$k_annot2} }
						else { $sc = 0 }
						
						if ( $max != 0 ) { $n_sc = $sc/$max }
						else { $n_sc = 0 }						
					}
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
			
			# save appris score and normalize score
			my ($m) = 'appris';
			my ($k_annot2) = $METHOD_LABELS->{$m}->[1];			
			$scores->{$transcript_id}->{$k_annot2} = $appris_score;
			push(@{$s_scores->{$appris_score}}, $transcript_id);			
			if ( !exists $max_scores->{$m} ) { $max_scores->{$m} = $appris_score }
			elsif ( $appris_score > $max_scores->{$m} ) { $max_scores->{$m} = $appris_score }
		}
	}
	
	# get normalize scores for appris
	my ($m) = 'appris';
	foreach my $transcript (@{$gene->transcripts}) {	
		my ($transcript_id) = $transcript->stable_id;
		if ( $transcript->translate and $transcript->translate->sequence ) {
			my ($n_sc) = 0;
			my ($max) = $max_scores->{$m};
			if ( $max != 0 ) {
				my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
				my ($sc) = $scores->{$transcript_id}->{$k_annot2};
				$n_sc = sprintf("%.2f",($sc/$max));
			}
			$nscores->{$transcript_id}->{$m} = $n_sc;			
		}
	}
		
	return ($scores, $s_scores, $nscores);
	
} # End get_final_scores

sub filter_by_class_codons($$\$\$)
{
	my ($gene, $i_s_scores, $ref_scores, $ref_nscores) = @_;
	my ($s_scores);
	my ($method) = 'appris';
	
	my (@sorted_scores) = keys (%{$i_s_scores});	
	for ( my $i = 0; $i < scalar(@sorted_scores); $i++ ) {
		my ($appris_score) = $sorted_scores[$i];
		foreach my $transcript_id (@{$i_s_scores->{$appris_score}}) {
			my ($index) = $gene->{'_index_transcripts'}->{$transcript_id};
			my ($transcript) = $gene->transcripts->[$index];
			my ($appris_nscore) = $$ref_nscores->{$transcript_id}->{$method};
			my ($frozen_score, $frozen_nscore) = ($appris_score, $appris_nscore);
			if ( $transcript->biotype and ($transcript->biotype eq 'nonsense_mediated_decay') ) {
				($frozen_score, $frozen_nscore) = (0,0);
			}
			if ( $transcript->translate->codons ) {
				my ($aux_codons) = '';
				foreach my $codon (@{$transcript->translate->codons}) {
					if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) {
						$aux_codons .= $codon->type.',';							
					}
				}
				unless ( ($aux_codons =~ /start/) and ($aux_codons =~ /stop/) ) {
					($frozen_score, $frozen_nscore) = (0,0);
				}
			}
			push(@{$s_scores->{$frozen_score}}, $transcript_id);
			$$ref_scores->{$transcript_id}->{$METHOD_LABELS->{$method}->[1]} = $frozen_score;
			$$ref_nscores->{$transcript_id}->{$method} = $frozen_nscore;
		}			
	}
	return $s_scores;
	
} # end filter_by_class_codons

# get the final annotation
sub get_annotations($$$)
{
	my ($gene, $scores, $sort_scores) = @_;
	
	my ($annotations);
		
	# scan transcripts sorted by appris score
	if ( scalar(keys(%{$sort_scores})) > 0 )
	{
		# 1. acquire the dominant transcripts from appris score.
		# They have to pass the cutoff to be added into "principal" list
		my ($princ_list, $isof_report) = step_appris($gene, $sort_scores);
		$logger->debug("PRINC_LIST_1: \n".Dumper($princ_list)."\n");
		$logger->debug("PRINC_REP_1: \n".Dumper($isof_report)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$annotations = step_tags(1, $scores, $princ_list, $isof_report);
			return $annotations;
		}

		# 2. from preserved transcript, we keep transcripts that they have got CCDS
		$princ_list = step_ccds($princ_list, $isof_report);
		$logger->debug("PRINC_LIST_2: \n".Dumper($princ_list)."\n");			
		if ( is_unique($princ_list, $isof_report) ) {
			$annotations = step_tags(2, $scores, $princ_list, $isof_report);
			return $annotations;
		}

		# 3. from preserved transcript, we keep transcripts that they have got TSL1
		$princ_list = step_tsl($princ_list, $isof_report);
		$logger->debug("PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$annotations = step_tags(3, $scores, $princ_list, $isof_report);
			return $annotations;
		}

		# 3_2. from preserved transcript, we keep transcripts that they have got eldest CCDS
		$princ_list = step_ccds_eldest($princ_list, $isof_report);
		$logger->debug("PRINC_LIST_3: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$annotations = step_tags(3, $scores, $princ_list, $isof_report);
			return $annotations;
		}

		# 4. from preserved transcript, we keep transcripts that they have got longest seq with CCDS
		$princ_list = step_ccds_longest($princ_list, $isof_report);
		$logger->debug("PRINC_LIST_4: \n".Dumper($princ_list)."\n");
		if ( is_unique($princ_list, $isof_report) ) {
			$annotations = step_tags(4, $scores, $princ_list, $isof_report);
			return $annotations;
		}

		# 5. from preserved transcript, we keep transcripts that they have got longest seq
		$princ_list = step_longest($princ_list, $isof_report);
		$logger->debug("PRINC_LIST_5: \n".Dumper($princ_list)."\n");
		$annotations = step_tags(5, $scores, $princ_list, $isof_report);
		
	}
	
	return $annotations;
	
} # End get_annotations

sub is_discarted_by_nmd_codons($)
{
	my ($transcript) = @_;
	my ($discarded) = 0;
	
	if ( $transcript->biotype and ($transcript->biotype eq 'nonsense_mediated_decay') ) {
		$discarded = 1;
	}
	if ( $transcript->translate->codons ) {
		my ($aux_codons) = '';
		foreach my $codon (@{$transcript->translate->codons}) {
			if ( ($codon->type eq 'start') or ($codon->type eq 'stop') ) {
				$aux_codons .= $codon->type.',';							
			}
		}
		unless ( ($aux_codons =~ /start/) and ($aux_codons =~ /stop/) ) {
			$discarded = 1;
		}
	}

	return $discarded;
	
} # end is_discarted_by_nmd_codons

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

sub step_appris($$)
{
	my ($gene, $appris_scores) = @_;
	my ($gene_id) = $gene->stable_id;
	my ($highest_score) = 0;
	my ($princ_list, $isof_report);
	my ($ccds_ids);
	
	# get TSL(1) annot
	my ($tsl_transcs) = get_tsl_annots($gene_id);

	# scan the transcripts from the sorted APPRIS scores
	my (@sorted_ap_scores) = sort { $b <=> $a } keys (%{$appris_scores});			
	for ( my $i = 0; $i < scalar(@sorted_ap_scores); $i++ ) {
		
		# save the highest score (the first one)
		if ( $i == 0 ) { $highest_score = $sorted_ap_scores[$i] }
		
		my ($ap_score) = $sorted_ap_scores[$i];		
		foreach my $transc_id (@{$appris_scores->{$ap_score}}) {
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
			}
			# save the APPRIS decision
			# say if the number total of transcripts is one			
			#my ($unique_transc) = ( scalar(keys(%{$appris_scores})) == 1 and scalar(@{$appris_scores->{$sorted_ap_scores[0]}}) == 1 ) ? 1 : 0;
			#if ( ( ($highest_score - $ap_score) <=  $APPRIS_CUTOFF) and ( is_discarted_by_nmd_codons($transcript) == 0 or $unique_transc == 1 ) ) {
			if ( ( ($highest_score - $ap_score) <=  $APPRIS_CUTOFF) ) {
				$transc_rep->{'principal'} = 1;
				$princ_list->{$transc_id} = 1;
			}
			#else { $transc_rep->{'principal'} = 0 }
			
			push(@{$isof_report}, $transc_rep);

		}			
	}
	return ($princ_list, $isof_report);
	
} # end step_appris

sub step_ccds($$)
{
	my ($i_princ_list, $isof_report) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			
			if ( exists $princ->{'ccds'} ) { push(@{$princ_isof}, $princ) }
			
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

sub step_tsl($$)
{
	my ($i_princ_list, $isof_report) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			
			if ( exists $princ->{'tsl'} ) { push(@{$princ_isof}, $princ) }
			
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

sub step_ccds_eldest($$)
{
	my ($i_princ_list, $isof_report) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			
			if ( exists $princ->{'ccds'} ) { push(@{$princ_isof}, $princ) }
			
		}
	}	

	# print princ isoforms with eldests CCDS
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
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

sub step_ccds_longest($$)
{
	my ($i_princ_list, $isof_report) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			
			if ( exists $princ->{'ccds'} ) { push(@{$princ_isof}, $princ) }
			
		}
	}	

	# print princ isoforms with eldests CCDS
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
		my (@sort_transc) = sort { $a->{'length'} <=> $b->{'length'} } @{$princ_isof};
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

sub step_longest($$)
{
	my ($i_princ_list, $isof_report) = @_;
	my ($report);
	
	# preliminar report
	my ($princ_isof);
	foreach my $princ ( @{$isof_report} ) {
		my ($transc_id) = $princ->{'id'};
		if ( exists $i_princ_list->{$transc_id} ) {
			
			if ( exists $princ->{'length'} ) { push(@{$princ_isof}, $princ) }
			
		}
	}	

	# print princ isoforms with eldests CCDS
	if ( defined $princ_isof and ( scalar(@{$princ_isof}) >= 1 ) ) {
		my (@sort_transc) = sort { $a->{'length'} <=> $b->{'length'} } @{$princ_isof};
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
		
} # end step_longest

sub step_tags($$$$)
{
	my ($step, $scores, $princ_list, $isof_report) = @_;
	my ($annotations);
	
	# add tags (reliability) and annotations for
	my ($m) = 'appris';
	my ($k_annot) = $METHOD_LABELS->{$m}->[0];
	my ($k_annot2) = $METHOD_LABELS->{$m}->[1];
	my ($k_annot3) = $METHOD_LABELS->{$m}->[2];	
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
		$annotations->{$transc_id}->{$k_annot} = $status;
		$annotations->{$transc_id}->{$k_annot3} = $annot;		
	}
		
	return $annotations;
	
} # step_tags

sub get_score_output($$$)
{
	my ($gene, $scores, $annots) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($content) = '';

	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		my ($status) = '-';
		my ($biotype) = '-';
		my ($translation) = 'TRANSLATION';
		my ($ccds_id) = '-';
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
		$status = $transcript->status if ($transcript->status);
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
						last;					
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
						$transcript_id."\t".
						$translation."\t".						
						$status."\t".
						$biotype."\t".
						$no_codons."\t".
						$ccds_id."\t".
						$transl_len."\t".
						$firestar_annot."\t".
						$matador3d_annot."\t".
						$corsair_annot."\t".
						$spade_annot."\t".
						$thump_annot."\t".
						$crash_sp_annot."\t".
						$crash_tp_annot."\t".
						$inertia_annot."\t".
						$proteo_annot."\t".
						$appris_annot."\n";
		}
		else {
			$translation = 'NO_TRANSLATION';
			$content .= $stable_id."\t".
						$transcript_id."\t".
						$translation."\t".						
						$status."\t".
						$biotype."\t".
						$no_codons."\t".
						$ccds_id."\n";						
		}
	}
	
	return $content;
	
} # End get_score_output

sub get_nscore_output($$$$)
{
	my ($gene, $scores, $nscores, $annots) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($content) = '';

	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		my ($biotype) = '-';
		my ($ccds_id) = '-';
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
			#$appris_annot = $scores->{$transcript_id}->{'score_principal_isoform'} if ( exists $scores->{$transcript_id}->{'score_principal_isoform'} );
			$appris_annot = $nscores->{$transcript_id}->{'appris'} if ( exists $nscores->{$transcript_id}->{'appris'} );
						
			$transl_len = length($transcript->translate->sequence);
			
			if ( $transcript->xref_identify ) {
				foreach my $xref (@{$transcript->xref_identify}) {
					if ( $xref->dbname eq 'CCDS') {
						$ccds_id = $xref->id;
						last;					
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
						$transcript_id."\t".
						$biotype."\t".
						$no_codons."\t".
						$ccds_id."\t".
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
		else {
			$content .= $stable_id."\t".
						$transcript_id."\t".
						$biotype."\t".
						$no_codons."\t".
						$ccds_id."\n";						
		}
	}
	
	return $content;
	
} # End get_nscore_output

sub get_label_output($$$)
{
	my ($gene, $scores, $annots) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($content) = '';

	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		my ($status) = '-';
		my ($biotype) = '-';
		my ($translation) = 'TRANSLATION';
		my ($ccds_id) = '-';
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
		$status = $transcript->status if ($transcript->status);
		$biotype = $transcript->biotype if ( defined $transcript->biotype);
		if ( $transcript->translate and $transcript->translate->sequence ) {
			
			$firestar_annot = $scores->{$transcript_id}->{'functional_residue'} if ( exists $scores->{$transcript_id}->{'functional_residue'} );
			$matador3d_annot = $scores->{$transcript_id}->{'conservation_structure'}if ( exists $scores->{$transcript_id}->{'conservation_structure'} );
			$corsair_annot = $scores->{$transcript_id}->{'vertebrate_signal'} if ( exists $scores->{$transcript_id}->{'vertebrate_signal'} );
			$spade_annot = $scores->{$transcript_id}->{'domain_signal'} if ( exists $scores->{$transcript_id}->{'domain_signal'} );
			$thump_annot = $scores->{$transcript_id}->{'transmembrane_signal'} if ( exists $scores->{$transcript_id}->{'transmembrane_signal'} );
			$crash_sp_annot = $scores->{$transcript_id}->{'peptide_signal'} if ( exists $scores->{$transcript_id}->{'peptide_signal'} );
			$crash_tp_annot = $scores->{$transcript_id}->{'mitochondrial_signal'} if ( exists $scores->{$transcript_id}->{'mitochondrial_signal'} );
			$inertia_annot = $scores->{$transcript_id}->{'unusual_evolution'} if ( exists $scores->{$transcript_id}->{'unusual_evolution'} );
			$proteo_annot = $scores->{$transcript_id}->{'peptide_evidence'} if ( exists $scores->{$transcript_id}->{'peptide_evidence'} );			
			
			$appris_annot = $annots->{$transcript_id}->{'principal_isoform'} if ( exists $annots->{$transcript_id}->{'principal_isoform'} );
			$appris_relia = $annots->{$transcript_id}->{'reliability'} if ( exists $annots->{$transcript_id}->{'reliability'} );
			
			$transl_len = length($transcript->translate->sequence);
			
			if ( $transcript->xref_identify ) {
				foreach my $xref (@{$transcript->xref_identify}) {
					if ( $xref->dbname eq 'CCDS') {
						$ccds_id = $xref->id;
						last;					
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
						$transcript_id."\t".
						$translation."\t".						
						$status."\t".
						$biotype."\t".
						$no_codons."\t".
						$ccds_id."\t".
						$transl_len."\t".				
						$firestar_annot."\t".
						$matador3d_annot."\t".
						$corsair_annot."\t".
						$spade_annot."\t".
						$thump_annot."\t".
						$crash_sp_annot."\t".
						$crash_tp_annot."\t".
						$inertia_annot."\t".
						$proteo_annot."\t".
						$appris_annot."\t".
						$appris_relia."\n";
		}
		else {
			$translation = 'NO_TRANSLATION';
			$content .= $stable_id."\t".
						$transcript_id."\t".
						$translation."\t".						
						$status."\t".
						$biotype."\t".
						$no_codons."\t".
						$ccds_id."\n";						
		}
	}
	
	return $content;
	
} # End get_label_output

# Get the content of strange cases
sub get_strange_output($$)
{
	my ($gene, $annots) = @_;
	my ($stable_id) = $gene->stable_id;
	my ($content) = '';
	
	foreach my $transcript (@{$gene->transcripts})
	{	
		my ($transcript_id) = $transcript->stable_id;
		if ( exists $annots->{$transcript_id}->{'strange'} ) {
			my ($chr) = $transcript->chromosome;
			my ($ccds_id) = '-';		
			if ( $transcript->xref_identify ) {
				foreach my $xref (@{$transcript->xref_identify}) {
					if ( $xref->dbname eq 'CCDS') {
						$ccds_id = $xref->id;
						last;					
					}
				}
			}			
			$content .= $chr."\t".
						$stable_id."\t".
						$transcript_id."\t".
						$ccds_id."\n";
		}
	}
	
	return $content;
	
} # End get_strange_output

# Get annotation of several methods
sub get_method_annot($$)
{
	my ($method, $result) = @_;
	
	my ($annot) = {
		'not_defined'	=> 0,
		'ok'			=> 0,
		'unknown'		=> 0,
		'rejected'		=> 0,
	};
		
	if ( $method eq 'inertia' ) { # get inertia
		if ( $result and $result->analysis and $result->analysis->inertia ) {
			my ($analysis) = $result->analysis->inertia;
			if ( $analysis->unusual_evolution and ($analysis->unusual_evolution eq $UNKNOWN_LABEL) ) {
				$annot->{'unknown'} = 1;
			}
			elsif ( $analysis->unusual_evolution and ($analysis->unusual_evolution eq $NO_LABEL) ) {
				$annot->{'rejected'} = 1;			
			}
			else {
				$annot->{'not_defined'} = 1;
			}
		}
		else {
			$annot->{'not_defined'} = 1;
		}
	}	
	elsif ( $method eq 'corsair' ) { # get corsair
		if ( $result and $result->analysis and $result->analysis->corsair ) {
			my ($analysis) = $result->analysis->corsair;
			if ( $analysis->vertebrate_signal and ($analysis->vertebrate_signal eq $OK_LABEL) ) {
				$annot->{'ok'} = 1;
			}
			elsif ( $analysis->vertebrate_signal and ($analysis->vertebrate_signal eq $UNKNOWN_LABEL) ) {
				$annot->{'unknown'} = 1;
			}
			elsif ( $analysis->vertebrate_signal and ($analysis->vertebrate_signal eq $NO_LABEL) ) {
				$annot->{'rejected'} = 1;
			}
			else {
				$annot->{'not_defined'} = 1;
			}
		}
		else {
			$annot->{'not_defined'} = 1;
		}
	}
	elsif ( $method eq 'crash_sp' ) { # get crash: signalp
		if ( $result and $result->analysis and $result->analysis->crash ) {
			my ($analysis) = $result->analysis->crash;
			if ( $analysis->peptide_signal and ($analysis->peptide_signal eq $NO_LABEL) ) {
				$annot->{'rejected'} = 1;
			}
			elsif ( $analysis->peptide_signal and ($analysis->peptide_signal eq $UNKNOWN_LABEL) ) {
				$annot->{'unknown'} = 1;
			}
			elsif ( $analysis->peptide_signal and ($analysis->peptide_signal eq $OK_LABEL) ) {
				$annot->{'ok'} = 1;
			}
			else {
				$annot->{'not_defined'} = 1;
			}
		} else {
			$annot->{'not_defined'} = 1;
		}
	}
	elsif ( $method eq 'crash_tp' ) { # get crash: targetp
		if ( $result and $result->analysis and $result->analysis->crash ) {
			my ($analysis) = $result->analysis->crash;
			if ( $analysis->mitochondrial_signal and ($analysis->mitochondrial_signal eq $NO_LABEL) ) {
				$annot->{'rejected'} = 1;
			}
			elsif ( $analysis->mitochondrial_signal and ($analysis->mitochondrial_signal eq $UNKNOWN_LABEL) ) {
				$annot->{'unknown'} = 1;
			}
			elsif ( $analysis->mitochondrial_signal and ($analysis->mitochondrial_signal eq $OK_LABEL) ) {
				$annot->{'ok'} = 1;
			}
			else {
				$annot->{'not_defined'} = 1;
			}
		} else {
			$annot->{'not_defined'} = 1;
		}
	}
	elsif ( $method eq 'thump' ) { # get thump
		if ( $result and $result->analysis and $result->analysis->thump ) {
			my ($analysis) = $result->analysis->thump;
			if ( $analysis->transmembrane_signal and ($analysis->transmembrane_signal eq $NO_LABEL) ) {
				$annot->{'rejected'} = 1;
			}
			elsif ( $analysis->transmembrane_signal and ($analysis->transmembrane_signal eq $UNKNOWN_LABEL) ) {
				$annot->{'unknown'} = 1;
			}
			elsif ( $analysis->transmembrane_signal and ($analysis->transmembrane_signal eq $OK_LABEL) ) {
				$annot->{'ok'} = 1;
			}
			else {
				$annot->{'not_defined'} = 1;
			}
		} else {
			$annot->{'not_defined'} = 1;
		}
	}
		
	return $annot;
}

# CRASH: get annotation from SignalP and TargetP
sub get_crash_annots($$$)
{
	my ($method, $transcript_id, $reports) = @_;
	my ($signal) = 0;
	my ($annot);

	# Scan every transcript if there is a pep-mit signal	
	foreach my $rst (@{$reports->transcripts}) {	
		my ($all_annots) = get_method_annot($method, $rst);
		if ( ($all_annots->{'ok'} == 1) ) {
				$signal = 1;			
		}
	}
	# Get consensus result for current transcript:
	# if there is signal for both methods, we take the consensus annotation. Otherwise, we don't know
	if ( $signal == 1 ) {
		my ($index) = $reports->{'_index_transcripts'}->{$transcript_id};
		my ($result) = $reports->transcripts->[$index];		
		$annot = get_crash_annot($method, $result);
	}
	else {
		$annot = $UNKNOWN_LABEL;
	}
	
	return ($annot,$signal);
		
} # End get_crash_annots

sub get_crash_annot($$)
{
	my ($method, $result) = @_;
	my ($annot);

	# get crash: signalp or targetp
	my ($m_annot) = get_method_annot($method, $result);
	
	# get consensus
	if ( $m_annot->{'not_defined'} == 1	) { # we don't results
		$annot = undef;		
	}
	elsif ( $m_annot->{'rejected'} == 1	) { # reject if all methods reject
		$annot = $NO_LABEL;
	}
	elsif ( $m_annot->{'ok'} == 1	) {
		$annot = $OK_LABEL;
	}
	else { # other cases, we don't know
		$annot = $UNKNOWN_LABEL;
	}
	
	return $annot;
		
} # End get_crash_annot

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