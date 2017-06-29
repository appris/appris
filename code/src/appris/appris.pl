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

use lib "$FindBin::Bin";
use appris;

###################
# Global variable #
###################

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($data_file) = undef;
my ($transcripts_file) = undef;
my ($translations_file) = undef;
my ($firestar_file) = undef;
my ($matador3d_file) = undef;
my ($matador3d2_file) = undef;
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
	'matador3d2=s'		=> \$matador3d2_file,
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
our $APPRIS_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'cutoff');
our $APPRIS_METHODS			= $cfg->val( 'APPRIS_VARS', 'methods');
our $FIRESTAR_MINRES		= $cfg->val( 'APPRIS_VARS', 'firestar_minres');
our $FIRESTAR_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'firestar_cutoff');
our $MATADOR3D_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'matador3d_cutoff');
our $MATADOR3D2_CUTOFF		= $cfg->val( 'APPRIS_VARS', 'matador3d2_cutoff');
our $SPADE_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'spade_cutoff');
our $CORSAIR_AA_LEN_CUTOFF	= $cfg->val( 'APPRIS_VARS', 'corsair_aa_cutoff');
our $CORSAIR_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'corsair_cutoff');
our $THUMP_CUTOFF			= $cfg->val( 'APPRIS_VARS', 'thump_cutoff');

# Required arguments
unless ( defined $config_file and 
		defined $translations_file and 
		defined $firestar_file and defined $matador3d_file and defined $matador3d2_file and defined $corsair_file and defined $spade_file and
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

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Local variables
	my ($firestar_result);
	my ($matador3d_result);
	my ($matador3d2_result);
	my ($corsair_result);
	my ($spade_result);
	my ($thump_result);
	my ($crash_result);
	my ($inertia_result);
	my ($proteo_result);

	# determine the methods involved in the final decision
	# By default, we use Matador3D2 but If we have genome coordinates, then we use Matador3D.
	my ($involved_methods) = $APPRIS_METHODS;
	if ( defined $data_file and defined $transcripts_file and defined $translations_file ) {
		$involved_methods =~ s/matador3d2/matador3d/g;
	}

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
	$logger->info("get matador3d2 result\n");
	if ( -e $matador3d2_file and (-s $matador3d2_file > 0) ) {
		$matador3d2_result = getStringFromFile($matador3d2_file);
		unless ( defined $matador3d2_result ) {
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
	my ($reports) = parse_appris_methods($gene, $firestar_result, $matador3d_result, $matador3d2_result, $spade_result, $corsair_result, $crash_result, $thump_result, $inertia_result, $proteo_result);
	$logger->debug("REPORTS:\n".Dumper($reports)."\n");

	# get scores of methods for each transcript
	$logger->info("-- get scores of methods for each variant\n");
	my ($scores,$s_scores) = appris::get_method_scores($gene, $reports, $involved_methods);
	$logger->debug("PRE_SCORES:\n".Dumper($scores)."\n");
	$logger->debug("PRE_S_SCORES:\n".Dumper($s_scores)."\n");
	
	# get annots of methods for each transcript
	$logger->info("-- get annots of methods for each variant\n");
	my ($annots) = appris::get_method_annots($gene, $s_scores, $involved_methods);
	$logger->debug("PRE_ANNOTS:\n".Dumper($annots)."\n");

	# get scores/annots of appris for each transcript
	$logger->info("--  get scores/annots of appris for each transcript\n");
	my ($nscores) = appris::get_final_scores($gene, $annots, $involved_methods, $scores, $s_scores);
	$logger->debug("PRE2_SCORES:\n".Dumper($scores)."\n");
	$logger->debug("PRE2_S_SCORES:\n".Dumper($s_scores)."\n");
	$logger->debug("PRE2_NSCORES:\n".Dumper($nscores)."\n");

	# get annotations indexing each transcript
	$logger->info("-- get final annotations\n");
	appris::get_final_annotations($gene, $scores, $s_scores, $nscores, $annots);
	$logger->debug("ANNOTS:\n".Dumper($annots)."\n");
	
	# print outputs
	my ($score_content) = appris::get_score_output($gene, $scores, $annots);
	my ($p_score_out) = APPRIS::Utils::File::printStringIntoFile($score_content, $output_main_file);
	unless( defined $p_score_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($nscore_content) = appris::get_nscore_output($gene, $nscores);
	my ($p_nscore_out) = APPRIS::Utils::File::printStringIntoFile($nscore_content, $output_nscore_file);
	unless( defined $p_nscore_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	my ($label_content) = appris::get_label_output($gene, $annots);
	my ($p_label_out) = APPRIS::Utils::File::printStringIntoFile($label_content, $output_label_file);
	unless( defined $p_label_out ) {
		$logger->error("Can not create output trace file: $!\n");
	}
	
	$logger->finish_log();
	
	exit 0;	
	
}

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