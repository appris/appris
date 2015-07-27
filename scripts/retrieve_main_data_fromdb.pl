#!/usr/bin/perl -W

use strict;
use warnings;
use Getopt::Long;
use FindBin;

use APPRIS::Registry;
use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;

use Data::Dumper;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_APPRIS_DB_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_APPRIS_DB_FILE	= $LOCAL_PWD.'/conf/apprisdb.ALL.ini';

# Input parameters
my ($species) = undef;
my ($chr) = undef;
my ($output_file) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'species=s'			=> \$species,
	'chr=s'				=> \$chr,
	'output=s'			=> \$output_file,
	'apprisdb-conf=s'	=> \$apprisdb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);
unless(defined $species and defined $output_file)
{
	print `perldoc $0`;
	exit 1;
}

# Optional arguments
# get vars of appris db
unless ( defined $apprisdb_conf_file ) {
	$apprisdb_conf_file = $CONFIG_INI_APPRIS_DB_FILE;
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
sub get_data_by_chr($$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# APPRIS Registry from given specie
	my ($cfg) = new Config::IniFiles( -file => $apprisdb_conf_file );
	my ($spe) = $species; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	my ($param) = {
			'-dbhost'       => $cfg->val($specie_db, 'host'),
			'-dbname'       => $cfg->val($specie_db, 'db'),
			'-dbuser'       => $cfg->val($specie_db, 'user'),
			'-dbpass'       => $cfg->val($specie_db, 'pass'),
			'-dbport'       => $cfg->val($specie_db, 'port'),
	};
	$logger->debug(Dumper($param)."\n");
	my ($registry) = APPRIS::Registry->new();
	$registry->load_registry_from_db(%{$param});
	
	# Get data by chr
	$logger->debug("-- get_data_by_chr -------\n");
	my ($output_content) = get_data_by_chr($registry, $chr);
	if ($output_content ne '') {
		my ($printing_file_log) = printStringIntoFile($output_content,$output_file);
		die ("ERROR: printing ") unless(defined $printing_file_log);		
	}

	$logger->finish_log();
	
	exit 0;	
	
}
sub get_data_by_chr($$)
{
	my ($registry, $ichr) = @_;
	my ($output) = '';
	my ($chromosome) = $registry->fetch_basic_by_region($ichr);

	foreach my $gene (@{$chromosome}) {
		my ($gene_id) = $gene->stable_id;
		my ($chr) = $gene->chromosome;
		my ($gene_name) = '-';
		if ($gene->external_name) {
			$gene_name = $gene->external_name;	
		}
		$logger->debug("-- $gene_id: ");
		
		my ($g_report);
		my ($len_report);
		my ($len_transcs);		
		foreach my $transcript (@{$gene->transcripts}) {		
			my ($transcript_id) = $transcript->stable_id;
			$logger->debug("$transcript_id ");
			
			my ($transcript_name) = '-';
			if ($transcript->external_name) {
				$transcript_name = $transcript->external_name;	
			}

			if ($transcript->translate) {
				my ($translation_length) = length($transcript->translate->sequence);
				my ($ccds_id) = '-';				
				if ( $transcript->xref_identify ) {
					foreach my $xref_identify (@{$transcript->xref_identify}) {								
						if ($xref_identify->dbname eq 'CCDS') {
							$ccds_id = $xref_identify->id;
						}
					}					
				}
				
				my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id,'appris');				
				if ( $analysis and $analysis->appris and $analysis->appris->principal_isoform_signal ) {
					my ($appris_annot) = $analysis->appris->principal_isoform_signal;
					if ( $appris_annot eq 'YES' ) {
#						$g_report->{'chr'}					= $chr;
#						$g_report->{'gene_id'}				= $gene_id;
#						$g_report->{'gene_name'}			= $gene_name;
#						$g_report->{'transcript_id'}		= $transcript_id;
#						$g_report->{'transcript_name'}		= $transcript_name;
#						$g_report->{'ccds_id'}				= $ccds_id;
#						$g_report->{'translation_length'}	= $translation_length;
#						$g_report->{'appris_annot'}			= 'appris_principal';
						my ($t_rep) = {
									'gene_id'				=> $gene_id,
									'gene_name'				=> $gene_name,
									'transcript_id'			=> $transcript_id,
									'transcript_name'		=> $transcript_name,
									'ccds_id'				=> $ccds_id,
									'translation_length'	=> $translation_length,
									'appris_annot'			=> 'appris_principal',				
						};
						$g_report->{$transcript_id}		= $t_rep;
					}
					elsif ( $appris_annot eq 'UNKNOWN' ) {
						my ($t_rep) = {
									'gene_id'				=> $gene_id,
									'gene_name'				=> $gene_name,
									'transcript_id'			=> $transcript_id,
									'transcript_name'		=> $transcript_name,
									'ccds_id'				=> $ccds_id,
									'translation_length'	=> $translation_length,
									'appris_annot'			=> 'appris_candidate',				
						};
						$g_report->{$transcript_id}		= $t_rep;
						push(@{$len_report->{$translation_length}}, $transcript_id);
						
#						unless ( defined $g_report ) {
#							$g_report->{'chr'}					= $chr;
#							$g_report->{'gene_id'}				= $gene_id;
#							$g_report->{'gene_name'}			= $gene_name;
#							$g_report->{'transcript_id'}		= $transcript_id;
#							$g_report->{'transcript_name'}		= $transcript_name;
#							$g_report->{'ccds_id'}				= $ccds_id;
#							$g_report->{'translation_length'}	= $translation_length;
#							$g_report->{'appris_annot'}			= 'LONGEST';
#						}
#						else {
#							if ( $translation_length > $g_report->{'translation_length'} ) { # maximum length
#								$g_report->{'chr'}					= $chr;
#								$g_report->{'gene_id'}				= $gene_id;
#								$g_report->{'gene_name'}			= $gene_name;
#								$g_report->{'transcript_id'}		= $transcript_id;
#								$g_report->{'transcript_name'}		= $transcript_name;
#								$g_report->{'ccds_id'}				= $ccds_id;
#								$g_report->{'translation_length'}	= $translation_length;
#								$g_report->{'appris_annot'}			= 'LONGEST';
#							}
#						}
					}
				}
			}
		}
		
		# get longest transc
		if ( defined $len_report ) {
			my (@long_transc) = sort {$b <=> $a} keys %{$len_report};
			my ($long_seq) = $long_transc[0];
			foreach my $transc_id (@{$len_report->{$long_seq}}) {
				$len_transcs->{$transc_id} = $long_seq;
			}
		}
		
		# print annots
		if ( defined $g_report ) {
			while (my ($transc_id, $t_rep) = each(%{$g_report}) ) {
				my ($annot) = $t_rep->{'appris_annot'};
				if ( defined $len_transcs and exists $len_transcs->{$transc_id} ) {
					$annot = 'appris_candidate_longest';
				}
				$output .=	$t_rep->{'gene_name'}."\t".
							$t_rep->{'gene_id'}."\t".
							$transc_id."\t".
							$t_rep->{'ccds_id'}."\t".
							$annot."\n";				
			}
		}		
		$logger->debug("\n");

		
#		if ( defined $g_report ) {
#			$output .=	$g_report->{'chr'}."\t".
#						$g_report->{'gene_id'}."\t".
#						$g_report->{'gene_name'}."\t".
#						$g_report->{'transcript_id'}."\t".
#						$g_report->{'transcript_name'}."\t".
#						$g_report->{'ccds_id'}."\t".
#						$g_report->{'appris_annot'}."\n";
#		}		
#		$logger->debug("\n");
	}
	return $output;
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

	--species= <Name of species -mammals->
	
	--output <Output file that has the main isoforms>	
	
=head2 Optional arguments:

	--chr  <Genomic region>

	--apprisdb-conf <Config file of APPRIS database>
				
=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_main_data.pl
	
		--species='Homo sapiens'
		
		--chr=21

		--output=../features/data/appris.results.rel7.v1.chr21.main.tsv
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
