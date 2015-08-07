#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Parser qw( parse_gencode );
use APPRIS::Utils::File qw( updateStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################
use vars qw(
	$SEQ_ID
	$SOURCE
	$FEATURE
	$CONFIG_INI_ENSEMBL_DB_FILE
);

$SEQ_ID						= 'SEQ';
$SOURCE						= 'Ensembl';
$FEATURE					= 'Protein';
#$CONFIG_INI_ENSEMBL_DB_FILE	= $FindBin::Bin.'/../conf/ensembldb.ini';
$CONFIG_INI_ENSEMBL_DB_FILE	= $ENV{APPRIS_CODE_CONF_DIR}.'/ensembl.ini';


# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($species) = undef;
my ($e_version) = undef;
my ($data_file) = undef;
my ($translations_file) = undef;
my ($output_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'species=s'			=> \$species,
	'e-version=s'		=> \$e_version,
	'data=s'			=> \$data_file,
	'translations=s'	=> \$translations_file,
	'output=s'			=> \$output_file,	
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined $species and
	defined $e_version and
	defined $data_file and
	defined $translations_file and
	defined $output_file
){
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

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Load Ensembl registry
	$logger->info("-- creates a new ensembl registry object\n");
	my ($cfg) = new Config::IniFiles( -file => $CONFIG_INI_ENSEMBL_DB_FILE );
	my ($e_core_param) = {
			'-host'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host'),
			'-user'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user'),
			'-verbose'    => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose'),
			'-db_version' => $e_version,
			'-species'    => $species,
	};
	$logger->debug(
		"\t-host		=> ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host')."\n".
		"\t-user		=> ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user')."\n".
		"\t-verbose	=> ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose')."\n".
		"\t-db_version	=> ".$e_version."\n".										
		"\t-species	=> ".$species."\n"
	);
	my ($registry) = 'Bio::EnsEMBL::Registry';
	eval {
		$registry->load_registry_from_db(%{$e_core_param});
	};
	$logger->error("can not load ensembl registry: $!\n") if $@;
		
	# Get gencode data
	$logger->info("##get:genconde:data ---------------\n");
	my ($gencode_data) = parse_gencode($data_file, $translations_file, $registry);
#$logger->debug("##DATA:\n".Dumper($gencode_data)."\n");	
	unless ( defined $gencode_data ) {
		$logger->error("can not parse gencode data: $!\n");
	}
	else {		
		$logger->info("get gencode\n");
		# Get transcript/transtlation sequence per gene
		foreach my $gene (@{$gencode_data}) {
			my ($output_transl_cont) = '';
			my ($output_cds_anot_cont) = '';
			my ($chr) = $gene->chromosome;
			my ($gene_id) = $gene->stable_id;
			$logger->info("$gene_id ---------------\n");
						
			# scan transctipt/translation/cds seq/cds info
			foreach my $transcript (@{$gene->transcripts}) {		
				my ($transcript_id) = $transcript->stable_id;
				$logger->info("$transcript_id ");
				
				if ( $transcript->translate ) {
					my ($translate) = $transcript->translate;

					if ( $transcript->exons and $translate->cds_sequence ) {
						my ($exons) = $transcript->exons;						
	
						for (my $icds = 0; $icds < scalar(@{$translate->cds_sequence}); $icds++) {
							my ($cds) = $translate->cds->[$icds];
							my ($exon) = $exons->[$icds];
							my ($exon_id) = $exon->stable_id; $exon_id = '-' unless (defined $exon_id);
							my ($pro_cds) = $translate->cds_sequence->[$icds];
		
							my ($cds_start) = $cds->start;
							my ($cds_end) = $cds->end;
							my ($cds_strand) = $cds->strand;
							my ($cds_phase) = $cds->phase;
							
							my ($pro_cds_start) = $pro_cds->start;
							my ($pro_cds_end) = $pro_cds->end;
							my ($pro_cds_end_phase) = $pro_cds->end_phase;
							my ($pro_cds_seq) = $pro_cds->sequence;
							
							# Delete the residue that is shared by two CDS						
							#if (defined $pro_cds_end_phase and $pro_cds_end_phase != 0) {
							#	$pro_cds_seq =~ s/\w{1}$//;
							#}
							
							# get CDS coordinates (GFF)
							if (defined $pro_cds_start and defined $pro_cds_end) {
								$output_cds_anot_cont .=	"$SEQ_ID\t".
															"$SOURCE\t".
															"$FEATURE\t".
															"$pro_cds_start\t".
															"$pro_cds_end\t".
															".\t".
															".\t".
															"$pro_cds_end_phase\t".
															"ID=$exon_id;Parent=$transcript_id;Gene=$gene_id\n";					
							}
						}
					}
				}			
			}
			# print individual file of cds annotations
			if ( $output_cds_anot_cont ne '' ) {
				my ($printing_file_log) = updateStringIntoFile($output_cds_anot_cont, $output_file);
				$logger->error("Printing output") unless ( defined $printing_file_log );
			}
	
			$logger->info("finished ---------------\n");
		}
	}
	
	$logger->finish_log();
	
	exit 0;	
}

main();


1;

__END__

=head1 NAME

getProCDS

=head1 DESCRIPTION

get protein coordinates (at genomic level) from GTF annotation, and protein sequence. 

=head2 Required arguments:

	--species= <Name of species -mammals->
	
	--e-version= <Number of Ensembl version of identifier>
	
	--data=  <GTF data file>
	
	--translations=  <Translations sequence file>
	
	--output= <Output file>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl getProCDS.pl
	
		--species='Homo sapiens'
		
		--e-version=69
			
		--data=../features/gencode.v7.annotation.gtf
		
		--translations=../features/gencode.v7.pc_translations.fa
		
		--output=../features/id/id.protein_cds.gff
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
