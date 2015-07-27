#!/usr/bin/perl -W

use strict;
use warnings;
use FindBin;
use Config::IniFiles;
use Getopt::Long;
use Data::Dumper;

use APPRIS::EnsEMBL;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD	
	$ENSEMBL_CONFIG_FILE
);

$LOCAL_PWD				= $FindBin::Bin;
$ENSEMBL_CONFIG_FILE	= $LOCAL_PWD.'/../../conf/ensembldb.ini';


# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($conf_file) = undef;
my ($id) = undef;
my ($species) = undef;
my ($e_version) = undef;
my ($out_data_file) = undef;
my ($out_pdata_file) = undef;
my ($out_transc_file) = undef;
my ($out_transl_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'id=s'					=> \$id,
	'species=s'				=> \$species,
	'e-version=s'			=> \$e_version,
	'out-data=s'			=> \$out_data_file,
	'out-pdata=s'			=> \$out_pdata_file,
	'out-transcripts=s'		=> \$out_transc_file,	
	'out-translations=s'	=> \$out_transl_file,	
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);

# Required arguments
unless (
	defined $id and
	defined $species and
	defined $e_version and
	defined $out_data_file and
	defined $out_pdata_file and
	defined $out_transc_file and
	defined $out_transl_file
){
	print `perldoc $0`;
	exit 1;
}

# Get conf file
my ($cfg) = new Config::IniFiles( -file =>  $ENSEMBL_CONFIG_FILE );
	
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
	my ($gtf_annot) = '';
	my ($gene_gtf_annot) = '';
	my ($transc_gtf_annot) = '';
	my ($transc_seq) = '';
	my ($transl_seq) = '';
	my ($gtf_pannot) = '';
	
	$logger->info("-- get given inputs\n");
	$logger->info("-- id: $id\n");
	$logger->info("-- species: $species\n");
		
	$logger->info("-- creates a new ensembl registry object\n");
	my ($ensembl) = APPRIS::EnsEMBL->new(
										-host		=> $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host'),
										-user		=> $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user'),
										-pass		=> $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'pass'),
										-verbose	=> $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose'),
										-db_version	=> $e_version,										
										-species	=> $species,										
	);
	$logger->debug(
		"\t-host		=> ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host')."\n".
		"\t-user		=> ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user')."\n".
		"\t-verbose	=> ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose')."\n".
		"\t-db_version	=> ".$e_version."\n".										
		"\t-species	=> ".$species."\n"
	);	
	#$logger->debug("\n".Dumper($ensembl)."\n");
	
	$logger->info("-- fetch by stable gene id\n");	
	my ($gene);
	if ( $id =~ /^ENS/ ) {
		my ($old_id) = $id; $old_id =~ s/\.\d*$//;
		$gene = $ensembl->fetch_by_stable_id($old_id);
	}
	else {
		 $gene = $ensembl->fetch_all_by_external_name($id);
	}
	#$logger->debug(Dumper($gene));

	if ( defined $gene and defined $gene->stable_id() ) {
		my ($gene_id) = $gene->stable_id();
		
		$logger->info("-- get gtf annot of gene: ".$gene_id."\n");
		$gene_gtf_annot = $ensembl->get_gtf_by_gene_adaptor($gene);
		$gtf_annot .= $gene_gtf_annot;
		
		my (@transcripts) = @{ $ensembl->get_all_transcripts($gene) };
		foreach my $transcript (@transcripts) {
			if ( defined $transcript and defined $transcript->stable_id() ) {
				my ($transcript_id) = $transcript->stable_id();
				
				$logger->info("-- get gtf annot of transcript: ".$transcript_id."\n");
				$transc_gtf_annot .= $ensembl->get_gtf_by_transc_adaptor($transcript);
				
				$logger->info("-- get gtf annot of exons transcript: ".$transcript_id."\n");
				$transc_gtf_annot .= $ensembl->get_gtf_exon_by_transc_adaptor($transcript);
				
				$logger->info("-- get transc seq of transcript: ".$transcript_id."\n");
				$transc_seq .= $ensembl->get_transc_seq_by_transc_adaptor($transcript);
				
				if ( $transcript->translation() ) {

					$logger->info("-- get gtf annot of start codon transcript: ".$transcript_id."\n");
					$transc_gtf_annot .= $ensembl->get_gtf_start_codon_by_transc_adaptor($transcript);
										
					$logger->info("-- get gtf annot of cds transcript: ".$transcript_id."\n");
					$transc_gtf_annot .= $ensembl->get_gtf_cds_by_transc_adaptor($transcript);
					
					$logger->info("-- get gtf annot of stop codon transcript: ".$transcript_id."\n");
					$transc_gtf_annot .= $ensembl->get_gtf_stop_codon_by_transc_adaptor($transcript);
					
					$logger->info("-- get transl seq of transcript: ".$transcript_id."\n");
					$transl_seq .= $ensembl->get_transl_seq_by_transc_adaptor($transcript);
										
					$logger->info("-- get gtf annot of protein: ".$transcript_id."\n");
					$gtf_pannot .= $ensembl->get_gtf_prot_by_transc_adaptor($transcript);
				}
				else {
					$logger->info("-- non-coding transcript\n");
				}

			}
			else {
				$logger->error("-- undefined transcript\n");
			}
		}
		$gtf_annot .= $transc_gtf_annot;
	}
	else {
		$logger->error("-- undefined gene\n");
	}	

	# print information
	if ( $gtf_annot ne '' ) {
		my ($printing_file_log) = printStringIntoFile($gtf_annot, $out_data_file);
		$logger->error("-- printing gtf annot") unless ( defined $printing_file_log );
	}
	if ( $transc_seq ne '' ) {
		my ($printing_file_log) = printStringIntoFile($transc_seq, $out_transc_file);
		$logger->error("-- printing transc seq") unless ( defined $printing_file_log );
	}
	if ( $transl_seq ne '' ) {
		my ($printing_file_log) = printStringIntoFile($transl_seq, $out_transl_file);
		$logger->error("-- printing transl seq") unless ( defined $printing_file_log );
	}
	if ( $gtf_pannot ne '' ) {
		my ($printing_file_log) = printStringIntoFile($gtf_pannot, $out_pdata_file);
		$logger->error("-- printing gtf annot") unless ( defined $printing_file_log );
	}
	

	$logger->finish_log();
	
	exit 0;	
}

main();


1;

__END__

=head1 NAME

getGTF

=head1 DESCRIPTION

get GTF annotations, transcripts, and translations from Ensembl gene identifier and species. 

=head2 Required arguments:

	--id= <Ensembl gene identifier>
	
	--species= <Name of species -mammals->
	
	--e-version= <Number of Ensembl version of identifier>
	
	--out-data= <Output data, GTF file>
	
	--out-pdata= <Output data pf protein, GTF file>
	
	--out-transcripts= <Transcript sequences, fasta format>
	
	--out-translations= <Protein sequences, fasta format>

=head2 Optional arguments:
	
	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl getGTF.pl
	
		--id=ENSG00000140416
		
		--species='Homo sapiens'
		
		--e-version=69
		
		--out-data=../../examples/ENSG00000140416/ENSG00000140416.annot.gtf
		
		--out-pdata=../../examples/ENSG00000140416/ENSG00000140416.pannot.gtf
		
		--out-transcripts=../../examples/ENSG00000140416/ENSG00000140416.transc.fa
		
		--out-translations=../../examples/ENSG00000140416/ENSG00000140416.transl.fa
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
