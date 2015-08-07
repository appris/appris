#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use APPRIS::Parser qw( parse_gencode );
use APPRIS::Utils::File qw( prepare_workspace printStringIntoFile updateStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################
use vars qw(
	$SEQ_ID
	$SOURCE
	$FEATURE
);

$SEQ_ID = 'SEQ';
$SOURCE = 'Ensembl';
$FEATURE = 'Protein';


# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($data_file) = undef;
my ($transcripts_file) = undef;
my ($translations_file) = undef;
my ($output_dir) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'data=s'			=> \$data_file,
	'transcripts=s'		=> \$transcripts_file,
	'translations=s'	=> \$translations_file,
	'output_dir=s'		=> \$output_dir,	
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $data_file and defined $transcripts_file and defined $translations_file and defined $output_dir )
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

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Get gencode data
	$logger->info("##get:genconde:data ---------------\n");
	my ($gencode_data) = parse_gencode($data_file, $transcripts_file, $translations_file);
#$logger->debug("##DATA:\n".Dumper($gencode_data)."\n");	
	unless ( defined $gencode_data ) {
		$logger->error("can not parse gencode data: $!\n");
	}
	else {		
		$logger->info("get gencode\n");
		# Get transcript/transtlation sequence per gene
		foreach my $gene (@{$gencode_data}) {
			my ($output_transc_cont) = '';
			my ($output_transl_cont) = '';
			my ($output_cds_seq_cont) = '';
			my ($output_cds_anot_cont) = '';
			my ($chr) = $gene->chromosome;
			my ($gene_id) = $gene->stable_id;
			$logger->info("$gene_id ---------------\n");
						
			# scan transctipt/translation/cds seq/cds info
			foreach my $transcript (@{$gene->transcripts}) {		
				my ($transcript_id) = $transcript->stable_id;
				$logger->info("$transcript_id ");
				
				# get transcript sequence				
				if ( $transcript->sequence ) {
					$output_transc_cont .= ">$transcript_id\n";
					$output_transc_cont .= $transcript->sequence."\n";					
				}
				
				if ( $transcript->translate ) {
					my ($translate) = $transcript->translate;

					# get translation sequence
					if ( $translate->sequence ) {
						my ($translate_seq) = $translate->sequence;						
						#if ( length($transcript->translate->sequence) <= 2 ) {
						#	$translate_seq .= 'X';
						#}
						$output_transl_cont .= ">$transcript_id\n";
						$output_transl_cont .= $translate_seq."\n";					
					}
					
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
							
							# get CDS sequence
							if (defined $pro_cds_seq and $pro_cds_seq ne '') {
								
								$output_cds_seq_cont .= ">$transcript_id|$gene_id|$exon_id|$chr|$cds_start|$cds_end|$cds_strand|$cds_phase\n";
								$output_cds_seq_cont .= $pro_cds_seq."\n";							
							}
							
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
			# print individual file of transcript sequences
			if ( $output_transc_cont ne '' ) {
				my ($workspace) = $output_dir.'/'.'chr'.$chr.'/transcript_sequences';
				$workspace = prepare_workspace($workspace);
				my ($output_file) = $workspace.'/'.$gene_id.'.faa';
				my ($printing_file_log) = printStringIntoFile($output_transc_cont, $output_file);
				$logger->error("Printing output") unless ( defined $printing_file_log );		
			}
			# print individual file of translation sequences
			if ( $output_transl_cont ne '' ) {
				my ($workspace) = $output_dir.'/'.'chr'.$chr.'/peptide_sequences';
				$workspace = prepare_workspace($workspace);
				my ($output_file) = $workspace.'/'.$gene_id.'.faa';
				my ($printing_file_log) = printStringIntoFile($output_transl_cont, $output_file);
				$logger->error("Printing output") unless ( defined $printing_file_log );		
			}
			# print individual file of cds sequences
			if ( $output_cds_seq_cont ne '' ) {
				my ($workspace) = $output_dir.'/'.'chr'.$chr;
				$workspace = prepare_workspace($workspace);				
				my ($output_file) = $workspace.'/'.'protein_cds.'."chr$chr".'.faa';
				my ($printing_file_log) = updateStringIntoFile($output_cds_seq_cont, $output_file);
				$logger->error("Printing output") unless ( defined $printing_file_log );
			}
			# print individual file of cds annotations
			if ( $output_cds_anot_cont ne '' ) {
				my ($workspace) = $output_dir.'/'.'chr'.$chr;
				$workspace = prepare_workspace($workspace);				
				my ($output_file2) = $workspace.'/'.'protein_cds.'."chr$chr".'.gff';
				my ($printing_file_log) = updateStringIntoFile($output_cds_anot_cont, $output_file2);
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

getGencodeEntity

=head1 DESCRIPTION

get GENCODE sequence/annotations

=head1 SYNOPSIS

inertia

=head2 Required arguments:

	--data=  <Gencode data file>
	
	--transcripts=  <Gencode transcript file>
	
	--translations=  <Gencode translations file>
	
	--output_dir= <Directory where output files will save>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl getGencodeEntity.pl
			
		--data=../features/gencode.v7.annotation.gtf
		
		--transcripts=../features/gencode.v7.pc_transcripts.fa
		
		--translations=../features/gencode.v7.pc_translations.fa
		
		--output_dir=../features/sequences/
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
