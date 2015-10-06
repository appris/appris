#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;

use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################

# Input parameters
my ($species) = undef;
my ($translations_file) = undef;
my ($output_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'species|s=s'		=> \$species,	
	'transl|i=s'		=> \$translations_file,
	'outfile|o=s'		=> \$output_file,
	'loglevel|l=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $species and defined $translations_file and defined $output_file )
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
$logger->init_log();

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->info("-- create uniprot report\n");
	my ($report) = _uniprot_gene_species_report($translations_file);
	$logger->debug("UNIPROT_REPORT:\n".Dumper($report)."\n");
	
	$logger->info("-- create isoforms sequences from given species\n");
	my ($output) = _isoforms_species_seqs($species,$report);
	$logger->debug("UNIPROT_REPORT:\n".Dumper($report)."\n");	
	if ( $output ne '' ) {
		my ($print_log) = printStringIntoFile($output, $output_file);
		$logger->error("printing _isoforms_species_seqs\n") unless(defined $print_log);			
	}
	
	$logger->finish_log();
	
	exit 0;	
}

sub _uniprot_gene_species_report($)
{
	my ($file) = @_;
	my ($report);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			my ($s_id) = $seq->id;
			my ($s_desc) = $seq->desc;
			my ($s_seq) = $seq->seq;
			if ( $s_id =~ /^(sp|tr)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences
				my ($isof_id) = $2;
				my ($gene_name) = $3;
				my ($gene_species);
				my ($gene_id) = $isof_id;
				if ( $isof_id =~ /([^\-]*)/ ) { $gene_id = $1 } 
				my (@desc) = split('OS=', $s_desc);
				if ( scalar(@desc) >= 2 ) {
					if ( $desc[1] =~ /([^\s]*\s[^\s]*)/ ) {
						$gene_species = $1;
					}					
				}
				if ( defined $gene_species ) {
					unless ( exists $report->{$gene_species}->{$gene_id} ) {
						$report->{$gene_species}->{$gene_id} = {
							'name'		=> $gene_name,
							'varsplic'	=> {}
						};
					}
					$report->{$gene_species}->{$gene_id}->{'varsplic'}->{$isof_id} = {
						'desc'		=> $s_desc,
						'seq' 		=> $s_seq
					}
				}
			}
		}		
	}
	return $report;
	
} # End _uniprot_gene_species_report

sub _isoforms_species_seqs($$)
{
	my ($species,$report) = @_;
	my ($output) = '';
	while (my ($aux_species,$spe_rep) = each(%{$report}) ) {
		if ( $species eq $aux_species ) {
			while (my ($gene_id,$gene_rep) = each(%{$spe_rep}) ) {
				my ($gene_name) = $gene_rep->{'name'};
				while (my ($isof_id,$isof_rep) = each(%{$gene_rep->{'varsplic'}}) ) {
					my ($isof_desc) = $isof_rep->{'desc'};
					my ($isof_seq) = $isof_rep->{'seq'};
					$output .= ">sp|$isof_id|$gene_name|$isof_desc\n$isof_seq\n";
				}
			}			
		}
	}
	return $output;
}

main();


1;

__END__

=head1 NAME

getSpeciesGenes

=head1 DESCRIPTION

get GENCODE sequence/annotations

=head1 SYNOPSIS

inertia

=head2 Required arguments:

	-s, --species='Homo sapiens'
	
	-i, --transl=  <UniProt isoforms>
	
	-o, --outfile= <Isoforms for unique species>

=head2 Optional arguments:

	-l, --loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl getSpeciesGenes.pl
			
		-s 'Homo sapiens'
		
		--transl=features/uniprot/uniprot_sprot.all.2015_09.fasta
		
		--outfile=features/uniprot/uniprot_sprot.all.Hsap.2015_09.fasta
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
