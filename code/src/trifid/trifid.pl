#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Getopt::Long;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile );


# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($input_file) = undef;
my ($trifid_pred_file) = undef;
my ($output_file) = undef;
my ($loglevel) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;


&GetOptions(
	'input=s'		=> \$input_file,
	'trifid=s'		=> \$trifid_pred_file,
	'output=s'		=> \$output_file,
	'loglevel=s'	=> \$loglevel,
	'logfile=s'		=> \$logfile,
	'logpath=s'		=> \$logpath,
	'logappend'		=> \$logappend,
);

# Required arguments
unless ( defined $input_file and defined $trifid_pred_file and defined $output_file )
{
	print `perldoc $0`;
	exit 1;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE	=> $logfile,
	-LOGPATH	=> $logpath,
	-LOGAPPEND	=> $logappend,
	-LOGLEVEL	=> $loglevel,
);
$logger->init_log($str_params);

# Main subroutine
sub main()
{
	# Declare vars
	my ($output_content) = "";

	my $fasta_object = Bio::SeqIO->new(
		-file => $input_file,
		-format => 'Fasta'
	);

	# obtain the gene id and a hash mapping transcript IDs to their translated sequence
	my ($gene_id);
	my (%transc_to_seq);
	while ( my $seq = $fasta_object->next_seq() )
	{
		if( $seq->id=~/^([^|]+)\|([^|]+)\|([^|]+)/ )
		{
			my ($transc_id) = $2;
			my ($seq_gene_id) = $3;
			if ( $transc_id =~ /^ENS/ ) { $transc_id =~ s/\.\d*$// }
			if ( $seq_gene_id =~ /^ENS/ ) { $seq_gene_id =~ s/\.\d*$// }
			$transc_to_seq{$transc_id} = $seq->seq;

			if ( ! defined($gene_id) ) {
				$gene_id = $seq_gene_id;
			} elsif ( $seq_gene_id ne $gene_id ) {
				$logger->error("gene ID mismatch: $seq_gene_id vs $gene_id\n");
			}
		}
	}

	# retrieve relevant Trifid predictions from tabix-indexed tab file
	if ( defined($gene_id) ) {
		$logger->info("-- retrieve Trifid predictions for gene $gene_id\n");
		my ($cmd) = "tabix -h $trifid_pred_file $gene_id";
		$logger->debug("$cmd\n");

		my ($header, @input_lines) = `$cmd`;
		if (@input_lines) {

			$header =~ s/^#//;
			my @col_names = split(/\t/, $header);
			my ($transc_col) = grep { $col_names[$_] eq "transcript_id" } (0 .. $#col_names);
			my ($seq_col) = grep { $col_names[$_] eq "sequence" } (0 .. $#col_names);

			if ( defined($transc_col) && defined($seq_col) ) {
				# match transcripts by identifier and sequence
				my (@output_lines);
				foreach my $line (@input_lines) {
					my ($transc_id, $trifid_seq) = (split(/\t/, $line))[$transc_col, $seq_col];
					if ( exists($transc_to_seq{$transc_id}) &&
							$trifid_seq eq $transc_to_seq{$transc_id} ) {
						push(@output_lines, $line);
					}
				}

				if (@output_lines) {
					$header =~ s/^/#/;
					unshift(@output_lines, $header);
					$output_content = join("", @output_lines);
				}
			}
		}
	} else {
		$logger->error("failed to retrieve gene_id\n");
	}

	if ( ! $output_content ) {
		$logger->warning("no relevant Trifid predictions found\n");
	}

	# print output
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}

	$logger->finish_log();

	exit 0;
}

main();


__END__

=head1 NAME

trifid

=head1 DESCRIPTION

Obtain Trifid results for the given gene

=head2 Required arguments:

	--input <Fasta sequence file>

	--trifid <Trifid predictions file>

	--output=FILE <Annotation output file>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>

	--logfile=FILE <Log to FILE (default: *STDOUT)>

	--logpath=PATH <Write logfile to PATH (default: .)>

	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl trifid.pl

	--input=examples/ENSG00000254647/transl.fa

	--trifid=examples/trifid_homo_sapiens_e99.tsv.gz

	--output=examples/ENSG00000254647/trifid

=head1 AUTHORS

Original APPRIS method script by Jose Manuel Rodriguez Carrasco,
adapted for TRIFID by Thomas Walsh.

For contact details see the L<APPRIS website|http://appris-tools.org>.

=cut
