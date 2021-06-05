#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Digest::SHA1 qw( sha1_hex );
use Getopt::Long;
use List::Util qw( all any );

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

	# Obtain the query gene ID and a hash mapping each
	# transcript ID to its relevant sequence metadata.
	my ($query_gene_id);
	my (%seq_meta);
	while ( my $seq = $fasta_object->next_seq() )
	{
		my (@match) = ( $seq->id =~ /^([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|([0-9]+)/ );

		if(@match)
		{
			my ($transl_id, $transc_id, $gene_id, $gene_name, $ccds_id, $length_aa) = @match;

			my ($gene_ver);  # TODO: check consistent gene version.
			if ( $gene_id =~ /^(ENS[^.]+)\.([0-9]+)$/ ) {
				($gene_id, $gene_ver) = ($1, $2);
			}

			my ($transc_ver);
			if ( $transc_id =~ /^(ENS[^.]+)\.([0-9]+)$/ ) {
				($transc_id, $transc_ver) = ($1, $2);
			}

			my ($transl_ver);
			if ( $transl_id =~ /^(ENS[^.]+)\.([0-9]+)$/ ) {
				($transl_id, $transl_ver) = ($1, $2);
			}

			$seq_meta{$transc_id} = {
				'transcript_id' => $transc_id,
				'translation_id' => $transl_id,
				'translation_seq_sha1' => sha1_hex($seq->seq),
				'length' => $length_aa
			};
			$seq_meta{$transc_id}{'transcript_version'} = $transc_ver if defined($transc_ver);
			$seq_meta{$transc_id}{'translation_version'} = $transl_ver if defined($transl_ver);

			if ( ! defined($query_gene_id) ) {
				$query_gene_id = $gene_id;
			} elsif ( $gene_id ne $query_gene_id ) {
				$logger->error("gene ID mismatch: $gene_id vs $query_gene_id\n");
			}
		}
	}

	# retrieve relevant Trifid predictions from tabix-indexed tab file
	if ( defined($query_gene_id) ) {
		$logger->info("-- retrieve Trifid predictions for gene $query_gene_id\n");
		my ($cmd) = "tabix -h $trifid_pred_file $query_gene_id";
		$logger->debug("$cmd\n");

		my ($header, @input_lines) = split(/\R/, `$cmd`);
		if (@input_lines) {

			$header =~ s/^#//;
			my @col_names = split(/\t/, $header);
			my (%col_name_set) = map { $_ => 1 } @col_names;

			my (@req_col_names) = ('gene_id', 'gene_name', 'transcript_index', 'transcript_id',
								   'translation_id', 'flags', 'ccdsid', 'appris', 'ann_type',
								   'length', 'trifid_score', 'norm_trifid_score');
			my ($all_req_cols_found) = all { exists($col_name_set{$_}) } @req_col_names;

			my ($seq_attr_col);
			if ( exists($col_name_set{'translation_seq_sha1'}) ) {  # TRIFID
				$seq_attr_col = 'translation_seq_sha1';
			} elsif ( exists($col_name_set{'sequence'}) ) {  # legacy TRIFID
				$seq_attr_col = 'sequence';
			}

			# Output any matching transcripts.
			if ( $all_req_cols_found && defined($seq_attr_col) ) {

				my (@exp_cmp_col_names) = ('transcript_id', 'translation_id', 'translation_seq_sha1',
										   'length', 'transcript_version', 'translation_version');
				my (@cmp_col_names) = grep { exists($col_name_set{$_}) } @exp_cmp_col_names;

				my (@output_lines);
				foreach my $line (@input_lines) {
					my (%rec);
					@rec{@col_names} = split(/\t/, $line);
					my ($transc_id) = $rec{'transcript_id'};

					if ( $seq_attr_col eq 'sequence' ) {
						$rec{'translation_seq_sha1'} = sha1_hex($rec{'sequence'})
					}

					if ( all { ! exists($seq_meta{$transc_id}{$_}) ||
							$rec{$_} eq $seq_meta{$transc_id}{$_} } @cmp_col_names ) {
						push(@output_lines, $line);
					}
				}

				if (@output_lines) {
					$header =~ s/^/#/;
					$output_content = join("\n", ($header, @output_lines)) . "\n";
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

For contact details see the L<APPRIS website|https://appris.bioinfo.cnio.es/>.

=cut
