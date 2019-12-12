#/usr/bin/env perl -w

use strict;

use File::stat;
use Set::Scalar;


my ($ref_dir_file, $alt_dir_file, $output_file) = @ARGV;

# Cutoff time to distinguish newly created output files.
my $cutoff_time = 1575936000;  # midnight on 2019-12-10

my %dir_files = (
	'ref' => $ref_dir_file,
	'alt' => $alt_dir_file,
);

my %cmp_info;
while (my($cmp_key, $dir_file) = each %dir_files) {
	open DIRFILE, $dir_file
		or die("failed to open file '${dir_file}'");
	while (<DIRFILE>) {
		
		chomp $_;
		my $gene_dir = $_;
		my $appris_file = $gene_dir.'/appris';

		if ( -f $appris_file ) {
			
			if (stat($appris_file)->mtime < $cutoff_time) {
				warn("skipping old APPRIS output file '${appris_file}'");
				next;
			}
		
			open my $fh, $appris_file
				or die("failed to open file '${appris_file}'");
			while (<$fh>) {
				my $line = $_;
				chomp $line;
				my @fields = split("\t", $line);
				if ( scalar(@fields) == 20 ) {  # i.e. full APPRIS output record
					my $gene_id = $fields[0];
					my $transc_id = $fields[2];
					my %transc_info = (
						'gene_name' => $fields[1],
						'pi_label' => $fields[19],
						'ccds_id' => $fields[7]
					);
					$cmp_info{$gene_id}{$cmp_key}{$transc_id} = \%transc_info;
				}
			}
			close $fh;
		}
	}
	close DIRFILE;
}

my %diff_info;
while (my($gene_id, $gene_info) = each %cmp_info) {

	my @transc_ids = keys %{$cmp_info{$gene_id}{'ref'}};
	my %p1_set_info = (
		'ref' => Set::Scalar->new(),
		'alt' => Set::Scalar->new()
	);
	
	foreach my $transc_id (@transc_ids) {
		foreach my $cmp_key (('ref', 'alt')) {
			my $pi_label = $cmp_info{$gene_id}{$cmp_key}{$transc_id}{'pi_label'};
			if ($pi_label eq 'PRINCIPAL:1') {
				$p1_set_info{$cmp_key}->insert($transc_id);
			}
		}
	}
	
	if ( ! $p1_set_info{'ref'}->is_equal($p1_set_info{'alt'}) ) {
		foreach my $transc_id (@transc_ids) {
			my $pi_label_ref = $cmp_info{$gene_id}{'ref'}{$transc_id}{'pi_label'};
			my $pi_label_alt = $cmp_info{$gene_id}{'alt'}{$transc_id}{'pi_label'};
			if ( grep(/^PRINCIPAL:1$/, ($pi_label_ref, $pi_label_alt)) ) {
				my %transc_info = (
					'gene_id' => $gene_id,
					'gene_name' => $cmp_info{$gene_id}{'ref'}{$transc_id}{'gene_name'},
					'transc_id' => $transc_id,
					'pi_label_ref' => $pi_label_ref,
					'pi_label_alt' => $pi_label_alt,
					'ccds_id' => $cmp_info{$gene_id}{'ref'}{$transc_id}{'ccds_id'}
				);
				$diff_info{$gene_id}{$transc_id} = \%transc_info;
			}
		}
	}
}

open(OUTFILE, '>', $output_file)
	or die("failed to open file '${output_file}'");
my @field_names = ('gene_id', 'gene_name', 'transc_id', 'pi_label_ref',
				   'pi_label_alt', 'ccds_id');
print OUTFILE join("\t", @field_names)."\n";
foreach my $gene_info (values %diff_info) {
	foreach my $transc_info (values %{$gene_info}) {
		my @row = (
			$transc_info->{'gene_id'},
			$transc_info->{'gene_name'},
			$transc_info->{'transc_id'},
			$transc_info->{'pi_label_ref'},
			$transc_info->{'pi_label_alt'},
			$transc_info->{'ccds_id'}
		); 
		my $line = join("\t", @row);
		print OUTFILE $line."\n";
	}
}
close OUTFILE;


__END__

=head1 NAME

cmp_scores.pl

=head1 DESCRIPTION

Compare annotations under different scoring methods, summarise differences.

=head2 Required arguments:
	
	ref_dir_file	Text file listing annotation directories with reference method scores.

	alt_dir_file	Text file listing annotation directories with alternative method scores.

	output_file		TSV file summarising transcripts with different annotations
					under the different scoring methods.

=head1 EXAMPLE

perl cmp_scores.pl dirs_ref.txt dirs_alt.txt score_cmp.tab
	
=cut

