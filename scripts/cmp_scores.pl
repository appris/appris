#/usr/bin/env perl -w

use strict;

use File::stat;
use Set::Scalar;


my ($reg_dir_file, $inv_dir_file, $output_file) = @ARGV;

# Cutoff time to distinguish newly created output files.
my $cutoff_time = 1575936000;  # midnight on 2019-12-10

my %dir_files = (
	'regular' => $reg_dir_file,
	'inverted' => $inv_dir_file,
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

	my @transc_ids = keys %{$cmp_info{$gene_id}{'regular'}};
	my %p1_set_info = (
		'regular' => Set::Scalar->new(),
		'inverted' => Set::Scalar->new()
	);
	
	foreach my $transc_id (@transc_ids) {
		foreach my $cmp_key (('regular', 'inverted')) {
			my $pi_label = $cmp_info{$gene_id}{$cmp_key}{$transc_id}{'pi_label'};
			if ($pi_label eq 'PRINCIPAL:1') {
				$p1_set_info{$cmp_key}->insert($transc_id);
			}
		}
	}
	
	if ( ! $p1_set_info{'regular'}->is_equal($p1_set_info{'inverted'}) ) {
		foreach my $transc_id (@transc_ids) {
			my $pi_label_regular = $cmp_info{$gene_id}{'regular'}{$transc_id}{'pi_label'};
			my $pi_label_inverted = $cmp_info{$gene_id}{'inverted'}{$transc_id}{'pi_label'};
			if ( grep(/^PRINCIPAL:1$/, ($pi_label_regular, $pi_label_inverted)) ) {
				my %transc_info = (
					'gene_id' => $gene_id,
					'gene_name' => $cmp_info{$gene_id}{'regular'}{$transc_id}{'gene_name'},
					'transc_id' => $transc_id,
					'pi_label_regular' => $pi_label_regular,
					'pi_label_inverted' => $pi_label_inverted,
					'ccds_id' => $cmp_info{$gene_id}{'regular'}{$transc_id}{'ccds_id'}
				);
				$diff_info{$gene_id}{$transc_id} = \%transc_info;
			}
		}
	}
}

open(OUTFILE, '>', $output_file)
	or die("failed to open file '${output_file}'");
my @field_names = ('gene_id', 'gene_name', 'transc_id', 'pi_label_regular',
				   'pi_label_inverted', 'ccds_id');
print OUTFILE join("\t", @field_names)."\n";
foreach my $gene_info (values %diff_info) {
	foreach my $transc_info (values %{$gene_info}) {
		my @row = (
			$transc_info->{'gene_id'},
			$transc_info->{'gene_name'},
			$transc_info->{'transc_id'},
			$transc_info->{'pi_label_regular'},
			$transc_info->{'pi_label_inverted'},
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
	
	reg_dir_file	Text file listing annotation directories with regular scores.

	inv_dir_file	Text file listing annotation directories with inverted scores.

	output_file		TSV file summarising transcripts with different annotations
					under the different scoring methods.

=head1 EXAMPLE

perl cmp_scores.pl dirs_regular.txt dirs_inverted.txt score_cmp.tab
	
=cut

