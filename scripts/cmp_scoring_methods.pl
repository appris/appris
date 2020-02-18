#/usr/bin/env perl

use 5.14.0;
use strict;
use warnings::register;

use Cwd;
use File::Find;
use Getopt::Long;
use Set::Scalar;

my ($a_file) = undef;
my ($b_file) = undef;
my ($ref_pi_label) = 'PRINCIPAL:1';
my ($out_file) = undef;

&GetOptions(
	'a-file|a=s' => \$a_file,
	'b-file|b=s' => \$b_file,
	'pi-label:s' => \$ref_pi_label,
	'out-file|o=s' => \$out_file
);

my %file_map = (
	'a' => $a_file,
	'b' => $b_file
);

my @known_pi_labels = (
	'PRINCIPAL:1',
	'PRINCIPAL:2',
	'PRINCIPAL:3',
	'PRINCIPAL:4',
	'PRINCIPAL:5',
	'ALTERNATIVE:1',
	'ALTERNATIVE:2',
	'MINOR'
);

if ( ! grep { $_ eq $ref_pi_label } @known_pi_labels ) {
	die("unknown APPRIS label: '$ref_pi_label'");
}

my %cmp_info;
my @known_opt_col_names = ('ccds_id', 'mane_select');
my %opt_col_map = map { $_ => Set::Scalar->new() } @known_opt_col_names;
while (my ($trial_key, $trial_file) = each(%file_map) ) {

	open my $fh, $trial_file
	  or die("failed to open file: '$trial_file'");

	my $first_line = <$fh>;
	chomp $first_line;
	my @field_names = split("\t", $first_line);
	my %name_to_idx = map { $field_names[$_],  $_ } 0 .. $#field_names;

	my @file_opt_col_names = grep { exists($name_to_idx{$_}) } @known_opt_col_names;
	foreach my $col_name (@file_opt_col_names) {
		$opt_col_map{$col_name}->insert($trial_key);
	}

	while (<$fh>) {
		my $line = $_;
		chomp $line;

		my @fields = split("\t", $line);
		my $gene_id = $fields[ $name_to_idx{'gene_id'} ];
		my $transc_id = $fields[ $name_to_idx{'transc_id'} ];
		$cmp_info{$gene_id}{$trial_key}{$transc_id} = {
			'pi_label' => $fields[ $name_to_idx{'pi_label'} ]
		};

		foreach my $col_name (@file_opt_col_names) {
			my $col_value = $fields[ $name_to_idx{$col_name} ];
			if ( $col_value && $col_value ne '-' ) {
				$cmp_info{$gene_id}{$trial_key}{$transc_id}{$col_name} = $col_value;
			}
		}
	}
	close $fh;
}

my %diff_info;
my @obs_opt_col_names = grep { $opt_col_map{$_}->size > 0 } @known_opt_col_names;
while (my($gene_id, $gene_info) = each %cmp_info) {

	next unless ( exists $cmp_info{$gene_id}{'a'} && exists $cmp_info{$gene_id}{'b'} );
	my $transc_set_a = Set::Scalar->new(keys %{$cmp_info{$gene_id}{'a'}});
	my $transc_set_b = Set::Scalar->new(keys %{$cmp_info{$gene_id}{'b'}});
	next unless ( $transc_set_a->is_equal($transc_set_b) );

	my @transc_ids = keys %{$cmp_info{$gene_id}{'a'}};
	my %p1_set_info = (
		'a' => Set::Scalar->new(),
		'b' => Set::Scalar->new()
	);

	foreach my $transc_id (@transc_ids) {
		foreach my $trial_key (('a', 'b')) {
			my $pi_label = $cmp_info{$gene_id}{$trial_key}{$transc_id}{'pi_label'};
			if ($pi_label eq $ref_pi_label) {
				$p1_set_info{$trial_key}->insert($transc_id);
			}
		}
	}

	if ( ! $p1_set_info{'a'}->is_equal($p1_set_info{'b'}) ) {
		foreach my $transc_id (@transc_ids) {
			my $pi_label_a = $cmp_info{$gene_id}{'a'}{$transc_id}{'pi_label'};
			my $pi_label_b = $cmp_info{$gene_id}{'b'}{$transc_id}{'pi_label'};
			if ( grep { $_ eq $ref_pi_label } ($pi_label_a, $pi_label_b) ) {
				my %transc_info = (
					'gene_id' => $gene_id,
					'transc_id' => $transc_id,
					'pi_label_a' => $pi_label_a,
					'pi_label_b' => $pi_label_b
				);

				foreach my $col_name (@obs_opt_col_names) {
					my @trial_keys = $opt_col_map{$col_name}->members;
					my $col_value;
					foreach my $trial_key (@trial_keys) {
						if ( ! (exists($cmp_info{$gene_id}{$trial_key}{$transc_id}{$col_name})
								&& defined($cmp_info{$gene_id}{$trial_key}{$transc_id}{$col_name})) ) {
							next;
						}
						my $trial_col_value = $cmp_info{$gene_id}{$trial_key}{$transc_id}{$col_name};
						if ( ! defined $col_value ) {
							$col_value = $trial_col_value;
						} elsif ( $trial_col_value ne $col_value ) {
							if ( $col_name eq 'ccds_id' and
								 ( $col_value =~ /\.\d+$/ xor $trial_col_value =~ /\.\d+$/ ) ) {
								my $bare_ccds_id = $col_value =~ s/\.\d+$//r;
								my $bare_trial_ccds_id = $trial_col_value =~ s/\.\d+$//r;
								if ( $bare_trial_ccds_id eq $bare_ccds_id ) {
									if ( $trial_col_value =~ /\.\d+$/ ) {
										$col_value = $trial_col_value;  # keep more specific CCDS ID
									}
									next;
								}
							}
							die("inconsistent values in '${col_name}' column (${trial_col_value} vs ${col_value})");
						}
					}
					$transc_info{$col_name} = $col_value;
				}

				$diff_info{$gene_id}{$transc_id} = \%transc_info;
			}
		}
	}
}

if ( scalar(keys %diff_info) > 0 ) {

	open(my $fh, '>', $out_file)
		or die("failed to open file '${out_file}'");
	my @field_names = ('gene_id', 'transc_id', 'pi_label_a', 'pi_label_b');
	push(@field_names, @obs_opt_col_names);
	print $fh join("\t", @field_names)."\n";
	foreach my $gene_info (values %diff_info) {
		foreach my $transc_info (values %{$gene_info}) {

			my @row = (
				$transc_info->{'gene_id'},
				$transc_info->{'transc_id'},
				$transc_info->{'pi_label_a'},
				$transc_info->{'pi_label_b'}
			);

			my @opt_col_values = map { $transc_info->{$_} } @obs_opt_col_names;
			push(@row, @opt_col_values);

			print $fh join("\t", @row)."\n";
		}
	}
	close $fh;
}
