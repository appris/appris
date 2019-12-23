#/usr/bin/env perl -w

use strict;

use Capture::Tiny 'capture_merged';
use Cwd;
use Data::Dumper;
use File::Find;
use Getopt::Long;
use Set::Scalar;


my ($baseline_dir) = undef;
my ($test_dir) = undef;
my ($baseline_features) = undef;
my ($test_features) = undef;
my ($out_file) = undef;

&GetOptions(
	'baseline-dir=s' => \$baseline_dir,
	'test-dir=s' => \$test_dir,
	'baseline-features=s' => \$baseline_features,
	'test-features=s' => \$test_features,
	'out-file=s' => \$out_file
);

sub _run_system_call($$) {
	my ($cmds, $log_file) = @_;
	my $exit_value;
	open (my $fh, '>', $log_file) or die("failed to open log file");
	print $fh capture_merged { 
		system(@{$cmds});
		if ($? == -1) {
			die("failed to execute: $!\n");
		}
		elsif ($? & 127) {
			die sprintf "child died with signal %d, %s coredump\n",
			            ($? & 127),  ($? & 128) ? 'with' : 'without';
		}
		else {
			$exit_value = $? >> 8;
		}
	};
	close $fh;
	if ( $exit_value != 0 ) {
		die("child exited with value ${exit_value}\n");
	}
}

sub match_ensembl_gene_dir ($) {
	my ($gene_dirs) = @_;
	push(@{$gene_dirs}, $File::Find::fullname) if ( -d && /^ENS(?:[A-Z]{3})?G\d+$/ );
}

my @baseline_gene_dirs;
my @baseline_run_dirs;
if ( ! -d $baseline_dir ) {
	die("baseline directory not found: ${baseline_dir}");
}
find({wanted => sub { match_ensembl_gene_dir(\@baseline_gene_dirs) }, follow => 1}, $baseline_dir);

my @test_gene_dirs;
my @test_run_dirs;
if ( ! -d $test_dir ) {
	die("test directory not found: ${test_dir}");
}
find({wanted => sub { match_ensembl_gene_dir(\@test_gene_dirs) }, follow => 1}, $test_dir);

my $feat_sets_info = {
	'baseline' => {
		'feat_set' => $baseline_features,
		'gene_dirs' => \@baseline_gene_dirs,
		'run_dirs' => \@baseline_run_dirs
	},
	'test' => {
		'feat_set' => $test_features,
		'gene_dirs' => \@test_gene_dirs,
		'run_dirs' => \@test_run_dirs
	}
};

while (my ($feat_set_key, $feat_set_info) = each(%{$feat_sets_info}) ) {

	my $feat_set = $feat_set_info->{'feat_set'};
	my $gene_dirs = $feat_set_info->{'gene_dirs'};
	my $run_dirs = $feat_set_info->{'run_dirs'};
	my $origin_dir = getcwd();
	foreach my $gene_dir (@{$gene_dirs}) {
		chdir $gene_dir;
		if ( -f 'annot.gtf' && -f 'pannot.gtf' && -f 'transc.fa' && -f 'transl.fa' ) {
			print "Running singĺe-gene APPRIS in: ${gene_dir}\n";
			my @cmds = ('appris_bin_g', '-s', 'Homo sapiens', '-x', "$feat_set", '-m', 'fm1sctra',
			            '-l', 'info');
			_run_system_call(\@cmds, 'log');
			push(@{$run_dirs}, $gene_dir);
		}
		chdir $origin_dir;
	}
}

my %cmp_info;
while (my ($feat_set_key, $feat_set_info) = each(%{$feat_sets_info}) ) {
	
	foreach my $run_dir (@{$feat_set_info->{'run_dirs'}}) {
		
		my $appris_file = $run_dir.'/appris';
		if ( ! -f $appris_file ) {
			die("file not found: '${appris_file}'");
		}
		
		open my $fh, $appris_file
			or die("failed to open file: '${appris_file}'");
		while (<$fh>) {
				my $line = $_;
				chomp $line;
				my @fields = split("\t", $line);
				if ( scalar(@fields) == 20 ) {  # i.e. full APPRIS output record
					my $gene_id = $fields[0];
					my $transc_id = $fields[2];
					$cmp_info{$gene_id}{$feat_set_key}{$transc_id} = {
						'gene_name' => $fields[1],
						'pi_label' => $fields[19],
						'ccds_id' => $fields[7]
					};
				}
		}
		close $fh;
	}
}

my %diff_info;
while (my($gene_id, $gene_info) = each %cmp_info) {

	my @transc_ids = keys %{$cmp_info{$gene_id}{'baseline'}};
	my %p1_set_info = (
		'baseline' => Set::Scalar->new(),
		'test' => Set::Scalar->new()
	);
	
	foreach my $transc_id (@transc_ids) {
		foreach my $feat_set_key (('baseline', 'test')) {
			my $pi_label = $cmp_info{$gene_id}{$feat_set_key}{$transc_id}{'pi_label'};
			if ($pi_label eq 'PRINCIPAL:1') {
				$p1_set_info{$feat_set_key}->insert($transc_id);
			}
		}
	}
	
	if ( ! $p1_set_info{'baseline'}->is_equal($p1_set_info{'test'}) ) {
		foreach my $transc_id (@transc_ids) {
			my $pi_label_baseline = $cmp_info{$gene_id}{'baseline'}{$transc_id}{'pi_label'};
			my $pi_label_test = $cmp_info{$gene_id}{'test'}{$transc_id}{'pi_label'};
			if ( grep(/^PRINCIPAL:1$/, ($pi_label_baseline, $pi_label_test)) ) {
				my %transc_info = (
					'gene_id' => $gene_id,
					'gene_name' => $cmp_info{$gene_id}{'baseline'}{$transc_id}{'gene_name'},
					'transc_id' => $transc_id,
					'pi_label_baseline' => $pi_label_baseline,
					'pi_label_test' => $pi_label_test,
					'ccds_id' => $cmp_info{$gene_id}{'baseline'}{$transc_id}{'ccds_id'}
				);
				$diff_info{$gene_id}{$transc_id} = \%transc_info;
			}
		}
	}
}

if ( scalar(keys %diff_info) > 0 ) {
	
	open(my $fh, '>', $out_file)
		or die("failed to open file '${out_file}'");
	my @field_names = ('gene_id', 'gene_name', 'transc_id', 'pi_label_baseline',
	                   'pi_label_test', 'ccds_id');
	print $fh join("\t", @field_names)."\n";
	foreach my $gene_info (values %diff_info) {
		foreach my $transc_info (values %{$gene_info}) {
			my @row = (
				$transc_info->{'gene_id'},
				$transc_info->{'gene_name'},
				$transc_info->{'transc_id'},
				$transc_info->{'pi_label_baseline'},
				$transc_info->{'pi_label_test'},
				$transc_info->{'ccds_id'}
			); 
			my $line = join("\t", @row);
			print $fh $line."\n";
		}
	}
	close $fh;
}
