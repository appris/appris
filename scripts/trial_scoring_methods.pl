#/usr/bin/env perl -w

use strict;

use Capture::Tiny 'capture_merged';
use Cwd qw(abs_path getcwd);
use File::Find;
use Getopt::Long;


my ($anno_dir) = undef;
my ($exp_features) = undef;
my ($out_file) = undef;

&GetOptions(
	'anno-dir=s' => \$anno_dir,
	'exp:s' => \$exp_features,
	'out-file=s' => \$out_file
);

my $curr_dir = getcwd();
$anno_dir = abs_path($curr_dir.'/'.$anno_dir);
$out_file = abs_path($curr_dir.'/'.$out_file);

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

my @gene_dirs;
my @run_dirs;
if ( ! -d $anno_dir ) {
	die("annotation directory not found: ${anno_dir}");
}
find({wanted => sub { match_ensembl_gene_dir(\@gene_dirs) }, follow => 1}, $anno_dir);


foreach my $gene_dir (@gene_dirs) {
	chdir $gene_dir;
	if ( -f 'annot.gtf' && -f 'pannot.gtf' && -f 'transc.fa' && -f 'transl.fa' && -f 'appris' ) {
		print "Running singĺe-gene APPRIS in: ${gene_dir}\n";
		my @cmds = ('appris_bin_g', '-s', 'Homo sapiens', '-x', "$exp_features", '-m', 'fm1scra',
		            '-l', 'info');
		_run_system_call(\@cmds, 'log');
		push(@run_dirs, $gene_dir);
	}
	chdir $curr_dir;
}

my @pi_recs;
foreach my $run_dir (@run_dirs) {

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
			my @pi_rec = (
				$fields[0],  # gene_id
				$fields[2],  # transc_id
				$fields[19],  # pi_label
				$fields[7]  # ccds_id
			);
			push(@pi_recs, \@pi_rec);
		}
	}
	close $fh;
}

open(OUTFILE, '>', $out_file) or die("failed to open file '${out_file}' for writing");
print OUTFILE join("\t", ('gene_id', 'transc_id', 'pi_label', 'ccds_id'))."\n";
foreach my $pi_rec (@pi_recs) {
	my $line = join("\t", @{$pi_rec})."\n";
	print OUTFILE $line;
}
