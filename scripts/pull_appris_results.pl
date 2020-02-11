#/usr/bin/env perl -w

use strict;
use warnings::register;

use Cwd qw(getcwd realpath);
use DateTime::Format::ISO8601;
use File::Find;
use File::Spec::Functions qw(catfile);
use Getopt::Long;

my ($anno_dir) = undef;
my ($after_dt_opt) = undef;
my ($out_file) = undef;

&GetOptions(
	'anno-dir=s' => \$anno_dir,
	'after:s' => \$after_dt_opt,
	'out-file=s' => \$out_file
);

print("Pulling APPRIS results from: $anno_dir\n");
if ( ! -d $anno_dir ) {
	die("annotation directory not found:' ${anno_dir}'");
}
$anno_dir = realpath($anno_dir);

my $after_dt;
if ( defined $after_dt_opt ) {
	$after_dt = DateTime::Format::ISO8601->parse_datetime($after_dt_opt);
	print("Pulling APPRIS results after: $after_dt\n");
}

$out_file = realpath($out_file);


sub _match_appris_run_dir ($) {
	my ($run_dirs) = @_;
	if ( -d && /^ENS(?:[A-Z]{3})?G\d+$/ ) {
		my $gene_dir = $File::Find::fullname;
		my $annot_file = catfile($gene_dir, 'annot.gtf');
		my $pannot_file = catfile($gene_dir, 'pannot.gtf');
		my $transc_file = catfile($gene_dir, 'transc.fa');
		my $transl_file = catfile($gene_dir, 'transl.fa');
		my $appris_file = catfile($gene_dir, 'appris');
		if ( -f $annot_file &&
             -f $pannot_file &&
			 -f $transc_file &&
			 -f $transl_file &&
             -f $appris_file ) {
			print("Found APPRIS run directory: '$gene_dir'\n");
			push(@{$run_dirs}, $gene_dir);
		}
	}
}

my @run_dirs;
if ( ! -d $anno_dir ) {
	die("annotation directory not found: ${anno_dir}");
}
find({
	preprocess => sub { return grep { -d } @_; },
	wanted => sub { _match_appris_run_dir(\@run_dirs) },
	follow => 1
}, $anno_dir);

my @pi_recs;
my $num_genes;
foreach my $run_dir (@run_dirs) {

	my $appris_file = catfile($run_dir, 'appris');
	
	if ( defined $after_dt ) {
		my $mod_dt = DateTime::Format::ISO8601->parse_datetime(
			POSIX::strftime("%Y-%m-%dT%H:00:00", localtime(( stat $appris_file )[9]))
		);
		if ( $mod_dt <= $after_dt ) {
			next;
		}
	}
	
	print("Pulling APPRIS result from: ${run_dir}\n");
	open(my $fh, $appris_file)
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
	close($fh);
	
	$num_genes++;
}

if (@pi_recs) {

	print("Pulled APPRIS results for ".$num_genes."/".scalar(@run_dirs)." genes.\n");
	open(my $fh, '>', $out_file) or die("failed to open file '${out_file}'");
	print $fh join("\t", ('gene_id', 'transc_id', 'pi_label', 'ccds_id'))."\n";
	foreach my $pi_rec (@pi_recs) {
		my $line = join("\t", @{$pi_rec})."\n";
		print $fh $line;
	}
	close($fh);

} else {
	print("No results found matching criteria.\n");
}


=head1 NAME

pull_appris_results

=head1 DESCRIPTION

Pull APPRIS results from annotation directory.

=head2 Required arguments:

  --anno-dir {path} <Path of annotations containing APPRIS results>

  --after {string} <Time (ISO-8601 format) after which results are taken>

  --out-file {path} <Path of APPRIS results summary file>

=head1 EXAMPLE

	perl scripts/pull_appris_results.pl --anno-dir=annotations/2019_12_hsapE98.v29/ \
	     --after 2020-02-10T12:00:00 --out-file=appris_results.tab

=cut

