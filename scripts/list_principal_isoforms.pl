#/usr/bin/env perl -w

use strict;

use Cwd;
use File::Find;
use Getopt::Long;

my ($anno_dir) = undef;
my ($out_file) = undef;

&GetOptions(
	'anno-dir=s' => \$anno_dir,
  'out-file=s' => \$out_file
);

sub match_ensembl_gene_dir ($) {
	my ($gene_dirs) = @_;
	push(@{$gene_dirs}, $File::Find::fullname) if ( -d && /^ENS(?:[A-Z]{3})?G\d+$/ );
}

my @gene_dirs;
if ( ! -d $anno_dir ) {
	die("annotation directory not found: ${anno_dir}");
}
find({wanted => sub { match_ensembl_gene_dir(\@gene_dirs) }, follow => 1}, $anno_dir);

my @pi_recs;
my $origin_dir = getcwd();
foreach my $gene_dir (@gene_dirs) {

	chdir $gene_dir;
  if ( -f 'annot.gtf' && -f 'pannot.gtf' && -f 'transc.fa' && -f 'transl.fa' ) {
		open my $fh, 'appris'
			or die("failed to open APPRIS file");
		while (<$fh>) {
				my $line = $_;
				chomp $line;
				my @fields = split("\t", $line);
				if ( scalar(@fields) == 20 ) {  # i.e. full APPRIS output record
					my $gene_id = $fields[0];
					my $transc_id = $fields[2];
					my $pi_label = $fields[19];
					if ( $pi_label =~ /PRINCIPAL:\d+/ ) {
						my @pi_rec = ($gene_id, $transc_id, $pi_label);
						push(@pi_recs, \@pi_rec);
					}
				}
		}
		close $fh;
	}
	chdir $origin_dir;
}

open(OUTFILE, '>', $out_file) or die("failed to open file '${out_file}'");
foreach my $pi_rec (@pi_recs) {
	my $line = join("\t", @{$pi_rec})."\n";
	print OUTFILE $line;
}
