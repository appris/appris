#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use List::Util qw(min max sum);
use APPRIS::Utils::File qw( getStringFromFile getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;
use Data::Dumper;

###################
# Global variable #
###################
use vars qw(
);

# Input parameters
my ($ccds_g_file) = undef;
my ($ccds_g_file2) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'ccdsfile1|1=s'		=> \$ccds_g_file,
	'ccdsfile2|2=s'		=> \$ccds_g_file2,
	'outfile|o=s'		=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless( defined $ccds_g_file and defined $ccds_g_file2 and defined $outfile )
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
	$logger->info("-- open files\n");
	my ($ccds_report) = ccds_report($ccds_g_file);
	my ($ccds_report2) = ccds_report($ccds_g_file2);
	
	$logger->info("-- cmp list of genes\n");
	my ($output) = cmp_ccds($ccds_report, $ccds_report2);
	if ($output ne '') {
		my ($printing_file_log) = printStringIntoFile($output,$outfile);
		$logger->error("printing") unless(defined $printing_file_log);		
	}

	$logger->finish_log();
	
	exit 0;
	
}

sub ccds_report($)
{
	my ($file) = @_;
	my ($report);
	my ($cont) = getTotalStringFromFile($file);
	if ( scalar(@{$cont}) > 1 ) {
		my (@head_cols) = split("\t", $cont->[0]);
		for (my $i=1; $i < scalar(@{$cont}); $i++) { # head (i=0) with the method name
			my ($line) = $cont->[$i];
			my (@line_cols) = split("\t", $line);
			if ( scalar(@line_cols) > 1 ) {
				my ($gene_id) = $line_cols[0];
				
				for (my $j=1; $j < scalar(@line_cols); $j++) { # col (j=0) is the gene_id 
					my ($col) = $line_cols[$j];
					if ( $col == 1 ) {
						my ($method) = $head_cols[$j]; $method =~ s/\s*//g;
						$report->{$method}->{$gene_id} = 1;
					}
				}				
			}
		}		
	}
	return $report;
}

sub cmp_ccds($$)
{
	my ($report1,$report2) = @_;
	my ($output) = '';
	while (my ($met1,$rep1) = each(%{$report1}) ) {
		$output .= "-- Method: $met1\n";
		if ( exists $report2->{$met1} ) {
			my ($rep2) = $report2->{$met1};
			while (my ($gene1,$disa) = each(%{$rep1}) ) {
				unless ( exists $rep2->{$gene1} ) { $output .= "< $gene1\n" }
			}
			$output .= "--\n";
			while (my ($gene2,$disa) = each(%{$rep2}) ) {
				unless ( exists $rep1->{$gene2} ) { $output .= "> $gene2\n" }
			}					
		}
		else {
			foreach my $gene1 (keys(%{$rep1})) { $output .= "< $gene1\n" }
		}
		$output .= "\n";
	}
	return $output;
}

main();

1;


__END__

=head1 NAME

cmp_stats_ccds.pl

=head1 DESCRIPTION

=head2 Required arguments:

	-1,--ccdsfile1= <List1 of genes with CCDS disagreement>
	
	-2,--ccdsfile2= <List2 of genes with CCDS disagreement>
	
	-o,--outfile= <Comparition Output file>

=head2 Optional arguments (logs):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl cmp_stats_ccds.pl
			
		-1 data/homo_sapiens/ens81.v9.17Jul2015/appris_stats.ccds_g.txt
		
		-2 data/homo_sapiens/ens80.v8.16Apr2015/appris_stats.ccds_g.txt
				
		--outfile=data/homo_sapiens/ens81.v9.17Jul2015/appris_cmp_ccds.e81VSe80.txt
		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut