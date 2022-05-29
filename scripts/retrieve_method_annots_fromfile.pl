#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Registry;
use APPRIS::Utils::File qw( prepare_workspace printStringIntoFile getStringFromFile );
use APPRIS::Parser qw(
	parse_firestar_rst
	parse_matador3d_rst
	parse_matador3d2_rst
	parse_spade_rst
	parse_corsair_rst
	parse_corsair_alt_rst
	parse_thump_rst
	parse_crash_rst
	parse_appris_rst
);
use APPRIS::Utils::Logger;


###################
# Global variable #
###################

# Input parameters
my ($inpath) = undef;
my ($methods) = undef;
my ($type) = undef;
my ($outpath) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'inpath=s'			=> \$inpath,
	'methods=s'			=> \$methods,
	'outpath=s'			=> \$outpath,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless( defined $inpath and defined $methods and defined $outpath )
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


#####################
# Method prototypes #
#####################
sub get_spade_annot($$);
sub get_corsair_alt_annot($$);


#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	# Output vars
	my ($tmp_outfiles);

	# Prepare input parameters
	$inpath=~s/\/$//;
	$outpath=~s/\/$//;

	# Get the temporal results of methods
	foreach my $method ( split(',',$methods) ) {
		if ( $method eq 'spade' or $method eq 'corsair' or $method eq 'corsair_alt' ) {
			# save tmp outfile
			my ($tmp_outfile) = "$outpath/tmp.appris_method.$method.txt";
			push(@{$tmp_outfiles}, [$method, $tmp_outfile] );
			# remove tmp file
			eval {
				my ($cmd) =	"rm $tmp_outfile";
				$logger->info("\n** script: $cmd\n");
				my (@cmd_out) = `$cmd`;
			};
			$logger->error("removing temporal file") if($@);
			# run
			eval {
				my ($cmd) =	"find $inpath -type f -name '$method' -exec cat {} \\; >> $tmp_outfile";
				$logger->info("\n** script: $cmd\n");
				my (@cmd_out) = `$cmd`;
			};
			$logger->error("finding the $method annotations") if($@);
		}
		else {
			$logger->error("$method method does not apply");
		}
	}

	# Read the temporal results of methods
	foreach my $tmp_outfile (@{$tmp_outfiles}) {
		# local vars
		my ($method,$tmp_file) = ($tmp_outfile->[0],$tmp_outfile->[1]);
		my ($rst) = '';
		my ($analysis);
		my ($output_cont) = '';

		# save the file content
		$logger->info("-- read $tmp_file\n");
		if ( -e $tmp_file and (-s $tmp_file > 0) ) {
			($rst) = getStringFromFile($tmp_file);
			$logger->error("can not open $tmp_file: $!\n") unless ( defined $rst );
		}
		else {
			$logger->error("tmp_file does not exit: $tmp_file\n");
		}
		# parse results
		$logger->info("-- parse $method report\n");
		if ( $method eq 'spade' ) {
			($analysis) = parse_spade_rst($rst);
		}
		elsif ( $method eq 'corsair' or $method eq 'corsair_alt' ) {
			($analysis) = parse_corsair_rst($rst);
		}
		

		# extract method report
		# print method output
		$logger->info("-- extract $method report\n");
		while ( my ($transcript_id,$report) = each(%{$analysis}) ) {
			if ( $method eq 'spade' ) {
				$output_cont .= get_spade_annot($transcript_id, $report);
			}
			elsif ( $method eq 'corsair' or $method eq 'corsair_alt' ) {
				$output_cont .= get_corsair_alt_annot($transcript_id, $report);
			}

			# print output
			if ( $output_cont ne '' ) {
				my ($output_file) = "$outpath/appris_method.$method.txt";
				my ($printing_file_log) = printStringIntoFile($output_cont, $output_file);
				$logger->error("Printing output") unless ( defined $printing_file_log );
			}
		}


	}

	$logger->info("retrieve_method_annots_fromfile:finished ---------------\n");
		
	$logger->finish_log();
	
	exit 0;
}


sub get_spade_annot($$)
{
	my ($transcript_id,$method) = @_;
	my ($report) = '';
	if ( ref $method eq ref {} and exists $method->{'domains'} ) {
		foreach my $region (@{$method->{'domains'}}) {
				my ($hmm_name) = $region->{'hmm_name'} || '';
				my ($pstart) = $region->{'alignment_start'} || '';
				my ($pend) = $region->{'alignment_end'} || '';
				my ($type_domain) = $region->{'type_domain'} || '';
				my ($bit_score) = $region->{'bit_score'} || '';
				$report .= $transcript_id."\t".
						$hmm_name."\t".
						$pstart."\t".
						$pend."\t".
						$bit_score."\n";
		}
	}
	return $report;
} # end get_spade_annot


sub get_corsair_alt_annot($$)
{
	my ($transcript_id,$method) = @_;
	my ($report) = '';
	if ( ref $method eq ref {} ) {
		my ($score) = $method->{'score'} || '0';
		$report .= $transcript_id."\t".
				$score."\n";
	}
	return $report;
} # end get_corsair_alt_annot


main();

1;


__END__

=head1 NAME

retrieve_method_annots_fromfile

=head1 DESCRIPTION

Exports the annotations for each method from the annotation files:

	* TODO: firestar: list of ligands/compounds.
	
	* TODO: Matador3D: list of pdb's, and coordinate alignments (at the moment, it's reported cds coordinates).
	
	* SPADE: list of domains, and alignment coordinates.
	
	* CORSAIR: list of transcripts with the score.

	* CORSAI_ALT: list of transcripts with the score.
	

=head1 SYNOPSIS

retrieve_method_annots_fromfile

=head2 Required arguments:

	--inpath= <Acquire input files from PATH>
	
	--method= <Name of Method(s)> (firestar,matador3d,spade,corsair,corsair_alt)	
	
	--outpath= <Output directory where the method results wil be saving>
	
=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_method_annots_fromfile.pl
		
		--inpath=annotations/2021_08.v45/homo_sapiens/g19v45

		--method=firestar,matador3d,spade,inertia
				
		--outpath=data/2021_08.v45/
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut