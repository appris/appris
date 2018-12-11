#!/usr/bin/perl -w
# _________________________________________________________________
# $Id$
# $Revision$
# Developed by:
#		Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es-
# _________________________________________________________________

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use JSON;
use LWP::UserAgent;

use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile prepare_workspace );

###################
# Global variable #
###################
		
# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($output_dir) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'conf=s'			=> \$config_file,
	'output=s'			=> \$output_dir,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,	
);

# Required arguments
unless ( defined $config_file and defined $output_dir )
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
$logger->init_log($str_params);

#####################
# Method prototypes #
#####################
sub parse_blast($$$);
sub check_alignment($$$$\$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# prepare tmp dir
	prepare_workspace($output_dir);

    # Scan all species
    my ($cfg) = JSON->new()->decode( getStringFromFile($config_file) );
	while ( my ($species_id, $species_rep) = each (%{$cfg}) ) {
		foreach my $spe_rep (@{$species_rep}) {
			my ($spe_name) = $spe_rep->{'name'};
			$spe_name =~ s/ /_/g;
			$logger->info("APPRISDATA:$spe_name\n".Dumper($spe_rep)."\n");
			my ($query) = '';
			eval
			{
				$logger->info("Running blast\n");
				my ($cmd) = "wget \"ftp://ftp.ncbi.nih.gov/genomes/$spe_name/protein/protein.fa.gz\" -O $output_dir/$spe_name.fasta.gz";
				$logger->debug("$cmd\n");
				system($cmd);
			};
			$logger->error("Downloading the species $spe_name: $!\n") if($@);
		}
	}


	
	$logger->finish_log();
	
	exit 0;	
}



main();

__END__

=head1 NAME

down_uniprot

=head1 DESCRIPTION

Download UniProt databases from the list of given species

=head1 VERSION

0.1

=head2 Required arguments:

	-c,--conf <Config file>

	-o,--output <Output directory>

=head2 Optional arguments:
	
	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
	

=head1 EXAMPLE

perl download_uniprot.pl

	--conf=/appris/conf/code/corsair_alt.diverge_time.human.json
	
	--output=/appris/db/uniprot_201810

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut