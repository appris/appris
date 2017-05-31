#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile printStringIntoTmpFile );
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################

# Input parameters
my ($l_dom_file) = undef;
my ($i_hmm_file) = undef;
my ($i_dat_file) = undef;
my ($loglevel) = undef;

&GetOptions(
	'list|l=s'			=> \$l_dom_file,
	'i_hmm|i1=s'		=> \$i_hmm_file,
	'i_dat|i2=s'		=> \$i_dat_file,
	'loglevel|l=s'		=> \$loglevel,
);

# Check required parameters
unless ( defined $l_dom_file and defined $i_hmm_file and defined $i_dat_file ) {
	print `perldoc $0`;
	print "\nCheck required inputs!!\n\n";
	exit 1;
}


#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	info("-- run pipeline for each gene datasets...");
	
}


main();


1;

__END__

=pod

=head1 NAME

selectGivenDomains

=head1 DESCRIPTION

Create Pfam database in HMM format with a given list of domains  

=head2 Required arguments:

  -l,  --list  {file} <List of Pfam domains>
		
  -i1, --i_hmm {file} <Pfam database  in HMM format>
  
  -i2, --i_dat {file} <Pfam data file in HMM format>
  
=head2 Optional arguments (log arguments):

  -l, --loglevel {vector} <Define log level: info,debug,none (default: NONE)>


=head1 EXAMPLE

	perl selectGivenDomains.pl -l scripts/spade/AllPfamDomainData.20170529.ID.txt -i1 pfam_201706/Pfam-A.hmm -i2 pfam_201706/Pfam-A.hmm.dat
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
