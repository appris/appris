#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile getStringFromFile printStringIntoFile );
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################

# Input parameters
my ($i_dom_file) = undef;
my ($i_hmm_file) = undef;
my ($i_dat_file) = undef;
my ($o_hmm_file) = undef;
my ($o_dat_file) = undef;
my ($loglevel) = undef;

&GetOptions(
	'i_dom|d=s'			=> \$i_dom_file,
	'i_hmm|i1=s'		=> \$i_hmm_file,
	'i_dat|i2=s'		=> \$i_dat_file,
	'o_hmm|o1=s'		=> \$o_hmm_file,
	'o_dat|o2=s'		=> \$o_dat_file,
	'loglevel|l=s'		=> \$loglevel,
);

# Check required parameters
unless ( defined $i_dom_file and defined $i_hmm_file and defined $i_dat_file and defined $o_hmm_file and defined $o_dat_file ) {
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
	my ($o_hmm_txt, $o_dat_txt) = ('','');
	
	info("-- create report with the list of domains...");
	my ($i_dom_rep) = APPRIS::Utils::File::getTotalStringFromFile($i_dom_file);
	my (%i_dom) = map { $_ =~ s/\.[0-9]*\n*//g; $_ => 1 } @{$i_dom_rep};
	
	info("-- extract domains from hmm db...");
	my ($i_hmm_txt) = APPRIS::Utils::File::getStringFromFile($i_hmm_file);	
	my (@i_hmm_rep) = split(/\/\/\n*/, $i_hmm_txt);
	foreach my $rep (@i_hmm_rep) {
		if ( $rep =~ /ACC\s+([^\n]*)\n+/ ) {
			my ($acc) = $1;
			$acc =~ s/\.[0-9]*//g;
			if ( exists $i_dom{$acc} ) {
				$o_hmm_txt .= $rep.'//'."\n";
			}
		}
	}

	info("-- extract domains from hmm db...");
	my ($i_dat_txt) = APPRIS::Utils::File::getStringFromFile($i_dat_file);	
	my (@i_dat_rep) = split(/\/\/\n*/, $i_dat_txt);
	foreach my $rep (@i_dat_rep) {
		if ( $rep =~ /AC\s+([^\n]*)\n+/ ) {
			my ($acc) = $1;
			$acc =~ s/\.[0-9]*//g;
			if ( exists $i_dom{$acc} ) {
				$o_dat_txt .= $rep.'//'."\n";
			}
		}
	}
	
	info("-- create new database files...");
	APPRIS::Utils::File::printStringIntoFile($o_hmm_txt, $o_hmm_file);
	APPRIS::Utils::File::printStringIntoFile($o_dat_txt, $o_dat_file);
	
	
#	print STDERR "REP:\n".Dumper(@i_hmm_rep);
	
#	open(FILE,$i_hmm_file) or return undef;
#	while ( <FILE> ) {
#		
#	}
#	my(@array)=<FILE>;
#	close(FILE);
	
	
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

  -d,  --i_dom {file} <List of Pfam domains>
		
  -i1, --i_hmm {file} <Pfam database  in HMM format>
  
  -i2, --i_dat {file} <Pfam data file in HMM format>
  
  -o1, --o_hmm {file} <Output file with the new Pfam database  in HMM format>
  
  -o2, --o_dat {file} <Output file with the new Pfam data file in HMM format>

=head2 Optional arguments (log arguments):

  -l, --loglevel {vector} <Define log level: info,debug,none (default: NONE)>


=head1 EXAMPLE

	perl selectGivenDomains.pl -d scripts/spade/AllPfamDomainData.20170529.ID.txt -i1 pfam_201706/raw/Pfam-A.hmm -i2 pfam_201706/raw/Pfam-A.hmm.dat -o1 pfam_201706/Pfam-A.hmm -o2 pfam_201706/Pfam-A.hmm.dat
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
