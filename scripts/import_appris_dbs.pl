#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use JSON;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Exception qw( info throw warning );

# Input parameters
my ($conf_file) = undef;
my ($indir) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf|c=s'			=> \$conf_file,
	'indir|d=s'			=> \$indir,
	'loglevel|l=s'		=> \$loglevel,
);

# Check required parameters
unless ( defined $conf_file and defined $indir ) {
	print `perldoc $0`;
	print "\nBad option combination of input data\n\n";
	exit 1;
}

###################
# Global variable #
###################
use vars qw(
	$SRV_DB_HOST
	$SRV_DATA_DIR
);
$SRV_DB_HOST = 'localhost';
$SRV_DATA_DIR = $indir;

# Extract Server config file
my ($config_json) = JSON->new();
my ($CONFIG) = $config_json->decode( getStringFromFile($conf_file) );
my (@CFG_SPECIES) = sort { $CONFIG->{'species'}->{$a}->{'order'} <=> $CONFIG->{'species'}->{$b}->{'order'} } keys(%{$CONFIG->{'species'}});
my ($CONFIG_VERSION) = $CONFIG->{'version'};


#################
# Method bodies #
#################


# Main subroutine
sub main()
{
	# create command
	info("-- create command...");	
	my ($cmd_imp) = "";
	foreach my $species_id ( @CFG_SPECIES ) {
		my ($cfg_species) = $CONFIG->{'species'}->{$species_id};
		
		foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
			foreach my $cfg_dataset (@{$cfg_assembly->{'datasets'}}) {				
				if ( exists $cfg_dataset->{'database'} and exists $cfg_dataset->{'database'}->{'name'} ) {
					my ($ds_id) = $cfg_dataset->{'id'};
					my ($ds_db) = $cfg_dataset->{'database'}->{'name'}.'_'.$ds_id;					
					my ($srv_relspe_dir) = $SRV_DATA_DIR.'/'.$species_id;
					my ($srv_reldat_dir) = $srv_relspe_dir.'/'.$ds_id;
					my ($srv_db_file) = $srv_reldat_dir.'/appris_db.dump.gz';
					$cmd_imp .= "appris_db_import -d $ds_db -h $SRV_DB_HOST -u root -i $srv_db_file && ";
				}
			}			
		}		
	}
	
	
	# import databases into server
	info("-- import databases into localhost...");
	if ( $cmd_imp ne '' ) {
		eval {
			$cmd_imp =~ s/\s*\&\&\s*$//g;
			my ($cmd) = "$cmd_imp";
			info($cmd);
			system($cmd);
		};
		throw("importing databases in localhost") if($@);		
	}
	
}


main();


1;

__END__

=pod

=head1 NAME

import_appris_db

=head1 DESCRIPTION

Import the databases of APPRIS annotations  

=head2 Arguments (data input):

  -c, --conf      {file}   <Config file for all gene datatasets (JSON format)>
  
  -d, --indir     {dir}    <Directory with the release directory with data from APPRIS>
  
=head1 EXAMPLE

	appristools_srv \
		-c ws/config.json \
		-d data/2017_08.v24

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
