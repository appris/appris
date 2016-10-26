#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use JSON;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Exception qw( info throw );

# Input parameters
my ($release) = undef;
my ($relnotes_file) = undef;
my ($server_file) = undef;
my ($loglevel) = undef;

&GetOptions(
	'release|r=s'		=> \$release,
	'notes|n=s'			=> \$relnotes_file,	
	'server|s=s'		=> \$server_file,
	'loglevel|l=s'		=> \$loglevel,
);

# Check required parameters
unless ( defined $release and defined $relnotes_file and defined $server_file ) {
	print `perldoc $0`;
	print "\nBad option combination of input data\n\n";
	exit 1;
}

###################
# Global variable #
###################
use vars qw(
	$SRV_NAME
	$SRV_PUB_RELEASE_DIR
	$SRV_DB_HOST
	$SRV_DB_USER
	$SRV_DB_PWD
	$SRV_RELEASE_DIR
	$SRV_RELDAT_DIR

	$APPRIS_DATA_DIR
	
	$LOC_WSDIR
	$LOC_RELEASE_DIR
	$LOC_RELDAT_DIR
		
);

$SRV_NAME				= 'appris@appris';
$SRV_PUB_RELEASE_DIR	= '/local/appris/pub/releases';
$SRV_DB_HOST			= 'localhost';
$SRV_DB_USER			= 'appris';
$SRV_DB_PWD				= 'appris.appris';
$SRV_RELEASE_DIR		= $SRV_PUB_RELEASE_DIR.'/'.$release;
#$SRV_RELDAT_DIR			= $SRV_RELEASE_DIR.'/datafiles';


#$APPRIS_DATA_DIR 		= $ENV{APPRIS_DATA_DIR};
$APPRIS_DATA_DIR		=  '/home/jmrodriguez/projects/APPRIS/data';
$LOC_WSDIR				= '/tmp';
$LOC_RELEASE_DIR 		= $LOC_WSDIR.'/'.$release;
#$LOC_RELDAT_DIR 		= $LOC_RELEASE_DIR.'/datafiles';


# Extract Server config file
my ($server_json) = JSON->new();
my ($SERVER) = $server_json->decode( getStringFromFile($server_file) );

#################
# Method bodies #
#################


# Main subroutine
sub main()
{
	# create local workspace
	info("create local workspace...");
	eval {
		my ($cmd) = "rm -rf $LOC_RELEASE_DIR && mkdir -p $LOC_RELEASE_DIR";
		info($cmd);
		system ($cmd);
	};
	throw("creating local workspace") if($@);
	
	# add into local workspace with datafiles, for each species
	info("add into local workspace with datafiles, for each species...");
	while (my ($species_id, $cfg_species) = each($SERVER->{'species'}) ) {

		foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
			#foreach my $cfg_dataset (@{$cfg_assembly->{'datasets'}}) {
			for ( my $i = 0; $i < scalar(@{$cfg_assembly->{'datasets'}}); $i++ ) {
				my ($cfg_dataset) = $cfg_assembly->{'datasets'}->[$i];
				if ( exists $cfg_dataset->{'database'} ) {
					my ($as_name) = $cfg_assembly->{'name'};
					my ($ds_id) = $cfg_dataset->{'id'};
					my ($ds_dir) = $APPRIS_DATA_DIR.'/'.$species_id.'/'.$ds_id;
					my ($relspe_dir) = $LOC_RELEASE_DIR.'/datafiles/'.$species_id;
					my ($reldat_dir) = $relspe_dir.'/'.$ds_id;
					
					eval {
						my ($cmd) = "mkdir -p $relspe_dir && cp -rp $ds_dir $relspe_dir/.";
						info($cmd);
						system ($cmd);
					};
					throw("creating local workspace") if($@);
					
					if ( exists $cfg_dataset->{'type'} and ($cfg_dataset->{'type'} eq 'current') and $i == 0 ) { # create assembly link to first current dataset
						eval {
							my ($cmd) = "cd $relspe_dir ; ln -s $ds_id $as_name";
							info($cmd);
							system ($cmd);
						};
						throw("creating assembly link") if($@);						
					}
				}
			}			
		}		
	}	
		
	# create release note for given release
	info("create release note for given release...");
	my ($release_notes) = getTotalStringFromFile($relnotes_file);
	my ($find) = 0;
	my ($rel_notes) = "";
	foreach my $nline (@{$release_notes}) {		
		if ( $nline =~ /^==APPRIS_RELEASE:[^\:]*:$release/ ) {
			$find = 1;
		}
		elsif ( $nline =~ /^==APPRIS_RELEASE:.*/ ) {
			$find = 0;
		}
		if ( $find == 1 ) {
			$rel_notes .= $nline;
		}
	}
	if ( $rel_notes ne '' ) {
		my ($relnotes_file) = $LOC_RELEASE_DIR.'/relnotes.txt';
		my ($plog) = printStringIntoFile($rel_notes, $relnotes_file);
		throw("Printing _get_ccds_stats_content\n") unless(defined $plog);		
	}

	# upload datafiles to server
	info("upload datafiles to server...");
	eval {
		my ($cmd) = "cd $LOC_RELEASE_DIR && tar -cf - . | ssh $SRV_NAME 'mkdir -p $SRV_RELEASE_DIR; cd $SRV_RELEASE_DIR/; tar -xf -'";
		info($cmd);
		system ($cmd);
	};
	throw("updating datafiles to server") if($@);
	
	# import databases into server
	info("import databases into server...");
	my ($cmd_imp) = "";
	while (my ($species_id, $cfg_species) = each($SERVER->{'species'}) ) {

		foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
			foreach my $cfg_dataset (@{$cfg_assembly->{'datasets'}}) {				
				if ( exists $cfg_dataset->{'database'} ) {
					my ($ds_id) = $cfg_dataset->{'id'};
					my ($ds_db) = $cfg_dataset->{'database'};					
					my ($srv_relspe_dir) = $SRV_RELEASE_DIR.'/datafiles/'.$species_id;
					my ($srv_reldat_dir) = $srv_relspe_dir.'/'.$ds_id;
					my ($srv_db_file) = $srv_reldat_dir.'/appris_db.dump.gz';
					$cmd_imp .= "appris_db_import -d $ds_db -h $SRV_DB_HOST -u root -i $srv_db_file && ";
				}
			}			
		}		
	}
	if ( $cmd_imp ne '' ) {
		eval {
			$cmd_imp =~ s/\s*\&\&\s*$//g;
			my ($cmd) = "ssh $SRV_NAME '$cmd_imp'";
			info($cmd);
			system ($cmd);
		};
		throw("importing databases in the server") if($@);		
	}
		
}

####################
# SubMethod bodies #
####################


main();


1;

__END__

=pod

=head1 NAME

appristools_srv

=head1 DESCRIPTION

Executes all APPRIS 'steps  

=head2 Required arguments (data input):

  -r, --release   {string} <Release identifier>
  
  -n, --notes     {file}   <Release Notes file - TXT format - >
		
  -s, --server   {file}    <Config file of Server - JSON format - >
  
=head2 Optional arguments (log arguments):

  -l, --loglevel {vector} <Define log level: info,debug,none (default: NONE)>


=head1 EXAMPLE

	appristools_srv \
		-r 2016_06.v17 \
		-n changelog.md \
		-s ws/server.json

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
