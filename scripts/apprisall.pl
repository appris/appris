#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile );
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;

# Input parameters
my ($steps) = undef;
my ($conf_species) = undef;
my ($methods) = undef;
my ($formats) = undef;
my ($email) = undef;
my ($loglevel) = undef;

&GetOptions(
	'steps|p=s'			=> \$steps,
	'conf|c=s'			=> \$conf_species,	
	'methods|m=s'		=> \$methods,
	'formats|f=s'		=> \$formats,
	'email|e=s'			=> \$email,
	'loglevel|l=s'		=> \$loglevel,
);

# Check required parameters
unless ( defined $steps and defined $conf_species ) {
	print `perldoc $0`;
	print "\nBad option combination of input data\n\n";
	exit 1;
}

# Get method's pipeline
unless ( defined $methods ) {
	$methods = $ENV{APPRIS_METHODS};
}

# get species name from config file
my ($APPRIS_SPECIES) = `. $conf_species; echo \$APPRIS_SPECIES`; $APPRIS_SPECIES =~ s/\s*$//g;
my ($APPRIS_SCRIPTS_DB_INI) = `. $conf_species; echo \$APPRIS_SCRIPTS_DB_INI`; $APPRIS_SCRIPTS_DB_INI =~ s/\s*$//g;
my ($APPRIS_DATA_DIR) = `. $conf_species; echo \$APPRIS_DATA_DIR`; $APPRIS_DATA_DIR =~ s/\s*//g;
#my ($APPRIS_WSERVER_DOWNLOAD_DATA_DIR) = $ENV{APPRIS_WSERVER_DOWNLOAD_DATA_DIR}; $APPRIS_WSERVER_DOWNLOAD_DATA_DIR =~ s/\s*//g;
#my ($APPRIS_WS_DATE) = `. $conf_species; echo \$APPRIS_WS_DATE`; $APPRIS_WS_DATE =~ s/\s*//g;


#################
# Method bodies #
#################
sub params_run_pipe();
sub param_check_files();
sub param_retrieve_data();
sub param_db_create();
sub param_db_insert();
sub param_db_backup();
sub param_retrieve_method();

# Main subroutine
sub main()
{
	# Step 1: execute APPRIS
	if ( $steps =~ /1/ )
	{		
		info("executing pipeline...");
		eval {
			my ($params) = params_run_pipe();
			my ($cmd) = "appris $params";
			system ($cmd);
		};
		throw("executing pipeline") if($@);

		info("checking results...");
		eval {
			my ($params) = param_check_files();
			my ($cmd) = "perl $ENV{APPRIS_SCRIPTS_DIR}/check_files.pl $params ";
			my (@output) = `$cmd`;
print STDOUT "OUTPUT: \n".Dumper(@output)."\n";
		};
		throw("checking results") if($@);
		
# TODO: Check if script retrives output (list of wrong genes).
# If is ERROR => Stop
# Otherwise => Keep going

# TODO: Send email with final decision of this step

	}
	
	# Step 2: retrieve the main data and stats comparison
	if ( $steps =~ /2/ )
	{
		info("retrieving main data...");
		eval {
			my ($params) = param_retrieve_data();
			my ($cmd) = "appris_retrieve_main_data $params";
		};
		throw("retrieving the main data") if($@);
		
# TODO: Compare numbers with older release: No. Principal Isoforms, CCDS comparison, etc.
# Eg.
# 		cut -f 2 appris_data.principal.txt | sort -u | wc -l

# TODO: Check if something is bad.
# If is ERROR => Stop
# Otherwise => Keep going

# TODO: Send email with final decision of this step
	}
	
	# Step 3: insert annotations into APPRIS database
	if ( $steps =~ /3/ )
	{
		# create database
		info("creating database...");
		eval {
			my ($params) = param_db_create();
			my ($cmd) = "appris_db_create $params";
			system ($cmd);
		};
		throw("creating database") if($@);
		
		info("inserting data into db...");
		eval {
			my ($params) = param_db_insert();
			my ($cmd) = "appris_insert_appris $params";
			system ($cmd);
		};
		throw("inserting data") if($@);
		
# TODO: Check if data has been inserted correctly (compare with the no. of gene dataset using g_annotation.sql file)
# If is ERROR => Stop
# Otherwise => Keep going

		info("creating db backup...");
		eval {
			my ($params) = param_db_backup();
			my ($cmd) = "appris_db_backup $params";
			system ($cmd);
		};
		throw("creating db backup") if($@);
		
# TODO: Send email with final decision of this step

	}


	# Step 4: retrieve data files for methods in GTF format and BED format (Download section of Web and TrackHUB)
	if ( $steps =~ /4/ )
	{
		my ($forks) = 0;
		foreach my $format ( split(',', $formats) ) {
			my ($params) = param_retrieve_method() . " -f $format ";			
			my ($pid) = fork();
			if (not defined $pid) {
				warn "Could not fork in $format";
				next;
			}
			if ($pid) {
				$forks++;
			} else {
			    #close STDOUT; close STDERR; # so parent can go on
				eval {
					my ($cmd) = "appris_retrieve_method_data $params";
					system ($cmd);
				};
				exit 1 if($@);
			    exit 0;
			}
		}
		for (1 .. $forks) {
			my ($pid) = wait();
		}
		

# TODO: Check if something is bad.
# If is ERROR => Stop
# Otherwise => Keep going
#
# TODO: Send email with final decision of this step

	}

## TODO: Step 5 Upload all data files to the server
#
#	# Step 5: Upload data files to the server
#	if ( $steps =~ /5/ )
#	{		
#
## TODO: Check if something is bad.
## If is ERROR => Stop
## Otherwise => Keep going
##
## TODO: Send email with final decision of this step
#
#	}
	
	
}

####################
# SubMethod bodies #
####################
sub params_run_pipe()
{
	my ($params) = '';
	
	if ( defined $conf_species ) { $params .= " -c $conf_species " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $email ) 	{ $params .= " -e $email " }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_check_files()
{
	my ($params) = '';
	
	if ( defined $conf_species ) { $params .= " -c $conf_species " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_retrieve_data()
{
	my ($params) = '';
	
	if ( defined $conf_species ) { $params .= " -c $conf_species " }
	else { throw("configuration is not provided") }
		
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_db_create()
{
	my ($params) = '';
	
	my ($cfg) = new Config::IniFiles( -file => $APPRIS_SCRIPTS_DB_INI );
	my ($spe) = $APPRIS_SPECIES; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	$params .= " -d ".$cfg->val($specie_db, 'db');
	$params .= " -h ".$cfg->val('APPRIS_DATABASES', 'host');
	$params .= " -u ".$cfg->val('APPRIS_DATABASES', 'user');
	$params .= " -p ".$cfg->val('APPRIS_DATABASES', 'pass');

	return $params;		
}
sub param_db_insert()
{
	my ($params) = '';
	
	if ( defined $conf_species ) { $params .= " -c $conf_species " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_db_backup()
{
	my ($params) = '';
	
	my ($cfg) = new Config::IniFiles( -file => $APPRIS_SCRIPTS_DB_INI );
	my ($spe) = $APPRIS_SPECIES; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	$params .= " -d ".$cfg->val($specie_db, 'db');
	$params .= " -h ".$cfg->val('APPRIS_DATABASES', 'host');
	$params .= " -u ".$cfg->val('APPRIS_DATABASES', 'user');
	$params .= " -p ".$cfg->val('APPRIS_DATABASES', 'pass');
	$params .= " -o ".$APPRIS_DATA_DIR.'/appris_db.dump.gz';

	return $params;		
}
sub param_retrieve_method()
{
	my ($params) = '';
		
	if ( defined $conf_species ) { $params .= " -c $conf_species " }
	else { throw("configuration is not provided") }
		
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}


main();


1;

__END__

=pod

=head1 NAME

apprisall

=head1 DESCRIPTION

Executes all APPRIS 'steps  

=head2 Required arguments (data input):

  -p, --steps {string} <Process steps>
	* 1 - Executes APPRIS pipeline -\n
	* 2 - Retrieves the list of principal isoforms and Compare stats between versions -\n
	* 3 - Inserts the annotations into database -\n
	* 4 - Retrieves the data files of methods -\n
		
  -c, --conf {file} <Config file for species>  

=head2 Optional arguments:
		
  -m, --methods {vector} <List of retrieved methods separated by commas: 'firestar,matador3d,spade,corsair,thump,crash,appris,proteo' (default: all)>
	* f  - Functionally important residues, firestar
	* m  - Protein structural information, Matador3D
	* s  - Presence of whole protein domains, SPADE
	* c  - Conservation against vertebrates, CORSAIR
	* t  - Presence of whole trans-membrane helices, THUMP
	* cr - Prediction of signal peptide and sub-cellular location, CRASH
	* p  - Proteomic evidence, PROTEO
  
  -f, --format {vector} <Output format: 'gtf,bed', or 'json' (default: JSON)>
  
  -e, --email {email} <E-mail address>

=head2 Optional arguments (log arguments):

  -l, --loglevel {vector} <Define log level: info,debug,none (default: NONE)>


=head1 EXAMPLE

	apprisall -p 1234 -c conf/scripts/apprisrc.Hsap -m 'firestar,matador3d,spade,corsair,thump,crash,appris,proteo' -e 'jmrodriguez@cnio.es' -f gtf

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
