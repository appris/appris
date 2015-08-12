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
my ($conf_db) = undef;
my ($methods) = undef;
my ($email) = undef;
my ($loglevel) = undef;

&GetOptions(
	'steps|p=s'			=> \$steps,
	'conf-sp|cs=s'		=> \$conf_species,	
	'conf-db|cd=s'		=> \$conf_db,
	'methods|m=s'		=> \$methods,
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
my ($species);
if ( getStringFromFile($conf_species) =~ /APPRIS_SPECIES="([^\"]*)"/) {
	$species = $1;
}


#################
# Method bodies #
#################
sub params_run_pipe();
sub param_check_files();
sub param_retrieve_data();
sub param_insert_db();

# Main subroutine
sub main()
{
	# Step 1: execute APPRIS
	if ( $steps =~ /1/ )
	{		
		info("executing pipeline...");
		my ($params_run_pipe) = params_run_pipe();
		eval {
			my ($cmd) = "appris $params_run_pipe";
			system ($cmd);
		};
		throw("executing pipeline") if($@);

		info("checking results...");
		my ($params_check) = param_check_files();
		eval {
			my ($cmd) = "perl $ENV{APPRIS_SCRIPTS_DIR}/check_files.pl $params_check ";
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
		my ($params_data_main) = param_retrieve_data();
		eval {
			my ($cmd) = "appris_retrieve_main_data $params_data_main";
			system ($cmd);
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
	if ( $steps =~ /3/ ) {
		
		# create database
		info("creating database...");
		my ($params_create_db) = param_create_db();
		eval {
			my ($cmd) = "appris_db_create $params_create_db";
			system ($cmd);
		};
		throw("creating database") if($@);
		
#		info("inserting APPRIS data into db...");
#		my ($ins_params) = param_insert_db();
#		eval {
#			my ($cmd) = "appris_insert_appris $ins_params";
#			system ($cmd);
#		};
#		throw("inserting data") if($@);
		
# TODO: Check if data has been inserted correctly (compare with the no. of gene dataset using g_annotation.sql file)
# If is ERROR => Stop
# Otherwise => Keep going

# TODO: Make a backup

# TODO: Send email with final decision of this step

	}


# TODO: Step 4 can be done in parallel

#	# Step 4: retrieve data files for methods in GTF format and BED format (Download section of Web and TrackHUB)
#	my ($dat2_params) = create_data_params();
#	if ( $steps =~ /4/ ) {	
#		$dat2_params .= " -f gtf ";
#		eval {
#			my ($cmd) = "appris_retrieve_method_data $dat2_params";
#			info("** script: $cmd\n");
#			system ($cmd);
#		};
#		throw("deleting log files of appris") if($@);
#
# TODO: Check if something is bad.
# If is ERROR => Stop
# Otherwise => Keep going
#
# TODO: Send email with final decision of this step
#
#	}


# TODO: Step 5 Upload all data files to the server


}
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
sub param_create_db()
{
	my ($params) = '';
	
	my ($cfg) = new Config::IniFiles( -file => $conf_db );
	my ($spe) = $species; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	$params .= " -d ".$cfg->val($specie_db, 'db');
	$params .= " -h ".$cfg->val('APPRIS_DATABASES', 'host');
	$params .= " -u ".$cfg->val('APPRIS_DATABASES', 'user');
	$params .= " -p ".$cfg->val('APPRIS_DATABASES', 'pass');

	return $params;		
}
sub param_insert_db()
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
		
  -c, --conf {string|file} <Species abbreviation | Config file for species>  
	* Hsap  - Homo sapiens -\n
	* Mmus  - Mus musculus -\n
	* Rnor  - Rattus norvegicus -\n
	* Drer  - Danio rerio -\n
	* Sscr  - Sus scrofa -\n
	* Ptro  - Pan troglodytes -\n
	* Dmel  - Drosophila melanogaster -\n
	* Cele  - Caenorhabditis elegans -\n
	* Lpar  - Lynx pardinus -\n  
 Or 
	<Config file name>

=head2 Optional arguments:
		
  -m, --methods {vector} <List of retrieved methods separated by commas: 'firestar,matador3d,spade,corsair,thump,crash,appris,proteo' (default: all)>
	* f  - Functionally important residues, firestar
	* m  - Protein structural information, Matador3D
	* s  - Presence of whole protein domains, SPADE
	* c  - Conservation against vertebrates, CORSAIR
	* t  - Presence of whole trans-membrane helices, THUMP
	* cr - Prediction of signal peptide and sub-cellular location, CRASH
	* p  - Proteomic evidence, PROTEO
  
  -f, --format {vector} <Output format: 'gtf', 'bed', or 'json' (default: JSON)>
  
  -e, --email {email} <E-mail address>

=head2 Optional arguments (log arguments):

  -l, --loglevel {vector} <Define log level: info,debug,none (default: NONE)>


=head1 EXAMPLE

	apprisall -p 1234 -s Hsap -m 'firestar,matador3d,spade,corsair,thump,crash,appris,proteo' -e 'jmrodriguez@cnio.es' -f gtf

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
