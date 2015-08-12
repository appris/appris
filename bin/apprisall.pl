#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use Data::Dumper;
use APPRIS::Utils::Exception qw( info verbose throw );

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
my ($outpath) = undef;
my ($email) = undef;
my ($loglevel) = undef;

&GetOptions(
	'steps|p=s'			=> \$steps,
	'conf|c=s'			=> \$conf_species,	
	'methods|m=s'		=> \$methods,
	'outpath|o=s'		=> \$outpath,
	'email|e=s'			=> \$email,
	'loglevel|l=s'		=> \$loglevel
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

#################
# Method bodies #
#################
sub params_appris();
sub create_data_params();
sub param_check_files();
sub param_data_main();

# Main subroutine
sub main()
{
	# Step 1: execute APPRIS
	if ( $steps =~ /1/ )
	{
#		info("executing APPRIS pipeline...");
#		my ($params_appris) = params_appris();
#		eval {
#			my ($cmd) = "appris $params_appris";
#			info("** script: $cmd\n");
#			system ($cmd);
#		};
#		throw("executing APPRIS pipeline") if($@);

		info("checking APPRIS results...");
		my ($params_check) = param_check_files();
		eval {
			my ($cmd) = "perl $ENV{APPRIS_SCRIPTS_DIR}/check_files.pl $params_check ";
			my (@output) = `$cmd`;
print STDOUT "OUTPUT: \n".Dumper(@output)."\n";
		};
		throw("checking APPRIS results") if($@);
	}
	
	# Step 2: retrieve the main data and stats comparison
	if ( $steps =~ /12/ )
	{
		info("retrieving APPRIS main data...");
		my ($params_data_main) = param_data_main();
		eval {
			my ($cmd) = "appris_retrieve_main_data $params_data_main";
			system ($cmd);
		};
		throw("retrieving the main data") if($@);
	}
	
#	# Step 3: insert annotations into APPRIS database
#	if ( $steps =~ /123/ ) {
#		my ($ins_params) = create_insert_params();
#		eval {
#			my ($cmd) = "appris_insert_appris $ins_params";
#			info("** script: $cmd\n");
#			system ($cmd);
#		};
#		throw("deleting log files of appris") if($@);
#	}
#
#	# Step 4: retrieve method data
#	my ($dat2_params) = create_data_params();
#	if ( $steps =~ /1234/ ) {	
#		$dat2_params .= " -f gtf ";
#		eval {
#			my ($cmd) = "appris_retrieve_method_data $dat2_params";
#			info("** script: $cmd\n");
#			system ($cmd);
#		};
#		throw("deleting log files of appris") if($@);
#	}

}
sub params_appris()
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
sub param_data_main()
{
	my ($params) = '';
	
	if ( defined $conf_species ) { $params .= " -c $conf_species " }
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
