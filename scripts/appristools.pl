#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use JSON;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile );
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$NUM_MAX_PROC
);

# Input parameters
my ($steps) = undef;
my ($conf_file) = undef;
my ($conf_data) = undef;
my ($methods) = undef;
my ($formats) = undef;
my ($num_process) = undef;
my ($email) = undef;
my ($loglevel) = undef;

&GetOptions(
	'steps|p=s'			=> \$steps,
	'conf|c=s'			=> \$conf_file,
	'data|d=s'			=> \$conf_data,	
	'methods|m=s'		=> \$methods,
	'formats|f=s'		=> \$formats,
	'nproc|n=s'			=> \$num_process,
	'email|e=s'			=> \$email,
	'loglevel|l=s'		=> \$loglevel,
);

# Check required parameters
unless ( defined $steps and ( defined $conf_file or defined $conf_data) ) {
	print `perldoc $0`;
	print "\nBad option combination of input data\n\n";
	exit 1;
}

# Get method's pipeline
unless ( defined $methods ) {
	$methods = $ENV{APPRIS_METHODS};
}
else {
	my ($met) = '';
	if ( $methods =~ /f/  ) { $met .= 'firestar,'   }
	if ( $methods =~ /m/  ) { $met .= 'matador3d,'  }
	if ( $methods =~ /m2/ ) { $met .= 'matador3d2,' }
	if ( $methods =~ /s/  ) { $met .= 'spade,'      }
	if ( $methods =~ /c/  ) { $met .= 'corsair,'    }
	if ( $methods =~ /t/  ) { $met .= 'thump,'      }
	if ( $methods =~ /r/  ) { $met .= 'crash,'      }
	if ( $methods =~ /p/  ) { $met .= 'proteo,'     }
	if ( $methods =~ /a/  ) { $met .= 'appris,'     }
	$met =~ s/\,$//g;
	$methods = $met;
}

# Extract Server config file
my ($server_json) = JSON->new();
my ($CONFIG) = $server_json->decode( getStringFromFile($conf_file) );
my (@CFG_SPECIES) = sort { $CONFIG->{'species'}->{$a}->{'order'} <=> $CONFIG->{'species'}->{$b}->{'order'} } keys(%{$CONFIG->{'species'}});

# get num. process (by default 1)
$NUM_MAX_PROC = ( defined $num_process ) ? $num_process : 1;
my ($NUM_PROC_CHILDS) = 0;
my ($PROC_CHILDS) = undef;



#################
# Method bodies #
#################
sub run_pipeline($);
sub params_run_pipe($);
sub param_check_files($);
sub param_retrieve_data($);
sub param_db_create($);
sub param_db_insert($);
sub param_db_backup($);
sub param_retrieve_method($);

# Main subroutine
sub main()
{
	# run pipeline for each gene datasets
	info("-- run pipeline for each gene datasets...");
	foreach my $species_id ( @CFG_SPECIES ) {
		my ($cfg_species) = $CONFIG->{'species'}->{$species_id};
		foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
			for ( my $i = 0; $i < scalar(@{$cfg_assembly->{'datasets'}}); $i++ ) {
				my ($cfg_dataset) = $cfg_assembly->{'datasets'}->[$i];
				if ( exists $cfg_dataset->{'pipeline'} and exists $cfg_dataset->{'pipeline'}->{'envfile'} ) {
					my ($conf_env_file) = ( $cfg_dataset->{'pipeline'}->{'envfile'} =~ /^\// ) ? $cfg_dataset->{'pipeline'}->{'envfile'} : $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/'.$cfg_dataset->{'pipeline'}->{'envfile'};
					
					my ($pid) = fork();
					
					# submit job
					if ( !defined $pid ) {
					    die "Cannot fork: $!";
					}
					if ( $pid ) {
						$NUM_PROC_CHILDS++;
					}	
					elsif ( $pid == 0 ) {
					    # client process
						info("\n** script($$): run_pipeline:$conf_env_file\n");
#					    close STDOUT; close STDERR; # so parent can go on
						eval {
							run_pipeline($conf_env_file);
						};
						exit 1 if($@);
					    exit 0;
					}
					
					# wait until at least one process finish
					while ( $NUM_PROC_CHILDS >= $NUM_MAX_PROC ) {
						my $pid = wait();
						$NUM_PROC_CHILDS--;						
						info("\n** script($pid) exits\n\n");
					}
					
				}
			}			
		}		
	}	

}

# run pipeline
sub run_pipeline($)
{
	my ($conf_data) = @_;
	
	# Step 1: execute APPRIS
	if ( $steps =~ /1/ )
	{		
		info("executing pipeline...");
		eval {
			my ($params) = params_run_pipe($conf_data);
			my ($cmd) = "apprisall $params";
			info($cmd);
			system ($cmd);
		};
		throw("executing pipeline") if($@);

#		info("checking results...");
#		eval {
#			my ($params) = param_check_files($conf_data);
#			my ($cmd) = "perl $ENV{APPRIS_SCRIPTS_DIR}/check_files.pl $params ";
#			info($cmd);
#			my (@output) = `$cmd`;
#print STDOUT "OUTPUT: \n".Dumper(@output)."\n";
#		};
#		throw("checking results") if($@);
		
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
			my ($params) = param_retrieve_data($conf_data);
			my ($cmd) = "appris_retrieve_main_data $params";
			info($cmd);
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
	if ( $steps =~ /3/ )
	{
		# create database
		info("creating database...");
		eval {
			my ($params) = param_db_create($conf_data);
			my ($cmd) = "appris_db_create $params";
			info($cmd);
			system ($cmd);
		};
		throw("creating database") if($@);
		
		info("inserting data into db...");
		eval {
			my ($params) = param_db_insert($conf_data);
			my ($cmd) = "appris_insert_appris $params";
			info($cmd);
			system ($cmd);
		};
		throw("inserting data") if($@);
		
# TODO: Check if data has been inserted correctly (compare with the no. of gene dataset using g_annotation.sql file)
# If is ERROR => Stop
# Otherwise => Keep going

		info("creating db backup...");
		eval {
			my ($params) = param_db_backup($conf_data);
			my ($cmd) = "appris_db_backup $params";
			info($cmd);
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
			my ($params) = param_retrieve_method($conf_data) . " -f $format ";			
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
					info($cmd);
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
sub params_run_pipe($)
{
	my ($conf_data) = @_;
	my ($params) = '';
	
	if ( defined $conf_data ) { $params .= " -c $conf_data " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $email ) 	{ $params .= " -e $email " }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_check_files($)
{
	my ($conf_data) = @_;
	my ($params) = '';
	
	if ( defined $conf_data ) { $params .= " -c $conf_data " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_retrieve_data($)
{
	my ($conf_data) = @_;
	my ($params) = '';
	
	if ( defined $conf_data ) { $params .= " -c $conf_data " }
	else { throw("configuration is not provided") }
		
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_db_create($)
{
	my ($conf_data) = @_;
	my ($params) = '';
	
	# get vars from config file
	my ($APPRIS_SPECIES) = `. $conf_data; echo \$APPRIS_SPECIES`; $APPRIS_SPECIES =~ s/\s*$//g;
	my ($APPRIS_SCRIPTS_DB_INI) = `. $conf_data; echo \$APPRIS_SCRIPTS_DB_INI`; $APPRIS_SCRIPTS_DB_INI =~ s/\s*$//g;
	my ($APPRIS_DATA_DIR) = `. $conf_data; echo \$APPRIS_DATA_DIR`; $APPRIS_DATA_DIR =~ s/\s*//g;
	
	my ($cfg) = new Config::IniFiles( -file => $APPRIS_SCRIPTS_DB_INI );
	my ($spe) = $APPRIS_SPECIES; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	$params .= " -d ".$cfg->val($specie_db, 'db');
	$params .= " -h ".$cfg->val('APPRIS_DATABASES', 'host');
	$params .= " -u ".$cfg->val('APPRIS_DATABASES', 'user');
	$params .= " -p ".$cfg->val('APPRIS_DATABASES', 'pass');

	return $params;		
}
sub param_db_insert($)
{
	my ($conf_data) = @_;
	my ($params) = '';
	
	if ( defined $conf_data ) { $params .= " -c $conf_data " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
	return $params;	
}
sub param_db_backup($)
{
	my ($conf_data) = @_;
	my ($params) = '';
	
	# get vars from config file
	my ($APPRIS_SPECIES) = `. $conf_data; echo \$APPRIS_SPECIES`; $APPRIS_SPECIES =~ s/\s*$//g;
	my ($APPRIS_SCRIPTS_DB_INI) = `. $conf_data; echo \$APPRIS_SCRIPTS_DB_INI`; $APPRIS_SCRIPTS_DB_INI =~ s/\s*$//g;
	my ($APPRIS_DATA_DIR) = `. $conf_data; echo \$APPRIS_DATA_DIR`; $APPRIS_DATA_DIR =~ s/\s*//g;
	
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
sub param_retrieve_method($)
{
	my ($conf_data) = @_;
	my ($params) = '';
		
	if ( defined $conf_data ) { $params .= " -c $conf_data " }
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

appristools

=head1 DESCRIPTION

Executes all APPRIS 'steps  

=head2 Required arguments (data input):

  -p, --steps {string} <Process steps>
	* 1 - Executes APPRIS pipeline -\n
	* 2 - Retrieves the list of principal isoforms and Compare stats between versions -\n
	* 3 - Inserts the annotations into database -\n
	* 4 - Retrieves the data files of methods -\n
		
  -c, --conf     {file} <Config file for all gene datatasets (JSON format)>
  	or  
  -d, --data     {file} <Config file for one gene datasets   (ENV file)>  

=head2 Optional arguments:
		
  -m, --methods {vector} <List of methods. Characters of methods: fmsctrpa (default: all)>
	* f  - Functionally important residues, firestar
	* m  - Protein structural information, Matador3D
	* s  - Presence of whole protein domains, SPADE
	* c  - Conservation against vertebrates, CORSAIR
	* t  - Presence of whole trans-membrane helices, THUMP
	* r  - Prediction of signal peptide and sub-cellular location, CRASH
	* p  - Proteomic evidence, PROTEO
	* a  - Principal Isoforms, APPRIS
  
  -f, --format {vector}   <Output format: 'gtf,bed,bed12', or 'json' (default: JSON)>
  
  -n, --nproc  {integer}  <Num. processes>
  
  -e, --email  {email}    <E-mail address>

=head2 Optional arguments (log arguments):

  -l, --loglevel {vector} <Define log level: info,debug,none (default: NONE)>


=head1 EXAMPLE

	appristools -p 1234 -c ws/config.json -m fmsctrpa -e 'jmrodriguez@cnio.es' -f gtf
	
	appristools -p 1234 -s conf/scripts/apprisrc.Hsap -m fmsctrpa -e 'jmrodriguez@cnio.es' -f gtf

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
