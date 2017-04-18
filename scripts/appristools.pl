#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Config::IniFiles;
use FindBin;
use JSON;
use Data::Dumper;
use APPRIS::Utils::File qw( getStringFromFile printStringIntoTmpFile );
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
my ($conf_db) = undef;
my ($methods) = undef;
my ($formats) = undef;
my ($num_process) = undef;
my ($email) = undef;
my ($loglevel) = undef;

&GetOptions(
	'steps|p=s'			=> \$steps,
	'conf|c=s'			=> \$conf_file,
	'data|d=s'			=> \$conf_data,
	'conf-db|b=s'		=> \$conf_db,
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

# Get database variables
unless ( defined $conf_db ) {
	$conf_db = $ENV{APPRIS_CONF_DB_FILE};	
}

# Extract Server config file
my ($CONFIG);
my ($CONFIG_VERSION);
my (@CFG_SPECIES);	
if ( defined $conf_file ) {
	my ($server_json) = JSON->new();
	$CONFIG = $server_json->decode( getStringFromFile($conf_file) );
	@CFG_SPECIES = sort { $CONFIG->{'species'}->{$a}->{'order'} <=> $CONFIG->{'species'}->{$b}->{'order'} } keys(%{$CONFIG->{'species'}});
	$CONFIG_VERSION = $CONFIG->{'version'};
}

# get num. process (by default 1)
$NUM_MAX_PROC = ( defined $num_process ) ? $num_process : 1;
my ($NUM_PROC_CHILDS) = 0;
my ($PROC_CHILDS) = undef;



#################
# Method bodies #
#################
sub create_tmp_db_ini($);
sub run_pipeline($;$);
sub params_run_pipe($);
sub param_check_files($);
sub param_retrieve_data($);
sub param_db_insert($;$);
sub param_retrieve_method($);

# Main subroutine
sub main()
{
	if ( defined $conf_file and defined $CONFIG and defined $CONFIG_VERSION and scalar(@CFG_SPECIES) > 0 ) {
		# run pipeline for each gene datasets
		info("-- run pipeline for each gene datasets...");
		foreach my $species_id ( @CFG_SPECIES ) {
			my ($cfg_species) = $CONFIG->{'species'}->{$species_id};
			foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
				for ( my $i = 0; $i < scalar(@{$cfg_assembly->{'datasets'}}); $i++ ) {
					my ($cfg_dataset) = $cfg_assembly->{'datasets'}->[$i];
					if ( exists $cfg_dataset->{'id'} and
						 exists $cfg_dataset->{'pipeline'} and exists $cfg_dataset->{'pipeline'}->{'envfile'} and
						 exists $cfg_dataset->{'database'} and exists $cfg_dataset->{'database'}->{'name'} and
						 exists $cfg_dataset->{'database'} and exists $cfg_dataset->{'database'}->{'inifile'}
					) {
						# get env files
						#	get env vars for dataset
						my ($conf_env_file) = ( $cfg_dataset->{'pipeline'}->{'envfile'} =~ /^\// ) ? $cfg_dataset->{'pipeline'}->{'envfile'} : $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/'.$cfg_dataset->{'pipeline'}->{'envfile'};
						my ($cfg_dataset_id) = $cfg_dataset->{'id'};
						my ($cfg_dataset_name) = $cfg_dataset_id; $cfg_dataset_name =~ s/v[0-9]+$//g;
						# 	create tmp ini file for db
						my ($conf_db_file) = create_tmp_db_ini($cfg_dataset);
												
						my ($config_dataset) = {
							'id'       => $cfg_dataset_id,
							'name'     => $cfg_dataset_name,
							'species'  => $species_id,
							'ws_name'  => $species_id.'/'.$cfg_dataset_name,
							'ws_date'  => $CONFIG_VERSION.'/'.$species_id.'/'.$cfg_dataset_id,
							'env_file' => $conf_env_file,
							'db_file'  => $conf_db_file,
						};
												
						# submit job
						my ($pid) = fork();						
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
								run_pipeline($conf_env_file, $config_dataset);
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
	elsif ( defined $conf_data ) {
		# run pipeline for unique gene datasets
		info("-- run pipeline for unique gene dataset...");		
		run_pipeline($conf_data);
	}
}

sub create_tmp_db_ini($) {
	my ($cfg_dataset) = @_;
	my ($cfg_dataset_id) = $cfg_dataset->{'id'};
	my ($db_inifile) = ( $cfg_dataset->{'database'}->{'inifile'} =~ /^\// ) ? $cfg_dataset->{'database'}->{'inifile'} : $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/'.$cfg_dataset->{'database'}->{'inifile'};
	my ($db_name) = $cfg_dataset->{'database'}->{'name'}.'_'.$cfg_dataset_id;
	
	my ($fname) = $db_name.'.ini';
	my ($config_cont) = getStringFromFile($db_inifile);
	$config_cont =~ s/__APPRIS__DATABASES__DBNAME__/$db_name/g;
	my ($outfile) = printStringIntoTmpFile($fname, $config_cont);
	
	return $outfile;
}

# create the environment variables for data analysis
sub create_env_vars($)
{
	my ($conf_ds) = @_;
	my ($output) = '';
	if ( defined $conf_ds ) {
		$output .= "export APPRIS_WS_NAME='".$conf_ds->{'ws_name'}."' && " if ( exists $conf_ds->{'ws_name'} );
		$output .= "export APPRIS_WS_DATE='".$conf_ds->{'ws_date'}."' && " if ( exists $conf_ds->{'ws_date'} );
	}	
	return $output;
}

# run pipeline
sub run_pipeline($;$)
{
	my ($conf_data, $conf_ds) = @_;
	
	# Step 0: export the environment variables
	my ($exp_env) = ( defined $conf_ds ) ? create_env_vars($conf_ds) : ''; 
	
	# Step 1: execute APPRIS
	if ( $steps =~ /1/ )
	{		
		info("executing pipeline...");
		eval {
			my ($params) = params_run_pipe($conf_data);
			my ($cmd) = "$exp_env appris_run_appris $params";
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
			my ($cmd) = "$exp_env appris_retrieve_main_data $params";
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
		info("inserting data into db...");
		eval {
			my ($params) = param_db_insert($conf_data, $conf_ds);
			my ($cmd) = "$exp_env appris_insert_appris $params";
			info($cmd);
			system ($cmd);
		};
		throw("inserting data") if($@);
		
# TODO: Check if data has been inserted correctly (compare with the no. of gene dataset using g_annotation.sql file)
# If is ERROR => Stop
# Otherwise => Keep going

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
					my ($cmd) = "$exp_env appris_retrieve_method_data $params";
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
sub param_db_insert($;$)
{
	my ($conf_data, $conf_ds) = @_;
	my ($params) = '';
	
	if ( defined $conf_ds and exists $conf_ds->{'db_file'} ) { $params .= " -d ".$conf_ds->{'db_file'}." " }
	else { throw("configuration is not provided") }

	if ( defined $conf_data ) { $params .= " -c $conf_data " }
	else { throw("configuration is not provided") }
	
	if ( defined $methods ) { $params .= " -m $methods " }
	else { throw("configuration is not provided") }
	
	if ( defined $loglevel ) { $params .= " -l $loglevel " }
	
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

  -b, --conf-db  {file} <INI file with database configuration>  
	
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
	
	appristools -p 1234 -d conf/scripts/apprisrc.Hsap -m fmsctrpa -e 'jmrodriguez@cnio.es' -f gtf

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
