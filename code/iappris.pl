#!/usr/bin/perl -w

use strict;
use warnings;
use threads;

use Getopt::Long;
use FindBin;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Parser qw(
	parse_infiles
	parse_firestar
	parse_matador3d
	parse_spade
	parse_corsair
	parse_thump
	parse_crash
	parse_proteo
	parse_appris
);
use APPRIS::Intruder qw(
	feed_by_gencode
	feed_by_analysis
);
use APPRIS::Utils::File qw( getStringFromFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_APPRIS_DB_FILE
	$APPRIS_METHODS
	
	$APPRIS_PROTEO_DATA_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_APPRIS_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/apprisdb.ini';
foreach my $method ( split(',',$ENV{APPRIS_METHODS}) ) {
	$APPRIS_METHODS->{$method} = 1;
}
# TEMPORAL SOLUTION!!!!
$APPRIS_PROTEO_DATA_FILE = $ENV{APPRIS_PROTEO_DATA_FILE};
# TEMPORAL SOLUTION!!!!

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($id) = undef;
my ($inpath) = undef;
my ($species) = undef;
my ($methods) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'id=s'				=> \$id,
	'inpath=s'			=> \$inpath,
	'species=s'			=> \$species,
	'methods=s'			=> \$methods,
	'apprisdb-conf=s'	=> \$apprisdb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined  $id and
	defined  $inpath and
	defined  $species
){
	print `perldoc $0`;
	exit 1;
}

# Get method's pipeline
#unless ( defined $methods ) {
#	$methods = join( ',',keys(%{$APPRIS_METHODS}) );
#}

# Optional arguments
# get vars of appris db
unless ( defined $apprisdb_conf_file ) {
	$apprisdb_conf_file = $CONFIG_INI_APPRIS_DB_FILE;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

#################
# Method bodies #
#################
sub create_entity($$);
sub create_reports($$$);
sub parse_reports($$);
sub feed_entity($$);
sub feed_reports($$);

# Main subroutine
sub main()
{	
	# Create gene report
	$logger->info("-- create input entity\n");
	my ($entity) = create_entity($id, $inpath);
	unless ( defined $entity ) {
		$logger->error("creating input entity");
	}

	# Create reports for pipeline
	$logger->info("-- create input reports for appris 'intruder'\n");
	my ($reports) = create_reports($entity, $methods, $inpath);
	unless ( defined $reports ) {
		$logger->error("creating input reports for appris 'intruder'");
	}
	
	# Parse reports for pipeline
	$logger->info("-- parse input reports for appris 'intruder'\n");
	my ($analyses) = parse_reports($entity, $reports);
	unless ( defined $analyses ) {
		$logger->error("parsing input reports for appris 'intruder'");
	}
		
	# Intruder entity into db
	$logger->info("-- feed input entity\n");
	my ($f_entity) = feed_entity($apprisdb_conf_file, $entity);
	unless ( defined $f_entity ) {
		$logger->error("feeding input entity");
	}

	# Intruder reports into db
	$logger->info("-- feed input reports for appris 'intruder'\n");
	my ($feeds) = feed_reports($apprisdb_conf_file, $analyses);
	unless ( defined $feeds ) {
		$logger->error("feeding input reports for appris 'intruder'");
	}	
	
	$logger->finish_log();
	
	exit 0;	
}

sub create_entity($$)
{
	my ($id, $inpath) = @_;
	my ($entity);	
	my ($inputs) = {
		'data'		=> $inpath.'/'.$id.'.annot.gtf',
		'pdata'		=> $inpath.'/'.$id.'.pannot.gtf',
		'transc'	=> $inpath.'/'.$id.'.transc.fa',
		'transl'	=> $inpath.'/'.$id.'.transl.fa',
	};
	if ( -e $inputs->{'data'} and (-s $inputs->{'data'} > 0) and
			-e $inputs->{'pdata'} and (-s $inputs->{'pdata'} > 0) and
			-e $inputs->{'transc'} and (-s $inputs->{'transc'} > 0) and
			-e $inputs->{'transl'} and (-s $inputs->{'transl'} > 0) 
	) {
		# get inputs
		my ($data_file) = $inputs->{'data'};
		my ($pdata_file) = $inputs->{'pdata'};
		my ($transc_file) = $inputs->{'transc'};
		my ($transl_file) = $inputs->{'transl'};
		
		# get entity (important, conserve parentesis embedding variable)
		my ($entities) = parse_infiles($data_file, $transc_file, $transl_file);
		if ( scalar($entities) > 0 ) {
			$entity = $entities->[0];
		} else {
			$logger->error("parsing gencode-ensembl entity: $!\n");
		}
	}
	
	return $entity;
	
} # end create_entity

sub create_reports($$$)
{
	my ($entity, $methods, $inpath) = @_;
	my ($inputs);
	my ($reports);
	
	# acquire the inputs depending on given methods
	foreach my $met ( split(',',$methods) ) {
		$met = lc($met);
		if ( $APPRIS_METHODS->{$met} ) {
			if ( $met eq 'ensembl' ) {
				$inputs->{$met} = $entity;
			}
			elsif ( $met eq 'appris' ) {
				$inputs->{$met} = $inpath.'/'.$id.'.'.$met;
				$inputs->{$met.'_label'} = $inpath.'/'.$id.'.'.$met.'.label';
			}
			else {
				$inputs->{$met} = $inpath.'/'.$id.'.'.$met;				
			}
		}
		else {
			return $reports;
		}
	}
			
	# create reports
	while ( my ($i,$file) = each (%{$inputs}) ) {
		# works different ensembl-gencode data
		if ( $i eq 'ensembl' ) {
			$reports->{$i} = $file; # save entity
		}
		else {
			# save the file content
			if ( -e $file and (-s $file > 0) ) {			
				my ($rst) = getStringFromFile($file);
				$logger->error("can not open $file: $!\n") unless ( defined $rst );
				$reports->{$i} = $rst;		
			}
			else {
				$logger->debug("file does not exit: $file\n");
			}
		}		
	}
	
	return $reports;
	
} # end create_reports

sub parse_reports($$)
{
	my ($entity, $reports) = @_;
	my ($analyses);
	my ($method) = '';
	
	$method = 'ensembl';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		$analyses->{$method} = $entity;
	}
	$method = 'firestar';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_firestar($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'matador3d';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_matador3d($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'spade';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_spade($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'corsair';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_corsair($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'thump';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_thump($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'crash';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_crash($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'proteo';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};		
		my ($analysis) = parse_proteo($entity, $report);
		$analyses->{$method} = $analysis;
	}
	$method = 'appris';
	if ( exists $reports->{$method} ) {
		$logger->info("-- parse $method report\n");
		my ($report) = $reports->{$method};
		my ($report_label) = $reports->{$method.'_label'};
		my ($analysis) = parse_appris($entity, $report, $report_label);
		$analyses->{$method} = $analysis;
	}
	
	return $analyses;
	
} # end parse_reports

sub feed_entity($$)
{
	my ($conf_file, $analysis) = @_;
	my ($feed);
	
	# insert entity from given specie
	my ($cfg) = new Config::IniFiles( -file => $conf_file );
	my ($spe) = $species; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	my ($param) = {
			'-dbhost'       => $cfg->val('APPRIS_DATABASES', 'host'),
			'-dbuser'       => $cfg->val('APPRIS_DATABASES', 'user'),
			'-dbpass'       => $cfg->val('APPRIS_DATABASES', 'pass'),
			'-dbport'       => $cfg->val('APPRIS_DATABASES', 'port'),
			'-dbname'       => $cfg->val($specie_db, 'db'),
	};
	$logger->debug(Dumper($param)."\n");
	my ($intruder) = APPRIS::Intruder->new();
	$intruder->load_registry_from_db(%{$param});
	$intruder->feed_by_gencode($analysis);
	$feed = 1;

	return $feed;

} # end feed_entity

sub feed_reports($$)
{
	my ($conf_file, $analyses) = @_;
	my ($feed);
	
	# insert entity from given specie
	my ($cfg) = new Config::IniFiles( -file => $conf_file );
	my ($spe) = $species; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	my ($param) = {
			'-dbhost'       => $cfg->val('APPRIS_DATABASES', 'host'),
			'-dbuser'       => $cfg->val('APPRIS_DATABASES', 'user'),
			'-dbpass'       => $cfg->val('APPRIS_DATABASES', 'pass'),
			'-dbport'       => $cfg->val('APPRIS_DATABASES', 'port'),
			'-dbname'       => $cfg->val($specie_db, 'db'),
	};
	$logger->debug(Dumper($param)."\n");
	my ($intruder) = APPRIS::Intruder->new();
	$intruder->load_registry_from_db(%{$param});
	while ( my ($method, $analysis) = each(%{$analyses}) ) {		
		$logger->info("-- feed $method report\n");
		unless ( $method eq 'ensembl' ) {
			$intruder->feed_by_analysis($analysis, $method);
			$feed = 1;			
		}
	}
	
	return $feed;

} # end feed_reports

main();


1;

__END__

=head1 NAME

iappris

=head1 DESCRIPTION

script that insert results of APPRIS into database 

=head1 SYNOPSIS

run_appris

=head2 Input arguments:

	--id= <Ensembl gene identifier>
	
	--species= <Name of species -mammals->
	
	--inpath= <Acquire input files from PATH>
	
=head2 Optional arguments:

		--methods= <List of APPRIS's methods ('ensembl,firestar,matador3d,spade,corsair,thump,crash,appris'. Default: ALL)>
		
		--ensembldb-conf <Config file of Ensembl database (default: 'conf/ensembldb.ini' file)>
		
		--apprisdb-conf <Config file of APPRIS database (default: 'conf/apprisdb.ini' file)>
				
=head2 Optional arguments (log arguments):
	
		--loglevel=LEVEL <define log level (default: NONE)>	
	
		--logfile=FILE <Log to FILE (default: *STDOUT)>
		
		--logpath=PATH <Write logfile to PATH (default: .)>
		
		--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl iappris.pl

	--id=ENSMUSG00000017167
	
	--species='Mus musculus'
	
	--methods=firestar,matador3d
	
	--inpath=../features/ENSMUSG00000017167_e69/

=head1 EXAMPLE

perl iappris.pl

	--id=ENSG00000093072.11
	
	--species='Homo sapiens'
	
	--inpath=/home/jmrodriguez/projects/Encode/gencode15/annotations/chr22/ENSG00000093072.11
	
	--loglevel=INFO
	
	--logappend
	
	--logpath=/home/jmrodriguez/projects/Encode/gencode15/logs/
			
	--logfile=iappris.log


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
