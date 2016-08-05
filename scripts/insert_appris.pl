#!/usr/bin/perl -w

use strict;
use warnings;
use threads;

use Getopt::Long;
use FindBin;
use Config::IniFiles;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use appris qw( create_gene_list );
use APPRIS::Parser qw( parse_gencode );
use APPRIS::Utils::File qw( prepare_workspace printStringIntoFile getTotalStringFromFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_ENSEMBL_DB_FILE
	$CONFIG_INI_APPRIS_DB_FILE
	$LOGLEVEL
	$LOGAPPEND
);

$LOCAL_PWD					= $FindBin::Bin;
$CONFIG_INI_ENSEMBL_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/ensembldb.ini';
$CONFIG_INI_APPRIS_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/apprisdb.ini';
$LOGLEVEL					= 'INFO';
$LOGAPPEND					= '';

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($id) = undef;
my ($position) = undef;
my ($gene_list_file) = undef;
my ($data_file) = undef;
my ($translations_file) = undef;
my ($species) = undef;
my ($e_version) = undef;
my ($inpath) = undef;
my ($methods) = undef;
my ($type_of_input) = undef;
my ($ensembldb_conf_file) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'species=s'			=> \$species,
	'id=s'				=> \$id,
	'position=s'		=> \$position,
	'gene-list=s'		=> \$gene_list_file,
	'data=s'			=> \$data_file,
	'transl=s'			=> \$translations_file,
	'e-version=s'		=> \$e_version,
	'inpath=s'			=> \$inpath,
	'methods=s'			=> \$methods,
	'ensembldb-conf=s'	=> \$ensembldb_conf_file,	
	'apprisdb-conf=s'	=> \$apprisdb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Get the type of input (the order of conditions is important)
# GENCODE mode
if ( defined $species and defined $data_file and defined $inpath ) {
	$type_of_input = 'datafile';
	if		( defined $gene_list_file ) { $type_of_input = 'datafile-list'; }
	elsif	( defined $position ) { $type_of_input = 'datafile-position'; }
	$inpath .= '/';
}
# SEQUENCE mode
elsif ( defined $translations_file and defined $species and defined $inpath ) {
	$type_of_input = 'sequence';
	$inpath .= '/';	
}
# ENSEMBL mode
elsif ( defined $species and defined $e_version and defined $inpath ) {
	$type_of_input = 'ensembl';
	if		( defined $id ) { $type_of_input = 'ensembl-id'; }
	elsif	( defined $gene_list_file ) { $type_of_input = 'ensembl-list'; }
	elsif	( defined $position ) { $type_of_input = 'ensembl-position'; }
	$inpath .= '/';
}
else {
	print `perldoc $0`;
	print "\nCheck required inputs!!\n\n";
	exit 1;
}

# Optional arguments

# get vars of ensembl db
unless ( defined $ensembldb_conf_file ) {
	$ensembldb_conf_file = $CONFIG_INI_ENSEMBL_DB_FILE;
}
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
$LOGLEVEL	= $loglevel if ( defined $loglevel );
$LOGAPPEND	= "--logappend" if ( defined $logappend );
$logger->init_log();

#################
# Method bodies #
#################
sub insert_ensembl($$);
sub insert_gencode($$;$);
sub insert_pipeline($$;$);
sub insert_appris($$$);

# Main subroutine
sub main()
{
	# run appris pipeline for each gene depending on input
	$logger->info("-- from given input...");
	if ( $type_of_input =~ /datafile/ ) {
		$logger->info(" $type_of_input type\n");
		
		# create gene list
		$logger->info("-- create gene list\n");
		my ($gene_list, $gencode_data) = appris::create_gene_list(
															-type	=> $type_of_input,
															-id		=> $id,
															-list	=> $gene_list_file,
															-pos	=> $position,
															-gdata	=> $data_file,
		);
		
		$logger->info("-- run gencode data files\n");
		my ($runtimes) = insert_gencode($gencode_data, $inpath, $gene_list);
	}
	elsif ( $type_of_input =~ /ensembl/ ) {
		$logger->info(" $type_of_input type\n");	

		# create gene list
		$logger->info("-- create gene list\n");
		my ($gene_list) = appris::create_gene_list(
											-type	=> $type_of_input,		
											-id		=> $id,
											-list	=> $gene_list_file,
											-pos	=> $position,
											-econf	=> $ensembldb_conf_file,
											-ever	=> $e_version,
											-spe	=> $species,
		);
		
		$logger->info("-- insert ensembl ids\n");
		my ($runtimes) = insert_ensembl($gene_list, $inpath);
	}
	elsif ( $type_of_input =~ /sequence/ ) {
		$logger->info(" $type_of_input type\n");
		
		# create gene list
		$logger->info("-- create gene list\n");
		my ($gene_list,$gene_data) = appris::create_gene_list(
											-type		=> $type_of_input,
											-gtransl	=> $translations_file,
		);
				
		$logger->info("-- insert sequence\n");
		my ($runtimes) = insert_sequence($gene_data, $inpath);	
	}	
	else {
		$logger->error("analying input parameters");
	}	
	
	$logger->finish_log();
	
	exit 0;	
}

sub insert_ensembl($$)
{
	my ($data, $inpath) = @_;
	my ($runtimes) = undef;
			
	foreach my $gene_id (keys(%{$data})) {
		my ($runtime) = insert_pipeline($gene_id, $inpath);
		push(@{$runtimes}, $runtime);		
	}
	
	return $runtimes;
	
} # end insert_ensembl

sub insert_gencode($$;$)
{
	my ($data, $inpath, $gene_list) = @_;
	my ($runtimes) = undef;
		
	foreach my $gene (@{$data}) {
		my ($gene_id) = $gene->stable_id;
		my ($gene_eid) = $gene_id;
		#if ( $gene->version ) {
		#	$gene_eid = $gene_id.'.'.$gene->version;
		#}
		if ( defined $gene_list ) { # if there is a gene list, we run appris for the list
			if ( exists $gene_list->{$gene_id} ) {
				my ($runtime) = insert_pipeline($gene_eid, $inpath, $gene);
				push(@{$runtimes}, $runtime);				
			}
		}
		else { # if there is not gene list, we run appris for all of them
			my ($runtime) = insert_pipeline($gene_eid, $inpath, $gene);
			push(@{$runtimes}, $runtime);			
		}
	}
	
	return $runtimes;
	
} # end insert_gencode

sub insert_sequence($$)
{
	my ($data, $inpath) = @_;
	my ($runtimes) = undef;
	
	if ( defined $data ) {
		while ( my ($gene_id,$gene) = each(%{$data}) ) {
			my ($runtime) = insert_pipeline($gene_id, $inpath, $gene);
			push(@{$runtimes}, $runtime);
		}		
	}
				
	return $runtimes;
} # end insert_sequence

sub insert_pipeline($$;$)
{
	my ($gene_id, $workspace, $gene) = @_;
	my ($runtime) = undef;
	
	# create parameters
	my ($params);
	$logger->info("\t-- create parameters ");

	# data from gencode type
	if ( $type_of_input =~ /datafile/ and defined $gene ) {
		$logger->info("from $type_of_input\n");
		if ( $type_of_input eq 'gencode-list' ) {
			my ($chr) = $gene->chromosome;
			$workspace .= "chr$chr".'/'.$gene_id;
		}
		else {
			$workspace .= $gene_id;			
		}
		
		$logger->info("\t-- prepare params\n");	
		$params = {
			'id'				=> $gene_id,
			'species'			=> "'$species'",
			'inpath'			=> $workspace,
		};
	}
	# data from ensembl type || sequence type
	elsif ( $type_of_input =~ /ensembl/ or $type_of_input =~ /sequence/ ) {
		$logger->info("from $type_of_input\n");
		$workspace .= $gene_id;

		$logger->info("\t-- prepare params\n");		
		$params = {
			'id'				=> $gene_id,
			'species'			=> "'$species'",
			'inpath'			=> $workspace,
		};
	}	
	else {
		$logger->info("\t-- do not run\n");
		return $runtime;
	}
	$params->{'methods'} = $methods if ( defined $methods );	
	$params->{'apprisdb-conf'} = $apprisdb_conf_file if ( defined $apprisdb_conf_file );
	
	# run pipe
	$logger->info("\t-- run pipe\n");
	my ($inserttime) = insert_appris($gene_id, $workspace, $params);
	$logger->error("inserting results") unless ( defined $inserttime );

	return $runtime;
	
} # end run_pipe_iappris

sub insert_appris($$$)
{
	my ($id, $workspace, $params) = @_;
	
	# get inputs
	my ($parameters) = '';
	if ( defined $params ) {
		while ( my ($k,$v) = each(%{$params}) ) {
			$parameters .= " --$k=$v ";
		}
	}
	
	# create appris job script for cluster
	my ($c_wspace) = $workspace;
	my ($c_id) = $id;
	my ($c_logpath) = $c_wspace;
	my ($c_logfile) = 'log';

	# run
	eval {
		my ($cmd) =	" perl $ENV{APPRIS_CODE_DIR}/iappris.pl ".
					" $parameters ".
					" --loglevel=$LOGLEVEL --logpath=$c_logpath --logfile=$c_logfile $LOGAPPEND ";
		$logger->info("\n** script: $cmd\n");
		my (@cmd_out) = `$cmd`;
	};
	$logger->error("running appris") if($@);

} # end insert_appris


main();


1;

__END__

=head1 NAME

insert_appris

=head1 DESCRIPTION

global script that insert results of APPRIS into database 

=head1 SYNOPSIS

insert_appris

=head2 Required arguments (inputs):

=head3 GENCODE choice: executes appris from gencode data (http://www.gencodegenes.org/data.html)
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--data=  <Gene annotation file>
	
	--inpath= <Acquire input files from PATH>

=head4 Optional arguments (exclusived):
	
	--id= <Gene identifier>
						
	--gene-list= <File with a list of genes>
			
	--position= <Genome position (eg1: 21. eg2: 21,22)>
	
=head3 ENSEMBL choice: executes appris from ensembl data ( -api version- )
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--e-version= <Number of Ensembl version of identifier>
	
	--inpath= <Acquire input files from PATH>
	
=head4 Optional arguments (exclusived):
			
	--id= <Ensembl gene identifier>
						
	--gene-list= <File with a list of genes>
			
	--position= <Genome position (eg1: 21. eg2: 21,22)>
			
=head2 Optional arguments (methods):
		
  --methods= <List of APPRIS's methods ('firestar,matador3d,spade,corsair,thump,crash,appris'. Default: ALL)>
		
=head2 Optional arguments (config files):

  --ensembldb-conf= <Config file of Ensembl database (default: 'conf/ensembldb.ini' file)>
  
  --apprisdb-conf= <Config file of APPRIS database (default: 'conf/apprisdb.ini' file)>

=head2 Optional arguments (log arguments):
	
  --loglevel=LEVEL <define log level (default: NONE)>
  
  --logfile=FILE <Log to FILE (default: *STDOUT)>
  
  --logpath=PATH <Write logfile to PATH (default: .)>
  
  --logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE of GENCODE's type of input

=head2 1. Main execution:

perl insert_appris.pl

	--species='Homo sapiens'
	
	--data={ABSOLUTE_PATH}/gencode.v15.annotation.gtf
	
	--inpath={ABSOLUTE_PATH}/annotations/
	
	--apprisdb-conf={ABSOLUTE_PATH}/conf/apprisdb.mus70.ini
		
	--methods=appris
	
	--logpath={ABSOLUTE_PATH}/logs/
			
	--logfile=insert_appris.log
	
	--loglevel=INFO
	
	--logappend
	
	
=head1 EXAMPLE of ENSEMBL's type of 'input

=head2 1. With ensembl gene identifier:

perl insert_appris.pl

	--species='Mus musculus'
	
	--e-version=70
	
	--id=ENSMUSG00000017167.13

	--inpath={ABSOLUTE_PATH}/ENSMUSG00000017167.13/

	--methods=firestar,matador3d
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
