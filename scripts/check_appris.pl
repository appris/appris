#!/usr/bin/perl -W

use strict;
use warnings;
use threads;

use Getopt::Long;
use FindBin;
use Config::IniFiles;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use appris qw( get_gene_list create_ensembl_input create_gencode_input );
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
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_ENSEMBL_DB_FILE	= $LOCAL_PWD.'/conf/ensembldb.ini';
$LOGLEVEL					= 'INFO';

# Input parameters
my ($id) = undef;
my ($position) = undef;
my ($gene_list_file) = undef;
my ($data_file) = undef;
my ($species) = undef;
my ($e_version) = undef;
my ($inpath) = undef;
my ($outfile) = undef;
my ($type_of_input) = undef;
my ($ensembldb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'id=s'				=> \$id,
	'position=s'		=> \$position,	
	'gene-list=s'		=> \$gene_list_file,
	'data=s'			=> \$data_file,
	'species=s'			=> \$species,
	'e-version=s'		=> \$e_version,
	'inpath=s'			=> \$inpath,
	'outfile=s'			=> \$outfile,
	'ensembldb-conf=s'	=> \$ensembldb_conf_file,	
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
(
	# gencode-position choice
	(
		defined  $position and
		defined  $data_file
	) or 
	# gencode-list choice
	(
		defined  $gene_list_file and
		defined  $data_file
	) or
	# gencode choice
	(
		defined $data_file
	) or 
	# sequence choice
	(
		defined $id
	) or
	# ensembl-position choice
	defined $position or
	# ensembl-list choice
	defined $gene_list_file or
	# ensembl choice
	defined $id
) and 
	# required
	defined  $species and
	defined  $e_version and	
	defined  $inpath and
	defined  $outfile
){
	print `perldoc $0`;
	exit 1;
}

# Get the type of input (the order of conditions is important)
if ( defined $position and defined $data_file ) {
	$type_of_input = 'gencode-position';
}
elsif ( defined $gene_list_file and defined $data_file ) {	
	$type_of_input = 'gencode-list';
}
elsif ( defined $data_file ) {
	$type_of_input = 'gencode';
}
elsif ( defined $id ) {
	$type_of_input = 'sequence';
}
elsif ( defined $position ) {
	$type_of_input = 'ensembl-position';
}
elsif ( defined $gene_list_file ) {
	$type_of_input = 'ensembl-list';
}
elsif ( defined $id ) {
	$type_of_input = 'ensembl';
}

# Optional arguments
# get vars of ensembl db
unless ( defined $ensembldb_conf_file ) {
	$ensembldb_conf_file = $CONFIG_INI_ENSEMBL_DB_FILE;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();

#################
# Method bodies #
#################
sub check_ensembl($$$);
sub check_gencode($$$;$);
sub check_pipeline($$$;$);
sub check_appris($$$);

# Main subroutine
sub main()
{
	# using gene list
	my ($gene_list);
	if ( defined $gene_list_file ) {
		$logger->info("-- using gene list\n");
		my ($genes) = getTotalStringFromFile($gene_list_file);
		foreach my $gene_id (@{$genes}) {
			$gene_id =~ s/\s*//mg;			
			$gene_list->{$gene_id} = 1;
		}		
	}

	# run appris pipeline for each gene depending on input
	$logger->info("-- from given input...");
	if ( ($type_of_input eq 'ensembl') or ($type_of_input eq 'ensembl-position') or ($type_of_input eq 'ensembl-list') ) {
		$logger->info(" $type_of_input type\n");
		
		# get gene list
		my ($gene_list);
		if ( defined $id ) {
			$logger->info("-- using gene id\n");
			$gene_list->{$id} = 1;
		}
		elsif ( defined $position ) {
			$logger->info("-- using genome position\n");
			$gene_list = appris::create_ensembl_input($position, $ensembldb_conf_file, $e_version, $species);		
		}
		elsif ( defined $gene_list_file ) {
			$logger->info("-- using gene list\n");
			$gene_list = appris::get_gene_list($gene_list_file);
		}
		
		$logger->info("-- check ensembl ids\n");
		my ($runtimes) = check_ensembl($gene_list, $inpath, $outfile);

	}
	elsif ( ($type_of_input eq 'gencode') or ($type_of_input eq 'gencode-position') or ($type_of_input eq 'gencode-list') ) {
		$logger->info(" $type_of_input type\n");
		
		# get gene list
		my ($gene_list);
		my ($data_fh);
		if ( defined $id ) {
			$logger->info("-- using gene id\n");
			$gene_list->{$id} = 1;
		}
		elsif ( defined $position ) {
			$logger->info("-- using genome position\n");
			$data_fh = appris::create_gencode_input($data_file, $position);
			if ( UNIVERSAL::isa($data_fh,'File::Temp') ) {
				$data_file = $data_fh->filename;
			}
		}		
		elsif ( defined $gene_list_file ) {
			$logger->info("-- using gene list\n");
			$gene_list = appris::get_gene_list($gene_list_file);
		}
		
		$logger->info("-- create gencode data files\n");
		my ($gencode_data) = appris::create_gencode_data($data_file);		
		
		# delete tmp file
		if ( defined $data_fh and UNIVERSAL::isa($data_fh,'File::Temp') ) {
			$data_fh->unlink_on_destroy(1);
		}

		$logger->info("-- run gencode data files\n");
		my ($runtimes) = check_gencode($gencode_data, $inpath, $outfile, $gene_list);
				
		$logger->info("-- print times\n");
		foreach my $runtime (@{$runtimes}) {
			#$logger->info($runtime->{'gene_id'}."\t".$runtime->{'run'});
		}
	}
	else {
		$logger->error("analying input parameters");
	}	
	
	$logger->finish_log();
	
	exit 0;	
}

sub check_ensembl($$$)
{
	my ($data, $inpath, $outfile) = @_;
	my ($runtimes) = undef;
			
	foreach my $gene_id (keys(%{$data})) {
		my ($runtime) = check_pipeline($gene_id, $inpath, $outfile);
		push(@{$runtimes}, $runtime);		
	}
	return $runtimes;
} # end check_ensembl

sub check_gencode($$$;$)
{
	my ($data, $inpath, $outfile, $gene_list) = @_;
	my ($runtimes) = undef;
		
	foreach my $gene (@{$data}) {
		my ($gene_id) = $gene->stable_id;
		my ($gene_ver) = $gene->version;
		my ($gene_eid) = $gene_id.'.'.$gene_ver;		
		if ( defined $gene_list ) { # if there is a gene list, we run appris for the list
			if ( exists $gene_list->{$gene_id} ) {
				my ($runtime) = check_pipeline($gene_eid, $inpath, $outfile, $gene);
				push(@{$runtimes}, $runtime);				
			}
		}
		else { # if there is not gene list, we run appris for all of them
			my ($runtime) = check_pipeline($gene_eid, $inpath, $outfile, $gene);
			push(@{$runtimes}, $runtime);			
		}
	}
	return $runtimes;
} # end check_gencode

sub check_pipeline($$$;$)
{
	my ($gene_id, $inpath, $outfile, $gene) = @_;
	my ($runtime) = undef;
	my ($workspace) = $inpath.'/';

	$logger->info("-- $gene_id\n");	
	
	# create parameters
	$logger->info("\t-- create parameters ");
	my ($params);
	# data from gencode type
	if ( defined $gene and defined $gene->chromosome ) {
		$logger->info("from gencode data\n");
		#my ($chr) = $gene->chromosome;
		#$workspace .= "chr$chr". '/'.$gene_id;
		$workspace .= $gene_id;
				
		$logger->info("\t-- prepare params\n");
		$params = {
			'id'				=> $gene_id,
			'inpath'			=> $workspace,
			'outfile'			=> $outfile,
		};
	}
	# data from ensembl type
	elsif ( defined $gene_id ) {
		$logger->info("from ensembl data\n");
		$workspace .= '/'.$gene_id;

		$logger->info("\t-- prepare params\n");		
		$params = {
			'id'				=> $gene_id,
			'inpath'			=> $workspace,
			'outfile'			=> $outfile,
		};
	}
	else {
		return $runtime;
	}
		
	# run appris pipeline
	if ( defined $params ) {
		
		# check results of pipeline
		$logger->info("\t-- check results\n");
		my ($checktime) = check_appris($gene_id, $workspace, $params);
		throw("checking results") unless ( defined $checktime );
		return $runtime;		
		
	}
	else {
		$logger->info("\t-- do not run appris\n");
	}

	return $runtime;
	
} # end check_pipeline

sub check_appris($$$)
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
	my ($c_logfile) = $c_id.'.log';

	# run
	eval {
		my ($cmd) =	" perl $ENV{APPRIS_HOME}/scripts/cappris.pl ".
					" $parameters ".
					" --loglevel=$LOGLEVEL --logpath=$c_logpath --logfile=$c_logfile --logappend ";			
		$logger->info("\n** script: $cmd\n");
		my (@cmd_out) = `$cmd`;
	};
	throw("running appris") if($@);

} # end check_appris




main();


1;

__END__

=head1 NAME

check_appris

=head1 DESCRIPTION

global script that check the results of APPRIS 

=head1 SYNOPSIS

check_appris

=head2 Input arguments:

	* Gencode choice: executes appris for one or more genes (using data from Gencode -cds, exons, seqs, etc.- )

		--data=  <Gene annotation file>
		
		--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
		--e-version= <Number of Ensembl version of identifier>

		--inpath= <Acquire input files from PATH>
		
		--outfile= <Output file>
		
	* Gencode-position choice: executes appris for a genome region (using data from Gencode -cds, exons, seqs, etc.- )

		--position= <Genome position> (21 or 21,22)
	
		--data=  <Gene annotation file>
		
		--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
		--e-version= <Number of Ensembl version of identifier>

		--inpath= <Acquire input files from PATH>
		
		--outfile= <Output file>
		
	* Gencode-list choice: executes appris for list of genes (using data from Gencode -cds, exons, seqs, etc.- )

		--gene-list=  <File with a list of genes>
	
		--data=  <Gene annotation file>
		
		--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
		--e-version= <Number of Ensembl version of identifier>

		--inpath= <Acquire input files from PATH>
		
		--outfile= <Output file>
		
	* Ensembl choice: executes appris for one gene (using data from Ensembl -api version- )

		--id= <Ensembl gene identifier>
		
		--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
		--e-version= <Number of Ensembl version of identifier>

		--inpath= <Acquire input files from PATH>
		
		--outfile= <Output file>
		
	* Ensembl-position choice: executes appris for a genome region (using data from Ensembl -api version- )

		--position= <Genome position> (21 or 21,22)
		
		--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
		--e-version= <Number of Ensembl version of identifier>

		--inpath= <Acquire input files from PATH>
		
		--outfile= <Output file>
		
	* Ensembl-list choice: executes appris for list of genes (using data from Ensembl -api version- )

		--gene-list=  <File with a list of genes>
	
		--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
		--e-version= <Number of Ensembl version of identifier>

		--inpath= <Acquire input files from PATH>
		
		--outfile= <Output file>
								
=head2 Optional arguments (log arguments):
	
		--loglevel=LEVEL <define log level (default: NONE)>	
	
		--logfile=FILE <Log to FILE (default: *STDOUT)>
		
		--logpath=PATH <Write logfile to PATH (default: .)>
		
		--logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl check_appris.pl

	--id=ENSMUSG00000017167
	
	--species='Homo sapiens'
	
	--e-version=70

	--inpath=../features/ENSMUSG00000017167_e69/
	
	--outfile=../features/ENSMUSG00000017167_e69/ENSMUSG00000017167.runtime.out

=head1 EXAMPLE

perl check_appris.pl

	--data=/home/jmrodriguez/projects/Encode/gencode15/features/chr21.gencode.v15.annotation.gtf
	
	--species='Homo sapiens'
	
	--e-version=70
	
	--inpath=/home/jmrodriguez/projects/Encode/gencode15/annotations/
	
	--outfile=/home/jmrodriguez/projects/Encode/gencode15/annotations/chr21.runtime.out
	
	--loglevel=INFO
	
	--logappend
	
	--logpath=/home/jmrodriguez/projects/Encode/gencode15/logs/
			
	--logfile=check_appris.chr21.log
	
=head1 EXAMPLE

perl check_appris.pl

	--position=21
	
	--data=/home/jmrodriguez/projects/Encode/gencode15/features/gencode.v15.annotation.gtf
	
	--species='Homo sapiens'
	
	--e-version=70
	
	--inpath=/home/jmrodriguez/projects/Encode/gencode15/annotations/
	
	--outfile=/home/jmrodriguez/projects/Encode/gencode15/annotations/chr21.runtime.out
	
	--loglevel=INFO
	
	--logappend
	
	--logpath=/home/jmrodriguez/projects/Encode/gencode15/logs/
			
	--logfile=check_appris.chr21.log

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
