#!/usr/bin/perl -W

use strict;
use warnings;
use threads;

use Getopt::Long;
use FindBin;
use Config::IniFiles;
use File::Temp;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use appris qw( create_ensembl_input );
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
	$LOGGER_CONF
	@APPRIS_CHR_LIST
	$MAX_NUM_GENES
	$ENSEMBL_CHR_GENE_LIST
	$GENETYPE_PROTEIN_CODING
	$GENETYPE_LONG_NON_CODING
	$GENETYPE_SHORT_NON_CODING
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_ENSEMBL_DB_FILE	= $LOCAL_PWD.'/conf/ensembldb.ini';
$LOGGER_CONF				= '';
$ENSEMBL_CHR_GENE_LIST		= undef;
$GENETYPE_PROTEIN_CODING	= {
	'protein_coding' => undef,
};
$GENETYPE_LONG_NON_CODING	= {
	'3prime_overlapping_ncrna' => undef,
	'ambiguous_orf' => undef,
	'antisense' => undef,
	'antisense_RNA' => undef,
	'lincRNA' => undef,
	'ncrna_host' => undef,
	'non_coding' => undef,
	'non_stop_decay' => undef,
	'processed_transcript' => undef,
	'retained_intron' => undef,
	'sense_intronic' => undef,
	'sense_overlapping' => undef,
};
$GENETYPE_SHORT_NON_CODING	= {
	'miRNA' => undef,
	'miscRNA' => undef,
	'rRNA' => undef,
	'tRNA' => undef,
	'ncRNA' => undef,
	'scRNA' => undef,
	'snlRNA' => undef,
	'snoRNA' => undef,
	'snRNA' => undef,
	'tRNA' => undef,
	'pseudogenic' => undef,
};


# Input parameters
my ($genetype) = undef;
my ($random) = undef;
my ($position) = undef;
my ($species) = undef;
my ($e_version) = undef;
my ($outdir) = undef;
my ($outfile) = undef;
my ($outpath) = undef;
my ($methods) = undef;
my ($ensembldb_conf_file) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'genetype=s'		=> \$genetype,
	'random=s'			=> \$random,
	'position=s'		=> \$position,
	'species=s'			=> \$species,
	'e-version=s'		=> \$e_version,
	'outdir=s'			=> \$outdir,
	'outfile=s'			=> \$outfile,	
	'ensembldb-conf=s'	=> \$ensembldb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined  $genetype and
	defined  $species and
	defined  $e_version and
	defined  $outdir and
	defined  $outfile
){
	print `perldoc $0`;
	exit 1;
}

# Optional arguments
# get vars of ensembl db
unless ( defined $ensembldb_conf_file ) {
	$ensembldb_conf_file = $CONFIG_INI_ENSEMBL_DB_FILE;
}
if ( defined $random ) {
	$MAX_NUM_GENES = $random;
}


# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();
$LOGGER_CONF .= " --loglevel=$loglevel " if ( defined $loglevel );
$LOGGER_CONF .= " --logpath=$logpath " if ( defined $logpath );
$LOGGER_CONF .= " --logfile=$logfile " if ( defined $logfile );
$LOGGER_CONF .= " --logappend " if ( defined $logappend );

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	# get ensembl data from 
	my ($output) = '';
	if ( defined $random ) {
		$logger->info("-- using $random random\n");
		$output .= create_random_data($random);
	}
	elsif ( defined $position ) {
		$logger->info("-- using $position position\n");
$logger->info("-- POS: $position\n".UNIVERSAL::isa($position, 'SCALAR')."\n");		
		$output .= create_data($position);
	}
	else {
		$logger->info("-- using genome position\n");
		my (@chr_list) = split(',',$ENV{'APPRIS_CHR_LIST'});	
		$output .= create_data(\@chr_list);		
	}
		
	# print output
	$logger->info("-- print output\n");
	if ( $output ne '' ) {
		my ($printing_file_log) = printStringIntoFile($output, $outfile);
		throw("creating $outfile file") unless ( defined $printing_file_log );
	}
	
	$logger->finish_log();
	
	exit 0;	
}

sub create_random_data($)
{
	my ($max_num_genes) = @_;
	my ($output) = '';
	
	# acquire ensembl annot from gene list
	my ($num_genes) = 0;
	while ( $num_genes < $max_num_genes ) {	
		my ($gene_list) = create_ensembl_random($max_num_genes-$num_genes, $ensembldb_conf_file, $e_version, $species);
		$logger->info("-- acquire ensembl annot from gene list\n");
		while (my ($gene_id,$gene_report) = each(%{$gene_list}) ) {
			my ($input_files) = run_getedata($gene_id, $species, $e_version);
			if ( defined $input_files ) {
				$num_genes++;
				$output .= $gene_id."\n";
				last if ( $num_genes == $max_num_genes );
			}
		}
		$logger->info("\t: ".$num_genes."\n");		
	}
	return $output;
}

sub create_data($)
{
	my ($position) = @_;
	my ($output) = '';
	my ($gene_list);
	
	# get gene list
	if ( ref($position) eq 'ARRAY' ) {
		foreach my $chr (@{$position}) {
			$logger->info("-- using genome position\n");
			my ($g_list) = create_ensembl($chr, $ensembldb_conf_file, $e_version, $species);
			while (my ($g_id, $g_report) = each(%{$g_list}) ) {
				$gene_list->{$g_id} = $g_report;
			}			
		}
	}
	else {
		$gene_list = create_ensembl($position, $ensembldb_conf_file, $e_version, $species);
	}
	
	# acquire ensembl annot from gene list
	$logger->info("-- acquire ensembl annot from gene list\n");
	while (my ($gene_id,$gene_report) = each(%{$gene_list}) ) {
		my ($input_files) = run_getedata($gene_id, $species, $e_version);
		if ( defined $input_files ) {
			$output .= $gene_id."\n";
		}
	}
	return $output;
}

sub create_ensembl_random($$$$)
{
	my ($num_genes, $ensembldb_conf_file, $e_version, $species) = @_;
	my ($gene_list);
	
	# get gene list of ensembl from random chr
	my (@APPRIS_CHR_LIST) = split(',',$ENV{'APPRIS_CHR_LIST'});	
	my ($chr) = $APPRIS_CHR_LIST[int(rand(25))];
	$logger->info("-- get gene list of ensembl from $chr\n");
	unless ( exists $ENSEMBL_CHR_GENE_LIST->{$chr} ) {
		$ENSEMBL_CHR_GENE_LIST->{$chr} = appris::create_ensembl_input($chr, $ensembldb_conf_file, $e_version, $species);
	}
	my (@g_keys) = keys(%{$ENSEMBL_CHR_GENE_LIST->{$chr}});
	
	# get random num depending on the total number of genes of current chr
	my ($n_g_list) = scalar(@g_keys);
	my ($n_max) = $n_g_list;
	if ( $num_genes < $n_g_list ) { $n_max = $num_genes; }
	my ($n) = $n_max/int(rand(4)+1);
	my ($tmp_n_g) = int(rand($n));
	if ($tmp_n_g==0){$tmp_n_g=1};	
	my ($n_g) = 0;
	while ( $n_g < $tmp_n_g ) {
		my ($g_i) = int(rand($n_g_list)); # get random gene from the list
		my ($gene_id) = $g_keys[$g_i];
		unless ( exists $gene_list->{$gene_id} ) {
			my ($biotype) = $ENSEMBL_CHR_GENE_LIST->{$chr}->{$gene_id};
			#$gene_list->{$gene_id} = $chr;
			if ( $genetype eq 'coding' ) {
				$gene_list->{$gene_id} = $biotype if ( exists $GENETYPE_PROTEIN_CODING->{$biotype} );
				$n_g++;
			}
			elsif ( $genetype eq 'long' ) {
				$gene_list->{$gene_id} = $biotype if ( exists $GENETYPE_LONG_NON_CODING->{$biotype} );
				$n_g++;				
			}
			elsif ( $genetype eq 'short' ) {
				$gene_list->{$gene_id} = $biotype if ( exists $GENETYPE_SHORT_NON_CODING->{$biotype} );
				$n_g++;
			}
		}
	}
	return $gene_list;
}

sub create_ensembl($$$$)
{
	my ($chr, $ensembldb_conf_file, $e_version, $species) = @_;
	my ($gene_list);
	
	# get gene list of ensembl from random chr
	$logger->info("-- get gene list of ensembl from $chr\n");
	unless ( exists $ENSEMBL_CHR_GENE_LIST->{$chr} ) {
		$ENSEMBL_CHR_GENE_LIST->{$chr} = appris::create_ensembl_input($chr, $ensembldb_conf_file, $e_version, $species);
	}
	my (@g_keys) = keys(%{$ENSEMBL_CHR_GENE_LIST->{$chr}});
	foreach my $gene_id (@g_keys) {
		unless ( exists $gene_list->{$gene_id} ) {
			my ($biotype) = $ENSEMBL_CHR_GENE_LIST->{$chr}->{$gene_id};
			#$gene_list->{$gene_id} = $chr;
			if ( $genetype eq 'coding' ) {
				$gene_list->{$gene_id} = $biotype if ( exists $GENETYPE_PROTEIN_CODING->{$biotype} );
			}
			elsif ( $genetype eq 'long' ) {
				$gene_list->{$gene_id} = $biotype if ( exists $GENETYPE_LONG_NON_CODING->{$biotype} );
			}
			elsif ( $genetype eq 'short' ) {
				$gene_list->{$gene_id} = $biotype if ( exists $GENETYPE_SHORT_NON_CODING->{$biotype} );
			}
		}		
	}
	return $gene_list;
}

sub run_getedata($$$)
{
	my ($id, $species, $e_version) = @_;
	
	my ($wsdir) = $outdir.'/'.$id;
	mkdir($wsdir);
	my ($data_tmpfilename) = $wsdir.'/'.$id.'.annot.gtf';
	my ($pdata_tmpfilename) = $wsdir.'/'.$id.'.pannot.gtf';
	my ($transc_tmpfilename) = $wsdir.'/'.$id.'.transc.fa';
	my ($transl_tmpfilename) = $wsdir.'/'.$id.'.transl.fa';
	my ($input_files) = {
		'annot'			=> $data_tmpfilename,
		'pannot'		=> $pdata_tmpfilename,
		'transc'		=> $transc_tmpfilename,
		'transl'		=> $transl_tmpfilename
	};
	
	# check if gene data already exists
	if (
		-e $data_tmpfilename and (-s $data_tmpfilename > 0) and 
		-e $transc_tmpfilename and (-s $transc_tmpfilename > 0) 
	) {
		$logger->info("\ngene data exists: $id\n");		
		my (@ls_out) = `ls -l $wsdir/*.align.faa`;
		if ( scalar(@ls_out) > 0 ) {
			$logger->info("\ngene aligns exist: $id\n");
			return undef;
		}
	}
	
	# get annot data from ensembl
	eval {
		my ($cmd) = "perl $ENV{APPRIS_PROGRAMS_SRC_DIR}/ensembl/getGTF.pl ".
						"--id=$id ".
						"--species='$species' ".
						"--e-version=$e_version ".
						"--out-data=$data_tmpfilename ".
						"--out-pdata=$pdata_tmpfilename ".
						"--out-transcripts=$transc_tmpfilename ".
						"--out-translations=$transl_tmpfilename ".
						"$LOGGER_CONF ";
		$logger->debug("\n** script: $cmd\n");
		system ($cmd);
	};
	if($@) {
		eval { # rm rst if does not work
			my ($cmd) = "rm -rf $wsdir";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		return undef;		
	}
	unless (
		-e $data_tmpfilename and (-s $data_tmpfilename > 0) and 
		-e $transc_tmpfilename and (-s $transc_tmpfilename > 0) 				
	) {
		eval { # rm rst if does not work
			my ($cmd) = "rm -rf $wsdir";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		return undef;
	}
	
	# get ECompara alignments
	my ($ecomparameter) = '';
	if ( $genetype eq 'coding' ) {
		my ($atype) = 'cds';
		$ecomparameter .= "--atype=$atype ".
						"--data=$data_tmpfilename ".
						"--transcripts=$transc_tmpfilename ".
						"--translations=$transl_tmpfilename ";		
	}
	elsif ( ($genetype eq 'long') or ($genetype eq 'short') ) {
		my ($atype) = 'exon';
		$ecomparameter .= "--atype=$atype ".
						"--data=$data_tmpfilename ".
						"--transcripts=$transc_tmpfilename ";
	}	
	eval {
		my ($cmd) = "perl $ENV{APPRIS_PROGRAMS_SRC_DIR}/ensembl/getEComparaAlign.pl ".
						"--species='$species' ".
						"--e-version=$e_version ".
						$ecomparameter.						
						"--outpath=$wsdir ".
						"$LOGGER_CONF ";
		$logger->info("\n** script: $cmd\n");
		system ($cmd);
	};
	if($@) {
		eval { # rm rst if does not work
			my ($cmd) = "rm -rf $wsdir";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		return undef;		
	}
	my (@ls_out) = `ls -l $wsdir/*.align.faa`;
	if ( scalar(@ls_out) == 0 ) {
		eval { # rm rst if does not work
			my ($cmd) = "rm -rf $wsdir";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		return undef;		
	}
	return ($input_files);
}

main();


1;

__END__

=head1 NAME

random_ensembl_genes

=head1 DESCRIPTION

retrieve list of genes from ENSEMBL 

=head1 SYNOPSIS

apprisall

=head2 Arguments:

	--genetype= <Gene typ: [coding,long,short] (default: coding)>

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--e-version= <Number of Ensembl version of identifier>	

	--outdir= <Output dir>
	
	--outfile= <Output file>
		
=head2 Optional arguments:

	--random= <Num. random genes>
	
	--position= <Genome position: 21 (default: all genome)> 

	--ensembldb-conf= <Config file of Ensembl database (default: 'conf/ensembldb.ini' file)>
		
=head2 Optional arguments (log arguments):
	
		--loglevel=LEVEL <define log level (default: NONE)>	
	
		--logfile=FILE <Log to FILE (default: *STDOUT)>
		
		--logpath=PATH <Write logfile to PATH (default: .)>
		
		--logappend= <Append to logfile (default: truncate)>

=head1 EXAMPLE

random_ensembl_genes

	--species='Mus musculus'
	
	--e-version=70
	
	--outdir=rand_ensembl_genes/
	
	--outfile=genelist_musculus70.chr4.txt

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
