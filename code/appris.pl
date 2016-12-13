#!/usr/bin/perl -w
# _________________________________________________________________
# $Id$
# $Revision$
# Developed by:
#		Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es-
# _________________________________________________________________

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use Config::IniFiles;
use File::Basename;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::File qw( prepare_workspace rm_dir printStringIntoFile getStringFromFile );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$SRC_DIR
	$DEFAULT_CONFIG_FILE
	$DEFAULT_CFG
	$CFG
	$LOGGER_CONF
	$METHOD_STRUCT
);

$LOCAL_PWD				= $FindBin::Bin;
$SRC_DIR				= $LOCAL_PWD.'/src/';
$DEFAULT_CONFIG_FILE	= $ENV{APPRIS_CODE_CONF_DIR}.'/pipeline.ini';
$DEFAULT_CFG			= new Config::IniFiles( -file => $DEFAULT_CONFIG_FILE );
$CFG					= undef;
$METHOD_STRUCT			= [ split( ',', $DEFAULT_CFG->val('APPRIS_PIPELINE', 'structure') ) ];
$LOGGER_CONF			= '';

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($id) = undef;
my ($data_file) = undef;
my ($pdata_file) = undef;
my ($transc_file) = undef;
my ($transl_file) = undef;
my ($inpath) = undef;
my ($species) = undef;
my ($e_version) = undef;
my ($outpath) = undef;
my ($methods) = undef;
my ($type_of_input) = undef;
my ($type_of_align) = undef;
my ($cached_path) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'id=s'					=> \$id,
	'data=s'				=> \$data_file,
	'pdata=s'				=> \$pdata_file,
	'transc=s'				=> \$transc_file,
	'transl=s'				=> \$transl_file,
	'species=s'				=> \$species,
	'e-version=s'			=> \$e_version,		
	'inpath=s'				=> \$inpath,
	'outpath=s'				=> \$outpath,		
	'methods=s'				=> \$methods,
	't-align=s'				=> \$type_of_align,
	'cached-path=s'			=> \$cached_path,
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);

# Required arguments
unless (
	(
		# first choice
		defined $inpath
		or
		# second choice
		defined $id
		or
		# third choice
		(defined $data_file and
		defined  $pdata_file and
		defined  $transc_file and
		defined  $transl_file)
		or
		# fourth choice
		defined  $transl_file
	)
	and 	
	defined  $species and
	defined  $outpath
){
	print `perldoc $0`;
	exit 1;
}

# Get method's pipeline
#unless ( defined $methods ) {
#	$methods = $DEFAULT_CFG->val('APPRIS_PIPELINE', 'structure');
#}

# Get the type of input (the order of conditions is important)
if ( defined $inpath ) {
	$type_of_input = 'inpath';
}
elsif ( defined $data_file and defined $pdata_file and defined $transc_file and defined $transl_file ) {
	$type_of_input = 'gencode';
}
elsif ( defined $transl_file ) {
	$type_of_input = 'sequence';
}
elsif ( defined $id and defined  $e_version ) {
	$type_of_input = 'ensembl';
}

# Get the type of input
if ( defined $type_of_align ) {
	#if ( lc($type_of_align) =~ /^ucsc\d*/ ) {
	if ( lc($type_of_align) eq 'ucsc' ) {
		$type_of_align = lc($type_of_align);
	}
	#elsif ( lc($type_of_align) =~ /^compara\d*/ ) {
	elsif ( (lc($type_of_align) eq 'compara') and defined $e_version ) { # has to be defined e_version
		$type_of_align = lc($type_of_align);
	}
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);
$LOGGER_CONF .= " --loglevel=$loglevel " if ( defined $loglevel );
$LOGGER_CONF .= " --logpath=$logpath " if ( defined $logpath );
$LOGGER_CONF .= " --logfile=$logfile " if ( defined $logfile );
$LOGGER_CONF .= " --logappend " if ( defined $logappend );


#####################
# Method prototypes #
#####################
sub create_workspace();
sub run_getgtf($$$$);
sub create_ini($);
sub create_inputs($);
sub run_getmafucsc($$$$);
sub run_getecompara($$$$$);
sub run_pipeline($$$);

sub _subs_template($$$);
sub _rm_files($);
sub _cp_files($$);

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	# prepare workspace
	$logger->info("-- prepare workspace\n");
	$outpath = prepare_workspace($outpath);
	if ( ($type_of_input eq 'gencode') or ($type_of_input eq 'sequence') ) {
		my ($name, $inpath, $suffix) = fileparse($transl_file, qr/\.[^.]+/);
		my (@inpath_n) = split('/', $inpath);
		if ( scalar(@inpath_n) > 0 ) {
			my ($num) = scalar(@inpath_n);
			$id = $inpath_n[$num-1];
		}
		else {
			$logger->error("id not defined");
		}
	}
	elsif ( $type_of_input eq 'inpath' ) {
		my (@inpath_n) = split('/', $inpath);		
		if ( scalar(@inpath_n) > 0 ) {
			my ($num) = scalar(@inpath_n);
			$id = $inpath_n[$num-1];
		}
		else {
			$logger->error("id not defined");
		}
	}
	
	# Create workspace for pipeline taking into account if results are cached.
	$logger->info("-- get ticket id and prepare workspace\n");
	create_workspace();
	
	# Create ini file
	$logger->info("-- create ini files for appris pipeline\n");
	my ($config_file) = create_ini($methods);
	unless ( defined $config_file ) {
		$logger->error("create ini files for appris pipeline");
	}
	$CFG = new Config::IniFiles( -file =>  $config_file );
	
	# Create inputs for pipeline
	$logger->info("-- create input files for appris pipeline\n");
	my ($input_files) = create_inputs($id);
	unless ( defined $input_files ) {
		$logger->error("creating input files for appris pipeline");
	}

	# Run methods
	$logger->info("-- run pipeline\n");
	my ($output) = run_pipeline($config_file, $input_files, $methods);

	# Delete tmp dir if ...
	if ( !defined $loglevel or (defined $loglevel and $loglevel ne 'debug') ) {	
		$logger->info("-- delete workspace (tmp dir)\n");
		delete_workspace();
	}
	
		
	$logger->finish_log();
	
	exit 0;
}

sub _subs_template($$$)
{
	my ($conf_file, $old, $new) = @_;	
	$conf_file =~ s/$old/$new/g;		
	return $conf_file;		
}

sub _rm_files($)
{
	my ($files) = @_;
	my ($ok);
	my ($list) = '';
	foreach my $file ( @{$files} ) {
		if ( UNIVERSAL::isa($file,'HASH') ) {
			while ( my ($k,$v) = each(%{$file}) ) {
				$list .= " $v ";
			}		
		}
		else {
			$list .= " $file ";	
		}		
	}
	if ( defined $list and $list ne '' ) {
		eval {
			my ($cmd) = "rm $list";
			$logger->info("-- $cmd\n");
			system ($cmd);		
		};
		$logger->error("deleting files: ".$!) if($@);
		$ok = 1;		
	}	
	return $ok;
}

sub _cp_files($$)
{
	my ($org, $dst) = @_;
	my ($ok);
	my ($list) = '';
	foreach my $file ( @{$org} ) {
		$list .= " $file ";
	}
	if ( defined $list and $list ne '' ) {
		eval {
			my ($cmd) = "cp -rp $list $dst/.";
			$logger->info("-- $cmd\n");
			system ($cmd);		
		};
		$logger->error("copying results: ".$!) if($@);
		$ok = 1;		
	}	
	return $ok;
}

sub create_workspace()
{
	# create workspace for each variant
	# create tmp dir for each variant
	my ($in) = Bio::SeqIO->new(
						-file => $transl_file,
						-format => 'Fasta'
	);
	while ( my $seq = $in->next_seq() )
	{
		# get idx for seq
		my ($seq_id) = $seq->id;
		my ($seq_s) = $seq->seq;		
		my ($cache) = APPRIS::Utils::CacheMD5->new(
			-dat => $seq_s,
			-ws  => $ENV{APPRIS_PROGRAMS_CACHE_DIR}			
		);		
		my ($seq_idx) = $cache->idx;
		my ($seq_sidx) = $cache->sidx;
		my ($seq_dir) = $cache->dir;
				
		# prepare cache dir
		my ($ws_cache) = $cache->c_dir();
		# prepare tmp dir
		my ($ws_tmp) = $ENV{APPRIS_TMP_DIR}.'/'.$seq_idx;
		prepare_workspace($ws_tmp);		
		
		$logger->info("-- cache: $seq_id > $seq_sidx\n");
	}
}

sub delete_workspace()
{		
	# delete tmp dir for each variant
	my ($in) = Bio::SeqIO->new(
						-file => $transl_file,
						-format => 'Fasta'
	);
	while ( my $seq = $in->next_seq() )
	{
		# get idx for seq
		my ($seq_id) = $seq->id;
		my ($seq_s) = $seq->seq;
		my ($cache) = APPRIS::Utils::CacheMD5->new(
			-dat => $seq_s,
			-ws  => $ENV{APPRIS_PROGRAMS_CACHE_DIR}
		);		
		my ($seq_idx) = $cache->idx;
		my ($ws_tmp) = $ENV{APPRIS_TMP_DIR}.'/'.$seq_idx;
		rm_dir($ws_tmp);
	}		
}

sub run_getgtf($$$$)
{
	my ($id, $species, $e_version, $datadir) = @_;	
	my ($input_files) = {
		'annot'			=> $datadir.'/'.'annot.gtf',
		'pannot'		=> $datadir.'/'.'pannot.gtf',
		'transc'		=> $datadir.'/'.'transc.fa',
		'transl'		=> $datadir.'/'.'transl.fa'
	};
	my ($data_file, $pdata_file, $transc_file, $transl_file) = ( 
		$input_files->{'annot'},
		$input_files->{'pannot'},
		$input_files->{'transc'},
		$input_files->{'transl'},
	);
	
	# check if input exists	
	unless (
		-e $data_file and (-s $data_file > 0) and 
		-e $pdata_file and (-s $pdata_file > 0) and 
		-e $transc_file and (-s $transc_file > 0) and 
		-e $transl_file and (-s $transl_file > 0)
	) {	
		eval {
			my ($cmd) = "perl $SRC_DIR/ensembl/getGTF.pl ".
							"--id=$id ".
							"--species='$species' ".
							"--e-version=$e_version ".
							"--out-data=$data_file ".
							"--out-pdata=$pdata_file ".
							"--out-transcripts=$transc_file ".
							"--out-translations=$transl_file ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		if($@) {
			my ($a) = _rm_files([$input_files]);
			$logger->error("deleting files") unless ( defined $a );
			return undef;		
		}
		unless (
			-e $data_file and (-s $data_file > 0) and 
			-e $pdata_file and (-s $pdata_file > 0) and 
			-e $transc_file and (-s $transc_file > 0) and 
			-e $transl_file and (-s $transl_file > 0)
		) {
			my ($a) = _rm_files([$input_files]);
			$logger->error("deleting files") unless ( defined $a );
			return undef;				
		}
	}
	
	return $input_files;
}

sub run_getmafucsc($$$$)
{
	my ($species, $t_align, $input_files, $datadir) = @_;
	my ($outpath) = $datadir;
	my ($data_file) = $input_files->{'annot'};
	my ($transc_file) = $input_files->{'transc'};
	my ($transl_file) = $input_files->{'transl'};
	
	# check if input exists
	if ( (`grep -c '>' $transl_file` ne `ls -1 $datadir/*.$t_align.faa | wc -l`) or (`grep -c '>' $transl_file` ne `ls -1 $datadir/*.$t_align.nh | wc -l`) ) {
		eval {
			my ($cmd) = "perl $SRC_DIR/ucsc/getUCSCAlign.pl ".
							"--species='$species' ".
							"--data=$data_file ".
							"--translations=$transl_file ".
							"--outpath=$datadir ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		if($@) {
			my ($a) = _rm_files(["$datadir/*.$t_align.faa"]);
			$logger->error("deleting files") unless ( defined $a );
			my ($b) = _rm_files(["$datadir/*.$t_align.nh"]);
			$logger->error("deleting files") unless ( defined $b );
			return undef;		
		}
		my (@ls_out) = `ls -1 $datadir/*.$t_align.faa`;
		if ( scalar(@ls_out) == 0 ) {
			my ($a) = _rm_files(["$datadir/*.$t_align.faa"]);
			$logger->error("deleting files") unless ( defined $a );
			my ($b) = _rm_files(["$datadir/*.$t_align.nh"]);
			$logger->error("deleting files") unless ( defined $b );
			return undef;		
		}	
	}
		
	return $outpath;
}

sub run_getecompara($$$$$)
{
	my ($species, $t_align, $e_version, $input_files, $datadir) = @_;
	my ($outpath) = $datadir;
	my ($data_file) = $input_files->{'annot'};
	my ($transc_file) = $input_files->{'transc'};
	my ($transl_file) = $input_files->{'transl'};
	
	# check if input exists
	if ( (`grep -c '>' $transl_file` ne `ls -1 $datadir/*.$t_align.faa | wc -l`) or (`grep -c '>' $transl_file` ne `ls -1 $datadir/*.$t_align.nh | wc -l`) ) {
		eval {
			my ($cmd) = "perl $SRC_DIR/ensembl/getEComparaAlign.pl ".
							"--species='$species' ".
							"--e-version=$e_version ".
							"--data=$data_file ".
							"--transcripts=$transc_file ".
							"--translations=$transl_file ".
							"--outpath=$datadir ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		if($@) {
			#my ($a) = _rm_files(["$datadir/*.$t_align.faa"]);
			#$logger->error("deleting files") unless ( defined $a );
			#my ($b) = _rm_files(["$datadir/*.$t_align.nh"]);
			#$logger->error("deleting files") unless ( defined $b );
			return undef;
		}
		my (@ls_out) = `ls -1 $datadir/*.$t_align.faa`;
		if ( scalar(@ls_out) == 0 ) {
			#my ($a) = _rm_files(["$datadir/*.$t_align.faa"]);
			#$logger->error("deleting files") unless ( defined $a );
			#my ($b) = _rm_files(["$datadir/*.$t_align.nh"]);
			#$logger->error("deleting files") unless ( defined $b );
			return undef;
		}	
	}
	
	return $outpath;
}

sub create_ini($)
{
	my ($methods) = @_;
	
	my ($config_cont) = getStringFromFile($DEFAULT_CONFIG_FILE);
	$config_cont = _subs_template($config_cont, 'APPRIS__PIPELINE__METHODS', $methods);
	$config_cont = _subs_template($config_cont, 'APPRIS__SPECIES', $species);
	
	my ($config_file) = $outpath.'/pipeline.ini';
	my ($a) = printStringIntoFile($config_cont, $config_file);
	$logger->error("-- printing config annot") unless ( defined $a );
	
	return ($config_file);
}

sub create_inputs($)
{
	my ($id) = @_;
	my ($input_files);
	#my ($datadir) =  $wspace.'/'.$CFG->val( 'INPUT_VARS', 'name');
		
	# Obtain inputs for pipeline
	$logger->info("-- obtain gene annotations, transcript seq, and translate seq...");
	if ( $type_of_input eq 'inpath' ) {
		$logger->info("from $type_of_input\n");
		$input_files = {
			'id'			=> $id,
			'species'		=> "'$species'"
		};
		my ($ifiles) = {
			'annot'			=> $inpath.'/'.'annot.gtf',
			'pannot'		=> $inpath.'/'.'pannot.gtf',
			'transc'		=> $inpath.'/'.'transc.fa',
			'transl'		=> $inpath.'/'.'transl.fa'
		};
		while ( my ($i, $ifil) = each(%{$ifiles}) ) {
			if ( -e $ifil and (-s $ifil > 0) ) {
				$input_files->{$i} = $ifil;
			}
		}
	}
	elsif ( $type_of_input eq 'gencode' ) {
		$logger->info("from $type_of_input\n");
		$input_files = {
			'id'			=> $id,
			'species'		=> "'$species'",
			'annot'			=> $data_file,
			'pannot'		=> $pdata_file,
			'transc'		=> $transc_file,
			'transl'		=> $transl_file,
		};		
	}	
	elsif ( $type_of_input eq 'sequence' ) {
		$logger->info("from $type_of_input\n");
		$input_files = {
			'id'			=> $id,
			'species'		=> "'$species'",
			'transl'		=> $transl_file,
		};		
	}	
	elsif ( $type_of_input eq 'ensembl' ) {
		$logger->info("from $type_of_input\n");
		$input_files = run_getgtf($id, $species, $e_version, $outpath);
		$logger->error("creating input files for appris pipeline") unless ( defined $input_files );
	}
	else {
		$logger->error("analying input parameter");
	}
	
	# create alignments for pipeline
	if ( defined $type_of_align ) {	
		if ( $type_of_align eq 'ucsc' ) {
			my ($t_align) = 'ucsc';
			$logger->info("-- create alignments...from $t_align\n");
			#my ($alignpath) = run_getmafucsc($species, $t_align, $input_files, $outpath);
			#$logger->error("creating alignments for appris pipeline") unless ( defined $alignpath );
		}
		#if ( $type_of_align =~ /^compara(\d*)/ ) {
		if ( ($type_of_align eq 'compara') and defined $e_version ) {
			#my ($t_align) = 'compara';
			#my ($c_version) = $1;
			my ($t_align) = $type_of_align;			
			$logger->info("-- create alignments...from $t_align\n");
			my ($alignpath) = run_getecompara($species, $t_align, $e_version, $input_files, $outpath);
			$logger->error("creating alignments for appris pipeline") unless ( defined $alignpath );
		}
		$input_files->{'alignpath'} = $outpath;
	}
	# TEMPORAL: For the moment we include always the alignpath
	$input_files->{'alignpath'} = $outpath;
	# TEMPORAL: For the moment we include always the alignpath
		
	return $input_files;
}

sub run_pipeline($$$)
{
	my ($config_file, $files, $methods_list) = @_;
	
	# acquire the outputs for each method
	foreach my $method ( split(',',$methods_list) ) {		
		if ( $method eq 'appris' ) {
			foreach my $met ( @{$METHOD_STRUCT} ) {
				$files->{$method}->{$met} = $outpath.'/'.$met;
			}
		}
		else {
			$files->{$method} = $outpath.'/'.$method;
		}
	}
	# appris is allways executed
	foreach my $met ( @{$METHOD_STRUCT} ) {
		$files->{'appris'}->{$met} = $outpath.'/'.$met;
	}
	
	# execute sequentially
	if ( exists $files->{'firestar'} ) {
		my ($m) = 'firestar';
		eval {
			my ($cmd) = "perl $SRC_DIR/firestar/firestar.pl ".
							"--conf=".$config_file." ".
							"--input=".$files->{'transl'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);
	}
#	if ( exists $files->{'matador3d'} ) {
#		my ($m) = 'matador3d';
#		eval {
#			my ($cmd) = "perl $SRC_DIR/matador3d/matador3d.pl ".
#							"--conf=".$config_file." ".
#							"--gff=".$files->{'annot'}." ".
#							"--input=".$files->{'transl'}." ".
#							"--output=".$files->{$m}." ".
#							"$LOGGER_CONF ";
#			$logger->info("\n** script: $cmd\n");
#			system ($cmd);
#		};
#		$logger->error("runing $m: ".$!) if($@);
#	}
	if ( exists $files->{'matador3d'} ) {
		my ($m) = 'matador3d';
		eval {
			my ($cmd) = "perl $SRC_DIR/matador3d_v2/matador3d.pl ".
							"--conf=".$config_file." ".
							"--input=".$files->{'transl'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);
	}
	if ( exists $files->{'spade'} ) {
		my ($m) = 'spade';
		eval {
			my ($cmd) = "perl $SRC_DIR/spade/spade.pl ".
							"--conf=".$config_file." ".
							"--input=".$files->{'transl'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);
	}
	if ( exists $files->{'corsair'} ) {
		my ($m) = 'corsair';
		eval {
			my ($cmd) = "perl $SRC_DIR/corsair/corsair.pl ".
							"--conf=".$config_file." ".
							"--input=".$files->{'transl'}." ".
							"--gff=".$files->{'pannot'}." ".
							"--output=".$files->{$m}." ".							
							"$LOGGER_CONF ";
							#if ( exists $files->{'pannot'} ) {
							#	$cmd .= "--gff=".$files->{'pannot'}." ";								
							#}							
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);
	}
	if ( exists $files->{'thump'} ) {
		my ($m) = 'thump';
		eval {
			my ($cmd) = "perl $SRC_DIR/thump/thump.pl ".
							"--conf=".$config_file." ".
							"--input=".$files->{'transl'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);
	}
	if ( exists $files->{'crash'} ) {
		my ($m) = 'crash';
		eval {
			my ($cmd) = "perl $SRC_DIR/crash/crash.pl ".
							"--conf=".$config_file." ".
							"--input=".$files->{'transl'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);
	}
	if ( exists $files->{'inertia'} ) {
		my ($m) = 'inertia';
		eval {
			my ($cmd) = "perl $SRC_DIR/inertia/inertia.pl ".
							"--conf=".$config_file." ".
							"--gff=".$files->{'annot'}." ".
							"--inpath=".$files->{'alignpath'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);	
	}
	if ( exists $files->{'proteo'} ) {
		my ($m) = 'proteo';
		eval {
			my ($cmd) = "perl $SRC_DIR/proteo/proteo.pl ".
							"--conf=".$config_file." ".
							"--data=".$files->{'annot'}." ".
							"--output=".$files->{$m}." ".
							"$LOGGER_CONF ";
			$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		$logger->error("runing $m: ".$!) if($@);	
	}
	# appris is allways executed
	my ($m) = 'appris';
	eval {
		my ($cmd) = "perl $SRC_DIR/appris/appris.pl ".
						"--conf=".$config_file." ".
						"--data=".$files->{'annot'}." ".
						"--transcripts=".$files->{'transc'}." ".
						"--translations=".$files->{'transl'}." ".

						"--firestar=".$files->{$m}->{'firestar'}." ".
						"--matador3d=".$files->{$m}->{'matador3d'}." ".
						"--spade=".$files->{$m}->{'spade'}." ".
						"--corsair=".$files->{$m}->{'corsair'}." ".
						"--thump=".$files->{$m}->{'thump'}." ".
						"--crash=".$files->{$m}->{'crash'}." ".
						"--inertia=".$files->{$m}->{'inertia'}." ".
						"--proteo=".$files->{$m}->{'proteo'}." ".

						"--output=".$files->{$m}->{'appris'}." ".
						"--output_nscore=".$files->{$m}->{'appris'}.".nscore ".
						"--output_label=".$files->{$m}->{'appris'}.".label ".
						"$LOGGER_CONF ";
		$logger->info("\n** script: $cmd\n");
		system ($cmd);
	};
	$logger->error("runing $m: ".$!) if($@);
}


main();


__END__

=head1 NAME

appris

=head1 DESCRIPTION

Script that execute APPRIS

=head1 SYNOPSIS

appris

=head2 Required arguments (inputs):

=head3 GENCODE choice: executes appris from gencode data (http://www.gencodegenes.org/data.html)
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--data=  <Gene annotation file>
	
	--transc= <Transcript sequences file>
	
	--transl= <Translation sequences file>
		
=head3 ENSEMBL choice: executes appris from ensembl data ( -api version- )
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--e-version= <Number of Ensembl version of identifier>
	
=head4 Optional arguments (exclusived):
			
	--id= <Ensembl gene identifier>
						
=head3 SEQUENCE choice: executes appris from sequence data (fasta format)
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--transl= <Translation sequences file>

=head4 Optional arguments (exclusived):
	
	--id= <Sequence identifier>
	
=head3 CLUSTER choice: executes appris from inpath

	--inpath= <Acquire input files from PATH>
	
=head2 Optional arguments (methods):
		
  --methods= <List of APPRIS's methods ('firestar,matador3d,spade,corsair,thump,crash,appris'. Default: ALL)>

=head2 Required arguments (outputs):
	
  --outpath= <Output directory>
		
=head2 Optional arguments (type of align):

  --t-align= <Type of alignment (['compara','ucsc'] default: NONE)>

=head2 Optional arguments (config files):

  --cached-path= <Cached path>
  	
=head2 Optional arguments (log arguments):
	
  --loglevel=LEVEL <define log level (default: NONE)>
  
  --logfile=FILE <Log to FILE (default: *STDOUT)>
  
  --logpath=PATH <Write logfile to PATH (default: .)>
  
  --logappend= <Append to logfile (default: truncate)>
    


=head1 EXAMPLE of GENCODE's type of input

=head2 1. Main execution:

perl appris.pl

	--species='Homo sapiens'
	
	--data={ABSOLUTE_PATH_WITH_ID_WITH_ID}/ENSG00000140416/annot.gtf
	
	--data={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/pannot.gtf
	
	--transc={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/transc.fa
	
	--transl={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/transl.fa
	
	--outpath={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/
	
	--logpath={ABSOLUTE_PATH_WITH_ID}/logs/
			
	--logfile=apprisall.log
	
	--loglevel=INFO
	
	--logappend


=head2 2. With given position:

perl appris.pl

	--species='Homo sapiens'
	
	--data={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/annot.gtf
	
	--data={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/pannot.gtf
	
	--transc={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/transc.fa
	
	--transl={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/transl.fa
	
	--position=22,10,M
	
	--outpath={ABSOLUTE_PATH_WITH_ID}/ENSG00000140416/
	
	--logpath={ABSOLUTE_PATH_WITH_ID}/logs/
			
	--logfile=apprisall.log
	
	--loglevel=INFO
	
	--logappend


=head1 EXAMPLE of ENSEMBL's type of 'input

=head2 1. With ensembl gene identifier:

perl appris.pl

	--species='Mus musculus'
	
	--e-version=70
	
	--id=ENSMUSG00000017167.13

	--outpath={ABSOLUTE_PATH_WITH_ID}

	--methods=firestar,matador3d
	
	--t-align=compara
	
	
=head1 EXAMPLE of SEQUENCE's type of input

=head2 1. Main execution:

perl appris.pl

	--species='Homo sapiens'
	
	--id=ENSG00000160404.13

	--transl={ABSOLUTE_PATH_WITH_ID}/ENSG00000160404.13/transl.fa

	--outpath={ABSOLUTE_PATH_WITH_ID}/ENSG00000160404.13/
	
	--logfile={ABSOLUTE_PATH_WITH_ID}/ENSG00000160404.13/log
	
	--loglevel=debug
	
	--logappend
	
	
=head1 EXAMPLE that run appris within CLUSTER

=head2 1. With a protein sequence as input:

perl appris.pl

	--species='Homo sapiens'
	
	--id=WSERVER_ENSG00000160404_XXXX

	--transl={ABSOLUTE_PATH_WITH_ID}/WSERVER_ENSG00000160404_XXXX/transl.fa

	--outpath={ABSOLUTE_PATH_WITH_ID}/WSERVER_ENSG00000160404_XXXX
	
	--cluster-conf={ABSOLUTE_PATH_WITH_ID}/appris/scripts/conf/cluster.ini.wserver
	
	--logfile={ABSOLUTE_PATH_WITH_ID}/WSERVER_ENSG00000160404_XXXX/log
	
	--loglevel=debug
	
	--logappend
	

=head2 1. With a ensembl version:

perl appris.pl

	--species='Homo sapiens'
	
	--id=ENSG00000160404.13
	
	--e-version=74
	
	--outpath={ABSOLUTE_PATH_WITH_ID}/WSERVER_ENSG00000160404_XXXX
	
	--cluster-conf={ABSOLUTE_PATH_WITH_ID}/appris/scripts/conf/cluster.ini.wserver
	
	--logfile={ABSOLUTE_PATH_WITH_ID}/WSERVER_ENSG00000160404_XXXX/ENSG00000160404.13/log
	
	--loglevel=debug
	
	--logappend
	
	
=head1 EXAMPLE of CACHED results

=head2 1. Using cached directory:

perl appris.pl

	--species='Homo sapiens'
	
	--outpath={ABSOLUTE_PATH_WITH_ID}/ENSG00000168556.5/
	
	--transl={ABSOLUTE_PATH_WITH_ID}/ENSG00000168556.5/transl.fa
	
	--transc={ABSOLUTE_PATH_WITH_ID}/ENSG00000168556.5/transc.fa
	
	--pdata={ABSOLUTE_PATH_WITH_ID}/ENSG00000168556.5/pannot.gtf
	
	--data={ABSOLUTE_PATH_WITH_ID}/ENSG00000168556.5/annot.gtf
	
	--methods=firestar
	
	--t-align=compara
	
	--cached-path={ABSOLUTE_PATH_WITH_ID}/appris/code/cached/homo_sapiens/e_70
	
	--logfile={ABSOLUTE_PATH_WITH_ID}/ENSG00000168556.5/log
	
	--loglevel=debug
	
	--logappend
	
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut




=head1 EXAMPLE

perl appris.pl


	
=head1 EXAMPLE

perl appris.pl

	--id=ENSG00000099999	

	--species='Homo sapiens'
	
	--e-version=70
		
	--outpath=examples/ENSG00000099999
	
	--t-align=compara
	
	--loglevel=DEBUG
	
	--logappend
	
	--logpath=examples/ENSG00000099999/
	
	--logfile=log
	
=head1 EXAMPLE

perl appris.pl


	
=head1 EXAMPLE

perl appris.pl

	--species='Homo sapiens'
	
	--e-version=74

	--transl=/local/jmrodriguez/gencode19/annotations/chr4/ENSG00000168556.5/transl.fa
	
	--transc=/local/jmrodriguez/gencode19/annotations/chr4/ENSG00000168556.5/transc.fa
	
	--pdata=/local/jmrodriguez/gencode19/annotations/chr4/ENSG00000168556.5/pannot.gtf
	
	--data=/local/jmrodriguez/gencode19/annotations/chr4/ENSG00000168556.5/annot.gtf
	
	--outpath=/local/jmrodriguez/gencode19/annotations/chr4/ENSG00000168556.5/
	
	--methods=firestar
	
	--t-align=compara
	
	--cached-path=/local/jmrodriguez/appris/code/cached/homo_sapiens/e_70
	
	--loglevel=debug
	
	--logappend
	
	--logpath=/local/jmrodriguez/gencode19/annotations/chr4/ENSG00000168556.5/
	
	--logfile=log
	

=head1 AUTHOR

Created and Developed by

	Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut