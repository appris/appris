#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use appris qw(
	create_gene_list
	create_appris_input
	create_appris_seqinput
	send_email
);
use APPRIS::Utils::File qw(
	prepare_workspace
	printStringIntoFile
	updateStringIntoFile
	getTotalStringFromFile
);
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );
use APPRIS::Utils::Clusters;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_ENSEMBL_DB_FILE
	$CONFIG_INI_EMAIL_FILE
	$LOGLEVEL
	$LOGAPPEND
	
	$NUM_MAX_PROC
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_ENSEMBL_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/ensembldb.ini';
$CONFIG_INI_EMAIL_FILE		= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/email.ini';
$LOGLEVEL					= 'INFO';
$LOGAPPEND					= '';


# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($id) = undef;
my ($position) = undef;
my ($data_file) = undef;
my ($transcripts_file) = undef;
my ($translations_file) = undef;
my ($trifid_pred_file) = undef;
my ($gene_list_file) = undef;
my ($species) = undef;
my ($e_version) = undef;
my ($outpath) = undef;
my ($methods) = undef;
my ($type_of_input) = undef;
my ($type_of_align) = undef;
my ($num_process) = undef;
my ($c_conf_file) = undef;
my ($ensembldb_conf_file) = undef;
my ($cached_path) = undef;
my ($wserver) = undef;
my ($email) = undef;
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
	'transc=s'			=> \$transcripts_file,
	'transl=s'			=> \$translations_file,
	'trifid=s'			=> \$trifid_pred_file,
	'e-version=s'		=> \$e_version,		
	't-align=s'			=> \$type_of_align,
	'methods=s'			=> \$methods,
	'outpath=s'			=> \$outpath,
	'num-process|p=s'	=> \$num_process,
	'cluster-conf=s'	=> \$c_conf_file,	
	'ensembldb-conf=s'	=> \$ensembldb_conf_file,
	'cached-path=s'		=> \$cached_path,
	'wserver=s'			=> \$wserver,
	'email=s'			=> \$email,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments for the different modes of inputs: (the order of conditions is IMPORTANT). Gencode is priority!!
# GENCODE mode
if ( defined $data_file and defined $transcripts_file and defined $translations_file and defined $species and defined $outpath ) {
	$type_of_input = 'datafile';
	$outpath .= '/';
}
elsif ( defined $gene_list_file and defined $data_file and defined $transcripts_file and defined $translations_file and defined $species and defined $outpath ) {
	$type_of_input = 'datafile-list';
	$outpath .= '/';
}
elsif ( defined $position and defined $data_file and defined $transcripts_file and defined $translations_file and defined $species and defined $outpath ) {
	$type_of_input = 'datafile-position';
	$outpath .= '/';
}
# SEQUENCE mode
elsif ( defined $translations_file and defined $species and defined $outpath ) {
	$type_of_input = 'sequence';
	$outpath .= '/';	
}
# ENSEMBL mode
elsif ( defined $id and defined $e_version and defined $species and defined $outpath ) {
	$type_of_input = 'ensembl';
	$outpath .= '/';
}
elsif ( defined $gene_list_file and defined $e_version and defined $species and defined $outpath ) {
	$type_of_input = 'ensembl-list';
	$outpath .= '/';
}
elsif ( defined $position and defined $e_version and defined $species and defined $outpath ) {
	$type_of_input = 'ensembl-position';
	$outpath .= '/';
}
else {
	print `perldoc $0`;
	print "\nCheck required inputs!!\n\n";
	exit 1;
}

# Optional arguments

# get the type of aligns
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

# get num. process
$NUM_MAX_PROC = $ENV{APPRIS_NUM_PARALLEL_PROC};
if ( defined $num_process ) { $NUM_MAX_PROC = $num_process }
my ($NUM_PROC_CHILDS) = 0;

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
$logger->init_log($str_params);
$LOGLEVEL	= $loglevel if ( defined $loglevel );
$LOGAPPEND	= "--logappend" if ( defined $logappend );


#################
# Method bodies #
#################
sub run_datafile($$);
sub run_sequence($$);
sub run_ensembl($);
sub run_pipeline($$;$;$);
sub prepare_appris_seq_params($$);
sub prepare_appris_params($$);
sub create_appris_input($$);
sub run_cluster_appris($$$);
sub run_appris($$$);
sub create_wsrunnerctr($);

# Main subroutine
sub main()
{
	# run appris pipeline for each gene depending on input
	$logger->info("-- from given input...");
	if ( $type_of_input =~ /datafile/ ) {
		$logger->info(" $type_of_input type\n");
		
		# create gene list
		$logger->info("-- create gene list\n");
		my ($gene_list, $gene_data) = appris::create_gene_list(
															-type		=> $type_of_input,
															-gdata		=> $data_file,
															-gtransc	=> $transcripts_file,
															-gtransl	=> $translations_file,
															-list		=> $gene_list_file,
															-pos		=> $position,
		);
		
		$logger->info("-- run data files\n");
		my ($runtimes) = run_datafile($gene_data, $gene_list);
						
	}
	elsif ( $type_of_input =~ /ensembl/ ) {
		$logger->info(" $type_of_input type\n");	

		# create gene list
		$logger->info("-- create gene list\n");
		my ($gene_list) = appris::create_gene_list(
											-type	=> $type_of_input,		
											-econf	=> $ensembldb_conf_file,
											-id		=> $id,
											-ever	=> $e_version,
											-spe	=> $species,
											-list	=> $gene_list_file,
											-pos	=> $position,
		);
		
		$logger->info("-- run ensembl ids\n");
		my ($runtimes) = run_ensembl($gene_list);
				
	}
	elsif ( $type_of_input =~ /sequence/ ) {
		$logger->info(" $type_of_input type\n");
		
		# create gene list
		$logger->info("-- create gene list\n");
		my ($gene_list,$gene_data) = appris::create_gene_list(
											-type		=> $type_of_input,
											-gtransl	=> $translations_file,
		);
				
		$logger->info("-- run sequence\n");
		my ($runtimes) = run_sequence($translations_file, $gene_data);
				
	}
	else {
		$logger->error("analying input parameters");
	}	
	
	# TODO: CHECK RESULTS!!
	# At least, check if the number of result files is correct.
	
	# send email
	if ( defined $email ) {
		$logger->info("-- send email\n");
		my ($methods) = $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};		
		appris::send_email($CONFIG_INI_EMAIL_FILE, $email, $wserver, $species, $methods, $outpath);		
	}
	
	$logger->finish_log();
	
	exit 0;	
}

sub run_datafile($$)
{
	my ($data, $gene_list) = @_;
	my ($runtimes) = undef;
		
	foreach my $gene (@{$data}) {
		my ($gene_id) = $gene->stable_id;
		my ($gene_name) = $gene_id;
		if ( $gene->external_name ) {
			$gene_name = $gene->external_name;
		}
		if ( defined $gene_list ) { # if there is a gene list, we run appris for the list
			if ( exists $gene_list->{$gene_id} or exists $gene_list->{$gene_name} ) {
				my ($runtime) = run_pipeline($gene_id, $outpath, $gene);
				push(@{$runtimes}, $runtime);				
			}
		}
		else { # if there is not gene list, we run appris for all of them
			my ($runtime) = run_pipeline($gene_id, $outpath, $gene);
			push(@{$runtimes}, $runtime);			
		}
	}
	return $runtimes;
} # end run_datafile

sub run_ensembl($)
{
	my ($gene_list) = @_;
	my ($runtimes) = undef;
			
	foreach my $gene_id (keys(%{$gene_list})) {
		my ($runtime) = run_pipeline($gene_id, $outpath);
		push(@{$runtimes}, $runtime);		
	}
	return $runtimes;
} # end run_ensembl

#sub run_sequence($)
#{
#	my ($transl_file) = @_;
#	my ($runtimes) = undef;
#	
#	my ($runtime) = run_pipeline('', $outpath, undef, $transl_file);
#	push(@{$runtimes}, $runtime);
#				
#	return $runtimes;
#} # end run_sequence
sub run_sequence($$)
{
	my ($transl_file, $data) = @_;
	my ($runtimes) = undef;
	
	if ( defined $data ) {
		while ( my ($gene_id,$gene) = each(%{$data}) ) {
			my ($runtime) = run_pipeline($gene_id, $outpath, $gene);
			push(@{$runtimes}, $runtime);
		}		
	}
	else {
		my ($runtime) = run_pipeline('', $outpath, undef, $transl_file);
		push(@{$runtimes}, $runtime);		
	}
				
	return $runtimes;
} # end run_sequence

sub run_pipeline($$;$;$)
{
	my ($gene_id, $workspace, $gene, $transl_file) = @_;	
	my ($runtime) = undef;
	
	$logger->info("\t>> $gene_id\n");	
	
	# create parameters
	$logger->info("\t-- create parameters ");
	my ($params) = {
		'species'		=> $species
	};
	# datafile mode
	if ( $type_of_input =~ /datafile/ ) {
		$logger->info("from $type_of_input\n");

		$logger->info("\t-- prepare workspace\n");
		my ($g_workspace) = $workspace;
		if ( $type_of_input eq 'datafile-position' ) {
			my ($chr) = $gene->chromosome;
			$g_workspace .= $chr.'/'.$gene_id;
		}
		else {
			$g_workspace .= $gene_id;
		}
		$g_workspace = prepare_workspace($g_workspace);
		$params->{'outpath'} = $g_workspace;
				
		$logger->info("\t-- prepare params\n");
		$params = prepare_appris_params($gene, $g_workspace);
		$params->{'e-version'} = $e_version if ( defined $e_version );
	}
	# ensembl mode
	elsif ( $type_of_input =~ /ensembl/ ) {
		$logger->info("from $type_of_input\n");
		if ( defined $trifid_pred_file ) {
			$logger->error("TRIFID not supported for input of type '${type_of_input}'\n");
		}

		$logger->info("\t-- prepare workspace\n");
		my ($g_workspace) = $workspace;
		$g_workspace .= $gene_id;
		$g_workspace = prepare_workspace($g_workspace);
		$params->{'outpath'} = $g_workspace;
		
		$logger->info("\t-- prepare params\n");
		$params->{'id'} = $gene_id;
		$params->{'e-version'} = $e_version;
	}	
	# sequence mode
	elsif ( $type_of_input =~ /sequence/ ) {
		$logger->info("from $type_of_input\n");

		$logger->info("\t-- prepare workspace\n");
		my ($g_workspace) = $workspace;
		if ( $gene_id ne '' ) {
			$g_workspace .= $gene_id;
#			# split id to avoid a directory with a lot of subdirectories (UniProt case)
#			my ($len) = length($gene_id);
#			my ($half) = int($len / 2);
#			my ($g1) = substr($gene_id, 0, $half);
#			my ($g2) = substr($gene_id, $half, $len);
#			my ($split_id) = "$g1/$g2";
#			$g_workspace .= $split_id;
		}
		$g_workspace = prepare_workspace($g_workspace);
		$params->{'outpath'} = $g_workspace;
				
		$logger->info("\t-- prepare params\n");
		$params->{'transl'} = $transl_file;
		# UniProt/neXtProt cases
		if ( $gene_id ne '' ) {
			$params = prepare_appris_seq_params($gene, $g_workspace);
			$params->{'id'} = $gene_id;
		}
	}
	else {
		$logger->info("\t-- do not run appris\n");
		return $runtime;
	}	
	# data for aligns
	if ( defined $type_of_align ) {
		if ( lc($type_of_align) eq 'ucsc' ) {
			$params->{'t-align'}	= $type_of_align;
		}
		elsif ( (lc($type_of_align) eq 'compara') and defined $e_version ) { # compara needs ensembl version
			$params->{'t-align'}	= $type_of_align;
			$params->{'e-version'}	= $e_version;
		}	
	}
	# the rest of parameters
	$params->{'methods'} = $methods if ( defined $methods );	
	$params->{'cached-path'} = $cached_path if ( defined $cached_path );	
	$params->{'cluster-conf'} = $c_conf_file if ( defined $c_conf_file );

	# run pipe (only if methods are defined)
	if ( defined $methods ) {
		if ( defined $c_conf_file ) {
			if ( defined $trifid_pred_file ) {
				$logger->error("TRIFID not supported for cluster pipeline\n");
			}

			$logger->info("\t-- run cluster pipe\n");
			my ($runtime) = run_cluster_appris($gene_id, $workspace, $params);
		}
		else {
			$logger->info("\t-- run local pipe\n");
			my ($runtime) = run_appris($gene_id, $workspace, $params);
		}		
	}

	return $runtime;
	
} # end run_pipeline

sub prepare_appris_seq_params($$)
{
	my ($gene, $workspace) = @_;
	my ($gene_id) = $gene->{'id'};
	
	# require parameter
	my ($params) = {
		'species'		=> $species,
		'outpath'		=> $workspace,
	};

	$params->{'trifid'} = $trifid_pred_file if ( defined $trifid_pred_file );
	
	# input files
	$logger->info("\t-- create datafile input\n");
	my ($in_files) = {
			'transl'		=> $workspace.'/'.'transl.fa',
	};		
	my ($create) = appris::create_appris_seqinput($gene, $in_files);
	if ( defined $create ) {
		if ( exists $in_files->{'transl'} and -e ($in_files->{'transl'}) ) {
			$params->{'transl'} = $in_files->{'transl'};
		}			
	}

	return $params;
	
} # end prepare_appris_seq_params

sub prepare_appris_params($$)
{
	my ($gene, $workspace) = @_;
	my ($gene_id) = $gene->stable_id;
	
	# require parameter
	my ($params) = {
		'species'		=> $species,
		'outpath'		=> $workspace,
	};

	$params->{'trifid'} = $trifid_pred_file if ( defined $trifid_pred_file );
	
	# input files
	my ($in_files) = {
			'data'			=> $workspace.'/'.'annot.gtf',
			'pdata'			=> $workspace.'/'.'pannot.gtf',
			'transc'		=> $workspace.'/'.'transc.fa',
			'transl'		=> $workspace.'/'.'transl.fa',
			'cdsseq'		=> $workspace.'/'.'cdsseq.fa',
			'logfile'		=> $workspace.'/'.'log',
			'errfile'		=> $workspace.'/'.'err',
	};
	
	# delete any log file for the new execution
	if ( -e $in_files->{'logfile'} and -e $in_files->{'errfile'} ) {
		eval {
			my ($cmd) = "rm -rf ".$in_files->{'logfile'}." ".$in_files->{'errfile'}." ";
			#$logger->info("\n** script: $cmd\n");
			system ($cmd);
		};
		throw("deleting log files of appris") if($@);
	}
		
	# If methods var is defined (we run appris), then 
	# check if files exists	
	if (
		defined $methods and
		-e $in_files->{'data'} and (-s $in_files->{'data'} > 0) and
		-e $in_files->{'pdata'} and (-s $in_files->{'pdata'} > 0) and
		-e $in_files->{'transc'} and (-s $in_files->{'transc'} > 0) and
		-e $in_files->{'transl'} and (-s $in_files->{'transl'} > 0) 
	) {
		$logger->info("\t-- use input path\n");
		$params->{'data'} = $in_files->{'data'};
		$params->{'pdata'} = $in_files->{'pdata'};
		$params->{'transc'} = $in_files->{'transc'};
		$params->{'transl'} = $in_files->{'transl'};
	}
	else {
		$logger->info("\t-- create datafile input\n");
		#my ($create) = create_appris_input($gene, $in_files);
		my ($create) = appris::create_appris_input($gene, $in_files);
		if ( defined $create ) {
			if ( exists $in_files->{'data'} and -e ($in_files->{'data'}) ) {
				$params->{'data'} = $in_files->{'data'};	
			}			
			if ( exists $in_files->{'pdata'} and -e ($in_files->{'pdata'}) ) {
				$params->{'pdata'} = $in_files->{'pdata'};
			}
			if ( exists $in_files->{'transc'} and -e ($in_files->{'transc'}) ) {
				$params->{'transc'} = $in_files->{'transc'};
			}
			if ( exists $in_files->{'transl'} and -e ($in_files->{'transl'}) ) {
				$params->{'transl'} = $in_files->{'transl'};
			}			
		}
	}	

	return $params;
	
} # end prepare_appris_params

#sub create_appris_input($$)
#{
#	my ($gene, $in_files) = @_;
#	my ($create) = undef;
#	
#	# gene vars
#	my ($chr) = $gene->chromosome;
#	my ($gene_id) = $gene->stable_id;
#	my ($data_cont) = '';
#	my ($pdata_cont) = '';
#	my ($transc_cont) = '';
#	my ($transl_cont) = '';
#	my ($cdsseq_cont) = '';
#	my ($datafile_file);
#	
#	# get gene annots
#	my ($cmd) = "grep -E 'gene_id \"$gene_id\[\.\\d+]{0,1}|GeneID:$gene_id' $data_file";
#	my (@global_data_cont);
#	eval { @global_data_cont = `$cmd`; };
#	throw("getting gene annots\n") if($@);
#	if ( scalar(@global_data_cont) == 0 ) {
#    	throw("empty gene annots\n");			
#	}		
#
#	# scan transctipt/translation/cds seq/cds info
#	if ( $gene->transcripts ) {
#		foreach my $transcript (@{$gene->transcripts}) {		
#			my ($transcript_id) = $transcript->stable_id;
#			my ($transcript_eid) = $transcript_id;
#			#if ( $transcript->version ) {
#			#	$transcript_eid = $transcript_id.'.'.$transcript->version;
#			#}
#			my ($transcript_name) = $transcript_id;
#			if ( $transcript->external_name ) {
#				$transcript_name = $transcript->external_name;
#			}
#					
#			# get transcript sequence				
#			if ( $transcript->sequence ) {
#				my ($seq) = $transcript->sequence;
#				my ($len) = length($transcript->sequence);
#				$transc_cont .= ">$transcript_eid|$gene_id|$transcript_name|$len\n";
#				$transc_cont .= $seq."\n";
#			}
#			
#			if ( $transcript->translate ) {
#				my ($translate) = $transcript->translate;
#							
#				# get translation sequence
#				if ( $translate->sequence ) {
#					my ($seq) = $translate->sequence;
#					my ($len) = length($translate->sequence);
#					my ($translate_id) = $transcript_eid;
#					$translate_id = $translate->protein_id if ( defined $translate->protein_id );
#					# mask short sequences
#					#if ( $len <= 2 ) { $seq .= 'X'; }
#					$transl_cont .= ">$transcript_eid|$translate_id|$gene_id|$transcript_name|$len\n";
#					$transl_cont .= $seq."\n";					
#				}
#									
#				if ( $transcript->exons and $translate->cds_sequence ) {
#					my ($exons) = $transcript->exons;
#					
#					for (my $icds = 0; $icds < scalar(@{$translate->cds_sequence}); $icds++) {
#						my ($cds) = $translate->cds->[$icds];
#						my ($exon) = $exons->[$icds];
#						my ($exon_id) = $exon->stable_id; $exon_id = '-' unless (defined $exon_id);
#						my ($pro_cds) = $translate->cds_sequence->[$icds];
#	
#						my ($cds_start) = $cds->start;
#						my ($cds_end) = $cds->end;
#						my ($cds_strand) = $cds->strand;
#						my ($cds_phase) = $cds->phase;
#						
#						my ($pro_cds_start) = $pro_cds->start;
#						my ($pro_cds_end) = $pro_cds->end;
#						my ($pro_cds_end_phase) = $pro_cds->end_phase;
#						my ($pro_cds_seq) = $pro_cds->sequence;
#							
#						# delete the residue that is shared by two CDS						
#						#if (defined $pro_cds_end_phase and $pro_cds_end_phase != 0) { $pro_cds_seq =~ s/\w{1}$//; }
#							
#						# retrieve cds sequence
#						if (defined $pro_cds_seq and $pro_cds_seq ne '') {
#							my ($len) = length($pro_cds_seq);
#							$cdsseq_cont .= ">$transcript_eid|$gene_id|$transcript_name|$len|$exon_id|$chr|$cds_start|$cds_end|$cds_strand|$cds_phase\n";
#							$cdsseq_cont .= $pro_cds_seq."\n";							
#						}
#													
#						# acquire protein coordinates
#						if (defined $pro_cds_start and defined $pro_cds_end) {
#							$pdata_cont .=	'SEQ'."\t".
#														'Ensembl'."\t".
#														'Protein'."\t".
#														$pro_cds_start."\t".
#														$pro_cds_end."\t".
#														'.'."\t".
#														'.'."\t".
#														$pro_cds_end_phase."\t".
#														"ID=$exon_id;Parent=$transcript_eid;Gene=$gene_id;Note=cds_coord>$cds_start-$cds_end:$cds_strand\n";					
#						}
#					}
#				}
#			}			
#		}
#	}
#
#	# create files
#	if ( scalar(@global_data_cont) > 0 ) {
#		$data_cont = join('',@global_data_cont);
#		my ($output_file) = $in_files->{'data'};		
#		my ($printing_file_log) = printStringIntoFile($data_cont, $output_file);
#		throw("creating $output_file file") unless ( defined $printing_file_log );
#	}
#	if ( $pdata_cont ne '' ) {
#		my ($output_file) = $in_files->{'pdata'};
#		my ($printing_file_log) = printStringIntoFile($pdata_cont, $output_file);
#		throw("creating $output_file file") unless ( defined $printing_file_log );
#	}
#	if ( $transc_cont ne '' ) {
#		my ($output_file) = $in_files->{'transc'};
#		my ($printing_file_log) = printStringIntoFile($transc_cont, $output_file);
#		throw("creating $output_file file") unless ( defined $printing_file_log );
#	}
#	if ( $transl_cont ne '' ) {
#		my ($output_file) = $in_files->{'transl'};
#		my ($printing_file_log) = printStringIntoFile($transl_cont, $output_file);
#		throw("creating $output_file file") unless ( defined $printing_file_log );
#	}
#	if ( $cdsseq_cont ne '' ) {
#		my ($output_file) = $in_files->{'cdsseq'};
#		my ($printing_file_log) = printStringIntoFile($cdsseq_cont, $output_file);
#		throw("creating $output_file file") unless ( defined $printing_file_log );
#	}
#	
#	# determine if appris has to run
#	if ( (scalar(@global_data_cont) > 0) and ($pdata_cont ne '') and ($transl_cont ne '') and ($cdsseq_cont ne '') ) {
#		$create = 1;
#	}
#	
#	return $create;
#	
#} # end create_appris_input

sub run_cluster_appris($$$)
{
	my ($gene_id, $workspace, $params) = @_;
	
	# prepare clusters	
	$logger->info("-- prepare clusters\n");
	my ($clusters) = new APPRIS::Utils::Clusters( -conf => $c_conf_file );
	$logger->error("preparing clusters") unless (defined $clusters);
	
	# get local vars
	my ($l_host) = $clusters->host;
	my ($l_user) = $clusters->user;
	
	# get free cluster
	$logger->info("-- get free cluster\n");
	my ($c_free) = $clusters->free();
	while ( !$c_free ) {
		$logger->info(".");
		sleep(5);
		$c_free = $clusters->free();
	}
	$logger->info("\n");
	$logger->info("-- free cluster: $c_free\n");
	
	# if the server which calls to cluster is equal or not
	# create workspace in cluster depending on gene_id, webserver jobid or input data (by default is gene_id)	
	my ($c_id) = $gene_id;
	if ( defined $wserver ) { $c_id = $wserver }	
	my ($c_jobhome);
	my ($c_wspace);
	my ($c_parameters) = '';
	my ($c_cmd2) = '';		
	unless ( $l_host eq $c_free )
	{
		$logger->info("-- calling outside of cluster\n");
		
		# create workspace in cluster depending on c_id
		$logger->info("-- create workspace within cluster\n");
		my ($c_wsbase) = $clusters->wspace($c_free);
		$c_wspace = $c_wsbase.'/'.$c_id;
		$c_wspace = $clusters->srmdir($c_free, $c_wspace);		
		#$c_wspace = $clusters->smkdir($c_free, $c_wspace, 1);
		sleep(1);
		$c_wspace = $clusters->wspace($c_free, $c_wspace);
		$logger->error("preparing cluster workspace") unless (defined $c_wspace);
		sleep(1);
		
		# create inputs within cluster
		$logger->info("-- create inputs within cluster\n");
		#my ($orig_files) = $workspace; #$orig_files =~ s/\/*$//; $orig_files .= '/*';
		#my ($c_file) = $clusters->dscopy($c_free, $c_wspace, $orig_files);
		my ($c_file) = $clusters->dscopy($c_free, $c_wsbase, $workspace);
		$logger->error("copying $c_wsbase into cluster") unless (defined $c_file);
		$c_file =~ s/\/\*$//;
		sleep(2);
		
		# obtain appris home of code
		$c_jobhome = $clusters->cluster($c_free)->j_home;
		
		# create parameters for job script
		$logger->info("-- create job script for cluster\n");	
		if ( defined $params ) {
			while ( my ($k,$v) = each(%{$params}) ) {
				if ( $k eq 'cluster-conf' ) {
					# do nothing => delete parameters for running on a cluster
				}
				# modify the outpaths because in a cluster is different
				elsif ( ($k eq 'data') or ($k eq 'transc') or ($k eq 'transl') or ($k eq 'outpath') ) {
					my ($v2) = $v;
					$v2 =~ s/$ENV{APPRIS_WORKSPACE}/$c_wsbase/;
					$c_parameters .= " --$k='$v2' ";
				}
				else {
					$c_parameters .= " --$k='$v' ";
				}
			}
			# DEPRECATED: we always use a path as input when execute in a cluster!!! Except we want to run ENSEMBL choice
			#unless ( defined $e_version ) {
			#	$c_parameters .= " --inpath=$c_wspace ";	
			#}
		}
		
		# copy results from cluster to localhost
		$c_cmd2	.=	"cd $c_wspace && tar -cf - * | ssh $l_user\@$l_host 'cd $workspace; tar -xf - ' \n";
		
		# clean cluster workspace, after execution
		#c_cmd2	.= "rm -rf $c_wspace \n";
	}
	else # localhost is the cluster host
	{
		$logger->info("-- running within the cluster\n");
		
		# obtain appris home of code
		$c_jobhome = "$ENV{APPRIS_HOME}";
		
		# get inputs
		if ( defined $params ) {
			while ( my ($k,$v) = each(%{$params}) ) {
				if ( $k eq 'cluster-conf' ) {
					# do nothing => delete parameters for running on a cluster
				}
				else {
					$c_parameters .= " --$k='$v' ";
				}
			}
		}
		
		# create workspace in cluster
		$c_wspace = $workspace.'/'.$c_id;
	}
	
	# user environment
	my ($c_env) = '';
	if ( $clusters->cluster($c_free)->u_env ) {
		my ($u_env) = $clusters->cluster($c_free)->u_env;
		$c_env .= "source $u_env \n";
	}

	# job environment
	$c_env .= "source $c_jobhome/code/conf/apprisrc \n";
	if ( defined $wserver ) {
		$c_env .= "source $c_jobhome/code/conf/apprisrc.WS \n";	
	}
	if ( defined $ENV{APPRIS_WS_NAME} ) {
		$c_env .= "export APPRIS_WS_NAME=$ENV{APPRIS_WS_NAME} \n";
	}
	if ( defined $e_version ) {
		$c_env .= "module load ensembl$e_version \n";	
	}
		
	# create main script
	my ($c_logpath) = $c_wspace;
	my ($c_logfile) = 'log';
	my ($c_cmd) = "\n".
				"perl $c_jobhome/code/appris.pl ".
				" $c_parameters ".
				" --loglevel='$LOGLEVEL' --logpath='$c_logpath' --logfile='$c_logfile' $LOGAPPEND \n";
	my ($c_stderr) = $c_id.'.err';
				
	# create script for cluster
	my ($c_script) = 	$c_env."\n".
						$c_cmd."\n".
						$c_cmd2."\n";
	my ($script) = $clusters->script($c_free, { 
												'wdir'		=> $c_wspace,
												'stdout'	=> $c_stderr,
												'stderr'	=> $c_stderr,
												'script'	=> $c_script,
	});
	$logger->error("creating job script for cluster") unless (defined $script);
	$logger->debug("\n** script:\n $script\n");	
	
	# create local script that will execute into cluster
	my ($l_script) = $c_wspace.'/'.$c_id.'.sh';
	$l_script = printStringIntoFile($script, $l_script);
	$logger->error("creating job script into local server") unless (defined $l_script);
	sleep(1);

	# submit job
	$logger->info("-- submit job: ");
	my ($job_id) = $clusters->submit($c_free, $l_script);
	$logger->error("can not submit the job for $c_id: $!\n") unless (defined $job_id);
	
	if ( defined $wserver ) {
		my ($pid) = $$;
		$c_wspace =~ s/\/*$//g;
		$c_logpath =~ s/\/*$//g;		
		my ($wsrunnerctr) = {
			'file'		=> $c_wspace.'/wsrunner.ctr',
			'host'		=> $c_free,
			'pid'		=> $pid,
			'jobid'		=> $job_id,
			'wspace'	=> $c_wspace,
			'script'	=> $c_cmd,
			'params'	=> $params,
			'logfile'	=> $c_logpath.'/'.$c_logfile,
		};
		create_wsrunnerctr($wsrunnerctr);		
	}
		
} # end run_cluster_appris

sub run_appris($$$)
{
	my ($id, $workspace, $params) = @_;
	
	# get inputs
	my ($parameters) = '';
	if ( defined $params ) {
		while ( my ($k,$v) = each(%{$params}) ) {
			$parameters .= " --$k='$v' ";
		}
	}
	
	# create workspace for cluster depending on gene_id, webserver jobid or input data (by default is gene_id)
	# create command line and FORK
	my ($c_id) = $id; if ( defined $wserver and $type_of_input =~ /sequence/ ) { $c_id = $wserver }
	my ($c_wspace) = $workspace.'/'.$id;
	my ($c_logpath) = $c_wspace;
#	# split id to avoid a directory with a lot of subdirectories (UniProt case)	
#	if ( $type_of_input =~ /sequence/ ) {
#		my ($len) = length($id);
#		my ($half) = int($len / 2);
#		my ($g1) = substr($id, 0, $half);
#		my ($g2) = substr($id, $half, $len);
#		my ($split_id) = "$g1/$g2";
#		$c_logpath = $workspace.'/'.$split_id;
#	}	
	my ($c_logfile) = 'log';
	my ($cmd) =	" perl $ENV{APPRIS_CODE_DIR}/appris.pl ".
				" $parameters ".
				" --loglevel='$LOGLEVEL' --logpath='$c_logpath' --logfile='$c_logfile' $LOGAPPEND ";
	my ($pid) = fork();	
	if ( defined $wserver ) {
		$c_wspace =~ s/\/*$//g;
		$c_logpath =~ s/\/*$//g;
		my ($wsrunnerctr) = {
			'file'		=> $workspace.'/wsrunner.ctr',
			'host'		=> 'localhost',
			'pid'		=> $pid,
			'jobid'		=> $c_id,
			'wspace'	=> $c_wspace,
			'script'	=> $cmd,
			'params'	=> $params,
			'logfile'	=> $c_logpath.'/'.$c_logfile,
		};
		create_wsrunnerctr($wsrunnerctr);
	}
		
	# submit job
	if ( !defined $pid ) {
	    die "Cannot fork: $!";
	}
	if ( $pid ) {
		$NUM_PROC_CHILDS++;
	}	
	elsif ( $pid == 0 ) {
	    # client process
		$logger->info("\n** script($$): $cmd\n");
	    close STDOUT; close STDERR; # so parent can go on
		eval {
			system($cmd);
		};
		exit 1 if($@);
	    exit 0;
	}
	
	# wait until at least one process finish
	while ( $NUM_PROC_CHILDS >= $NUM_MAX_PROC ) {
		my $pid = wait();
		$NUM_PROC_CHILDS--;
		$logger->info("\n** script($pid) exits\n\n");
	}
				
} # end run_appris

sub create_wsrunnerctr($)
{
	my ($ctr) = @_;
	
	# control file
	my ($ctrcont) = '';
	my ($file) = $ctr->{'file'};
	
	# control parameters
	$ctrcont = "RUNNER_HOST:".$ctr->{'host'}."\n";
	$file = updateStringIntoFile($ctrcont, $file);
	
	$ctrcont = "RUNNER_PID:".$ctr->{'pid'}."\n";
	$file = updateStringIntoFile($ctrcont, $file);

	$ctrcont = "RUNNER_JOBID:".$ctr->{'jobid'}."\n";
	$file = updateStringIntoFile($ctrcont, $file);

	$ctrcont = "RUNNER_WSPACE:".$ctr->{'wspace'}."\n";
	$file = updateStringIntoFile($ctrcont, $file);

	$ctrcont = "RUNNER_SCRIPT:".$ctr->{'script'}."\n";
	$file = updateStringIntoFile($ctrcont, $file);
	
	my (@outputs);
	if ( exists $ctr->{'params'}->{'data'} and exists $ctr->{'params'}->{'pdata'} and exists $ctr->{'params'}->{'transc'} and exists $ctr->{'params'}->{'transl'} ) {
		push(@outputs, $ctr->{'params'}->{'data'});
		push(@outputs, $ctr->{'params'}->{'pdata'});
		push(@outputs, $ctr->{'params'}->{'transc'});
		push(@outputs, $ctr->{'params'}->{'transl'});
	}
	elsif ( exists $ctr->{'params'}->{'transl'} ) {
		push(@outputs, $ctr->{'params'}->{'transl'});
	}
	else {
		my ($ctrlog) = 'Gene without Coding Region (CDS)';
		$ctrcont = "RUNNER_STATUS:".'ERROR'."\n";
		$file = updateStringIntoFile($ctrcont, $file);
		$ctrcont = "RUNNER_LOG:".$ctrlog."\n";
		$file = updateStringIntoFile($ctrcont, $file);
		
	}
		
	if ( exists $ctr->{'params'}->{'methods'} ) {
		my (@meths) = split(',', $ctr->{'params'}->{'methods'});
		foreach my $met (@meths) {
			push(@outputs, $ctr->{'wspace'}.'/'.$met);
			if ( $met eq 'appris' ) {
				push(@outputs, $ctr->{'wspace'}.'/'.$met);
				push(@outputs, $ctr->{'wspace'}.'/'.$met.'.label');
				push(@outputs, $ctr->{'wspace'}.'/'.$met.'.nscore');
			}
		}		
	}
	if ( exists $ctr->{'logfile'} ) {
		push(@outputs, $ctr->{'logfile'});
	}	
	my ($outfiles) = join(',',@outputs);

	$ctrcont = "RUNNER_OUTFILES:".$outfiles."\n";
	$file = updateStringIntoFile($ctrcont, $file);
	
} # end print_wsrunner_control

main();


1;

__END__

=pod

=head1 NAME

run_appris

=head1 DESCRIPTION

global script that runs APPRIS for GENCODE and ENSEMBL 

=head1 SYNOPSIS

run_appris

=head2 Required arguments (inputs):

=head3 DATAFILE choice: executes appris from datafile data (http://www.gencodegenes.org/data.html)
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--data=  <Gene annotation file>
	
	--transc= <Transcript sequences file>
	
	--transl= <Translation sequences file>

=head4 Optional arguments (exclusived):
	
	--id= <Gene identifier>
						
	--gene-list= <File with a list of genes>
			
	--position= <Genome position (eg1: 21. eg2: 21,22)>
	
=head3 ENSEMBL choice: executes appris from ensembl data ( -api version- )
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--e-version= <Number of Ensembl version of identifier>
	
=head4 Optional arguments (exclusived):
			
	--id= <Ensembl gene identifier / External gene name>
						
	--gene-list= <File with a list of genes>
			
	--position= <Genome position (eg1: 21. eg2: 21,22)>
			
=head3 SEQUENCE choice: executes appris from sequence data (fasta format)
	
=head4 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--transl= <Translation sequences file>

=head4 Optional arguments (exclusived):
	
	--id= <Sequence identifier>
						
=head2 Optional arguments (methods):
		
  --methods= <List of APPRIS's methods ('firestar,matador3d,spade,corsair,corsair_alt,thump,crash,appris')>
  Note: If it is not defined, we only create input data files.
  
  --trifid= <Trifid predictions file>

=head2 Required arguments (outputs):
	
  --outpath= <Output directory>
		
=head2 Optional arguments (type of align):

  --t-align= <Type of alignment with version (['compara74','ucsc'] default: NONE)>
	
=head2 Optional arguments (parallel):

  -p/--num-process= <Number of process to run in parallel>

  --cluster-conf= <Config file of cluster execution (default: none)>
  
=head2 Optional arguments (config files):

  --ensembldb-conf= <Config file of Ensembl database (default: 'conf/ensembldb.ini' file)>
  
  --apprisdb-conf= <Config file of APPRIS database (default: 'conf/apprisdb.ini' file)>
  
  --cached-path= <Cached path>
  
  --wserver= <JobId of web service execution>
  
  --email= <Send results by email>  

=head2 Optional arguments (log arguments):
	
  --loglevel=LEVEL <define log level (default: NONE)>
  
  --logfile=FILE <Log to FILE (default: *STDOUT)>
  
  --logpath=PATH <Write logfile to PATH (default: .)>
  
  --logappend= <Append to logfile (default: truncate)>


=head1 EXAMPLE of GENCODE's type of input

=head2 1. Main execution:

run_appris

	--species='Homo sapiens'
	
	--data={ABSOLUTE_PATH}/gencode.v15.annotation.gtf
	
	--transc={ABSOLUTE_PATH}/gencode.v15.pc_transcripts.fa
	
	--transl={ABSOLUTE_PATH}/gencode.v15.pc_translations.fa
	
	--methods=appris
	
	--outpath={ABSOLUTE_PATH}/annotations/
	
	--logpath={ABSOLUTE_PATH}/logs/
			
	--logfile=run_appris.log
	
	--loglevel=INFO
	
	--logappend


=head2 2. With given position:

run_appris

	--species='Homo sapiens'
	
	--data={ABSOLUTE_PATH}/gencode.v15.annotation.gtf
	
	--transc={ABSOLUTE_PATH}/gencode.v15.pc_transcripts.fa
	
	--transl={ABSOLUTE_PATH}/gencode.v15.pc_translations.fa
	
	--position=22,10,M

	--outpath={ABSOLUTE_PATH}/annotations/
	
	--logpath={ABSOLUTE_PATH}/logs/
			
	--logfile=run_appris.log
	
	--loglevel=INFO
	
	--logappend
	

=head2 3. With given gene list:

run_appris

	--species='Homo sapiens'
	
	--data={ABSOLUTE_PATH}/gencode.v15.annotation.gtf
	
	--transc={ABSOLUTE_PATH}/gencode.v15.pc_transcripts.fa
	
	--transl={ABSOLUTE_PATH}/gencode.v15.pc_translations.fa
	
	--gene-list={ABSOLUTE_PATH}/gene_list.txt

	--outpath={ABSOLUTE_PATH}/annotations/
	
	--logpath={ABSOLUTE_PATH}/logs/
			
	--logfile=run_appris.log
	
	--loglevel=INFO
	
	--logappend


=head1 EXAMPLE of ENSEMBL's type of 'input

=head2 1. With ensembl gene identifier:

run_appris

	--species='Mus musculus'
	
	--e-version=70
	
	--id=ENSMUSG00000017167.13

	--outpath={ABSOLUTE_PATH}

	--methods=firestar,matador3d
	
	--t-align=compara70
	
=head2 2. With external gene name:

run_appris

	--species='Mus musculus'
	
	--e-version=70
	
	--id=RNF215

	--outpath={ABSOLUTE_PATH}

	--methods=firestar,matador3d
	
		
=head1 EXAMPLE of SEQUENCE's type of input

=head2 1. Main execution:

run_appris

	--species='Homo sapiens'
	
	--id=ENSG00000160404.13

	--transl={ABSOLUTE_PATH}/ENSG00000160404.13/transl.fa

	--outpath={ABSOLUTE_PATH}/ENSG00000160404.13/
	
	--logfile={ABSOLUTE_PATH}/ENSG00000160404.13/log
	
	--loglevel=debug
	
	--logappend
	
	
=head1 EXAMPLE that run appris within CLUSTER

=head2 1. With a protein sequence as input:

run_appris

	--species='Homo sapiens'
	
	--id=WSERVER_ENSG00000160404_XXXX

	--transl={ABSOLUTE_PATH}/WSERVER_ENSG00000160404_XXXX/transl.fa

	--outpath={ABSOLUTE_PATH}/WSERVER_ENSG00000160404_XXXX
	
	--cluster-conf={ABSOLUTE_PATH}/appris/scripts/conf/ws/cluster.ini
	
	--logfile={ABSOLUTE_PATH}/WSERVER_ENSG00000160404_XXXX/log
	
	--loglevel=debug
	
	--logappend
	

=head2 1. With a ensembl version:

run_appris

	--species='Homo sapiens'
	
	--id=ENSG00000160404.13

	--e-version=74

	--outpath={ABSOLUTE_PATH}/WSERVER_ENSG00000160404_XXXX
	
	--cluster-conf={ABSOLUTE_PATH}/appris/scripts/conf/scripts/cluster.ini
	
	--logfile={ABSOLUTE_PATH}/WSERVER_ENSG00000160404_XXXX/ENSG00000160404.13/log
	
	--loglevel=debug
	
	--logappend
	
	
=head1 EXAMPLE of CACHED results

=head2 1. Using cached directory:

run_appris

	--species='Homo sapiens'
	
	--outpath={ABSOLUTE_PATH}/ENSG00000168556.5/
	
	--transl={ABSOLUTE_PATH}/ENSG00000168556.5/transl.fa
	
	--transc={ABSOLUTE_PATH}/ENSG00000168556.5/transc.fa
	
	--pdata={ABSOLUTE_PATH}/ENSG00000168556.5/pannot.gtf
	
	--data={ABSOLUTE_PATH}/ENSG00000168556.5/annot.gtf
	
	--methods=firestar
	
	--t-align=compara74
	
	--cached-path={ABSOLUTE_PATH}/appris/code/cached/homo_sapiens/e_70
	
	--logfile={ABSOLUTE_PATH}/ENSG00000168556.5/log
	
	--loglevel=debug
	
	--logappend
	
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
