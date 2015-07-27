#!/usr/bin/perl -W
# _________________________________________________________________
# $Id$
# $Revision$
# Developed by:
#		Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es-
# _________________________________________________________________

use strict;
use FindBin;
use Getopt::Long;
use Config::IniFiles;
use Statistics::Descriptive;
use POSIX;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile open_dir );
use lib "$FindBin::Bin/lib";
use Exon;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$BIN_DIR
	$SRC_DIR
	$LOGGER_CONF
	
	$WSPACE_BASE
	$WSPACE_CACHE
	$WSPACE_DATA
	$SPECIES
	$RUN_PROGRAM
	$PROG_IN_ALIGNS
	$PROG_SLR_TYPES
	$PROG_PVALUE
	$PROG_OMEGA_CUTOFF
	
	$IN_ALIGN_SUFFIX
	$IN_TREE_SUFFIX
	
	$UNKNOWN_LABEL
	$NO_LABEL
	
);

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($data_file) = undef;
my ($inpath) = undef;
my ($output_file) = undef;
my ($appris) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'conf=s'				=> \$config_file,
	'gff=s'					=> \$data_file,	
	'inpath=s'				=> \$inpath,
	'output=s'				=> \$output_file,
	'appris'				=> \$appris,	
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);
	
# Required arguments
unless (
	defined  $config_file and
	defined  $data_file and
	defined  $inpath and
	defined  $output_file
){
    print `perldoc $0`;
    exit 1;
}

# Get conf vars
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$BIN_DIR			= $LOCAL_PWD.'/bin/';
$SRC_DIR			= $LOCAL_PWD.'/../';
$LOGGER_CONF		= '';

$WSPACE_BASE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('INERTIA_VARS', 'name').'/';
$WSPACE_CACHE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('CACHE_VARS', 'name').'/';
$WSPACE_DATA		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('INPUT_VARS', 'name').'/';
$SPECIES			= $cfg->val( 'APPRIS_PIPELINE', 'species');
$RUN_PROGRAM		= $cfg->val( 'INERTIA_VARS', 'program');
$PROG_IN_ALIGNS		= [ split( ',', $cfg->val('INERTIA_VARS', 'aligns') ) ];
$PROG_SLR_TYPES		= [ split( ',', $cfg->val('INERTIA_VARS', 'stypes') ) ];
$PROG_PVALUE		= $cfg->val('INERTIA_VARS', 'pvalue');
$PROG_OMEGA_CUTOFF	= $cfg->val('INERTIA_VARS', 'omega');

$IN_ALIGN_SUFFIX	= '.faa';
$IN_TREE_SUFFIX		= '.nh';

$UNKNOWN_LABEL		= 'UNKNOWN';
$NO_LABEL			= 'NO';


my ($gawk_bin) = "/usr/bin/awk";
my ($cl) = "";
my ($rl) = "";


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
sub create_aligns($$);
sub inertia($$$);
sub omega($$$$$$$);
sub get_exon_annotation($$$);
sub print_omega_report($$$$$$);
sub _realign($$$$$);
sub _get_id_version($);
sub _get_gene_list_from_gff($$);
sub _filter_transcripts_from_gff($$);
sub _get_transcript_list_from_gff($$);
sub _get_appris_annotations($$$$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{	
	# Create aligns from the main aligns
	$logger->info("-- re-aligning for each align file\n");
	my ($align_files);
	foreach my $align ( @{$PROG_IN_ALIGNS} ) {
		$align_files = create_aligns($data_file, $align);		
	}
			
	# Execute slr
	$logger->info("-- execute slr for each align of every transc\n");
	while ( my ($transc_id, $t_report) = each(%{$align_files}) ) {
		while ( my ($align, $a_report) = each(%{$t_report}) ) {
			my ($align_file) = $a_report->{'align'};
			my ($tree_file) = $a_report->{'tree'};
			my ($local_align_file) = "$align_file\.slr";
			unless ( -e $local_align_file and (-s $local_align_file > 0) ) {
				if ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {
					# first try with the oficial version
					eval {
						my ($cmd) = "time $RUN_PROGRAM ".
										"-seqfile $align_file ".
										"-treefile $tree_file ".
										"-outfile $local_align_file ".
										" 1>&2 2> /dev/null ";
						$logger->info("** script: $cmd\n");
						system ($cmd);
					};
					if ( $? != 0 or !(-e $local_align_file) or (-s $local_align_file == 0) ) {
						# second try with the version 1.2
						eval {
							my ($cmd) = "time Slr_1.2 ".
											"-seqfile $align_file ".
											"-treefile $tree_file ".
											"-outfile $local_align_file ".
											" 1>&2 2> /dev/null ";
							$logger->info("** script: $cmd\n");
							system ($cmd);
						};
						$logger->error("runing $RUN_PROGRAM: ".$!) if($@);						
					}
				}			
			}
			else { $logger->info("\t** cached\n"); }		
		}
	}

	# Parse slr results
	$logger->info("-- parsing slr results\n");
	my ($inertia_out_file) = inertia($data_file, $WSPACE_DATA, $WSPACE_BASE);

	# Getting result
	if ( defined $inertia_out_file ) {
		eval {
			my ($cmd) = "cp -rp $inertia_out_file $output_file";
			$logger->info("-- $cmd\n");
			system ($cmd);		
		};
		$logger->error("copying results: ".$!) if($@);
	}

	$logger->finish_log();
	
	exit 0;	
}

sub create_aligns($$)
{
	my ($data_file, $inalign) = @_;	
	my ($align_files);
	my ($inalign_suffix) = '.'.$inalign.$IN_ALIGN_SUFFIX;
	my ($intree_suffix) = '.'.$inalign.$IN_TREE_SUFFIX;	
	
	# Get the main align files
	$logger->info("-- -- getting aligns from $inalign\n");
	my ($main_align_files) = open_dir($WSPACE_DATA, $inalign_suffix.'$');
	
	foreach my $main_align_filename (@{$main_align_files}) {
		my ($main_align_file) = $WSPACE_DATA.'/'.$main_align_filename;
		if ( ($main_align_filename =~ /^(.*)$inalign_suffix/) and -e $main_align_file and (-s $main_align_file > 0) ) {
			my ($id) = $1;
			
			# check if transcript is protein coding (and cds start-end found)
			if ( _filter_transcripts_from_gff($data_file,$id) ) {
				my ($main_tree_file) = $WSPACE_DATA.'/'.$id.$intree_suffix;
				
				# if not exits the tree and the alignment is ucsc, the copy the main phy tree
				if ( -e $main_tree_file and (-s $main_tree_file > 0) ) {
					# create aligns for each type
					foreach my $stype ( @{$PROG_SLR_TYPES} ) {
						$align_files->{$id}->{$stype} = _realign($stype, $main_align_file, $main_tree_file, $id, $SPECIES);					
					}
				}
				else {
					$logger->error("no phylogenetic tree");
				}
			}
			else {
				$logger->warning("$id is not protein-coding or cds start-end dont found");
			}
		}
	}
	return $align_files;
}

sub _realign($$$$$)
{
	my ($type, $main_align_file, $main_tree_file, $id, $species) = @_;
	my ($main_specie) = lc($species); $main_specie =~ s/^\s*//; $main_specie =~ s/\s*$//; $main_specie =~ s/\s/\_/;
	my ($align_file) = $WSPACE_BASE.'/'.$id.'.'.$type;
	my ($tree_file) = $WSPACE_BASE.'/'.$id.'.'.$type.'.tree';
	my ($align_files) = {
				'align'		=> $align_file,
				'tree'		=> $tree_file, 		
	};
	if ( $type eq 'compara' )
	{
		# prepare the align for slr
		$logger->info("-- -- preparing the alignment and tree for slr\n");		
		unless ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {
			eval {
				my ($cmd) = "python $BIN_DIR/pslr.py $main_specie $main_align_file $main_tree_file $align_file $tree_file 1>&2 2> /dev/null";				
				$logger->info("\t** script: $cmd\n");
				system($cmd);
			};
			$logger->error("filtering alignment and tree") if $@;			
		}
		else { $logger->info("\t** cached\n"); }
	}
	elsif ( $type eq 'ucsc' )
	{
		# prepare the align for slr
		$logger->info("-- -- preparing the alignment and tree for slr\n");		
		unless ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {
			eval {
				my ($cmd) = "python $BIN_DIR/pslr.py $main_specie $main_align_file $main_tree_file $align_file $tree_file 1>&2 2> /dev/null";				
				$logger->info("\t** script: $cmd\n");
				system($cmd);
			};
			$logger->error("filtering alignment and tree") if $@;			
		}
		else { $logger->info("\t** cached\n"); }
	}
	elsif ( $type =~ /prank/ )
	{
		my ($local_align_file) = "$align_file\.best\.fas";
		
		# re-align the main alignment
		$logger->info("-- -- re-align the main alignment using $type\n");
		# IMPORTANT: THIS COMMENT IS TEMPORAL
		#unless ( -e $local_align_file and (-s $local_align_file > 0) ) {
		unless ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {		
		# IMPORTANT: THIS COMMENT IS TEMPORAL
			eval {
				my ($cmd) = "prank -prunetree -d=$main_align_file -o=$align_file 1>&2 2> /dev/null";
				$logger->info("\t** script: $cmd\n");
				system($cmd);
			};
			$logger->error("filtering alignment and tree") if $@;
		}
		else { $logger->info("\t** cached\n"); }

		# prepare the align for slr
		$logger->info("-- -- preparing the alignment and tree for slr\n");
		unless ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {		
			eval {
				my ($cmd) = "python $BIN_DIR/pslr.py $main_specie $local_align_file $main_tree_file $align_file $tree_file 1>&2 2> /dev/null";				
				$logger->info("\t** script: $cmd\n");
				system($cmd);
			};
			$logger->error("filtering alignment and tree") if $@;
		}
		else { $logger->info("\t** cached\n"); }
	}
	elsif ( $type =~ /kalign/ )
	{
		my ($local_align_file) = "$align_file\.fas";
		
		# re-align the main alignment
		# IMPORTANT: THIS COMMENT IS TEMPORAL
		#unless ( -e $local_align_file and (-s $local_align_file > 0) ) {
		unless ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {
		# IMPORTANT: THIS COMMENT IS TEMPORAL
			$logger->info("-- -- re-align the main alignment using $type\n");
			eval {
				my ($cmd) = "kalign -f fasta -c input -d wu -i $main_align_file -o $local_align_file 1>&2 2> /dev/null";
				$logger->info("\t** script: $cmd\n");
				system($cmd);
			};
			$logger->error("filtering alignment and tree") if $@;
		}
		else { $logger->info("\t** cached\n"); }

		# prepare the align for slr
		$logger->info("-- -- preparing the alignment and tree for slr\n");
		unless ( -e $align_file and (-s $align_file > 0) and -e $tree_file and (-s $tree_file > 0) ) {
			eval {
				my ($cmd) = "python $BIN_DIR/pslr.py $main_specie $local_align_file $main_tree_file $align_file $tree_file 1>&2 2> /dev/null";				
				$logger->info("\t** script: $cmd\n");
				system($cmd);
			};
			$logger->error("filtering alignment and tree") if $@;
		}
		else { $logger->info("\t** cached\n"); }
	}
	return $align_files;
}

sub inertia($$$)
{
	my ($gff_file, $data_dir, $slr_dir) = @_;
	my ($inertia_out_file);
	
	# Get gene list
	$logger->info("-- -- -- getting genes from data\n");
	my ($all_g) = _get_gene_list_from_gff($gff_file, $data_dir);
	#$logger->debug(Dumper($all_g)."\n");

	foreach my $gene (@{$all_g})
	{
		# Get transcript list
		$logger->info("-- -- -- getting transcripts from data\n");		
		my ($all_ts) = _get_transcript_list_from_gff($gff_file, $gene);
		#$logger->debug(Dumper($all_ts)."\n");
				
		if(defined $all_ts and scalar(@{$all_ts})>0)
		{
			my($omega_out_file_list);
			my($omega_ts_list);
			
			# Check if every transcript has got slr result (compara,prank,kalign)
			$logger->info("-- -- -- getting slr for each transcript\n");			
			foreach my $type (@{$PROG_SLR_TYPES})
			{
				my($omega_type_out_file)=$slr_dir.'/'.$gene.'.'.$type.'.omega';
				$omega_out_file_list->{$type}=$omega_type_out_file;
						
				my($ts_list);
				my($slr_file_list);
				foreach my $ts (@{$all_ts})
				{

					my($ts_slr_file)=$slr_dir.'/'."$ts.$type.slr";
					if(-e $ts_slr_file and (-s $ts_slr_file > 0))
					{
						#$logger->debug("$ts_slr_file\n");
						push(@{$ts_list}, $ts);
						$slr_file_list->{$ts}=$ts_slr_file;
						unless (exists $omega_ts_list->{$ts}) { $omega_ts_list->{$ts} = undef };
					}
				}
				# get the annotations (Omega) for one type of alignment
				$logger->info("-- -- -- acquire omega values for given transcript\n");
				unless ( -e $omega_type_out_file and (-s $omega_type_out_file > 0) ) {
					if(defined $ts_list and scalar(@{$ts_list}) and defined $slr_file_list) {
						omega($appris, $gene,$ts_list,$gff_file,$slr_file_list,$type,$omega_type_out_file);
					}
				}
				else { $logger->info("\t** cached\n"); }
			}
						
			# Get the annotations for the main isoform /* APPRIS */ ----------------
			if ( defined $appris )
			{				
				# get the annotations (of INERTIA) the whole types of SLR alignments
				$logger->info("-- -- -- acquire inertia values\n");
				$inertia_out_file=$slr_dir.'/'.$gene.'.inertia';
				_get_appris_annotations($gene,$gff_file,$omega_ts_list,$omega_out_file_list,$inertia_out_file);
			}
		}
	}
	return $inertia_out_file;
}

# Get the id and the version from Ensembl identifiers of gencode7
sub _get_id_version($)
{
	my ($i_id) = @_;
	my ($id, $version) = (undef,undef);
	
	if ( $i_id =~ /^([^\.]*)\.(\d)*$/ ) {
		($id, $version) = ($1, $2);
	}
	return ($id, $version);
		
} # End _get_id_version

# Get the gene list from input files of alignments
sub _get_gene_list_from_gff($$)
{
	my ($gff_file, $dir) = @_;

	my ($all_t);	
	my ($aux_all_g);
	my ($all_g);
	my (@all_g_content);
	
	# get transcript that has alignment file
	unless ( opendir(DIR, $dir) ) {
		$logger->error("can't opendir $dir: $!");
		return undef;	
	}
	my (@files) = grep { !/^\./ && /\.faa$/ } readdir(DIR);
	closedir DIR;
	
	# get unique list of transcripts
	foreach my $file (@files) {
		if ( $file =~ /^([^\.]*)\./ ) {
			$all_t->{$1} = undef;
		}
	}
	
	# seek the genes for every transcript
	my ($transcript_codition) = '';	
	foreach my $t ( keys(%{$all_t}) ) {
		$transcript_codition .= ' $12 ~ /'.$t.'/ ||'; # for rel7 version		
	}	
    if ($transcript_codition ne '') {
    	$transcript_codition =~ s/\|\|$//;
    	$transcript_codition = '('.$transcript_codition.')';
		my (@all_g_content) = `awk '{if( \$3=="transcript" && $transcript_codition ){print \$10}}' $gff_file`;
		if ( scalar(@all_g_content) == 0 ) {
    		$logger->error("gff has not transcript information");			
		}
		else {
		   	foreach my $g_content (@all_g_content)
		   	{
		   		$g_content =~ s/\;$//;
		   		$g_content =~ s/\"//g;
		   		$g_content =~ s/\s*//g;
				$aux_all_g->{$g_content} = undef;
		   	}
		}		
    }
    else {
    	$logger->error("gff has not transcript information");
    }
    
    # save the gene within array
    foreach my $g ( keys(%{$aux_all_g}) ) {
    	my ($id, $version) = _get_id_version($g); # for rel7 version
    	push(@{$all_g}, $id); # for rel7 version
    	#push(@{$all_g}, $g);
	}

	return $all_g;	
}

# Get transcript list (from a gene)
sub _get_transcript_list_from_gff($$)
{
	my ($gff_file, $gene) = @_;
	
	my ($all_ts);
   	my (@all_ts_content) = `grep $gene $gff_file | awk '{if(\$3==\"transcript\"){print \$12}}'`;
   	foreach my $t_content (@all_ts_content)
	{	
   		$t_content =~ s/\;$//;
   		$t_content =~ s/\"//g;
   		$t_content =~ s/\s*//g;
   		my ($id, $version) = _get_id_version($t_content); # for rel7 version
		push(@{$all_ts}, $id); # for rel7 version
		#push(@{$all_ts}, $t_content);
	}
	return $all_ts;
}

# Say if a transcript is protein-coding (has cds), cds start-end have been founding and the length of protein
sub _filter_transcripts_from_gff($$)
{
	my ($gff_file, $trans_id) = @_;
	
	my ($count);	
   	my (@all_ts_content) = `grep $trans_id $gff_file | awk '{if(\$3==\"transcript\" || \$3==\"CDS\" || \$3==\"start_codon\" || \$3==\"stop_codon\"){print \$0}}'`;
   	foreach my $t_content (@all_ts_content)
	{
		my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$attributes) = split("\t", $t_content);
		if (defined $trans_id and ($type eq 'CDS') ) { 
			$count->{'cds'} = 1;
			unless ( exists $count->{'cds_length'} ) {
				$count->{'cds_length'} = 0;
			}
			$count->{'cds_length'} = abs($end - $start) + $count->{'cds_length'};
		}
		elsif (defined $trans_id and ($type eq 'start_codon') ) {
			$count->{'start_codon'} = 1;
		}
		elsif (defined $trans_id and ($type eq 'stop_codon') ) {
			$count->{'stop_codon'} = 1;
		}
	}
	my ($proceed) = 0;
	if ( exists $count->{'cds'} and exists $count->{'start_codon'} and exists $count->{'stop_codon'} and ($count->{'cds_length'} <= 27000) ) {
		$proceed = 1;
	}
	return $proceed;
}

sub omega($$$$$$$)
{
	my ($appris, $gene, $ts_list, $gff_file, $slr_file_list, $type, $output_file)=@_;
	
	my($strand) = "";
	my($i) = 0;
	my($j) = 0;
	my($k) = 0;
	my($l) = 0;
	my($z) = 0;
	my(@gff) = ();
	my(@pos) = ();
	my(@codon) = ();
	my(@fields) = ();
	my(@omega) = ();
	my($total) = 0;
	my($tot_exon) = 0;
	my($pos_exon) = 0;
	my($tot_pos) = 0;
	my(@exons) = ();
	my(@ts_exons) = ();

    my($gene_stat) = Statistics::Descriptive::Full->new(); #deprecated, i think :)
	 	
	#get strand
    $cl = "grep $gene $gff_file | head -1 | cut -f7 | ";
    open(ST,$cl) or die "$!";
    $strand = <ST>;
    close(ST) or die "$!";
    chomp($strand);
	
    # get all exons of gene and  sorted
    if($strand eq "-")
	{
		$cl = "grep $gene $gff_file | $gawk_bin \'\$3==\"CDS\"\' | cut -f4,5 | sort -ru | ";
	}	
    elsif($strand eq "+")
    {
		$cl = "grep $gene $gff_file | $gawk_bin \'\$3==\"CDS\"\' | cut -f4,5 | sort -u | ";
    }
    #$cl = "grep $gene $gff_file | $gawk_bin \'\$3==\"CDS\"\' | cut -f4,5 | sort -un | ";
    open(EX,$cl) or die "$!";
    while($rl = <EX>)
    {
		my($ex) = Exon::new(); #instatiate exon object
		chomp($rl);
		my(@fields) = split(/\s+/,$rl);
		$ex->{'start'} = $fields[0];
		$ex->{'end'} = $fields[1];
		$ex->{'strand'} = $strand;
		$ex->{'om_average'} = 0; #will hold average omega, initiate as 0
		$ex->{'std'} = 0;# will hold stdandard deviation of omega, initiate as 0
		$ex->{'d_value'}=0; 
		$ex->{'p_value'}=0;

		push(@exons,$ex); #build array of exons for all exons of gene
	}


	##############################################################################
	# for each transcript
	##############################################################################
	for($l = 0; $l < scalar(@{$ts_list}); $l++)
	{ 
		my($slr_file) = $slr_file_list->{$ts_list->[$l]}; #slr input file

		#get gff coords
		#$cl = "grep $ts_list->[$l] $gff_file | $gawk_bin \'\$3==\"CDS\"\' | sort > $WSPACE_BASE/$ts_list->[$l].gff";
		$cl = "grep $ts_list->[$l] $gff_file | $gawk_bin \'\$3==\"CDS\"\' | sort -nk 4,5 > $WSPACE_BASE/$ts_list->[$l].gff";
		system($cl);
		#field 4 and 5 of gff file are the exon coordinates
		if($strand eq "+")
		{
		    $cl = "cut -f4,5 $WSPACE_BASE/$ts_list->[$l].gff | sort -u | "; 
		}
		elsif($strand eq "-")
		{
		    $cl = "cut -f4,5 $WSPACE_BASE/$ts_list->[$l].gff | sort -ru | "; 
		}
		#$cl = "cut -f4,5 $WSPACE_BASE/$ts_list->[$l].gff | sort -un | ";
	
		#reverse sort for negative strand
		#TODO maybe other gff files are sorted differently, keep this in mind
		@ts_exons = (); # array with just the coordinates of each exon of the transcript
	    open(TE,$cl) or die "$!";
		while($rl = <TE>)
		{
		    chomp($rl);
		    push(@ts_exons,[split(/\s+/,$rl)]);
		}
		for($j = 0; $j <= $#ts_exons; $j++)
		{
		    for($k = 0; $k <= $#exons; $k++)
		    {
				if(($ts_exons[$j][0] == $exons[$k]->start()) && ($ts_exons[$j][1] == $exons[$k]->end()))
				{
					$ts_exons[$j][2] = $k; #find out the number of the exon in the gene by comparing coordinates
					push(@{$exons[$k]->{'ts'}},$ts_list->[$l]); #add transcript id to array
				}
	    	}
		}	

		#put coordinates rel. to transcript
		if($strand eq "-")
		{
		    $cl ="perl ".$BIN_DIR.'/'."get_locus_mrna.pl $WSPACE_BASE/$ts_list->[$l].gff > $WSPACE_BASE/$ts_list->[$l].mrna.gff";
		    system($cl);
		    $cl ="perl ".$BIN_DIR.'/'."adj-gff.pl $WSPACE_BASE/$ts_list->[$l].gff $WSPACE_BASE/$ts_list->[$l].mrna.gff > $WSPACE_BASE/$ts_list->[$l].plus.gff";
		    system($cl);
		    $cl = "perl ".$BIN_DIR.'/'."adj-mrna.pl $WSPACE_BASE/$ts_list->[$l].plus.gff > $WSPACE_BASE/$ts_list->[$l].relts.gff";
		    system($cl);
		}
		else
		{
			#strand == "+"
		    $cl = "perl ".$BIN_DIR.'/'."adj-mrna.pl $WSPACE_BASE/$ts_list->[$l].gff > $WSPACE_BASE/$ts_list->[$l].relts.gff";
		    system($cl);
		}		

		@gff = ();
		open(GFF,"$WSPACE_BASE/$ts_list->[$l].relts.gff") or die "$!";
		while($rl = <GFF>)
		{
		    chomp($rl);
		    push(@gff,[split(/\t/,$rl)]);
		}
		close(GFF) or die "$!";

		@pos = ();	
		push(@pos,[1,floor($gff[0][4]/3)]);#put coordinates relative to peptide, 3 codons = 1 aminoacid position
	
		for($j = 1; $j <= $#gff; $j++)
		{
			push(@pos,[$pos[$j-1][1] + 1,floor($gff[$j][4]/3)]);
		}
		$pos[$#pos][1] -= 1;#last codon is stopcodon and slr doesn't consider this (probably obvious ;)

#		$cl = "grep $ts_list->[$l] $p_gff_file > $WSPACE_BASE/$ts_list->[$l].relts.gff";
#		system($cl);
#
#		@pos = ();
#		$cl = "cut -f4,5 $WSPACE_BASE/$ts_list->[$l].relts.gff | ";		
#	    open(TE,$cl) or die "$!";
#		while($rl = <TE>)
#		{
#		    chomp($rl);
#		    push(@pos,[split(/\s+/,$rl)]);
#		}

	
		@omega = ();
		open(SLR,$slr_file) or die "$!";
		<SLR>;#first line contains field names etc
		while($rl = <SLR>)
		{
			$rl =~ s/^\s*//; #starts with whitespace
			@fields = split(/\s+/,$rl);
			push(@omega,$fields[3]);
		}
		close(SLR) or die "$!";
		
		for($j = 0; $j <= $#ts_exons; $j++)
		{
			#exon has not been seen before
			if($exons[$ts_exons[$j][2]]->om_average() == 0)
			{
				my($exon_stat) = Statistics::Descriptive::Full->new(); #statistics object for each exon
				#exon coordinates relative to peptide
				for($k = $pos[$j][0]; $k <= $pos[$j][1]; $k++)
				{
					$exon_stat->add_data($omega[$k-1]);#add omega values for each position of the exon to statistics object
				    $gene_stat->add_data($omega[$k-1]); #add omega values for each position of the exon to gene statistics object
				    push(@{$exons[$ts_exons[$j][2]]->{'omega'}},$omega[$k-1]); #add omega values to "omega" array of exon object
			    }
				$exons[$ts_exons[$j][2]]->{'om_average'} = $exon_stat->mean() if(defined $exon_stat->mean()); #average omega 
				$exons[$ts_exons[$j][2]]->{'std'} = $exon_stat->standard_deviation() if(defined $exon_stat->standard_deviation());
		    }
		}
	} # End loop for each trancript
	
	
		
	##############################################################################
	# for each exon of gene	
	##############################################################################
	for($j = 0; $j <= $#exons; $j++)
	{
		my(@ex_w) = ();
		#array to take omega values of exon with average omega > 0.25 (this cut-off was used by Goldman group for ENCODE paper)
		my(@ts_w) = ();
		#array to take omega values of all remaining of exons
		@ex_w = @{$exons[$j]->omega()};
		for($k = 0; $k <= $#exons; $k++)
		{
		    if($k == $j)
		    {
				next; #omit exon that has average omega > 0.25, see above if statement
			}
			push(@ts_w,@{$exons[$k]->omega()});
		}

		my($ks_table_file)="$WSPACE_BASE/ks.input";
		open(OU,">$ks_table_file") or die "$!"; #write data for R input, values of both arrays, tab separated
		if($#ex_w < $#ts_w)
		{
		    for($z=0;$z<=$#ex_w;$z++)
		    {
				print OU $ex_w[$z],"\t",$ts_w[$z],"\n";
		    }
		    for($z=scalar(@ex_w);$z<=$#ts_w;$z++)
		    {
				print OU "NA","\t",$ts_w[$z],"\n"; #have to put "NA" if one array is bigger than the other one, otherwise R doesn't like it
		    }
		}
		else
		{
		    for($z=0;$z<=$#ts_w;$z++)
		    {
				print OU $ex_w[$z],"\t",$ts_w[$z],"\n";
		    }
		    for($z=scalar(@ts_w);$z<=$#ex_w;$z++)
		    {
				print OU $ex_w[$z],"\t","NA","\n";
		    }
		}

		close(OU) or die "$!";
		# If run inertia in CATON cluster
		#$cl = "source /home/inb/.bashrc && source /etc/profile && source /etc/bash.bashrc && R --no-save --args $ks_table_file < ".$BIN_DIR.'/'."ks.R | "; # R command line for Kolmogorov-Smirnov test, see ks.R script		
		$cl = "R --no-save --args $ks_table_file < ".$BIN_DIR.'/'."ks.R | "; # R command line for Kolmogorov-Smirnov test, see ks.R script
		my(@R) = ();
		open(R,$cl) or die "$!"; #read R input
		@R=<R>;
		close(R);
		chomp(@R);
		for($z = 0; $z<=$#R;$z++)
		{
			#add KS test statistic and p-value to exon object (D = 0.2693, p-value = 0.001961)
		    if($R[$z] =~ m/^D = ([^\,]+), p-value = ([^\n]+)/)
		    {
		    	my($d_val)=$1;
		    	my($p_val)=$2;
				$exons[$j]->{'d_value'}=$d_val if(defined $d_val); 
				$exons[$j]->{'p_value'}=$p_val if(defined $p_val);
		    }
		}
	} # End loop for each exon of gene


	##############################################################################
	# print results
	##############################################################################

	print_omega_report($appris, $type, $gene, $ts_list, \@exons, $output_file);
	
}

# Get annotation per exon
sub get_exon_annotation($$$)
{
	my ($exon_report, $ts_list, $ts_exon_list) = @_;

	my ($exon_annot)=$UNKNOWN_LABEL;
	
	# this cut-offa were used by Goldman group for ENCODE paper
	if ( $exon_report->{'om_average'} < $PROG_OMEGA_CUTOFF ) { 
		$exon_annot=$UNKNOWN_LABEL;
	}
	else {
		if ( $exon_report->{'p_value'} < $PROG_PVALUE ) {
			$exon_annot=$NO_LABEL;
		}
		else {
			$exon_annot=$UNKNOWN_LABEL;
		}
	}
	# but if one exon is in more than the half of transcrips then the exon can not be rejected.
	if (defined $ts_list and (scalar(@{$ts_list}) > 0 ) and
		defined $ts_exon_list and (scalar(@{$ts_exon_list}) > 0) ) {
		if ( (scalar(@{$ts_exon_list}) / scalar(@{$ts_list})) > 0.5 ) {
			$exon_annot = $UNKNOWN_LABEL;
		}
	}	
	return $exon_annot;
}

# Get printed report of omega
sub print_omega_report($$$$$$)
{
	my ($appris, $type, $gene, $ts_list, $exons, $output_file) = @_;
	
	my($output_content)='';
	my($transcript_report);
	my(@sorted_exons)=reverse sort { $a->om_average <=> $b->om_average } @{$exons};

	# Print by sorted omega mean
	$output_content .=	"# omega_average\t".
						"omega_exon_id\t".
						"start_exon\t".
						"end_exon\t".
						"strand_exon\t".
						"difference_value\t".
						"p_value\t".
						"st_desviation\t".
						"exon_annotation\t".
						"transcript_list\n";

	for(my $j = 0; $j < scalar(@{$exons}); $j++)	
	{
		my ($ts_exon_list) = $sorted_exons[$j]->{'ts'};
		
		# get the exon annotation from slr scores
		my ($omega_exon_annot) = get_exon_annotation($sorted_exons[$j], $ts_list, $ts_exon_list);

		my($transcript_list)='NULL';
		if(scalar(@{$sorted_exons[$j]->{'ts'}})>0)
		{
			$transcript_list=join(";",@{$sorted_exons[$j]->{'ts'}});
			
			# we keep the annotatios per transcript.
			foreach my $ts (@{$sorted_exons[$j]->{'ts'}}) {
				
				# we keep exon annotations
				my ($e_index) = $sorted_exons[$j]->{'start'};
				$transcript_report->{$ts}->{'exons'}->{$e_index} = {
					'annot'			=> $omega_exon_annot,
					'start'			=> $sorted_exons[$j]->{'start'},
					'end'			=> $sorted_exons[$j]->{'end'},
					'strand'		=> $sorted_exons[$j]->{'strand'},
					'om_average'	=> $sorted_exons[$j]->{'om_average'},
					'p_value'		=> $sorted_exons[$j]->{'p_value'}					
				};
				
				# if one exon was rejected then the transcript is rejected
				if (exists $transcript_report->{$ts}->{'annot'} and defined $transcript_report->{$ts}->{'annot'}) {
					if ($transcript_report->{$ts}->{'annot'} ne $NO_LABEL) { # we keep the global value of NO					
						$transcript_report->{$ts}->{'annot'} = $omega_exon_annot;
					}					
				}
				else {
					$transcript_report->{$ts}->{'annot'} = $omega_exon_annot;
				}
			}
		}

		my($omega_report_content)=
								$sorted_exons[$j]->{'om_average'}."\t".
								($j+1)."\t".
								$sorted_exons[$j]->{'start'}."\t".
								$sorted_exons[$j]->{'end'}."\t".
								$sorted_exons[$j]->{'strand'}."\t".
								$sorted_exons[$j]->{'d_value'}."\t".
								$sorted_exons[$j]->{'p_value'}."\t".
								$sorted_exons[$j]->{'std'}."\t".
								$omega_exon_annot."\t".
								$transcript_list."\n";

		$output_content .= $omega_report_content;
	}
	
	# Get the annotations for the main isoform /* APPRIS */ ----------------
	if ( defined $appris )
	{	
		# print annotations of transcripts and theirs sorted exons
		if (defined $transcript_report) {		
			$output_content .= "----------------------------------------------------------------------\n";
			while (my ($ts,$ts_report) = each(%{$transcript_report})) {
				$output_content .= ">".$ts."\t".$ts_report->{'annot'}."\n";
								
				if (exists $ts_report->{'exons'} and defined $ts_report->{'exons'}) {
					my @sort_e_index = sort {$a <=> $b} keys %{$ts_report->{'exons'}};
					 foreach my $sort_i (@sort_e_index) {
					 	$output_content .= "\t".
					 						$ts_report->{'exons'}->{$sort_i}->{'start'}.":".
					 						$ts_report->{'exons'}->{$sort_i}->{'end'}.":".
					 						$ts_report->{'exons'}->{$sort_i}->{'strand'}."\t".
					 						$ts_report->{'exons'}->{$sort_i}->{'annot'}."\t".
					 						$ts_report->{'exons'}->{$sort_i}->{'om_average'}."\t".
					 						$ts_report->{'exons'}->{$sort_i}->{'p_value'}."\n";
					 }					
				}
			}		
		}
		else {
			# if there is not results. we don't known
			foreach my $ts (@{$ts_list}) { 
				$output_content .= ">".$ts."\t".$UNKNOWN_LABEL."\n"; 	
			}				
		}
	}
	
	if (defined $output_content and $output_content ne '') {
		local(*OUTPUT_FILE);
		open(OUTPUT_FILE, ">$output_file");
		print OUTPUT_FILE $output_content;
		close(OUTPUT_FILE);		
	}	
}

sub _get_appris_annotations($$$$$)
{
	my ($gene, $gff_file, $omega_ts_list, $omega_out_file_list, $inertia_out_file) = @_;
	
	my($output_content)='';	
	my ($omega_trans_report);
	my ($inertia_trans_report);
	
	# get the annotations per alignment
	while( my ($type,$omega_file) = each(%{$omega_out_file_list}) )
	{
		if( -e $omega_file and (-s $omega_file > 0) )
		{
			my ($ts_id);
			my ($ts_annot);
			open(OMEGA_FILE, $omega_file);
			while (my $line = <OMEGA_FILE>) {
				if ($line =~ /\>([^\t+]*)\t+([^\n]*)\n/) {
					$ts_id = $1;
					$ts_annot = $2;
					$omega_trans_report->{$ts_id}->{$type}->{'annot'} = $ts_annot;				
				}
				elsif ( $line =~ /^\t([^\:]*)\:([^\:]*)\:([^\t]*)\t+([^\t]*)\t+([^\t]*)\t+([^\n]*)\n+$/ ) {
					my ($e_index) = $1;
					$omega_trans_report->{$ts_id}->{$type}->{'exons'}->{$e_index} = {
						'start'			=> $1,
						'end'			=> $2,
						'strand'		=> $3,
						'annot'			=> $4,
						'om_average'	=> $5,
						'p_value'		=> $6,
					};
				}
				
			}
			close(OMEGA_FILE);			
		}
	}
		
	# Pick the exon and transcript info of Inertia from each Omega results
	foreach my $ts (keys (%{$omega_ts_list}) )
	{
		my ($ts_inertia_annot) = $NO_LABEL;

		# Flag that control the number of results from Omega (filter,prank,kalign SLR results)
		my ($num_omega_results) = 0;

		# IMPORTANT: THIS COMMENT IS TEMPORAL
		my ($ts2) = $ts;
		if ( $ts2 =~ /^ENS/ ) { $ts2 =~ s/\.\d*$// }
		my ($omega_trans_rep); 
		if ( exists $omega_trans_report->{$ts} ) {
			$omega_trans_rep = $omega_trans_report->{$ts};
		}
		elsif ( exists $omega_trans_report->{$ts2} ) {
			$omega_trans_rep = $omega_trans_report->{$ts2};
		}
		# IMPORTANT: THIS COMMENT IS TEMPORAL		
		while ( my ($type,$o_t_report) = each(%{$omega_trans_rep}) )
		{
			my ($omega_annot) = $o_t_report->{'annot'};			
			if($omega_annot eq $UNKNOWN_LABEL) {
				$ts_inertia_annot = $UNKNOWN_LABEL;
			}			
			$num_omega_results++; # count the number of Omega results
			# initialize the exons with NO
			foreach my $e_index (keys %{$o_t_report->{'exons'}}) {
				$inertia_trans_report->{$ts}->{'exons'}->{$e_index} = {
					'annot'			=> $NO_LABEL,
					'start'			=> $o_t_report->{'exons'}->{$e_index}->{'start'},
					'end'			=> $o_t_report->{'exons'}->{$e_index}->{'end'},
					'strand'		=> $o_t_report->{'exons'}->{$e_index}->{'strand'},
					'om_average'	=> $o_t_report->{'exons'}->{$e_index}->{'om_average'},
					'p_value'		=> $o_t_report->{'exons'}->{$e_index}->{'p_value'}
				};
			}			
		}
		
		# We don't reject a transcript if at least has got 2 Omega results
		if ( $num_omega_results < 2 ) { $ts_inertia_annot = $UNKNOWN_LABEL; }
		
		# Inertia transcript
		$inertia_trans_report->{$ts}->{'annot'} = $ts_inertia_annot;
		
		# Inertia exons: exon is rejected if all method of omega reject itself
		while ( my ($type,$o_t_report) = each(%{$omega_trans_rep}) )
		{
		   	foreach my $e_index (keys(%{$inertia_trans_report->{$ts}->{'exons'}})) {				
				if ( exists $o_t_report->{'exons'}->{$e_index} and defined $o_t_report->{'exons'}->{$e_index} and
					 exists $o_t_report->{'exons'}->{$e_index}->{'annot'} and defined $o_t_report->{'exons'}->{$e_index}->{'annot'}
				) {
					if ( $o_t_report->{'exons'}->{$e_index}->{'annot'} eq $UNKNOWN_LABEL ) {
						$inertia_trans_report->{$ts}->{'exons'}->{$e_index}->{'annot'} = $UNKNOWN_LABEL;
					}
				}
			}
		}
	}
	
	# We reject a transcript if at least one consensus exon has been rejected
	my ($inertia_rejected_exon);
	my ($inertia_length_trans);
	while ( my ($ts,$i_t_report) = each(%{$inertia_trans_report}) )
   	{
		if (exists $i_t_report->{'exons'} and defined $i_t_report->{'exons'}) {
   			my ($ts_inertia_annot) = $UNKNOWN_LABEL; # init
			
			# scan consensus exons by order
			my @sort_e_index = sort {$a <=> $b} keys %{$i_t_report->{'exons'}};			
			for ( my $e_i = 0; $e_i < scalar(@sort_e_index); $e_i++ ) {
				my ($e_index) = $sort_e_index[$e_i];
				
				if ( exists $i_t_report->{'exons'}->{$e_index} and
					 exists $i_t_report->{'exons'}->{$e_index}->{'annot'} and 
					 defined $i_t_report->{'exons'}->{$e_index}->{'annot'} )
				{
					# if one consensus exon is rejected, then the transcript is rejected
					if ( $i_t_report->{'exons'}->{$e_index}->{'annot'} eq $NO_LABEL ) {
						$ts_inertia_annot = $NO_LABEL;
						
						# get the rejected exons sorting by omega value
						if ( exists $i_t_report->{'exons'}->{$e_index}->{'om_average'} ) {
							my ($om_avg) = $i_t_report->{'exons'}->{$e_index}->{'om_average'};
							my ($rej_exon) = $i_t_report->{'exons'}->{$e_index};
							$inertia_rejected_exon->{$om_avg}->{$ts}->{'exons'}->{$e_index} = $rej_exon;
						}
					}
					# get the transcript length from exons
					if ( $e_i == 0 ) {
						$inertia_length_trans->{$ts}->{'start'} = $i_t_report->{'exons'}->{$e_index}->{'start'};
					}
					if ( $e_i == (scalar(@sort_e_index)-1) ) {
						$inertia_length_trans->{$ts}->{'end'} = $i_t_report->{'exons'}->{$e_index}->{'end'};
						$inertia_length_trans->{$ts}->{'strand'} = $i_t_report->{'exons'}->{$e_index}->{'strand'};
					}
				}
			}
			$i_t_report->{'annot'} = $ts_inertia_annot;
   		}
   	}
	
	# Get the genome area of the whole transcripts
	my ($area_transcripts);
	while ( my ($ts, $ts_rep) = each(%{$inertia_length_trans}) ) {
		if ( defined $area_transcripts ) {
			if ( $ts_rep->{'start'} < $area_transcripts->{'start'} ) {
				$area_transcripts->{'start'} = $ts_rep->{'start'};
			}
			if ( $ts_rep->{'end'} > $area_transcripts->{'end'} ) {
				$area_transcripts->{'end'} = $ts_rep->{'end'};
			}
		}
		else {
			$area_transcripts->{'start'} = $ts_rep->{'start'};
			$area_transcripts->{'end'} = $ts_rep->{'end'};
			$area_transcripts->{'strand'} = $ts_rep->{'strand'};
		}		
	}
		
   	# The consensus exon that rejects a transcript must fall over other transcript
   	# If this is not the case, we ignore that consensus exon.
   	# We scan the exons sorting by omega value
   	if ( defined $inertia_rejected_exon )
   	{
   		my @sort_o_index = sort {$b <=> $a} keys %{$inertia_rejected_exon};
   		foreach my $o_index (@sort_o_index) {
   			my ($rej_exons) = $inertia_rejected_exon->{$o_index};
   			my ($possibly_reject_exon_trans);
			while ( my ($ts,$rej_e_report) = each(%{$rej_exons}) )
			{
				if (exists $rej_e_report->{'exons'} and defined $rej_e_report->{'exons'}) {
					while ( my ($ref_e_index, $rej_e_rep) = each(%{$rej_e_report->{'exons'}}) )
					{
						my ($rej_e_start) = $rej_e_rep->{'start'};
						my ($rej_e_end) = $rej_e_rep->{'end'};
						my ($num_overlapping) = 0;
						
						while ( my ($ts2, $length_trans_rep) = each(%{$inertia_length_trans}) ) {
							if ( $ts ne $ts2 ) {
								my ($e_start) = $length_trans_rep->{'start'};
								my ($e_end) = $length_trans_rep->{'end'};
								if ( ( $rej_e_start < $e_end ) and ( $rej_e_end > $e_start ) ) {
									$num_overlapping++;
									$possibly_reject_exon_trans->{$ts} = $ref_e_index;
								}
							}
						}
					   	# if there is not overlapping then we discard this rejection
					   	if ( $num_overlapping == 0 ) {
					   		#$inertia_trans_report->{$ts}->{'exons'}->{$ref_e_index}->{'annot'} = $UNKNOWN_LABEL;
					   	} 
					}
				}		
			}
			# we get the transcripts that are not involved in the possibly rejected transcripts
			my (@the_other_transcripts);
			foreach my $ts ( keys(%{$inertia_length_trans}) ) {
				unless ( defined $possibly_reject_exon_trans and exists $possibly_reject_exon_trans->{$ts} ) {
					push(@the_other_transcripts, $ts);
				}
			}
			# we can not reject an exon coming from several transcripts if modify the genome area of all transcripts
			my ($new_area_transcripts);
			foreach my $ts ( @the_other_transcripts ) {
				if ( defined $new_area_transcripts ) {
					if ( $inertia_length_trans->{$ts}->{'start'} < $new_area_transcripts->{'start'} ) {
						$new_area_transcripts->{'start'} = $inertia_length_trans->{$ts}->{'start'};
					}
					if ( $inertia_length_trans->{$ts}->{'end'} > $new_area_transcripts->{'end'} ) {
						$new_area_transcripts->{'end'} = $inertia_length_trans->{$ts}->{'end'};
					}
				}
				else {
					$new_area_transcripts->{'start'} = $inertia_length_trans->{$ts}->{'start'};
					$new_area_transcripts->{'end'} = $inertia_length_trans->{$ts}->{'end'};
					$new_area_transcripts->{'strand'} = $inertia_length_trans->{$ts}->{'strand'};
				}
			}
			if (($new_area_transcripts->{'start'} > $area_transcripts->{'start'}) or
				($new_area_transcripts->{'end'} < $area_transcripts->{'end'})
			) {
				while ( my ($ts, $ref_e_index) = each (%{$possibly_reject_exon_trans}) ) {
					$inertia_trans_report->{$ts}->{'exons'}->{$ref_e_index}->{'annot'} = $UNKNOWN_LABEL;
				}
			}			 	
   		}
   	}   	
   	
   	# Again, We check if at least one consensus exon has been rejected, then a transcript is rejected
	while ( my ($ts,$i_t_report) = each(%{$inertia_trans_report}) )
   	{
		if (exists $i_t_report->{'exons'} and defined $i_t_report->{'exons'}) {
   			my ($ts_inertia_annot) = $UNKNOWN_LABEL; # init
			
			# scan consensus exons
			while ( my ($e_index, $e_report) = each(%{$i_t_report->{'exons'}}) )
			{
				# if one consensus exon is rejected, then the transcript is rejected				
				if ( exists $e_report->{'annot'} and defined $e_report->{'annot'} ) {					
					if ( $e_report->{'annot'} eq $NO_LABEL ) {
						$ts_inertia_annot = $NO_LABEL;						
					}
				}
			}
			$i_t_report->{'annot'} = $ts_inertia_annot;
   		}
   	}   	

	# Discard the transcripts whose CDS start and CDS end have not found   	
   	foreach my $ts (keys(%{$inertia_trans_report}))
   	{   	
    	my($gff_trans_content)=`grep $ts $gff_file | awk '{if(\$3==\"transcript\"){print \$0}}'`; # one line
   		if ( ($gff_trans_content =~ /CDS start not found/) or ($gff_trans_content =~ /CDS end not found/) ) {
	   		$inertia_trans_report->{$ts}->{'annot'} = $NO_LABEL;
	   	}
   	}

	# We don't admit all transcript has got NO annotation
	my ($global_inertia_rejections) = 0;	
   	foreach my $ts (keys(%{$inertia_trans_report}))
   	{
   		if ( $inertia_trans_report->{$ts}->{'annot'} eq $NO_LABEL ) {
   			$global_inertia_rejections++;
   		}
   	}
	if ( $global_inertia_rejections == scalar( keys(%{$inertia_trans_report})) ) {
	   	foreach my $ts (keys(%{$inertia_trans_report}))
	   	{
	   		$inertia_trans_report->{$ts}->{'annot'} = $UNKNOWN_LABEL;
	   	}		
	}
	
	# Print annotations of transcripts and theirs sorted exons
	$output_content .= "### inertia 1.0 prediction results ##################################\n";	
	while ( my ($ts,$i_t_report) = each(%{$inertia_trans_report}) )
	{
		$output_content .= ">".$ts."\t".$i_t_report->{'annot'}."\n";
		if (exists $i_t_report->{'exons'} and defined $i_t_report->{'exons'}) {
			my @sort_e_index = sort {$a <=> $b} keys %{$i_t_report->{'exons'}};
			 foreach my $sort_i (@sort_e_index) {
			 	$output_content .= "\t".
			 						$i_t_report->{'exons'}->{$sort_i}->{'start'}.":".
			 						$i_t_report->{'exons'}->{$sort_i}->{'end'}.":".
			 						$i_t_report->{'exons'}->{$sort_i}->{'strand'}."\t".
			 						$i_t_report->{'exons'}->{$sort_i}->{'annot'}."\n";
			 }					
		}		
	}

	if (defined $output_content and $output_content ne '') {
		local(*OUTPUT_FILE);
		open(OUTPUT_FILE, ">$inertia_out_file");
		print OUTPUT_FILE $output_content;
		close(OUTPUT_FILE);		
	}			
}
	
main();

1;


__END__

=head1 NAME

inertia

=head1 DESCRIPTION

Script to calculate omega average in order to flag "strange" exons based on gff coordinates and slr output.
Also, it reports the INERTIA annotations

=head1 SYNOPSIS

inertia

=head2 Required arguments:

    --conf <Config file>
    
    --gff <GFF file that contains CDS information>

	--inpath <Input path where alignments are located>

    --outfile <Annotation output file>

=head2 Optional arguments:

	--appris <Flag that enables the output for APPRIS (default: NONE)>

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl inertia.pl
	
	--conf=../conf/pipeline.ini	
		
	--gff=examples/ENSG00000140416/ENSG00000140416.annot.gtf
	
	--inpath=examples/ENSG00000140416/
	
	--output=examples/ENSG00000140416/ENSG00000140416.inertia
	
	--appris	
	
	--loglevel=DEBUG
	
	--logappend
	
	--logpath=examples/ENSG00000140416/
	
	--logfile=ENSG00000140416.log

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)
		
=cut
