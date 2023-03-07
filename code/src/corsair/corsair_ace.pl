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
use Bio::SeqIO;
use JSON;
use Config::IniFiles;
use Data::Dumper;
use File::Copy qw( move );

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile prepare_workspace );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$GIVEN_SPECIES
	$WSPACE_TMP
	$WSPACE_CACHE
	$CODE_CONF_DIR
	$CACHE_FLAG
	$RUN_PROGRAM
	$PROG_DB_PREFIX
	$PROG_DB
	$PROG_DB_V
	$PROG_DB_INV
	$PROG_DB_UID
	$PROG_EVALUE
	$PROG_IDEPER
	$PROG_CUTOFF
	$PROG_GAPLEN
	$PROG_GAPTOT
	$PROG_CDSLEN
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
	$META_CONFIG
	$DIVERGE_TIME
	$NONSCORING_CLADE
);
		
# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($input_file) = undef;
my ($output_file) = undef;
my ($appris) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'conf=s'			=> \$config_file,
	'input=s'			=> \$input_file,
	'output=s'			=> \$output_file,
	'appris'			=> \$appris,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,	
);

# Required arguments
unless ( defined $config_file and defined $input_file and defined $output_file )
{
	print `perldoc $0`;
	exit 1;
}

# Get conf vars
my ($cfg)			= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$GIVEN_SPECIES		= $cfg->val('APPRIS_PIPELINE', 'species');

# $WSPACE_TMP		= $ENV{APPRIS_TMP_DIR};
# $WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$WSPACE_TMP			= "/local/ljmrodriguez/CORSAIR_ACE/tmp";
$WSPACE_CACHE		= "/local/ljmrodriguez/CORSAIR_ACE/cache";

$CODE_CONF_DIR		= $ENV{APPRIS_CODE_CONF_DIR};
$CACHE_FLAG			= $cfg->val('CORSAIR_ACE_VARS', 'cache');
$RUN_PROGRAM		= $cfg->val( 'CORSAIR_ACE_VARS', 'program');

# $PROG_DB_PREFIX		= $ENV{APPRIS_PROGRAMS_DB_DIR};
$PROG_DB_PREFIX		= "/local/ljmrodriguez/CORSAIR_ACE/db";

$PROG_DB_V			= $cfg->val('CORSAIR_ACE_VARS', 'db_v');
$PROG_DB_INV		= $cfg->val('CORSAIR_ACE_VARS', 'db_inv');
$PROG_DB			= undef;
$PROG_DB_UID		= undef;
$PROG_EVALUE		= $cfg->val('CORSAIR_ACE_VARS', 'evalue');
$PROG_IDEPER		= int($cfg->val('CORSAIR_ACE_VARS', 'ideper'));
$PROG_GAPLEN		= int($cfg->val('CORSAIR_ACE_VARS', 'gaplen'));
$PROG_GAPTOT		= int($cfg->val('CORSAIR_ACE_VARS', 'gaptot'));
$PROG_CDSLEN		= int($cfg->val('CORSAIR_ACE_VARS', 'cdslen'));
$PROG_CUTOFF		= $cfg->val('APPRIS_VARS', 'corsair_alt_cutoff');
$OK_LABEL			= 'YES';
$UNKNOWN_LABEL		= 'UNKNOWN';
$NO_LABEL			= 'NO';

my ($meta_config_file) = $CODE_CONF_DIR.'/corsair_alt.meta.json';
$META_CONFIG = JSON->new()->decode( getStringFromFile($meta_config_file) );

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

#####################
# Method prototypes #
#####################
sub _get_id($);
sub extend_exon_until_min_length($$);
sub calculate_score($$);
sub parse_blast($$);
sub check_alignment($$$$\$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Declare variables
	my (@seq_list);
	my ($seq_report);
	my ($vert_score);

	# Get the variables for given species
	if ( exists $META_CONFIG->{$GIVEN_SPECIES} and exists $META_CONFIG->{$GIVEN_SPECIES}->{'animal'} ) {
		my ($animal) = $META_CONFIG->{$GIVEN_SPECIES}->{'animal'};

		if ( $animal eq 'vertebrates' ) {
			$PROG_DB = $PROG_DB_V;
			$PROG_DB_UID = (split('/', $PROG_DB))[0];
		}
		elsif ( $animal eq 'invertebrates' ) {
			$PROG_DB = $PROG_DB_INV;
			$PROG_DB_UID = join('_', (split('/', $PROG_DB))[0], 'invert');
		}
		else {
			$logger->error("Unknown/invalid animal group: '$animal'");
		}

		my ($diverge_time_file) = $CODE_CONF_DIR.'/'.$META_CONFIG->{$GIVEN_SPECIES}->{'diverge_time_file'};
		$DIVERGE_TIME = JSON->new()->decode( getStringFromFile($diverge_time_file) );

		my ($clade_file) = $CODE_CONF_DIR.'/'.$META_CONFIG->{$GIVEN_SPECIES}->{'clade_file'};
		$NONSCORING_CLADE = JSON->new()->decode( getStringFromFile($clade_file) );
	}
	else {
		$logger->error("Species not supported: '$GIVEN_SPECIES'");
	}

	# Handle sequence file
	local(*INFILE);
	open(INFILE, $input_file) or $logger->error("Can not open alignment directory: $!\n");
	my (@infile) = <INFILE>;
	close(INFILE);
	# $logger->debug("INFILE ---------------\n".Dumper(@infile)."\n");

	# Create the input report
	# {
    #       'ENST00000291700' => [
    #                              {
    #                                'seq' => 'MSELEKAMVALIDVFHQYSGREGDKHKLKKSELKELINNELSHFLE',
    #                                'id' => '>ENST00000291700|ENSG00000160307|S100B|46|ENSE00003846634|21|46602278|46602415|-|0',
    #                                'len' => 46
    #                              },
    #                              {
    #                                'len' => 46,
    #                                'id' => '>ENST00000291700|ENSG00000160307|S100B|46|ENSE00001051347|21|46599366|46599503|-|0',
    #                                'seq' => 'EIKEQEVVDKVMETLDNDGDGECDFQEFMAFVAMVTTACHEFFEHE'
    #                              }
    #                            ]
	# }
	my ($infasta);
	my ($s_id,$s) = (undef,undef);
	for (my $i = 0; $i < scalar(@infile); $i++) {
		# get line and remove last newline
		my ($l) = $infile[$i];
		chomp($l);
		# comment line
		if ( $l =~ /^>([^\|]*)\|/ ) {
			# create the new sequence report
			$s_id = $1;
			$s = {
				'id'  => $l,
				'seq' => ''
			};
		}
		# sequence line
		elsif ( $l ne '' ) { $s->{'seq'} .= $l; }
		# if defined the last sequence report with id and seq
		if ( defined $s_id and defined $s and $s->{'seq'} ne '' ) {
			# add the last report into output array
			$s->{'len'} = length($s->{'seq'});
			push(@{$infasta->{$s_id}}, $s);
			($s_id,$s) = (undef,undef);
		}
	}
	# $logger->debug("INFASTA ---------------\n".Dumper($infasta)."\n");


	# Go through every transcript
	while ( my ($transc_id, $t_report) = each (%{$infasta}) ) {
		for (my $i = 0; $i < scalar(@{$t_report}); $i++) {

			# Extend the current exon using the flanking cds's...
			# my ($sequence_id, $sequence, $seq_length, $ext_sequence, $ext_seq_length) = extend_exon($i, $t_report);

			# Add the flanking cds's if the cds length is minus than...
			my ($sequence_id, $sequence, $seq_length, $ext_sequence, $ext_seq_length) = extend_exon_until_min_length($i, $t_report);

			# Calculate the scores for the current CDS and the extesion
			$logger->info("$sequence_id ---------------\n");
			# $logger->info("$seq_length\n");
			# $logger->info("$sequence\n");
			my ($species_score, $species_report) = calculate_score($sequence, $seq_length);
			# $logger->info("***\n");
			# $logger->info("$ext_seq_length\n");
			# $logger->info("$ext_sequence\n");
			my ($ext_species_score, $ext_species_report) = calculate_score($ext_sequence, $ext_seq_length);

			# Save specie report
			$seq_report->{$sequence_id} = {
				'score'			=> $species_score,
				'report'		=> $species_report,
				'ext_score'		=> $ext_species_score,
				'ext_report'	=> $ext_species_report,
			};
						
			# Save score data
			push(@{$vert_score->{$species_score}}, $sequence_id);
			
			# Save seq ids
			push(@seq_list, $sequence_id);
			
		}
	}
	
	# Print records per sequence
	my ($output_content) = "";
	foreach my $seq_id (@seq_list) {
		# scores of species for global aligns
		if ( exists $seq_report->{$seq_id} ) {
			my ($seq_rep) = $seq_report->{$seq_id};
			my ($seq_score) = $seq_rep->{'score'};
			my ($sp_rep) = $seq_rep->{'report'};
			my ($ext_seq_score) = $seq_rep->{'ext_score'};
			my ($ext_sp_rep) = $seq_rep->{'ext_report'};
			$output_content .= ">".$seq_id."\t".$seq_score."\t".$ext_seq_score."\n";
			foreach my $sp_found (@{$sp_rep}) {
				$output_content .= $sp_found."\n";
			}
			$output_content .= "***\n";
			foreach my $sp_found (@{$ext_sp_rep}) {
				$output_content .= $sp_found."\n";
			}
		}
	}	
	
	# Get the annotations for the main isoform /* APPRIS */ ----------------
	if ( defined $appris )
	{		
		$output_content .= get_appris_annotations($vert_score);
	}
	
	# Print records by transcript ---------------
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}	
	

	$logger->finish_log();
	
	exit 0;	
}

sub _get_id($)
{
	my ($s_id) = @_;
	my ($id) = '';
	my (@s_ids) = split(/\|/, $s_id);
	if ( scalar(@s_ids) >= 10 ) { $id = $s_ids[3]."|".$s_ids[4]."|".$s_ids[5]."|".$s_ids[6]."|".$s_ids[7]."|".$s_ids[8]."|".$s_ids[9] }
	return $id;
}

sub extend_exon_until_min_length($$)
{
	my ($i, $t_report) = @_;
	
	# get variables from transcript report
	my ($sequence_id) = $t_report->[$i]->{'id'};
	$sequence_id =~ s/^>//;
	my ($sequence) = $t_report->[$i]->{'seq'};
	my ($seq_length) = $t_report->[$i]->{'len'};
	my ($ext_sequence) = $sequence;
	my ($ext_seq_length) = $seq_length;

	# Add the flanking cds's if the cds length is minus than...
	if ( $i == 0) { # N-Term
		my ($k) = $i+1;
		do {
			my ($s) = $t_report->[$k];
			# get id
			my ($s_id) = _get_id($s->{'id'});
			# concatenate id and seq
			$sequence_id = $sequence_id ."__C-Term|". $s_id;
			$ext_sequence = $ext_sequence . $s->{'seq'};
			$ext_seq_length = length($ext_sequence);
			$k++;
		} while ( $k <= scalar(@{$t_report})-1 and $ext_seq_length < $PROG_CDSLEN );
	}
	elsif ( $i == scalar(@{$t_report})-1 ) { # C-Term
		my ($j) = $i-1;
		do {
			my ($s) = $t_report->[$j];
			# get id
			my ($s_id) = _get_id($s->{'id'});
			# concatenate id, seq, and seq_len
			$sequence_id = $sequence_id .'__N-Term|'. $s_id;
			$ext_sequence = $s->{'seq'} . $ext_sequence;
			$ext_seq_length = length($ext_sequence);
			$j--;
		} while ( $j >= 0 and $ext_seq_length < $PROG_CDSLEN );
	}
	else { # Internal CDS's
		my ($j) = $i-1;
		my ($k) = $i+1;
		do {
			if ( $j >= 0 ) {
				my ($s) = $t_report->[$j];
				# get id
				my ($s_id) = _get_id($s->{'id'});
				# concatenate id, seq, and seq_len
				$sequence_id = $sequence_id .'__N-Term|'. $s_id;
				$ext_sequence = $s->{'seq'} . $ext_sequence;
				$ext_seq_length = length($ext_sequence);
				$j--;
			}
			if ( $k <=  scalar(@{$t_report})-1 ) {
				my ($s) = $t_report->[$k];
				# get id
				my ($s_id) = _get_id($s->{'id'});
				# concatenate id, seq, and seq_len
				$sequence_id = $sequence_id .'__C-Term|'. $s_id;
				$ext_sequence = $ext_sequence . $s->{'seq'};
				$ext_seq_length = length($ext_sequence);
				$k++;
			}
		} while ( ($j >= 0 or $k <=  scalar(@{$t_report})-1) and $ext_seq_length < $PROG_CDSLEN );
	}
	return ($sequence_id, $sequence, $seq_length, $ext_sequence, $ext_seq_length);
} # ebd extend_exon_until_min_length

sub calculate_score($$)
{
	my ($sequence, $seq_length) = @_;

	# Create cache obj
	my ($cache) = APPRIS::Utils::CacheMD5->new(
		-dat => $sequence,
		-ws  => $WSPACE_CACHE			
	);		
	my ($seq_idx) = $cache->idx;
	my ($seq_sidx) = $cache->sidx;
	my ($seq_dir) = $cache->dir;
			
	# prepare cache dir
	my ($ws_cache) = $cache->c_dir();
	# prepare tmp dir
	my ($ws_tmp) = $WSPACE_TMP.'/'.$seq_idx;
	prepare_workspace($ws_tmp);			

	# Create temporal file
	my ($fasta_sequence_file) = $ws_cache.'/seq.faa';
	unless(-e $fasta_sequence_file and (-s $fasta_sequence_file > 0) ) # Cached fasta
	{			
		my ($fasta_sequence_cont) = ">Query\n$sequence";
		my ($print_out) = printStringIntoFile($fasta_sequence_cont, $fasta_sequence_file);
		unless( defined $print_out ) {
			$logger->error("Can not create tmp file: $!\n");
		}
	}

	# Get the db file	
	my ($PROG_DB_PATH) = $PROG_DB_PREFIX.'/'.$PROG_DB;

	# Run blast
	my ($tmp_blast_file) = $ws_tmp.'/seq.'.$PROG_DB_UID;
	my ($blast_sequence_file) = $ws_cache.'/seq.'.$PROG_DB_UID;
	unless (-e $blast_sequence_file and (-s $blast_sequence_file > 0) and ($CACHE_FLAG eq 'yes')) # Blast Cache
	{
		eval
		{
			$logger->info("Running blast\n");
			my ($cmd) = "$RUN_PROGRAM -d $PROG_DB_PATH -i $fasta_sequence_file -e$PROG_EVALUE -o $tmp_blast_file";
			$logger->debug("$cmd\n");						
			system($cmd) == 0 or $logger->error("system call exit code: $?");
		};
		$logger->error("Running blast: $!\n") if($@);

		# Given that no errors were detected, cache the BLAST result file.
		if ( -e $tmp_blast_file && -s $tmp_blast_file > 0 ) {  # avoid caching empty BLAST result
			move($tmp_blast_file, $blast_sequence_file) or $logger->error("Caching blast result\n");
		}
	}
	
	# Parse blast
	$logger->info("Parsing blast\n");
	$logger->debug("Parsing blast file: $blast_sequence_file\n");
	my ($species_score, $species_report) = parse_blast($blast_sequence_file, $seq_length);

	return ($species_score, $species_report);
} # end calculate_score

sub parse_blast($$)				# reads headers for each alignment in the blast output
{
	my ($blast_file,$seq_length) = @_;
	my ($species_score, $species_report) = (0, undef);
	my ($species_found);
	my ($aln_report);

	# TODO: use BioPerl BLAST parser
	open (BLASTFILE, $blast_file) or die "on the spot";

	my $string = "";
	my $length = 0;
	my $faalen = 0;
	my $species = "";

	while (<BLASTFILE>)
	{
		if ($_ eq "\n")
			{$string = "";}
		chomp;
		$_ =~ s/^\s+//;
		if(/letters+\)/)			# gets query seqlen once
		{
			my @temp = split " ";
			$faalen = $temp[0];
			$faalen =~ s/(^\(|,)//g;
			if ($faalen != $seq_length)
				{ $logger->warning("Length of Blast's query is different than input sequence: $faalen != $seq_length\n"); }
		}
		if ( $_ =~ /^>/ ) { $string = $_ }
		else { $string = join " ", $string, $_ } # cos BLAST has bad habit of printing species name over two lines			
		# gets subject sequence length
		# Extract the species name
		if(/Length =/)						
		{
			# get the name from RefSeq db and UniProt db
			if ( $string =~ /\[([^\]]*)\]\s*Length = ([0-9]*)/ or $string =~ /OS=([^\s]*\s+[^\s]*)\s*.*Length = ([0-9]*)\s*$/ ) {
				$species = $1;
				$length = $2;
			}
		}
		if(/Identities/)
		{
			my @data = split " ";
			my @iden = split /\//, $data[2];
			my @identities = split " ", $iden[0];
			my $identity = $iden[0]/$iden[1]*100;
			my $species_firstname = ( $species =~ /^([^\s]*)/ ) ? $1 : ""; #get the name from RefSeq db and UniProt db

			if ($identity < $PROG_IDEPER) # no points when the identity is less than...
				{ }
			elsif (exists $species_found->{$species}) # no points when the species has been already recorded
				{ }
			elsif ( exists $NONSCORING_CLADE->{$species_firstname} ) # skip closely related species
				{ }
			else
			{
				# check if the alignment keep all requirements
				my ($aln_score,$aln_sms) = check_alignment($length,$species,$faalen,$seq_length,$aln_report);
				if ( defined $aln_score and ($aln_score > 0) ) {
					# get identity score
					my ($aln_iden) = sprintf '%.2f', $identity;
					my ($iden_score) = $aln_iden/100;
					my ($divtime_score) = divtime_score($species);
					# save global score for transc
					if ( $aln_score > 0 ) {
						push(@{$species_report}, "$species\t$iden_score\t$divtime_score"); # Record species and score
						$species_found->{$species} = 1;							
					}
					$species_score = $divtime_score if ( $divtime_score > $species_score );
				}
			}
		}
	}

	return ($species_score, $species_report);
} # end parse_blast

# (Normalized) Score based on diverge time
sub divtime_score($)
{
	my ($species) = shift;
	my ($species_firstname) = ( $species =~ /^([^\s]*)/ ) ? $1 : "";	
	my ($score) = 0;
	if ( exists $DIVERGE_TIME->{$species_firstname} ) {
		foreach my $div_time (@{$DIVERGE_TIME->{$species_firstname}}) {
			if ( $div_time->{'name'} eq $species ) {
				$score = $div_time->{'norm'};
			}
		}
	}
	return $score;
}

# check if the alignment keep all requirements
sub check_alignment($$$$\$)
{
	my $candlength = shift;
	my $specie = shift;	
	my $targlength = shift;
	my $seq_length = shift;
	my $ref_report = shift;
	my $oldinput = "dummy";
	my @target = ();
	my @startq = ();
	my @endq = ();
	my @aln_lines = ();
	my @candidate = ();
	my @startc = ();
	my @endc = ();
	my ($aln_score) = 0;
	my ($aln_sms) = '';

	my ($aln_offset) = undef;
	while (<BLASTFILE>) {
		chomp;
		if (/^Query:\s+([0-9]+)\s+(\S+)\s+([0-9]+)/) {		# read in query line
			if ( ! defined $aln_offset ) {
				$aln_offset = $-[2];
			}
			elsif ( $-[2] != $aln_offset ) {
				$logger->error("cannot parse query line - inconsistent BLAST alignment offset\n");
			}
			push @target, $2;
			push @startq, $1;
			push @endq, $3;
		}
		elsif (/^Sbjct:\s+([0-9]+)\s+(\S+)\s+([0-9]+)/) {		# read in subject line
			if ( ! defined $aln_offset ) {
				$logger->error("cannot parse subject line - undefined BLAST alignment offset\n");
			}
			elsif ( $-[2] != $aln_offset ) {
				$logger->error("cannot parse subject line - inconsistent BLAST alignment offset\n");
			}
			push @candidate, $2;
			push @startc, $1;
			push @endc, $3;
		}
		elsif (/^\s/) { 									# read in match line
			if ( ! defined $aln_offset ) {
				$logger->error("cannot parse match line - undefined BLAST alignment offset\n");
			}
			push @aln_lines, substr($_, $aln_offset);
		}
		elsif ($_ eq '' && $oldinput eq '')			# two carriage returns in a row mean alignment has ended
			{last}
		$oldinput = $_;
		}
	
	# process alignment ...

	# process the target and candidate lines
	# get the residues in array for the target and candidate
	my $target = join "", @target;
	my $candidate = join "", @candidate;
	@target = split "", $target;
	@candidate = split "", $candidate;

	# $logger->debug("---------------------------------\n"); 
	# $logger->debug("TARGET:$target|\n");
	# $logger->debug("CANDID:$candidate|\n");
	# $logger->debug("ALIGN_:@aln_lines|\n");
	
	# filter: the align line can not contain more than N gaps
	my $aln = join "", @aln_lines;
	# reject if Alignment has longer unmatches together
	if ( $aln =~ /[ ]{$PROG_GAPLEN,}/ ) {
		# $logger->debug("Alignment has longer unmatches\n");
		return (0,"Alignment has longer unmatches");
	}
	my $unmatchres = $aln =~ tr/ //;
	# $logger->debug("COUNT TOTAL UNMATCHES: $unmatchres\n");
	# reject if Alignment has too many unmatches in total
	if ( $unmatchres > $PROG_GAPTOT ) {
		# $logger->debug("Alignment has too many unmatches in total\n");
		return (0,"Alignment has too many unmatches in total");
	}
	
	# get the start and end indexes for the target and candidate
	my $targstart = $startq[0];
	my $targend = $endq[$#endq];
	my $candstart = $startc[0];
	my $candend = $endc[$#endc];	
	
	# get the start and end for the query
	my @boundaries = (1, $seq_length);
	if ($boundaries[0] < $targstart)
		{ $boundaries[0] = $targstart; }
	
	# count the identity resifues
	# count the gaps for target and candidate
	my $identities = 0;
	my $gapres = 0;
	my ($gapresconttarg, $gaprescontcand) = (0,0);
	my $totalres = $boundaries[1] - $boundaries[0] + 1;
	for ( my $res=0, my $j=$boundaries[0]; $j <= $boundaries[1]; $j++, $res++ ) {
		if(defined $target[$res] and $candidate[$res]) {
			if ($target[$res] eq $candidate[$res])
				{ $identities++ }
			if ($target[$res] eq "-") {
				$gapres++; $j--;
				$gapresconttarg++;
				if ( $gapresconttarg > $PROG_GAPLEN ) {
					# $logger->debug("Long gap in target\n");
					return (0,"Long gap in target"); }
			}
			if ($candidate[$res] eq "-") {
				$gapres++;
				$gaprescontcand++;
				if ( $gaprescontcand > $PROG_GAPLEN ) {
					# $logger->debug("Long gap in candidate\n");
					return (0,"Long gap in candidate"); }
			}
			# reset counters for target and candidate gaps together
			if ($target[$res] ne "-") { $gapresconttarg=0 }
			if ($candidate[$res] ne "-") { $gaprescontcand=0 }
		}
	}
	# $logger->debug("GAP_RES_TARGET: $gapresconttarg\n");
	# $logger->debug("GAP_RES_CANDIDATE: $gaprescontcand\n");
	# $logger->debug("TOTAL_RES: $totalres\n");
	# $logger->debug("TOTAL_UNMATCHES: $unmatchres\n");
	# $logger->debug("GAP_RES: $gapres\n");
	# $logger->debug("IDE_RES: $identities\n");
	
	# reject if Alignment has too many gaps in total
	if ( $gapres > $PROG_GAPTOT ) {
		# $logger->debug("Alignment has too many gaps in total\n");
		return (0,"Alignment has too many gaps in total");
	}

	# save the values
	my $aln_flag = 1;
	my $identity = 0;
	my $gaps = 0;
	my ($specie_point) = 1; # in the case of species has weight but it is not the case
	if ($totalres > 0) {
		$identity = $identities/$totalres*100;
		$gaps = $gapres/$totalres*100;
	}
	$$ref_report->{$specie} = {
		'iden'		=> $identity,
		'gaps'		=> $gaps,
		'totalres'	=> $totalres
	};

	# score of global align
	$aln_score = $aln_flag*$specie_point;
	
	return ($aln_score, $aln_sms);							# if sequence passes all tests
	
} # end check_alignment

# This method selects the best isoform. It works taking into account the corsair score, or using the num. of difference species 
sub get_appris_annotations($)
{
	my ($vertebrate_score) = @_;
	my ($cutoffs);
	my ($output_content) = '';

	$output_content .= "\n";
	$output_content .= "# ================================ #\n";
	$output_content .= "# Conservation against vertebrates #\n";
	$output_content .= "# ================================ #\n";
	
	if(defined $vertebrate_score)
	{
		my(@vertebrate_score_list) = sort { $a <= $b } keys(%{$vertebrate_score} );
		
		if(scalar(@vertebrate_score_list)>0)
		{
			# We tag the transcript as UNKOWN whose num domains are biggest
			my(@trans_biggest_score);
			my($biggest_vertebrate_score)=$vertebrate_score_list[0];
			foreach my $trans_id (@{$vertebrate_score->{$biggest_vertebrate_score}})
			{
				push(@trans_biggest_score,$trans_id);
			}
			my($unique)=1;
			for(my $i = 1; $i < scalar(@vertebrate_score_list); $i++)
			{
				my ($current_score) = $vertebrate_score_list[$i];
				
				# If the biggest score is bigger or equal than $PROG_CUTOFF => the transcripts are rejected
				if ( ($biggest_vertebrate_score - $current_score) >= $PROG_CUTOFF )
				{
					foreach my $trans_id (@{$vertebrate_score->{$current_score}})
					{
						$cutoffs->{$trans_id} = $NO_LABEL;
					}
				}
				else
				{
					$unique=0;
					foreach my $trans_id (@{$vertebrate_score->{$current_score}})
					{
						$cutoffs->{$trans_id} = $UNKNOWN_LABEL;
					}					
				}
			}
			# There is one transcript with the biggest score
			if (($unique == 1) and scalar(@trans_biggest_score) == 1)
			{
				foreach my $trans_id (@trans_biggest_score)
				{
					$cutoffs->{$trans_id} = $OK_LABEL;				
				}
			}
			else
			{
				foreach my $trans_id (@trans_biggest_score)
				{
					$cutoffs->{$trans_id} = $UNKNOWN_LABEL;				
				}
			}
		}
	}
	
	# Get appris output
	while ( my ($trans_id,$annot) = each (%{$cutoffs}) ) {
		$output_content .= ">".$trans_id."\t".$annot."\n";
	}
	return $output_content;
}


main();

__END__

=head1 NAME

CORSAIR_ACE

=head1 DESCRIPTION

Run CORSAIR_ACE program

The CORSAIR_ACE is like Alt-CORSAIR but the scores are for each CDS.
The alignment is with the CDS and the flaking CDS until to have a minimum of 50aa in total

=head1 VERSION

0.1

=head2 Required arguments:

	--conf <Config file>

	--input <Fasta sequence file with the CDS sequences>

	--output <Annotation output file>

=head2 Optional arguments:

	--appris <Flag that enables the output for APPRIS (default: NONE)>
	
	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
	

=head1 EXAMPLE

perl corsair.pl

	--conf=../conf/pipeline.ini

	--input=example/OTTHUMG00000020713/cdsseq.faa
	
	--output=example/OTTHUMG00000020713/corsair_ace

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
