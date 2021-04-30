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
	$CACHE_FLAG
	$RUN_PROGRAM
	$PROG_DB_PREFIX
	$PROG_DB
	$PROG_DB_V
	$PROG_DB_INV
	$PROG_EVALUE
	$PROG_MINLEN
	$PROG_CUTOFF
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
	$DEFALULT_CORSAIR_SPECIES_FILE
	$SPECIES
	$UNMATCHES_THRESHOLD
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
	'input=s'		=> \$input_file,
	'output=s'		=> \$output_file,
	'appris'			=> \$appris,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'		=> \$logfile,
	'logpath=s'		=> \$logpath,
	'logappend'		=> \$logappend,	
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
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$CACHE_FLAG			= $cfg->val('CORSAIR_VARS', 'cache');
$RUN_PROGRAM		= $cfg->val('CORSAIR_VARS', 'program');
$PROG_DB_PREFIX		= $ENV{APPRIS_PROGRAMS_DB_DIR};
$PROG_DB			= undef;
$PROG_DB_V			= $cfg->val('CORSAIR_VARS', 'db_v');
$PROG_DB_INV		= $cfg->val('CORSAIR_VARS', 'db_inv');
$PROG_EVALUE		= $cfg->val('CORSAIR_VARS', 'evalue');
$PROG_MINLEN		= $cfg->val('CORSAIR_VARS', 'minlen');
$PROG_CUTOFF		= $cfg->val('CORSAIR_VARS', 'cutoff');
$OK_LABEL			= 'YES';
$UNKNOWN_LABEL		= 'UNKNOWN';
$NO_LABEL			= 'NO';
$DEFALULT_CORSAIR_SPECIES_FILE	= $ENV{APPRIS_CODE_CONF_DIR}.'/corsair_species.json';
$SPECIES = JSON->new()->decode( getStringFromFile($DEFALULT_CORSAIR_SPECIES_FILE) );
$UNMATCHES_THRESHOLD	= 10;

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
sub parse_blast($$$);
sub check_alignment($$$$\$$);
sub _get_specie_score($);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Declare variables
	my (@seq_list);
	my ($seq_report);
	my ($seq_e_report);
	my ($vert_score);
	
	# Get the variables from determined specie
	my ($PROG_DB_UID);
	if ( exists $SPECIES->{$GIVEN_SPECIES} and exists $SPECIES->{$GIVEN_SPECIES}->{'animal'} ) {
		if ( $SPECIES->{$GIVEN_SPECIES}->{'animal'} eq 'vertebrates' ) {
			$PROG_DB = $PROG_DB_V;
			$PROG_DB_UID = (split('/', $PROG_DB))[0];
		}
		elsif ( $SPECIES->{$GIVEN_SPECIES}->{'animal'} eq 'invertebrates' ) {
			$PROG_DB = $PROG_DB_INV;
			$PROG_DB_UID = join('_', (split('/', $PROG_DB))[0], 'invert');
		}
		else {
			$logger->error("Animal category does not exit");			
		}
	}
	else {
		$logger->error("Species does not exit");
	}

	my ($PROG_DB_PATH) = $PROG_DB_PREFIX.'/'.$PROG_DB;

	# Handle sequence file
	my $in = Bio::SeqIO->new(
						-file => $input_file,
						-format => 'Fasta'
	);
	
	# Scan every fasta sequence
	while ( my $seq = $in->next_seq() )
	{
		if ( $seq->id=~/^([^|]*)\|([^|]*)/ )
		{
			my ($sequence_id) = $2;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
			my ($sequence) = $seq->seq;
			my ($sequence_length) = $seq->length;			
			$logger->info("$sequence_id ---------------\n");
			
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
			
			# If apply, use exon data
			my ($exons);
			my ($edges) = join ":", 1, $sequence_length;
			$exons->{$edges} =  undef;
			
			# Run blast
			my ($tmp_blast_file) = $ws_tmp.'/seq.'.$PROG_DB_UID;
			my ($blast_sequence_file) = $ws_cache.'/seq.'.$PROG_DB_UID;
			unless (-e $blast_sequence_file and (-s $blast_sequence_file > 0) and ($CACHE_FLAG eq 'yes')) # Blast Cache
			{
				eval
				{
					$logger->info("Running blast\n");
					my ($cmd) = "$RUN_PROGRAM -d $PROG_DB_PATH -i $fasta_sequence_file -e $PROG_EVALUE -o $tmp_blast_file";
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
			my ($species_score, $species_report, $exon_species_report) = parse_blast($blast_sequence_file, $exons, $sequence_length);
			
			# Save score data
			push(@{$vert_score->{$species_score}},$sequence_id);
			
			# Save specie report
			$seq_report->{$sequence_id} = {
				'score'		=> $species_score,
				'report'	=> $species_report
			};
			
			# Save specie report exon per exon
			$seq_e_report->{$sequence_id} = $exon_species_report;
			
			# Save seq ids
			push(@seq_list, $sequence_id);
			
		}
	}
	
	# Get the best score among exon sequences (only for exon reports)
	my ($best_seq_e_report);
	while ( my ($seq_id, $seq_rep) = each (%{$seq_e_report}) ) {
		while ( my ($e_coord, $e_rep) = each (%{$seq_rep}) ) {
			if ( exists $e_rep->{'score'} ) {
				my ($e_score) = $e_rep->{'score'};
				unless ( exists $best_seq_e_report->{$e_coord} ) {
					$best_seq_e_report->{$e_coord} = { $e_score => "$seq_id|" };
				}
				else { 
					if ( exists $best_seq_e_report->{$e_coord}->{$e_score} ) {
						$best_seq_e_report->{$e_coord}->{$e_score} .= "$seq_id|";
					}
					else {
						my (@e_scores) = keys(%{$best_seq_e_report->{$e_coord}});
						my ($old_e_score) = $e_scores[0];
						if ( $e_score > $old_e_score ) {
							$best_seq_e_report->{$e_coord} = { $e_score => "$seq_id|" };
						}
						elsif ( $e_score == $old_e_score ) {
							$best_seq_e_report->{$e_coord}->{$e_score} .= "$seq_id|";
						}						
					}
				}				
			}			
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
			$output_content .= ">".$seq_id."\t".$seq_score."\n";
			foreach my $sp_found (@{$sp_rep}) {
				$output_content .= $sp_found."\n";
			}			
		}
		# scores of species for aligns exons by exon (sorted exons)
		if ( exists $seq_e_report->{$seq_id} ) {
			my ($seq_rep) = $seq_e_report->{$seq_id};
			my (@sort_exon_species_report) = sort { # sort by transc start coord
							my ($a2,$b2);
							$seq_rep->{$a}->{'pep_index'} =~ /^([^\:]*)\:/; $a2=$1;
							$seq_rep->{$b}->{'pep_index'} =~ /^([^\:]*)\:/; $b2=$1;
							$a2 <=> $b2
						} keys %{$seq_rep};			
			foreach my $exon_coords_txt (@sort_exon_species_report) {
				my ($exon_sp_rep) = $seq_rep->{$exon_coords_txt};				
				if ( exists $exon_sp_rep->{'pep_index'} and exists $exon_sp_rep->{'score'} and exists $exon_sp_rep->{'species'} and (scalar($exon_sp_rep->{'species'}) > 0) ) {
					my ($exon_pep_index) = $exon_sp_rep->{'pep_index'};
					my ($exon_score) = $exon_sp_rep->{'score'};
					my ($exon_score_txt) = $exon_score;										
					if ( exists $best_seq_e_report->{$exon_coords_txt} and defined $best_seq_e_report->{$exon_coords_txt} ) { # get the best score from other aligned seq
						my ($best_e_rep) = $best_seq_e_report->{$exon_coords_txt};
						my (@best_e_scores) = keys(%{$best_e_rep});
						if ( scalar(@best_e_scores) > 0 ) {
							my ($best_e_score) = $best_e_scores[0];
							my ($best_e_score_int) = sprintf '%.1f', $best_e_scores[0];
							my ($exon_score_int) = sprintf '%.1f', $exon_score;
							if ( $best_e_score_int > $exon_score_int ) {
								my ($best_seq_txt) = '';;
								if ( exists $best_e_rep->{$best_e_score} ) {
									my ($seq_list) = $best_e_rep->{$best_e_score};
									$seq_list =~ s/\|$//;									
									$best_seq_txt = "{".$best_e_score."-".$seq_list."}";
								}
								$exon_score_txt = $exon_score." ".$best_seq_txt;
							}							
						}						
					}					
					$output_content .= "\t- ".$exon_coords_txt."[".$exon_pep_index."]"."\t".$exon_score_txt."\n";
					while (my ($sp,$sp_rep) = each(%{$exon_sp_rep->{'species'}}) ) {
						if ( exists $sp_rep->{'score'} and exists $sp_rep->{'iden'} ) {
							my ($sc) = $sp_rep->{'score'};
							my ($si) = sprintf '%.2f', $sp_rep->{'iden'};
							$output_content .= "\t\t".$sp."\t".$si."\t".$sc."\n";							
						}
					}
				}
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

sub parse_blast($$$)				# reads headers for each alignment in the blast output
{
	my ($blast_file,$exons,$sequence_length) = @_;
	my ($species_score, $species_report, $exon_species_report) = (0, undef, undef);
	my ($species_found);
	my ($aln_report);

	# TODO: use BioPerl BLAST parser
	open (BLASTFILE, $blast_file) or die "on the spot";

	my $string = "";
	my $length = 0;
	my $length_diff = 0;
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
			if ($faalen != $sequence_length)
				{ $logger->warning("Length of Blast's query is different than input sequence\n"); }
		}
		$string = join " ", $string, $_;	# cos BLAST has bad habit of printing species name over two lines
		if($string =~ /[a-z]+\]/)		# ... has species name
		{
			my @data = split /\[/, $string;
			$string = "";
			$species = $data[$#data];
			$species =~ s/\]//;
		}
		if(/Length =/)						# gets subject sequence length
		{
			my @data = split " ";
			$length = $data[2];
		}
		if(/Identities/)
		{
			my @data = split " ";
			my @iden = split /\//, $data[2];
			$length_diff = abs($length - $faalen); # difference in length between query and subject
			my @identities = split " ", $iden[0];
			my $identity = $iden[0]/$iden[1]*100;
			if ($identity < 50)				# gets no points for this sequence
				{ }
			elsif ($length_diff > $PROG_MINLEN) # gets no points for this sequence
				{ }
			elsif (exists $species_found->{$species})				# gets no points for this sequence
				{ }
			elsif ( exists $SPECIES->{$GIVEN_SPECIES} and exists $SPECIES->{$GIVEN_SPECIES}->{'scores'} and !(exists $SPECIES->{$GIVEN_SPECIES}->{'scores'}->{$species}) ) # only we accept determined species
				{ }
			#elsif ( !(exists $SPECIES->{'REST_OF_SPECIES'}->{'scores'}->{$species}) ) # only we accept determined species
			#	{ }
			else
			{
					my ($aln_score,$aln_sms) = check_alignment($length,$species,$faalen,$exons,$aln_report, $_);
					if ( defined $aln_score and ($aln_score >= 0) ) {
						# save global score for transc
						if ( $aln_score > 0 ) {
							my ($aln_iden) = sprintf '%.2f', $identity;
							push(@{$species_report}, "$species\t$aln_iden\t$aln_score"); # Record species and score
							$species_found->{$species} = 1;							
						}
						# save scores per exons (genomics region), if apply
						if ( defined $aln_report ) {
							while ( my ($pep_index, $aln_rep) = each(%{$aln_report}) ) {
								if ( defined $aln_rep ) {
									if ( defined $exons->{$pep_index} ) { # apply
										my ($transc_index) = $exons->{$pep_index};
										my ($exon_species_score) = 0;
										my ($exon_species);
										while (my ($spe, $spe_rep) = each(%{$aln_rep}) ) {
											if ( exists $spe_rep->{'score'} and exists $spe_rep->{'iden'} ) {
												my ($spe_sc) = $spe_rep->{'score'};
												my ($spe_iden) = $spe_rep->{'iden'};
												$exon_species_score += $spe_sc;
												$exon_species->{$spe} = {
													'score'	=> $spe_sc,
													'iden'	=> $spe_iden
												};										
											}
										}
										$exon_species_report->{$transc_index} = {
											'pep_index'	=> $pep_index,
											'score'		=> $exon_species_score,
											'species'	=> $exon_species
										};
									}
								}
							}
						}
						$species_score += $aln_score;						
					}
			}
		}
	}

	return ($species_score, $species_report, $exon_species_report);
}


sub check_alignment($$$$\$$) #parses BLAST alignments exon by exon
{
	my $candlength = shift;	
	my $specie = shift;	
	my $targlength = shift;
	my $exons = shift;
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
	my ($specie_point) = _get_specie_score($specie);

	my ($aln_offset) = undef;
	while (<BLASTFILE>)
		{
		chomp;
		if (/^Query:\s+([0-9]+)\s+(\S+)\s+([0-9]+)/)		# read in query line
			{
			if ( ! defined $aln_offset )
				{
				$aln_offset = $-[2];
				}
			elsif ( $-[2] != $aln_offset )
				{
				$logger->error("cannot parse query line - inconsistent BLAST alignment offset\n");
				}
			push @target, $2;
			push @startq, $1;
			push @endq, $3;
			}
		elsif (/^Sbjct:\s+([0-9]+)\s+(\S+)\s+([0-9]+)/)		# read in subject line
			{
			if ( ! defined $aln_offset )
				{
				$logger->error("cannot parse subject line - undefined BLAST alignment offset\n");
				}
			elsif ( $-[2] != $aln_offset )
				{
				$logger->error("cannot parse subject line - inconsistent BLAST alignment offset\n");
				}
			push @candidate, $2;
			push @startc, $1;
			push @endc, $3;
			}
		elsif (/^\s/)						# read in match line
		{
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
	
	my $target = join "", @target;
	my $candidate = join "", @candidate;
	
	my $targstart = $startq[0];
	my $targend = $endq[$#endq];
	my $candstart = $startc[0];
	my $candend = $endc[$#endc];	
	
	if ($targstart > 4 or $candstart > 4)		# reject if different N-terminal
		{return (0,"It has different N-terminal")}
	
	if ( (abs($candlength - $candend) > 4) or (abs($targlength - $targend) > 4) ) # reject if subject has longer C-terminal
		{return (0,"Subject has longer C-terminal")}

	my $aln = join "", @aln_lines;
	if ($aln =~ /([+\s]{$UNMATCHES_THRESHOLD,})/)  # reject if subject has longer unmatches
		{return (0,"Subject has longer unmatches")}

	@target = split "", $target;
	@candidate = split "", $candidate;
	
	my $aln_flag = 1;
	my $loopstart = 0;
	my (@sort_pep_exons) = sort { # sort pep coord from start value
							my ($a2,$b2);
							$a =~ /([^\:]*)\:/; $a2=$1;
							$b =~ /([^\:]*)\:/; $b2=$1;
							$a2 <=> $b2
						} keys %{$exons};
	for (my $i=0; $i < scalar(@sort_pep_exons); $i++)
	{
		my $pep_index = $sort_pep_exons[$i];
		my @boundaries = split ":", $pep_index;
		if ($boundaries[0] < $targstart)
			{ $boundaries[0] = $targstart; }
		my $identities = 0;
		my $gapres = 0;
		my ($gapresconttarg, $gaprescontcand) = (0,0);
		my ($gapconttarg, $gapcontcand) = ('false','false');		
		my $totalres = $boundaries[1] - $boundaries[0] + 1;
		my $res = 0;
		my $j = 0;
		for ($res=$loopstart,$j=$boundaries[0];$j<=$boundaries[1];$j++,$res++)
			{
				if(defined $target[$res] and $candidate[$res])
				{				
					if ($target[$res] eq $candidate[$res])
						{$identities++}
					if ($target[$res] eq "-") {
						$gapres++;$j--;$gapconttarg = 'true';
						if ( $gapconttarg eq 'true' ) {
							$gapresconttarg++;
							if ( $gapresconttarg >= $UNMATCHES_THRESHOLD ) { return (0,"Long gap in target") }
						}
					}
					if ($candidate[$res] eq "-") {
						$gapres++; $gapcontcand = 'true';
						if ( $gapcontcand eq 'true' ) {
							$gaprescontcand++;
							if ( $gaprescontcand >= $UNMATCHES_THRESHOLD ) { return (0,"Long gap in candidate") }
						}
					}
					if ($target[$res] ne "-") {$gapconttarg = 'false';$gapresconttarg=0}
					if ($candidate[$res] ne "-") {$gapcontcand = 'false';$gaprescontcand=0}
				}
			}
		$loopstart = $res;
		
		my $cds_flag = 1;
		my $identity = 0;
		my $gaps = 0;
		if ($totalres > 0) {
			$identity = $identities/$totalres*100;
			$gaps = $gapres/$totalres*100;
		}			
		if ($identity < 40 && $totalres > 8) { # reject if two exons are too different
			$aln_sms = "Two exons are too different";
			$cds_flag = 0;
			$aln_flag = 0;
		}
		if ($gapres >= $UNMATCHES_THRESHOLD) { # reject if exons have substantial gaps
			$aln_sms = "Exons have substantial gaps";
			$cds_flag = 0;
			$aln_flag = 0;
		}

		my ($cds_score) = $cds_flag*$specie_point;				
		$$ref_report->{$pep_index}->{$specie} = {
			'score'		=> $cds_score,
			'iden'		=> $identity,
			'gaps'		=> $gaps,
			'totalres'	=> $totalres
		};		
	}
	
	# score of global align
	$aln_score = $aln_flag*$specie_point;
	
	return ($aln_score,$aln_sms);							# if sequence passes all tests
	
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

sub _get_specie_score($)
{
	my ($specie) = @_;
	my ($point) = 0;
	if ( exists $SPECIES->{$GIVEN_SPECIES} and exists $SPECIES->{$GIVEN_SPECIES}->{'scores'} and exists $SPECIES->{$GIVEN_SPECIES}->{'scores'}->{$specie} ) { # Punish the specie that is equal than given as input.
		$point = $SPECIES->{$GIVEN_SPECIES}->{'scores'}->{$specie};
	}
	elsif ( exists $SPECIES->{'REST_OF_SPECIES'} and exists $SPECIES->{'REST_OF_SPECIES'}->{'scores'} and exists $SPECIES->{'REST_OF_SPECIES'}->{'scores'}->{$specie} ) { # the rest of species
		$point = $SPECIES->{'REST_OF_SPECIES'}->{'scores'}->{$specie};
	}
	return $point;						
}

main();

__END__

=head1 NAME

corsair

=head1 DESCRIPTION

Run CORSAIR program

=head1 VERSION

0.1

=head2 Required arguments:

	--conf <Config file>

	--input <Fasta sequence file>
	
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

	--input=example/OTTHUMG00000020713.faa
	
	--output=example/OTTHUMG00000020713.output

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
