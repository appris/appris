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

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( parse_file printStringIntoFile getStringFromFile prepare_workspace );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$GIVEN_SPECIES
	$WSPACE_TMP
	$WSPACE_CACHE
	$RUN_PROGRAM
	$PROG_DB
	$PROG_DB_V
	$PROG_DB_INV
	$PROG_EVALUE
	$PROG_MINLEN
	$PROG_CUTOFF
	$PROG_MAXPRO
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
	$DEFALULT_CORSAIR_PRIMATES_FILE
	$DEFALULT_CORSAIR_DIVERGE_TIME_FILE
	$PRIMATES
	$DIVERGE_TIME
);
		
# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($config_file) = undef;
my ($gff_file) = undef;
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
	'gff=s'				=> \$gff_file,
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
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
#$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$WSPACE_CACHE		= $cfg->val( 'CORSAIR2_VARS', 'cache_dir');
$RUN_PROGRAM		= $cfg->val( 'CORSAIR2_VARS', 'program');
$PROG_DB_V			= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('CORSAIR2_VARS', 'db_v');
$PROG_DB_INV		= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('CORSAIR2_VARS', 'db_inv');
# $PROG_DB			= $PROG_DB_V; # vertebrate by default
$PROG_DB			= $cfg->val('CORSAIR2_VARS', 'db');
$PROG_EVALUE		= $cfg->val('CORSAIR2_VARS', 'evalue');
$PROG_MINLEN		= $cfg->val('CORSAIR2_VARS', 'minlen');
$PROG_CUTOFF		= $cfg->val('CORSAIR2_VARS', 'cutoff');
$PROG_MAXPRO		= $cfg->val('CORSAIR2_VARS', 'maxpro');
$OK_LABEL			= 'YES';
$UNKNOWN_LABEL		= 'UNKNOWN';
$NO_LABEL			= 'NO';
$DEFALULT_CORSAIR_PRIMATES_FILE		= $ENV{APPRIS_CODE_CONF_DIR}.'/corsair_alt.primates.json';
$DEFALULT_CORSAIR_DIVERGE_TIME_FILE	= $ENV{APPRIS_CODE_CONF_DIR}.'/corsair_alt.diverge_time.human.json';
$PRIMATES = JSON->new()->decode( getStringFromFile($DEFALULT_CORSAIR_PRIMATES_FILE) );
$DIVERGE_TIME = JSON->new()->decode( getStringFromFile($DEFALULT_CORSAIR_DIVERGE_TIME_FILE) );

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
sub _get_cds_coordinates_from_gff($);
sub _create_block_cdsseq_file($$$);
sub parse_blast($$$);
sub check_alignment($$$$\$$);

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
	my ($vert_score); #DEPRECATED
	
	# Get the CDS coordinates info from gff file ---------------
	$logger->info("Get cds coordinates from gff ---------------\n");
	my ($trans_cds_coords, $cdsblock_coords) = _get_cds_coordinates_from_gff($gff_file);
    $logger->debug("CDS coordenates ---------------\n".Dumper($trans_cds_coords));
    $logger->debug("Sorted CDS coordenates ---------------\n".Dumper($cdsblock_coords));

	# Create fasta file with the joined CDS for the analysis by exons
	my ($cds_file) = _create_block_cdsseq_file($cdsblock_coords, $input_file, $output_file);


	# Handle sequence file
    my $in = Bio::SeqIO->new(
                        -file   => $cds_file,
                        -format => 'Fasta'
    );
	
	# Scan every fasta sequence
	while ( my $seq = $in->next_seq() )
	{
		if ( $seq->id=~/([^|]*)/ )
		{
			my ($cds_id) = $1;
			my ($sequence) = $seq->seq;
			my ($sequence_length) = $seq->length;			
			$logger->info("$cds_id ---------------\n");
			
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
			if ( defined $gff_file and (-e $gff_file) and (-s $gff_file > 0) ) {			
				my(@global_sequence_exon_info)=`grep $cds_id $gff_file`;
				foreach my $sequence_exon_info (@global_sequence_exon_info)
				{
					my(@exon_info) = split /\t/, $sequence_exon_info;				
					my($edges) = join ":", $exon_info[3], $exon_info[4];
					my ($attributes) = $exon_info[8];
					my ($trans_edges);
					if ( $attributes =~ /cds\_coord\>([^\-]*)\-([^\:]*)\:([\-|\+])/ ) {
						$trans_edges = $1.':'.$2;
					}
					elsif ( ($attributes =~ /start\_cds \"([^\"]*)\"/) and ($attributes =~ /end\_cds \"([^\"]*)\"/) ) {
						if ( $attributes =~ /start\_cds \"([^\"]*)\"/ ) { $trans_edges = $1; }
						if ( $attributes =~ /end\_cds \"([^\"]*)\"/ ) {	$trans_edges = $trans_edges.":".$1; }						
					}
					if ( defined $edges and defined $trans_edges ) {
						$exons->{$edges} =  $trans_edges;			
					}
					else {
						$logger->error("getting protein coords\n");	
					}
				}
			}
			else {
				my ($edges) = join ":", 1, $sequence_length;
				$exons->{$edges} =  undef;
			}
			
			# Run blast
			my ($blast_sequence_file) = $ws_cache.'/seq.refseq';
			$logger->debug("Blast file: $blast_sequence_file\n");
			unless (-e $blast_sequence_file and (-s $blast_sequence_file > 0) ) # Blast Cache
			{
				eval
				{
					$logger->info("Running blast\n");
					my ($cmd) = "$RUN_PROGRAM -a $PROG_MAXPRO -d $PROG_DB -i $fasta_sequence_file -e0.001 -G 13 -o $blast_sequence_file";
					$logger->debug("$cmd\n");						
					system($cmd);
				};
				$logger->error("Running blast: $!\n") if($@);
			}
			
			# Parse blast
			$logger->info("Parsing blast\n");			
			my ($species_score, $species_report, $exon_species_report) = parse_blast($blast_sequence_file, $exons, $sequence_length);
			
			# Save score data
			push(@{$vert_score->{$species_score}},$cds_id); # DEPRECATED
					
			# Save specie report exon per exon
			# $seq_e_report->{$sequence_id} = $exon_species_report;
			$seq_e_report->{$cds_id} = {
				'score'		=> $species_score,
				'report'	=> $species_report
			};

			# Get the Maximum value of Blast of the middle exons which that have been created from 3 combinations
			my ($sequence_id) = $cds_id;
			$sequence_id =~ s/\_\_\d//;
			if ( exists $seq_report->{$sequence_id} ) {
				if ( $species_score > $seq_report->{$sequence_id}->{'score'} ) {
					$logger->debug("MAXIMIM in ".$sequence_id." for CDS_ID: ".$cds_id."\n");
					$seq_report->{$sequence_id} = {
						'score'		=> $species_score,
						'report'	=> $species_report
					};
				}
			}
			else {
				$seq_report->{$sequence_id} = {
					'score'		=> $species_score,
					'report'	=> $species_report
				};
			}
			
			# Save seq ids
			push(@seq_list, $sequence_id);
		}
	}
	$logger->debug("Seq Exon report ---------------\n".Dumper($seq_e_report));
	$logger->debug("Seq report ---------------\n".Dumper($seq_report));
	
	# Print records per sequence
	my ($output_content) = "";
	while ( my ($seq_id, $seq_rep) = each (%{$seq_report}) ) {
		my ($seq_score) = $seq_rep->{'score'};
		my ($sp_rep) = $seq_rep->{'report'};
		$output_content .= ">".$seq_id."\t".$seq_score."\n";
		foreach my $sp_found (@{$sp_rep}) {
			$output_content .= $sp_found."\n";
		}			
	}
	
	# Get the annotations for the main isoform /* APPRIS */ ----------------
	# BEGIN: DEPRECATED
	if ( defined $appris )
	{		
		$output_content .= get_appris_annotations($vert_score);
	}
	# END: DEPRECATED
	
	# Print records by transcript ---------------
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}	
	

	$logger->finish_log();
	
	exit 0;	
}

# Get the CDS coordinates info from gff file
sub _get_cds_coordinates_from_gff($)
{
	my ($gff_file) = @_;		
	
	# Get the CDS info from gff file
	my (@global_trans_cds_info)=`awk -F "\t" '{if( \$3=="CDS" ){print \$0}}' $gff_file`;
	$logger->debug('awk -F "\t" \'{if( $3=="CDS" ){print $0}}\' '.$gff_file."\n");
	if ( scalar(@global_trans_cds_info) == 0 ) {
		$logger->error("gff has not transcript information");			
	}

	# Save CDS coordinates within structure
	# my ($total_trans_cds_coords);
	my ($trans_cds_coords);
	foreach my $trans_cds_info (@global_trans_cds_info)
	{
		if (defined $trans_cds_info and $trans_cds_info ne "")
		{			
			my (@cds_info) = split /\t/, $trans_cds_info;
			if(defined $cds_info[8] and ($cds_info[8] =~ /transcript_id "([^\"]*)"/))
			{
				my ($sequence_id) = $1;				
				my ($strand);
				if(defined $cds_info[6])
				{
					$strand=$cds_info[6];
					$trans_cds_coords->{$sequence_id}->{'strand'} = $strand
				}
				if(defined $cds_info[0] and defined $cds_info[3] and defined $cds_info[4])
				{					
					my ($chr) = $cds_info[0];
					my ($start) = $cds_info[3];
					my ($end) = $cds_info[4];
					my ($edges);
					# Get the CDS coord for each transcript (with frame)
					if(defined $cds_info[7]) # save the frame
					{
						my ($frame) = $cds_info[7];
						$edges = $chr."-". join ":", $start, $end, $frame;
					}
					else {
						$edges = $chr."-". join ":", $start, $end;
					}
					push(@{$trans_cds_coords->{$sequence_id}->{'coord'}}, $edges);
				}
			}
		}
	}
	# Only applied for transcripts with more than one CDS!!
    # Extract the Block coordinates 3 by 3 for the middle CDS
	# Extract the Block coordinates 2 by 2 when they are N-term and C-term
    my ($cdsblock_coords);
    while(my ($sequence_id,$cds_coords) = each(%{$trans_cds_coords}))
    {
		if ( exists $cds_coords->{'coord'} and scalar(@{$cds_coords->{'coord'}}) > 1 )
		{
			my ($coords) = $cds_coords->{'coord'};
			my ($strand) = $cds_coords->{'strand'};
			for( my $i=0; $i<scalar(@{$coords}); $i+=1 ) {
				my ($cdsblock_comp);
				if ($i == 0) {
					$cdsblock_comp = $coords->[$i]."_".$sequence_id   ."+". $coords->[$i+1]."_".$sequence_id;
					$cdsblock_coords->{$cdsblock_comp}->{$coords->[$i]} = $strand;
				} elsif ($i == (scalar(@{$coords})-1) ) {
					$cdsblock_comp = $coords->[$i-1]."_".$sequence_id ."+". $coords->[$i]."_".$sequence_id;
					$cdsblock_coords->{$cdsblock_comp}->{$coords->[$i]} = $strand;
				} else {
					my ($c) = $coords->[$i]."__0";
					$cdsblock_comp = $coords->[$i-1]."_".$sequence_id ."+". $coords->[$i]."_".$sequence_id ."+". $coords->[$i+1]."_".$sequence_id;
					$cdsblock_coords->{$cdsblock_comp}->{$c} = $strand;
					my ($c) = $coords->[$i]."__1";
					$cdsblock_comp = $coords->[$i-1]."_".$sequence_id ."+". $coords->[$i]."_".$sequence_id;
					$cdsblock_coords->{$cdsblock_comp}->{$c} = $strand;
					my ($c) = $coords->[$i]."__2";
					$cdsblock_comp = $coords->[$i]."_".$sequence_id ."+". $coords->[$i+1]."_".$sequence_id;
					$cdsblock_coords->{$cdsblock_comp}->{$c} = $strand;
				}
			}		
		}
		elsif ( exists $cds_coords->{'coord'} and scalar(@{$cds_coords->{'coord'}}) == 1 )
		{
			my ($coords) = $cds_coords->{'coord'};
			my ($strand) = $cds_coords->{'strand'};
			my ($cdsblock_comp) = $coords->[0]."_".$sequence_id;
			$cdsblock_coords->{$cdsblock_comp}->{$coords->[0]} = $strand;
		}
    }

    return ($trans_cds_coords, $cdsblock_coords);	
}

# Create fasta file with the joined CDS for the analysis by exons
sub _create_block_cdsseq_file($$$)
{
	my ($cdsblock_coords, $input_file, $output_file) = @_;
	my ($report);
	my ($cdsblock_report);
	my ($cdsblock_seqs) = '';
	my ($dirname,$basename) = parse_file($output_file);
	my ($outfile) = $dirname."/cdsblock.fa";

	# Get the sequence for CDS
	if (-e $input_file and (-s $input_file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $input_file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() )
		{
			my (@attr) = split(/\|/, $seq->id);
			if ( scalar(@attr) >= 10 )
			{
				my ($trans_id) = $attr[0];
				my ($cds_chr) = $attr[5];
				my ($cds_start) = $attr[6];
				my ($cds_end) = $attr[7];
				my ($cds_phase) = $attr[9];
				my ($cds_index) = $cds_chr."-".$cds_start.":".$cds_end.":".$cds_phase."_".$trans_id;
				$report->{$cds_index} = $seq->seq;
			}
		}	
	}
	$logger->debug("CDS Report ---------------\n".Dumper($report));	

	# Create a report with the CDSblocks
    while(my ($cdsblock_index, $cdsblock_comps) = each(%{$cdsblock_coords}))
    {
		my ($cdsblock_seq) = '';
		my ($trans_list);
		foreach my $cds_index ( split(/\+/, $cdsblock_index) )
		{				
			if ( exists $report->{$cds_index} ) {
				if ( $cds_index =~ /\_([^\$]*)$/ ) { $trans_list->{$1} = "" }
				$cdsblock_seq .= $report->{$cds_index};
			}
		}
		if ( $cdsblock_seq ne '' ) {
			my ($cdsblock_comp) = $cdsblock_index;
			$cdsblock_comp =~ s/\_(\w*)//g;
			while(my ($cdsblock_cds_coord, $cdsblock_cds_rep) = each(%{$cdsblock_comps}))
			{
				# my ($cdsblock_cds_coord) = join("", keys($cdsblock_comps));
				my ($cdsblock_trans) = join("", keys($trans_list));
				$cdsblock_report->{$cdsblock_seq}->{$cdsblock_comp}->{'coord'}->{$cdsblock_cds_coord} = "";
				$cdsblock_report->{$cdsblock_seq}->{$cdsblock_comp}->{'trans'}->{$cdsblock_trans} = "";
			}
		}
	}
	$logger->debug("CDSBLOCK Report ---------------\n".Dumper($cdsblock_report));	

	# Create the FASTA sequences with the CDSblocks
    while(my ($cdsblock_seq, $cdsblock_rep) = each(%{$cdsblock_report}))
    {
		while(my ($cdsblock_comp, $cdsblock_c) = each(%{$cdsblock_rep}))
		{
			my ($cdsblock_trns) = join("+", keys($cdsblock_c->{'trans'}));
			while(my ($cdsblock_cds_id, $cdsblock_cds_rep) = each(%{$cdsblock_c->{'coord'}}))
			{
				# my ($cdsblock_cds_id) = join("+", keys($cdsblock_c->{'coord'}));
				my ($block_cds_com) = $cdsblock_cds_id."_".$cdsblock_trns."|".$cdsblock_comp;
				$cdsblock_seqs .= ">".$block_cds_com."\n".$cdsblock_seq."\n";
			}
		}
	}

	# Print FASTA sequences into a file
	my ($print_out) = printStringIntoFile($cdsblock_seqs, $outfile);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}	
	return ($outfile);
}

sub parse_blast($$$)				# reads headers for each alignment in the blast output
{
	my ($blast_file,$exons,$sequence_length) = @_;
	my ($species_score, $species_report, $exon_species_report) = (0, undef, undef);
	my ($species_found);
	my ($aln_report);

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
			$faalen =~ s/\(//;
			if ($faalen != $sequence_length)
				{ $logger->warning("Length of Blast's query is different than input sequence\n"); }
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
			$length_diff = abs($length - $faalen); # difference in length between query and subject
			my @identities = split " ", $iden[0];
			my $identity = $iden[0]/$iden[1]*100;
			my $species_firstname = ( $species =~ /^([^\s]*)/ ) ? $1 : ""; #get the name from RefSeq db and UniProt db

			# if ($identity < 50)				# gets no points for this sequence
			# 	{ }
			# elsif (exists $species_found->{$species}) # gets no points for this sequence
			if (exists $species_found->{$species}) # gets no points for this sequence
				{ }
			elsif ( exists $PRIMATES->{$species_firstname} ) # only we accept species out of primates
				{ }
			else
			{
				my ($aln_score,$aln_sms) = check_alignment($length,$species,$faalen,$exons,$aln_report, $_);
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
					$species_score = $divtime_score if ( $divtime_score > $species_score );
				}
			}
		}
	}

	return ($species_score, $species_report, $exon_species_report);
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
	my @candidate = ();
	my @startc = ();
	my @endc = ();
	my $longestunmatches = 0;
	my ($aln_score) = 0;
	my ($aln_sms) = '';
	my ($specie_point) = 1;
	
	my $longest_elem = sub {
	    my $max = -1;
	    for (@_) {
	        if (length > $max) {  # no temp variable, length() twice is faster
	            $max = length;
	        }
	    }
	    return $max;
	};
	
	while (<BLASTFILE>)
		{
		chomp;
		if (/^Query:/)						# read in query line
			{
			my @query = split " ";
			push @target, $query[2];
			push @startq, $query[1];
			push @endq, $query[3];
			}
		if (/^Sbjct:/)						# read in subject line
			{
			my @subject = split " ";
			push @candidate, $subject[2];
			push @startc, $subject[1];
			push @endc, $subject[3];
			}
		if (/^\s+/)						# read in match line
		{
			my $aln = $_;
			$aln =~ s/^\s+//g;
			my @unmatches = $aln =~ /[\s\+]{4,}/g;
			if (scalar(@unmatches) > 0) {
				my $unmatches = $longest_elem->(@unmatches);
				if ( defined $unmatches and $unmatches > $longestunmatches )
				  { $longestunmatches = $unmatches } 
			}
		}
		if ($_ eq $oldinput)					# two carriage returns in a row mean alignment has ended
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
		
	if ($targstart > 4)		# reject if different N-terminal, only for the Query
		{return (0,"It has different N-terminal")}
	
	# if ( (abs($targlength - $targend) > 4) ) # reject if subject has longer C-terminal, only for the Query
	# 	{return (0,"Subject has longer C-terminal")}

	# if ( $longestunmatches > 4 ) # reject if subject has longer unmatches
	# 	{return (0,"Subject has longer unmatches")}

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
							if ( $gapresconttarg > 4 ) { return (0,"Long gap in target") }			
						}
					}
					if ($candidate[$res] eq "-") {
						$gapres++; $gapcontcand = 'true';
						if ( $gapcontcand eq 'true' ) {
							$gaprescontcand++;
							if ( $gaprescontcand > 4 ) { return (0,"Long gap in candidate") }			
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
		if ($gapres > 4) { # reject if exons have substantial gaps
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


main();

__END__

=head1 NAME

corsair_alt

=head1 DESCRIPTION

Run CORSAIR program

=head1 VERSION

0.1

=head2 Required arguments:

	--conf <Config file>

	--input <Fasta sequence file with the sequence separated by CDS coordinates>
	
	--output <Annotation output file>

=head2 Optional arguments:
	
	--appris <Flag that enables the output for APPRIS (default: NONE)>
	
	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
	

=head1 EXAMPLE

perl corsair.pl

	--conf=pipeline.ini

	--gff=OTTHUMG00000020713/annto.gtf
	
	--input=OTTHUMG00000020713/cdsseq.fa
	
	--output=OTTHUMG00000020713/corsair_alt

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut