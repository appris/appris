#!/usr/bin/perl

use strict;
use warnings::register;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Config::IniFiles;
use Data::Dumper;
use POSIX qw(ceil);

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile prepare_workspace );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$WSPACE_TMP
	$WSPACE_CACHE
	$CACHE_FLAG
	$RUN_PROGRAM
	$PROG_DB
	$PROG_DB_DIR
	$PROG_EVALUE
	$APPRIS_CUTOFF
	$MIN_LENGTH_CDS
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
	'gff=s'				=> \$gff_file,
	'input=s'			=> \$input_file,
	'output=s'			=> \$output_file,
	'appris'			=> \$appris,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $config_file and defined $gff_file and defined $input_file and defined $output_file )
{
    print `perldoc $0`;
    exit 1;
}

# Get conf vars
my ($cfg) 			= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$CACHE_FLAG			= $cfg->val('MATADOR3D_VARS', 'cache');
$RUN_PROGRAM		= $cfg->val('MATADOR3D_VARS', 'program');
$PROG_DB			= $cfg->val('MATADOR3D_VARS', 'db');
$PROG_DB_DIR		= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$PROG_DB;
$PROG_EVALUE		= $cfg->val('MATADOR3D_VARS', 'evalue');
$APPRIS_CUTOFF		= $cfg->val('MATADOR3D_VARS', 'cutoff');
$MIN_LENGTH_CDS		= 6;

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

my $EXP_CFG = new Config::IniFiles( -file => $ENV{APPRIS_EXP_CONF_FILE} );
my $calc_frames = $EXP_CFG->val( 'matador3d', 'calc_frames', 0 );
my $best_hsp_only = $EXP_CFG->val( 'matador3d', 'best_hsp_only', 0 );
my $count_terminal_gaps = $EXP_CFG->val( 'matador3d', 'count_terminal_gaps', 0 );
my $recal_align_scores = $EXP_CFG->val( 'matador3d', 'recal_align_scores', 0 );

#####################
# Method prototypes #
#####################
sub _parse_blast($$);
sub _check_alignment($$$);
sub _run_blastpgp($$);
sub _get_best_cds_score($);
sub _get_biggest_cds($$$);
sub _get_biggest_mini_cds($$$);
sub _get_record_annotations($);
sub _get_appris_annotations($);
sub _get_cds_coordinates_from_gff($$);
sub _get_peptide_coordinate_with_frames($$$$$$);
sub _get_peptide_coordinate_sans_frames($$$$$$);
sub _get_mini_peptide_coordinate_with_frames($$$$$$$);
sub _get_mini_peptide_coordinate_sans_frames($$$$$$$);
sub _get_init_report($$$$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Get the CDS coordinates info from gff file ---------------
	$logger->info("##Get cds coordinates from gff ---------------\n");
	my ($trans_cds_coords, $mini_cds_interval_pool) = _get_cds_coordinates_from_gff($input_file, $gff_file);
    $logger->debug("##CDS coordenates ---------------\n".Dumper($trans_cds_coords));
    $logger->debug("##Mini-CDS interval pool ---------------\n".Dumper($mini_cds_interval_pool));

	# For every sequence run the method ---------------
    my ($transcript_report);
    my $fasta_object = Bio::SeqIO->new(
                        -file => $input_file,
                        -format => 'Fasta'
    );
	while ( my $seq = $fasta_object->next_seq() )
	{
		if($seq->id=~/^([^|]*)\|([^|]*)/)
		{			
			my ($sequence_id) = $2;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
            my ($sequence) = $seq->seq;
            
            $logger->info("\n##$sequence_id ###############################\n");


			# get the intersection of CDS coordinates creating "mini-cds"
			my ($cds_reports) = _get_init_report($sequence_id, $sequence, $trans_cds_coords, $mini_cds_interval_pool);
			$logger->debug("##Init CDS report ---------------\n".Dumper($cds_reports));

			if(defined $cds_reports)
			{
                # Init transcript report from cds list
				$transcript_report->{$sequence_id}->{'score'} = 0;
				#$transcript_report->{$sequence_id}->{'seq'} = $sequence;
				$transcript_report->{$sequence_id}->{'cds'} = $cds_reports;
				
				# Run blast
				$logger->info("\n##Running Blast ---------------\n");
				my ($blast_sequence_file) = _run_blastpgp($sequence_id, $sequence);


                # Parse blast
				$logger->info("\n##Parsing Blast ---------------\n");
                my ($score,$pdb_cds_reports) = _parse_blast($blast_sequence_file, $cds_reports);
				$transcript_report->{$sequence_id}->{'score'} = $score;
				$transcript_report->{$sequence_id}->{'cds'} = $pdb_cds_reports if (defined $pdb_cds_reports);
			}
			else
			{
				$transcript_report->{$sequence_id}->{'score'} = 0;
				#$logger->error("peptide $sequence_id does not have exon info\n");
			}
		}
    }
	$logger->debug("##CDS scores ---------------\n".Dumper($transcript_report)."\n");

	# Get the score of the CDS from the best score for each transcript ---------------
	_get_best_cds_score($transcript_report);
	$logger->debug("##The Best CDS scores ---------------\n".Dumper($transcript_report)."\n");


	# Get records by transcript ---------------
	my ($output_content) = _get_record_annotations($transcript_report);


	# Get the annotations for the main isoform /* APPRIS */ ----------------
	if ( defined $appris )
	{		
		$output_content .= _get_appris_annotations($transcript_report);
	}
	
		
	# Print records by transcript ---------------
	my ($print_out) = printStringIntoFile($output_content, $output_file);
	unless( defined $print_out ) {
		$logger->error("Can not create output file: $!\n");
	}


	$logger->finish_log();
	
	exit 0;	
}

# Run blast
sub _run_blastpgp($$)
{
	my ($sequence_id, $sequence) = @_;

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
	
	
	# Create temporal file for blast
	my ($fasta_sequence_file) = $ws_cache.'/seq.faa';
	unless(-e $fasta_sequence_file and (-s $fasta_sequence_file > 0) ) # Cached fasta
	{	
		my ($fasta_sequence_content_file) = ">Query\n$sequence";
		my ($print_fasta) = printStringIntoFile($fasta_sequence_content_file, $fasta_sequence_file);
		unless( defined $print_fasta ) {
			$logger->error("Can not create temporal file:$fasta_sequence_file: $!\n");
		}
	}
	
	# Run blast
	my ($blast_sequence_file) = $ws_cache.'/seq.'.$PROG_DB;
	unless(-e $blast_sequence_file and (-s $blast_sequence_file >0) and ($CACHE_FLAG eq 'yes')) # Cached Blast
	{
		eval
		{
			# parameters to run blaspgp
			$logger->debug("$RUN_PROGRAM -d $PROG_DB_DIR -i $fasta_sequence_file -e$PROG_EVALUE -o $blast_sequence_file\n");                   
			system("$RUN_PROGRAM -d $PROG_DB_DIR -i $fasta_sequence_file -e$PROG_EVALUE -o $blast_sequence_file");			
		};
		$logger->error("Running blast\n") if($@);
	}
	return $blast_sequence_file;
} # end _run_blastpgp

# Reads headers for each alignment in the blast output
sub _parse_blast($$)
{
	my ($blast_file, $reports) = @_;
	my $pdbid = "";
	my $points = 0; 
	
	my ($in) = new Bio::SearchIO(
							-format => 'blast', 
							-file   => $blast_file
	);
	while( my $result = $in->next_result )
	{
		while( defined $result and my $hit = $result->next_hit )
		{
			$pdbid = $hit->name;
			$pdbid =~ s/PDB://;

			while( defined $hit and my $hsp = $hit->next_hsp )
			{
				# get values for current alignemnt
				my ($pdb_points, $pdb_results) = _check_alignment($pdbid, $hsp, $reports);

				# get the best values for each mini-cds
				for ( my $icds = 0; $icds < scalar(@{$reports}); $icds++ )
			    {
			    	my ($cds) = $reports->[$icds];			    	
			    	my ($cds_index) = $cds->{'index'};
			    	my ($mini_cds_list) = $cds->{'mini_cds'};

			    	if ( exists $pdb_results->{$cds_index} ) {
			    		my ($pdb_result) = $pdb_results->{$cds_index};
			    		my ($mini_pdb_results) = $pdb_result->{'mini_cds'};
						for ( my $iminicds = 0; $iminicds < scalar(@{$mini_cds_list}); $iminicds++ )
					    {
					    	my ($mini_cds) = $mini_cds_list->[$iminicds];
					    	my ($mini_cds_index) = $mini_cds->{'index'};
							if ( exists $mini_pdb_results->{$mini_cds_index} ) {
								my ($mini_pdb_result) = $mini_pdb_results->{$mini_cds_index};
								if ( $mini_pdb_result->{'score'} ne '-') {									
									if ( $reports->[$icds]->{'mini_cds'}->[$iminicds]->{'score'} eq '-') {
										$reports->[$icds]->{'mini_cds'}->[$iminicds] = $mini_pdb_result;
									}
									else {
										if ( $mini_pdb_result->{'score'} > $reports->[$icds]->{'mini_cds'}->[$iminicds]->{'score'} ) {
											$reports->[$icds]->{'mini_cds'}->[$iminicds] = $mini_pdb_result;
										}
									}								
								}
							}
							
						}			    		
			    	}
			    }
			    last if ($best_hsp_only);
			}
		}
	}
	
	# get the cds score: the sum of mini-cds scores
	for ( my $icds = 0; $icds < scalar(@{$reports}); $icds++ )
    {
    	my ($cds_report) = $reports->[$icds];
		my ($mini_cds_list) = $cds_report->{'mini_cds'};
		my ($cds_score) = 0;	
		for ( my $iminicds = 0; $iminicds < scalar(@{$mini_cds_list}); $iminicds++ )
	    {
	    	my ($mini_cds_report) = $mini_cds_list->[$iminicds];
	    	my ($mini_cds_score) = $mini_cds_report->{'score'};
	    	if ( $mini_cds_score ne '-' ) {
	    		$cds_score += $mini_cds_score;
	    	}	    	
	    }
	    $cds_report->{'score'} = $cds_score;
	    $points += $cds_score;
    }			    	
	
	return ($points, $reports);	
} # end _parse_blast

# Parses BLAST alignments exon by exon
sub _check_alignment($$$)
{
	my ($pdb_id, $hsp, $cds_list) = @_;
	my ($cds_reports);
	my ($score) = 0;
	my (@target) = ();
	my (@candidate) = ();

	# process alignment ...
	my ($pdb_identity) = sprintf("%.1f", $hsp->percent_identity);
	my ($target_str) = $hsp->query_string;
	my ($candidate_str) = $hsp->hit_string;

	my ($targstart) = $hsp->start('query');
	my ($targend) = $hsp->end('query');
	my ($candstart) = $hsp->start('hit');
	my ($candend) = $hsp->end('hit');
		
	@target = split "", $target_str;
	@candidate = split "", $candidate_str;
    
	$logger->debug("\n#PDB: $pdb_id ---------------\n");
	$logger->debug("#IDENTITY:$pdb_identity|TARGET:$targstart:$targend|CANDIDATE:$candstart:$candend\n");           

	my ($hsp_idx_min) = 0;  # 0-offset HSP index at start of each mini-CDS alignment
	foreach my $cds (@{$cds_list})
	{
    	my ($cds_index) = $cds->{'index'};
    	my ($cds_coord) = $cds->{'coord'};
		my ($cds_frame) = $cds->{'frame'};
    	my ($mini_cds_list) = $cds->{'mini_cds'};
    	my ($mini_cds_score_list);
		my ($mini_score) = 0;    	
		$logger->debug("#cds_trans:$cds_coord\[$cds_index\]\n");

		foreach my $mini_cds (@{$mini_cds_list})
		{
	    	my ($mini_cds_index) = $mini_cds->{'index'};
	    	my ($mini_cds_seq) = $mini_cds->{'seq'};	    	
	    	my ($mini_cds_coord) = $mini_cds->{'coord'};
			my ($mini_cds_frame) = $mini_cds->{'frame'};
			my @boundaries = split ":", $mini_cds_index;
			my $mini_cds_start = $boundaries[0];
			my $mini_cds_end = $boundaries[1];
			$logger->debug("#mini_cds_trans:$mini_cds_coord\[$mini_cds_index\]\n");

			# Finish when target-candidate don't cover exons or
			# Next until exon is within target-candidate
			if ( ($mini_cds_start > $targend) or ($targstart > $mini_cds_end) )
			{
				$mini_cds_score_list->{$mini_cds_index} = {
							'index'		=> $mini_cds_index,
							'coord'		=> $mini_cds_coord,
							'frame'		=> $mini_cds_frame,
							'score'		=> '-',
							'seq'		=> $mini_cds_seq 
				};
				next;
			}

			my $gapres = 0;			
			my $identities = 0;
			my $pep_pos;
			my $hsp_idx;
			my $aln_len = 0;
			my $targ_gapres = 0;
	        
			# Get the minimum 1-offset peptide position common to the HSP and mini-CDS.
			my ($pep_pos_min);
			if ($mini_cds_start < $targstart)
			{
				$pep_pos_min = $targstart;

				if ($count_terminal_gaps)
				{
					# Count leading gaps in a notional alignment with 100% query coverage.
					$gapres += ($targstart - $mini_cds_start) + 1;
				}
			}
			else
			{
				$pep_pos_min = $mini_cds_start;
			}

			# Get the maximum 1-offset peptide position common to the HSP and mini-CDS.
			my ($pep_pos_max);
			if ($targend < $mini_cds_end)
			{
				$pep_pos_max = $targend;

				if ($count_terminal_gaps)
				{
					# Count trailing gaps in a notional alignment with 100% query coverage.
					$gapres += ($mini_cds_end - $targend) + 1;
				}
			}
			else
			{
				$pep_pos_max = $mini_cds_end;
			}

			# Count and classify the residues
			for ($pep_pos=$pep_pos_min, $hsp_idx=$hsp_idx_min ; $pep_pos<=$pep_pos_max ; $pep_pos++, $hsp_idx++)
			{
				if(defined $target[$hsp_idx] and defined $candidate[$hsp_idx])
				{
					if ($target[$hsp_idx] eq $candidate[$hsp_idx])
					{
						$identities++;
					}
					elsif ($target[$hsp_idx] eq "-")
					{
						$gapres++;
						$pep_pos--;
						$targ_gapres++;
					}
					elsif ($candidate[$hsp_idx] eq "-")
					{
						$gapres++;
					}
					$aln_len++;
				}
				else
				{
					$logger->error("invalid HSP index: $hsp_idx");
				}
			}

			$hsp_idx_min = $hsp_idx;
			my ($align_start) = $pep_pos_min;
			my ($align_end) = $pep_pos_min + $aln_len - 1;

			$logger->debug("#hsp_idx_min: $hsp_idx_min pep_pos_min: $pep_pos_min pep_pos_max: $pep_pos_max hsp_idx: $hsp_idx pep_pos: $pep_pos aln_len: $aln_len\n");

			$logger->debug("#mini_aligns_coord\[$align_start:$align_end\]\n");

			# Skip mini-CDS shorter than the minimum length, but only after having
			# advanced the HSP index to the start of the subsequent mini-CDS.
			if ( (($mini_cds_end - $mini_cds_start) + 1) < $MIN_LENGTH_CDS)
			{
				$mini_cds_score_list->{$mini_cds_index} = {
					'index'		=> $mini_cds_index,
					'coord'		=> $mini_cds_coord,
					'frame'		=> $mini_cds_frame,
					'score'		=> '-',
					'seq'		=> $mini_cds_seq
				};
				next;
			}

			my $mini_cds_residues = $mini_cds_end - $mini_cds_start + 1;
			my $align_residues = $align_end - $align_start + 1;
	
			# Value of identity score: the sort of conditions is important!!!
			my $totalidentity = 0;
			my $identity = 0;
			$identity = $identities/$mini_cds_residues*100;
			$identity = sprintf("%.2f", $identity);
			my $align_identity = 0;			
			$align_identity = $identities/$align_residues*100 if ( $align_residues != 0);
			$align_identity = sprintf("%.2f", $align_identity);
			if ( $recal_align_scores ) {
				if ($identity >= 40)
					{$totalidentity = 1}
				elsif ($identity >= 30)
					{$totalidentity = 0.80}
				elsif ($identity >= 20)
					{$totalidentity = 0.60}
				elsif ($identity >= 15)
					{$totalidentity = 0.40}
				elsif ($identity >= 10)
					{$totalidentity = 0.20}
				else
					{$totalidentity = 0}
			} else {
				if ($identity >= 50)
					{$totalidentity = 1}
				elsif ($identity >= 40)
					{$totalidentity = 0.80}
				elsif ($identity >= 30)
					{$totalidentity = 0.60}
				elsif ($identity >= 25)
					{$totalidentity = 0.50}
				elsif ($identity >= 20)
					{$totalidentity = 0.40}
				elsif ($identity >= 15)
					{$totalidentity = 0.20}
				else
					{$totalidentity = 0}
			}
			$logger->debug("\tIdentity: $identity = $identities/$mini_cds_residues * 100\n");
			$logger->debug("\tAlign_identity: $align_identity = $identities/$align_residues * 100\n");           
			$logger->debug("\tTotal identity: $totalidentity\n");
	
			# Value of gap score: the sort of conditions is important!!!
			my $totalgaps = 0;
			my $gaps = 0;
			$gaps = ( $gapres / ($mini_cds_residues + $targ_gapres) ) *100;
			if ($gaps <= 3)
				{$totalgaps = 1}
			elsif ($gaps <= 6)
				{$totalgaps = 0.80}
			elsif ($gaps <= 10)
				{$totalgaps = 0.50}
			elsif ($gaps <= 15)
				{$totalgaps = 0.33}
			elsif ($gaps <= 20)
				{$totalgaps = 0.20}
			elsif ($gaps <= 25)
				{$totalgaps = 0}
			elsif ($gaps > 50)
				{$totalgaps = -1}
			elsif ($gaps > 40)
				{$totalgaps = -0.75}
			elsif ($gaps > 30)
				{$totalgaps = -0.5}
			elsif ($gaps > 25)
				{$totalgaps = -0.25}
			else
				{$totalgaps = 0}			
			$logger->debug("\tGaps: $gaps = $gapres/$mini_cds_residues*100\n");
			$logger->debug("\tTotal gaps: $totalgaps\n");
	
			# Value of total residues: the sort of conditions is important!!!
			my $totalres = 0;
			if ($mini_cds_residues >= 30)
				{$totalres = 1}
			elsif ($mini_cds_residues >= 20)
				{$totalres = 0.75}	
			elsif ($mini_cds_residues >= 12)
				{$totalres = 0.50}
			elsif ($mini_cds_residues >= 6)
				{$totalres = 0.33}
			else
				{$totalres = 0}
			$logger->debug("\tTotal residues: $totalres\n");
	
			# Exon score
			my ($escore) = $totalidentity*$totalgaps*$totalres;
			if ( $totalgaps < 0 ) { $escore = $totalgaps*$totalres }
			$logger->debug("\tScore: $escore\n");

			my ($mini_cds_report) = {
						'index'			=> $mini_cds_index,
						'coord'			=> $mini_cds_coord,
						'frame'			=> $mini_cds_frame,
						'seq'			=> $mini_cds_seq,
						'score'			=> $escore,
						'score_iden'	=> $totalidentity,
						'score_gaps'	=> $totalgaps,
						'score_res'		=> $totalres,
						'score_desc'	=> "$totalidentity*$totalgaps*$totalres",
						'pdb'			=> $pdb_id,
						'pdb_identity'	=> $identity,
						'alignment_start' => $align_start,
						'alignment_end'   => $align_end,
						'alignment_iden'   => $align_identity
			};
			$mini_cds_score_list->{$mini_cds_index} = $mini_cds_report;
			$mini_score += $escore;
		}
		
		my ($cds_report) = {
					'index'					=> $cds_index,
					'coord'					=> $cds_coord,
					'frame'					=> $cds_frame,
					'score'					=> $mini_score,
					'pdb_identity'			=> $pdb_identity,					
					'mini_cds'				=> $mini_cds_score_list
		};
		$cds_reports->{$cds_index} = $cds_report;
		$score += $mini_score;		
    }
	return ($score,$cds_reports);
		
} # end _check_alignment

# Get the better score for each CDS of transcript
sub _get_best_cds_score($)
{
	my ($transcript_report) = @_;

	my (%aux_transcript_report) = %{$transcript_report};
	while(my ($sequence_id,$trans_report) = each(%{$transcript_report}))
	{
		if(	exists $trans_report->{'score'} and defined $trans_report->{'score'} and
			exists $trans_report->{'cds'} and defined $trans_report->{'cds'}
		)
		{
			my ($pdb_cds_score) = 0;

			foreach my $pdb_cds_report (@{$trans_report->{'cds'}})
			{
				# cds must be bigger than 6 residues
				my (@cds_pep_coord) = split(':', $pdb_cds_report->{'index'});		
				if ( ( ($cds_pep_coord[1] - $cds_pep_coord[0]) + 1 ) >= $MIN_LENGTH_CDS )
				{
					$logger->debug("#_get_biggest_score: $sequence_id\n");
					$pdb_cds_score += _get_biggest_score($sequence_id, $pdb_cds_report, \%aux_transcript_report);
				}
			}
			$trans_report->{'score'} = $pdb_cds_score;
		}
	}		
} # end _get_best_cds_score

sub _get_biggest_score($$$)
{
	my ($sequence_id, $pdb_cds_report, $transcript_report) = @_;

	my ($pdb_cds_score) = 0;
	
	# if exits, get the biggest score from the "mini" cds
	if (exists $pdb_cds_report->{'mini_cds'} and defined $pdb_cds_report->{'mini_cds'})
	{
		my ($mini_pdb_cds_score) = 0;
	
		foreach my $mini_pdb_cds_report (@{$pdb_cds_report->{'mini_cds'}})
		{							
			# mini cds must be bigger than 6 residues
			my (@mini_cds_pep_coord) = split(':', $mini_pdb_cds_report->{'index'});
			if ( ( ($mini_cds_pep_coord[1] - $mini_cds_pep_coord[0])+1 ) >= $MIN_LENGTH_CDS )
			{
				$logger->debug("#_get_biggest_mini_cds:\n");
				my ($external_seq_id, $big_mini_pdb_cds_report) = _get_biggest_mini_cds($sequence_id, $mini_pdb_cds_report, $transcript_report);
				if(defined $big_mini_pdb_cds_report)
				{
					if (defined $external_seq_id)
						{ $mini_pdb_cds_report->{'external_id'} = $external_seq_id; }
	
					$mini_pdb_cds_report->{'score'} = $big_mini_pdb_cds_report->{'score'};
					$mini_pdb_cds_report->{'score_res'} = $big_mini_pdb_cds_report->{'score_res'};
					$mini_pdb_cds_report->{'score_gaps'} = $big_mini_pdb_cds_report->{'score_gaps'};
					$mini_pdb_cds_report->{'score_iden'} = $big_mini_pdb_cds_report->{'score_iden'};
					$mini_pdb_cds_report->{'score_desc'} = $big_mini_pdb_cds_report->{'score_desc'};
					$mini_pdb_cds_report->{'pdb'} = $big_mini_pdb_cds_report->{'pdb'};
					$mini_pdb_cds_report->{'pdb_identity'} = $big_mini_pdb_cds_report->{'pdb_identity'};
					$mini_pdb_cds_report->{'alignment_iden'} = $big_mini_pdb_cds_report->{'alignment_iden'};
					
					# In the case, a mini-exon had not alignment (no coord) and there is an external alignment => we add the mini-cds coordinates to the alignment
					unless ( exists $mini_pdb_cds_report->{'alignment_start'} or $mini_pdb_cds_report->{'alignment_end'} ) {
						$mini_pdb_cds_report->{'alignment_start'} = $mini_cds_pep_coord[0];
						$mini_pdb_cds_report->{'alignment_end'} = $mini_cds_pep_coord[1];
					}
					
					if ( defined $big_mini_pdb_cds_report->{'score'} and ($big_mini_pdb_cds_report->{'score'} ne '-') ) {
						$pdb_cds_score += $big_mini_pdb_cds_report->{'score'};
						$mini_pdb_cds_score += $big_mini_pdb_cds_report->{'score'};
					}					
				} else {
					if ( defined $mini_pdb_cds_report->{'score'} and ($mini_pdb_cds_report->{'score'} ne '-') ) {
						$pdb_cds_score += $mini_pdb_cds_report->{'score'};
						$mini_pdb_cds_score += $mini_pdb_cds_report->{'score'};						
					}
				}				
			}
		}
		$pdb_cds_report->{'score'} = $mini_pdb_cds_score;					
	}
	else
	{
		# otherwise, get the biggest score from the cds
		$logger->debug("#_get_biggest_cds:\n");
		my ($external_seq_id, $big_pdb_cds_report) = _get_biggest_cds($sequence_id, $pdb_cds_report, $transcript_report);
						
		# first, get the biggest score from the cds 
		if ( defined $big_pdb_cds_report )
		{
			if (defined $external_seq_id)
				{ $pdb_cds_report->{'external_id'} = $external_seq_id; }
			
			$pdb_cds_report->{'score'} = $big_pdb_cds_report->{'score'};
			$pdb_cds_score += $big_pdb_cds_report->{'score'};						
		}
	}
	return $pdb_cds_score;
} # end _get_biggest_score

# Get the biggest score for every CDS
sub _get_biggest_cds($$$)
{
	my ($main_seq_id, $main_pdb_cds_report, $transcript_report) = @_;

	my ($pdb_cds_report);
	my ($external_id);
	my ($biggest_pdb_cds_score) = $main_pdb_cds_report->{'score'};
	my ($biggest_pdb_cds_coord) = $main_pdb_cds_report->{'coord'};
	
	# If we found one CDS coming from the same exon and it had bigger score than 
	# the given cds, we get it
	while(my ($sequence_id,$trans_report) = each(%{$transcript_report}))
	{
		next if ($main_seq_id eq $sequence_id); # Jump for the same sequence

		if(exists $trans_report->{'cds'} and defined $trans_report->{'cds'})
		{
			my ($trans_cds_list) = $trans_report->{'cds'};
			foreach my $trans_pdb_cds_report (@{$trans_cds_list})
			{
				$logger->debug("(".$trans_pdb_cds_report->{'coord'}." eq ".$biggest_pdb_cds_coord.") and (".$trans_pdb_cds_report->{'score'}." > ".$biggest_pdb_cds_score.") :\n");
				if(
					($trans_pdb_cds_report->{'coord'} eq $biggest_pdb_cds_coord ) and
					($trans_pdb_cds_report->{'score'} > $biggest_pdb_cds_score) 
				){
					$logger->debug("Into:_get_biggest_cds\n");
					$external_id = $sequence_id;
					$biggest_pdb_cds_score = $trans_pdb_cds_report->{'score'};
					$pdb_cds_report = $trans_pdb_cds_report;
				}
			}
		}
	}	
	return ($external_id,$pdb_cds_report);
} # end _get_biggest_cds

# Get the biggest score for every mini CDS
sub _get_biggest_mini_cds($$$)
{
	my ($main_seq_id, $main_mini_pdb_cds_report, $transcript_report) = @_;

	my ($mini_pdb_cds_report);
	my ($external_id);
	my ($biggest_mini_pdb_cds_score) = $main_mini_pdb_cds_report->{'score'};
	my ($biggest_mini_pdb_cds_coord) = $main_mini_pdb_cds_report->{'coord'};
	my ($biggest_mini_pdb_cds_pdb);
	if ( $calc_frames ) {
		$biggest_mini_pdb_cds_pdb = ( exists $main_mini_pdb_cds_report->{'pdb'} and defined $main_mini_pdb_cds_report->{'pdb'} ) ? $main_mini_pdb_cds_report->{'pdb'} : '' ;
	} else {
		$biggest_mini_pdb_cds_pdb = defined $main_mini_pdb_cds_report->{'pdb'} ? $main_mini_pdb_cds_report->{'pdb'} : '' ;
	}
	my ($biggest_mini_pdb_cds_frame) = ( exists $main_mini_pdb_cds_report->{'frame'} ) ? $main_mini_pdb_cds_report->{'frame'} : "-";
	$biggest_mini_pdb_cds_score = 0 if ( !defined $biggest_mini_pdb_cds_score or ($biggest_mini_pdb_cds_score eq '-') );

	# If we found one CDS coming from the same exon and it had bigger score than 
	# the given cds, we get it
	while(my ($sequence_id,$trans_report) = each(%{$transcript_report}))
	{
		next if ($main_seq_id eq $sequence_id); # Jump for the same sequence
				
		if(exists $trans_report->{'cds'} and defined $trans_report->{'cds'})
		{
			foreach my $trans_pdb_cds_report (@{$trans_report->{'cds'}})
			{
				# Check if the frame between the compared isoforms are equal
				my ($trans_cds_frame) = $trans_pdb_cds_report->{'frame'};
				if ( $calc_frames || $biggest_mini_pdb_cds_frame eq $trans_cds_frame ) {
					my ($trans_mini_cds_list) = $trans_pdb_cds_report->{'mini_cds'};	
					foreach my $trans_mini_cds (@{$trans_mini_cds_list})			
					{
						my $trans_mini_cds_score;
						my $trans_mini_cds_pdb;
						if ( $calc_frames ) {
							if ( ! defined $trans_mini_cds->{'score'} || $trans_mini_cds->{'score'} eq '-' ) {
								next;
							}
							if ( ! exists $trans_mini_cds->{'frame'} || ! defined $trans_mini_cds->{'frame'}
									|| $trans_mini_cds->{'frame'} != $biggest_mini_pdb_cds_frame ) {
								next;
							}
							$trans_mini_cds_score = $trans_mini_cds->{'score'};
							$trans_mini_cds_pdb = ( exists $trans_mini_cds->{'pdb'} and defined $trans_mini_cds->{'pdb'}
																      and $trans_mini_cds->{'pdb'} ne '' ) ? $trans_mini_cds->{'pdb'} : '';
						} else {
						  $trans_mini_cds_score = ( $trans_mini_cds->{'score'} ne '-' ) ? $trans_mini_cds->{'score'} : 0;
							$trans_mini_cds_pdb = ( exists $trans_mini_cds->{'pdb'} and $trans_mini_cds->{'pdb'} ne '' ) ? $trans_mini_cds->{'pdb'} : '';
						}

						my ($trans_mini_cds_index) = $trans_mini_cds->{'index'};
						my ($trans_mini_cds_seq) = $trans_mini_cds->{'seq'};
						my (@trans_mini_cds) = split(':',$trans_mini_cds->{'coord'}); # Sticking paster for reverse strand
						my ($trans_mini_cds_coord_aux) = $trans_mini_cds[1].':'.$trans_mini_cds[0];
						if(
							(($trans_mini_cds->{'coord'} eq $biggest_mini_pdb_cds_coord) or ($trans_mini_cds_coord_aux eq $biggest_mini_pdb_cds_coord))
							and 
							($trans_mini_cds_score > $biggest_mini_pdb_cds_score)
						){
							$logger->debug("Into:_get_biggest_mini_cds: $sequence_id [".$trans_mini_cds->{'coord'}."] - ".$trans_mini_cds->{'index'}." > ".$trans_mini_cds->{'score'}."\n");
							$external_id = $sequence_id;
							$biggest_mini_pdb_cds_score = $trans_mini_cds->{'score'};
							$biggest_mini_pdb_cds_score = 0 if ( !defined $biggest_mini_pdb_cds_score or ($biggest_mini_pdb_cds_score eq '-') );
							$mini_pdb_cds_report = $trans_mini_cds;
						}
					}
				}					
			}
		}
	}
	return ($external_id,$mini_pdb_cds_report);
} # end _get_biggest_mini_cds

# Get records by transcript
sub _get_record_annotations($)
{
	my ($transcript_report) = @_;
	my ($output_content) = '';
	
	while (my ($sequence_id, $trans_report) = each(%{$transcript_report}) )
	{
		my ($score) = 0;
		if(exists $trans_report->{'score'} and defined $trans_report->{'score'})
			{ $score=$trans_report->{'score'}; }
		$output_content.=">".$sequence_id."\t".$score."\n";

		if(exists $trans_report->{'cds'} and defined $trans_report->{'cds'}
		)
		{
			my ($pdb_cds_reports) = $trans_report->{'cds'};
			for ( my $cds_order = 0; $cds_order < scalar(@{$pdb_cds_reports}); $cds_order++ )
			{
				my ($cds_index) = $pdb_cds_reports->[$cds_order]->{'index'}; # original peptide coordinates of CDS
				my ($cds_report) = $pdb_cds_reports->[$cds_order];
				
				$output_content.=	'- '.$cds_index.
										'['.($cds_order+1).']'."\t".
											$cds_report->{'score'};
				$output_content.= "\n";
				
				if(exists $cds_report->{'mini_cds'} and defined $cds_report->{'mini_cds'})
				{
					my ($mini_cds_score_list) = $cds_report->{'mini_cds'};
					for(my $index=0; $index<scalar(@{$mini_cds_score_list});$index++)
					{

						my ($mini_cds_index) = $pdb_cds_reports->[$cds_order]->{'mini_cds'}->[$index]->{'index'}; # original peptide coordinates of CDS
						my ($mini_cds_report) = $mini_cds_score_list->[$index];
						
						my ($mini_cds_align) = '';
						if (exists $mini_cds_report->{'alignment_start'} and defined $mini_cds_report->{'alignment_start'} and ($mini_cds_report->{'alignment_start'} != 0) and 
							exists $mini_cds_report->{'alignment_end'} and defined $mini_cds_report->{'alignment_end'} and ($mini_cds_report->{'alignment_end'} != 0))
							{ $mini_cds_align = "[".$mini_cds_report->{'alignment_start'}.":".$mini_cds_report->{'alignment_end'}."]"; }
																				
						my ($mini_cds_score_desc) = '';
						if (exists $mini_cds_report->{'score_desc'} and defined $mini_cds_report->{'score_desc'} and ($mini_cds_report->{'score_desc'} ne '') )
							{ $mini_cds_score_desc = "[".$mini_cds_report->{'score_desc'}."]"; }
													
						my ($mini_cds_pdb) = '';
						if (exists $mini_cds_report->{'pdb'} and defined $mini_cds_report->{'pdb'} and ($mini_cds_report->{'pdb'} ne '') )
							{ $mini_cds_pdb = "\t".$mini_cds_report->{'pdb'}; }
							
						#my ($mini_cds_pdb_iden) = '';
						#if (exists $mini_cds_report->{'pdb_identity'} and defined $mini_cds_report->{'pdb_identity'} and ($mini_cds_report->{'pdb_identity'} ne '') )
						#	{ $mini_cds_pdb_iden = "[".$mini_cds_report->{'pdb_identity'}."]"; }
							
						my ($mini_cds_align_iden) = '';
						if (exists $mini_cds_report->{'alignment_iden'} and defined $mini_cds_report->{'alignment_iden'} and ($mini_cds_report->{'alignment_iden'} ne '') )
							{ $mini_cds_align_iden = "[".$mini_cds_report->{'alignment_iden'}."]"; }
														
						my ($mini_external_id) = '';						
						if (exists $mini_cds_report->{'external_id'} and defined $mini_cds_report->{'external_id'} and ($mini_cds_report->{'external_id'} ne '') )
							{ $mini_external_id = "\t".$mini_cds_report->{'external_id'}; }
							
						$output_content.=	"\t".$mini_cds_index.$mini_cds_align.
											"\t".$mini_cds_report->{'score'}.$mini_cds_score_desc.
											$mini_cds_pdb.$mini_cds_align_iden.
											$mini_external_id;
						$output_content.= "\n";             					
					}							
				}
			}
		}		
	}
	return $output_content;
} # end _get_record_annotations

# Get the annotations for the main isoform /* APPRIS */ ----------------
sub _get_appris_annotations($)
{
	my ($transcript_report) = @_;
	
	my ($score_transcripts);
	my (@unknow_tag);
	my (@no_tag);
	my ($output_content) = '';

	$output_content .= "\n";
	$output_content .= "# ==================================== #\n";
	$output_content .= "# Conservation of homologous structure #\n";
	$output_content .= "# ==================================== #\n";

	# get the whole list of scores
	while (my ($sequence_id, $trans_report) = each(%{$transcript_report}) )
	{
		if(exists $trans_report->{'score'} and defined $trans_report->{'score'})
		{
			push(@{$score_transcripts->{$trans_report->{'score'}}},$sequence_id);
		}
	}
	
	if ( defined $score_transcripts )
	{
		# sort by descending order
		my (@sorted_score_list) = sort { $b <=> $a } keys (%{$score_transcripts}); 

		for ( my $i=0; $i < scalar(@sorted_score_list); $i++ )
		{
			if ( $i == 0 ) # get the biggest score
			{
				foreach my $trans_id (@{$score_transcripts->{$sorted_score_list[$i]}})
					{ push(@unknow_tag,$trans_id); }
			}
			else
			{
				if( $sorted_score_list[0] - $sorted_score_list[$i] <= $APPRIS_CUTOFF )
				{
					foreach my $trans_id (@{$score_transcripts->{$sorted_score_list[$i]}})
						{ push(@unknow_tag,$trans_id); }
				}
				else
				{
					foreach my $trans_id (@{$score_transcripts->{$sorted_score_list[$i]}})
						{ push(@no_tag,$trans_id); }
				}
			}
		}
		if(scalar(@unknow_tag)==1)
		{
			foreach my $trans_id (@unknow_tag)
			{
				$output_content .= '>' . $trans_id . "\t" . 'YES' . "\n";
			}
		}
		else
		{
			foreach my $trans_id (@unknow_tag)
			{
				$output_content .= '>' . $trans_id . "\t" . 'UNKNOWN' . "\n";
			}				
		}
		foreach my $trans_id (@no_tag)
		{
				$output_content .= '>' . $trans_id . "\t" . 'NO' . "\n";
		}
	}	
	return $output_content;	
}

# Get the CDS coordinates info from gff file
sub _get_cds_coordinates_from_gff($$)
{
	my ($input_file, $gff_file) = @_;		
	my (@global_trans_cds_info);
	
	# Get the CDS info from gff file
    my ($in) = Bio::SeqIO->new(
						-file => $input_file,
						-format => 'Fasta'
	);	
	my ($transcript_codition) = '';
    while ( my $seq = $in->next_seq() )
    {
        if($seq->id=~/^([^|]*)\|([^|]*)/)
        {
            my ($sequence_id) = $2;
            if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
			$transcript_codition .= ' $9 ~ /transcript_id \"'.$sequence_id.'\"/ ||';
        }
    }    
    if ($transcript_codition ne '') {
    	$transcript_codition =~ s/\|\|$//;
    	$transcript_codition = '('.$transcript_codition.')';
		@global_trans_cds_info=`awk -F "\t" '{if( \$3=="CDS" && $transcript_codition ){print \$0}}' $gff_file`;
		$logger->debug('awk -F "\t" \'{if( $3=="CDS" && '.$transcript_codition.' ){print $0}}\' '.$gff_file."\n");
		if ( scalar(@global_trans_cds_info) == 0 ) {
    		$logger->error("gff has not transcript information");			
		}
    }
    else {
    	$logger->error("gff has not transcript information");
    }

	# Save CDS coordinates within structure
	my (%total_trans_cds_coords);
	my ($trans_cds_coords);
    my ($in2) = Bio::SeqIO->new(
						-file => $input_file,
						-format => 'Fasta'
	);
    while ( my $seq = $in2->next_seq() )
    {

        if($seq->id=~/^([^|]*)\|([^|]*)/)
        {
            my ($sequence_id) = $2;
            if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }

			foreach my $trans_cds_info (@global_trans_cds_info)
			{
				if (defined $trans_cds_info and ($trans_cds_info =~ /$sequence_id/))
				{
					my (@cds_info) = split /\t/, $trans_cds_info;
					my ($strand);
					if(defined $cds_info[6])
					{
						$strand=$cds_info[6];
						$trans_cds_coords->{$sequence_id}->{'strand'} = $strand
					}
					if(defined $cds_info[3] and defined $cds_info[4])
					{					
						my ($start) = $cds_info[3];
						my ($end) = $cds_info[4];
						my ($edges);
						if(defined $strand) # Get the whole CDS coordinates by strand
						{
							$total_trans_cds_coords{$strand}{'starts'}{$start} = undef unless(exists $total_trans_cds_coords{$strand}{'starts'}{$start});
							$total_trans_cds_coords{$strand}{'ends'}{$end} = undef unless(exists $total_trans_cds_coords{$strand}{'ends'}{$end});
						}

						# Get the CDS coord for each transcript (with frame)
						if(defined $cds_info[7]) # save the frame
						{
							my ($frame) = $cds_info[7];
							$edges = join ":", $start, $end, $frame;
						}
						else {
							$edges = join ":", $start, $end;
						}
						push(@{$trans_cds_coords->{$sequence_id}->{'coord'}}, $edges);
					}					
				}
			}
		}
    }
    
	my $mini_cds_interval_pool;
	while(my ($strand, $trans_cds_coords) = each(%total_trans_cds_coords))
    {
		my @mini_cds_starts;
		my @mini_cds_ends;
    	if($strand eq '-')
    	{
			my @trans_cds_starts = sort { $b <=> $a } keys %{$trans_cds_coords->{'ends'}};  # GFF 'end' is reverse-strand CDS start
			my @trans_cds_ends = sort { $b <=> $a } keys %{$trans_cds_coords->{'starts'}};  # GFF 'start' is reverse-strand CDS end
			# TODO: verify assumption that max coord is a start, min coord is an end
			my %mini_cds_start_set = map { $_ => 1 } @trans_cds_starts;
			foreach my $end (@trans_cds_ends[ 0 .. ($#trans_cds_ends-1) ]) {
				$mini_cds_start_set{$end - 1} = 1;
			}
			@mini_cds_starts = sort { $b <=> $a } keys %mini_cds_start_set;
			my %mini_cds_end_set = map { $_ => 1 } @trans_cds_ends;
			foreach my $start (@trans_cds_starts[ 1 .. $#trans_cds_starts ]) {
				$mini_cds_end_set{$start + 1} = 1;
			}
			@mini_cds_ends = sort { $b <=> $a } keys %mini_cds_end_set;

    	} else {

			my @trans_cds_starts = sort { $a <=> $b } keys %{$trans_cds_coords->{'starts'}};
			my @trans_cds_ends = sort { $a <=> $b } keys %{$trans_cds_coords->{'ends'}};
			# TODO: verify assumption that min coord is a start, max coord is an end
			my %mini_cds_start_set = map { $_ => 1 } @trans_cds_starts;
			foreach my $end (@trans_cds_ends[ 0 .. ($#trans_cds_ends-1) ]) {
				$mini_cds_start_set{$end + 1} = 1;
			}
			@mini_cds_starts = sort { $a <=> $b } keys %mini_cds_start_set;
			my %mini_cds_end_set = map { $_ => 1 } @trans_cds_ends;
			foreach my $start (@trans_cds_starts[ 1 .. $#trans_cds_starts ]) {
				$mini_cds_end_set{$start - 1} = 1;
			}
			@mini_cds_ends = sort { $a <=> $b } keys %mini_cds_end_set;
    	}

		for my $i ( 0 .. $#mini_cds_starts ) {
			my @mini_cds_interval = ($mini_cds_starts[$i], $mini_cds_ends[$i]);
			push(@{$mini_cds_interval_pool->{$strand}}, \@mini_cds_interval);
		}
    }

    return ($trans_cds_coords, $mini_cds_interval_pool);
}

# Get relative coordinates of peptide coming from transcriptome coordinates
sub _get_peptide_coordinate_sans_frames($$$$$$)
{
	my ($trans_cds_start, $trans_cds_end, $pep_cds_end, $frame, $cds_order_id, $num_cds) = @_;

	my ($pep_cds_start) = $pep_cds_end+1;
	my ($pep_cds_end_div) = int(abs(($trans_cds_end + 1 + $frame) - $trans_cds_start) / 3);
	my ($pep_cds_end_mod) = abs(($trans_cds_end + 1 + $frame) - $trans_cds_start) % 3;
	my ($accumulate) = 0;
	if($pep_cds_end_mod == 1)
	{
		$frame=1;
	}
	elsif($pep_cds_end_mod == 2)
	{
		$frame=2;
		$accumulate=1;
	}
	else
	{
		$frame=0;
	}
	if ( $pep_cds_end_div == 0 ) { # cases when the CDS is smaller than 3 bp (1st cds of ENST00000372776 -rel7-)
		$pep_cds_end=$pep_cds_start;
		$frame=0;
	}
	if ( $num_cds == ($cds_order_id+1) ) { # we add the accumulate residues in the last CDS
		$pep_cds_end+=$pep_cds_end_div+$accumulate;
	}
	else {
		$pep_cds_end+=$pep_cds_end_div;
	}

	#return ($pep_cds_start,$pep_cds_end,$frame);
	return ($pep_cds_start,$pep_cds_end);
}

# Get relative coordinates of peptide coming from transcriptome coordinates
sub _get_peptide_coordinate_with_frames($$$$$$)
{
	my ($trans_cds_start, $trans_cds_end, $pep_cds_end, $frame, $cds_order_id, $num_cds) = @_;

	my ($pep_cds_start) = $pep_cds_end+1;
	my ($pep_cds_len) = abs($trans_cds_end - $trans_cds_start) + 1;
	my ($offset_pep_cds_len) = $pep_cds_len - $frame;
	my ($pep_cds_end_div) = ceil($offset_pep_cds_len / 3);  # include partial codon at end, if present
	my ($pep_cds_end_mod) = $offset_pep_cds_len % 3;

	my $next_frame;
	if($pep_cds_end_mod == 1)
	{
		$next_frame=2;
	}
	elsif($pep_cds_end_mod == 2)
	{
		$next_frame=1;
	}
	else
	{
		$next_frame=0;
	}

	if ( $cds_order_id == 0 && $frame != 0 ) {  # include partial codon at start, if present
		$pep_cds_end_div += 1;
	}

	if ( $pep_cds_end_div == 0 ) { # cases when the CDS is smaller than 3 bp (1st cds of ENST00000372776 -rel7-)
			$pep_cds_end=$pep_cds_start;
	} else {
		$pep_cds_end = $pep_cds_start + $pep_cds_end_div - 1;
	}
	if ( $cds_order_id == ($num_cds - 1) && $next_frame != 0 ) {
		warnings::warn("CDS has incomplete reading frame\n");
	}

	if ($calc_frames) {
		return ($pep_cds_start,$pep_cds_end,$next_frame);
	} else {
		return ($pep_cds_start,$pep_cds_end);
	}
}

# Get relative coordinates of peptide coming from transcriptome coordinates
sub _get_mini_peptide_coordinate_sans_frames($$$$$$$)
{
	my ($trans_cds_start, $trans_cds_end, $pep_cds_start, $pep_cds_end, $frame, $cds_order_id, $num_cds) = @_;

	#my ($pep_cds_start) = $pep_cds_end+1;
	my ($pep_cds_end_div) = int(abs(($trans_cds_end + 1 + $frame) - $trans_cds_start) / 3);
	my ($pep_cds_end_mod) = abs(($trans_cds_end + 1 + $frame) - $trans_cds_start) % 3;
	my ($accumulate) = 0;
	if($pep_cds_end_mod == 1)
	{
		$frame=1;
	}
	elsif($pep_cds_end_mod == 2)
	{
		$frame=2;
		$accumulate=1;
	}
	else
	{
		$frame=0;
	}
	if ( $pep_cds_end_div == 0 ) { # cases when the CDS is smaller than 3 bp (1st cds of ENST00000372776 -rel7-)
		$pep_cds_end=$pep_cds_start;
		$frame=0;
	}
	if ( $num_cds == ($cds_order_id+1) ) { # we add the accumulate residues in the last CDS
		$pep_cds_end+=$pep_cds_end_div+$accumulate;
	}
	else {
		$pep_cds_end+=$pep_cds_end_div;
	}

	#return ($pep_cds_start,$pep_cds_end,$frame);
	return ($pep_cds_start,$pep_cds_end);
}

# Get relative coordinates of peptide coming from transcriptome coordinates
sub _get_mini_peptide_coordinate_with_frames($$$$$$$)
{
	my ($trans_cds_start, $trans_cds_end, $pep_cds_start, $pep_cds_end, $frame, $cds_order_id, $num_cds) = @_;

	my ($pep_cds_len) = abs($trans_cds_end - $trans_cds_start) + 1;
	my ($offset_pep_cds_len) = $pep_cds_len - $frame;
	my ($pep_cds_end_div) = ceil($offset_pep_cds_len / 3);  # include partial codon at end, if present
	my ($pep_cds_end_mod) = $offset_pep_cds_len % 3;

	my $next_frame;
	if($pep_cds_end_mod == 1)
	{
		$next_frame=2;
	}
	elsif($pep_cds_end_mod == 2)
	{
		$next_frame=1;
	}
	else
	{
		$next_frame=0;
	}

	if ( $pep_cds_start == 1 && $frame != 0 ) {  # include partial codon at start, if present
		$pep_cds_end_div += 1;
	}

	if ( $pep_cds_end_div == 0 ) { # cases when the CDS is smaller than 3 bp (1st cds of ENST00000372776 -rel7-)
			$pep_cds_end=$pep_cds_start;
	} else {
			$pep_cds_end = $pep_cds_start + $pep_cds_end_div - 1;
	}

	if ($calc_frames) {
		return ($pep_cds_start,$pep_cds_end,$next_frame);
	} else {
		return ($pep_cds_start,$pep_cds_end);
	}
}

# Get the intersection of CDS coordinates creating "mini-cds"
sub _get_init_report($$$$)
{
	my ($sequence_id, $sequence, $trans_cds_coords, $mini_cds_interval_pool) = @_;

	# Get CDS coordinates relative to peptide
	my ($cds_list) = (); # CDS info
	my ($num_cds) = scalar(@{$trans_cds_coords->{$sequence_id}->{'coord'}});
	my ($pep_cds_start) = 0;
	my ($pep_cds_end) = 0;
	my $pep_cds_start_frame;
	my $pep_cds_end_frame;
	for ( my $cds_order_id = 0; $cds_order_id < $num_cds; $cds_order_id++ )
	{
		my ($trans_strand) = $trans_cds_coords->{$sequence_id}->{'strand'};
		my ($cds_coords) = $trans_cds_coords->{$sequence_id}->{'coord'}->[$cds_order_id];
		my (@boundaries) = split ":", $cds_coords;
		my ($trans_cds_start) = $boundaries[0];
		my ($trans_cds_end) = $boundaries[1];
		my ($trans_cds_frame) = $boundaries[2];

		if ($calc_frames) {

			if ( $cds_order_id == 0 ) {
				$pep_cds_end_frame = $trans_cds_frame;
			}
			$pep_cds_start_frame = $pep_cds_end_frame;

			($pep_cds_start,$pep_cds_end,$pep_cds_end_frame) = _get_peptide_coordinate_with_frames($trans_cds_start, $trans_cds_end, $pep_cds_end, $pep_cds_start_frame, $cds_order_id, $num_cds);
			if ( $cds_order_id + 1 < $num_cds ) {
				my ($pred_next_frame) = $pep_cds_end_frame;
				my ($next_coords) = $trans_cds_coords->{$sequence_id}->{'coord'}->[$cds_order_id + 1];
				my (@next_boundaries) = split(':', $next_coords);
				my ($actual_next_frame) = $next_boundaries[2];
				if ( $pred_next_frame ne $actual_next_frame ) {
					$logger->debug("sequence $sequence_id has incorrect frame prediction: $pred_next_frame vs $actual_next_frame\n");
				}
			}
		} else {
			($pep_cds_start,$pep_cds_end) = _get_peptide_coordinate_sans_frames($trans_cds_start, $trans_cds_end, $pep_cds_end, 0, $cds_order_id, $num_cds);
		}

		my ($cds_index) = "$pep_cds_start:$pep_cds_end";
		my ($cds_frame) = $trans_cds_frame;
		my ($exon_index) = "$trans_cds_start:$trans_cds_end";
		
		$logger->debug("#Trans_CDS:$trans_cds_start-$trans_cds_end:$trans_cds_frame\[$pep_cds_start-$pep_cds_end\]\n");
				
		# Get the list of mini-CDS intervals depending on strand
		my (@mini_cds_intervals);
		if ($trans_strand eq '-') {
			my @interval_pool = @{$mini_cds_interval_pool->{'-'}};
			my @match_start_idxs = grep { $interval_pool[$_][0] == $trans_cds_end } 0 .. $#interval_pool;  # GFF 'end' is reverse-strand CDS start
			my @match_end_idxs = grep { $interval_pool[$_][1] == $trans_cds_start } 0 .. $#interval_pool;  # GFF 'start' is reverse-strand CDS end
			if ( @match_start_idxs && @match_end_idxs ) {
				@mini_cds_intervals = @interval_pool[ $match_start_idxs[0] .. $match_end_idxs[-1] ]
			}
		} elsif ($trans_strand eq '+') {
			my @interval_pool = @{$mini_cds_interval_pool->{'+'}};
			my @match_start_idxs = grep { $interval_pool[$_][0] == $trans_cds_start } 0 .. $#interval_pool;
			my @match_end_idxs = grep { $interval_pool[$_][1] == $trans_cds_end } 0 .. $#interval_pool;
			if ( @match_start_idxs && @match_end_idxs ) {
				@mini_cds_intervals = @interval_pool[ $match_start_idxs[0] .. $match_end_idxs[-1] ]
			}
		}

		if(@mini_cds_intervals)
		{
			my ($mini_cds_list);
			if (scalar(@mini_cds_intervals) == 1) {
				# Theres is not mini cds => the same CDS coordinates
				my ($mini_cds_edges) = join ":", $pep_cds_start,$pep_cds_end;
				my ($mini_cds_seq) = substr($sequence, $pep_cds_start-1, ($pep_cds_end - $pep_cds_start) + 1 );
				my ($mini_exon_edges) = join ":", $trans_cds_start,$trans_cds_end;
				if ($calc_frames) {
					push(@{$mini_cds_list}, {
								'index'	=> $mini_cds_edges,
								'coord'	=> $mini_exon_edges,
								'frame'	=> $pep_cds_start_frame,
								'score'	=> '-',
								'seq'	=> $mini_cds_seq,
					});
				} else {
					push(@{$mini_cds_list}, {
								'index'	=> $mini_cds_edges,
								'coord'	=> $mini_exon_edges,
								'score'	=> '-',
								'seq'	=> $mini_cds_seq,
					});
				}
			}
			else
			{
				my ($mini_trans_cds_start) = 0;
				my ($mini_trans_cds_end) =  0;
				my ($mini_pep_cds_start) = $pep_cds_start;
				my ($mini_pep_cds_end) = $pep_cds_start;
				my $mini_pep_start_frame;
				my $mini_pep_end_frame;
				my ($num_mini_cds) = scalar(@mini_cds_intervals);
				for(my $k=0;$k<$num_mini_cds;$k++)
				{
					my ($mini_trans_cds_start, $mini_trans_cds_end) = @{$mini_cds_intervals[$k]};

					if ($calc_frames) {
						if ( $k == 0 ) {
							$mini_pep_end_frame = $pep_cds_start_frame;
						}
						$mini_pep_start_frame = $mini_pep_end_frame;

						($mini_pep_cds_start,$mini_pep_cds_end,$mini_pep_end_frame) = _get_mini_peptide_coordinate_with_frames($mini_trans_cds_start, $mini_trans_cds_end, $mini_pep_cds_start, $mini_pep_cds_end, $mini_pep_start_frame, $k, $num_mini_cds);
					} else {
						($mini_pep_cds_start,$mini_pep_cds_end) = _get_mini_peptide_coordinate_sans_frames($mini_trans_cds_start, $mini_trans_cds_end, $mini_pep_cds_start, $mini_pep_cds_end, 0, $k, $num_mini_cds);
					}

					my ($mini_cds_edges) = join ":", $mini_pep_cds_start,$mini_pep_cds_end;
					my ($mini_cds_seq) = substr($sequence, $mini_pep_cds_start-1, ($mini_pep_cds_end - $mini_pep_cds_start) + 1 );
					my ($mini_exon_edges);
					if ($trans_strand eq '-') {
						$mini_exon_edges = join ":", $mini_trans_cds_end,$mini_trans_cds_start;
					}
					elsif ($trans_strand eq '+') {
						$mini_exon_edges = join ":", $mini_trans_cds_start,$mini_trans_cds_end;
					}
					if (defined $mini_exon_edges) {
						if ($calc_frames) {
							push(@{$mini_cds_list}, {
										'index' => $mini_cds_edges,
										'coord' => $mini_exon_edges,
										'frame' => $mini_pep_start_frame,
										'score'	=> '-',
										'seq'	=> $mini_cds_seq,
							});
						} else {
							push(@{$mini_cds_list}, {
										'index' => $mini_cds_edges,
										'coord' => $mini_exon_edges,
										'score'	=> '-',
										'seq'	=> $mini_cds_seq,
							});
						}

						$logger->debug("\tMini_trans_CDS:$mini_exon_edges\[$mini_cds_edges\]\n");					
					}
					else {
						$logger->error("\tMini_trans_CDS_2:ERROR\n");
					}
					$mini_pep_cds_start=$mini_pep_cds_end+1;
				}

				# Until the end of the transcript (if not exceed the main coordinates)
				if(	($mini_pep_cds_start <= $pep_cds_end) and ($mini_pep_cds_end <= $pep_cds_end)) 
				{
					my ($mini_cds_edges) = join ":", $mini_pep_cds_start,$pep_cds_end;
					my ($mini_cds_seq) = substr($sequence, $mini_pep_cds_start-1, ($pep_cds_end - $mini_pep_cds_start) + 1 );
					my ($mini_exon_edges);
					if ($trans_strand eq '-') {
						$mini_exon_edges = join ":", $trans_cds_start,$mini_trans_cds_end;
					}
					elsif ($trans_strand eq '+') {
						$mini_exon_edges = join ":", $mini_trans_cds_end,$trans_cds_end;
					}
					if (defined $mini_exon_edges) {
						if ($calc_frames) {
							push(@{$mini_cds_list}, {
										'index' => $mini_cds_edges,
										'coord' => $mini_exon_edges,
										'frame' => $mini_pep_start_frame,
										'score'	=> '-',
										'seq'	=> $mini_cds_seq,
							});
						} else {
							push(@{$mini_cds_list}, {
										'index' => $mini_cds_edges,
										'coord' => $mini_exon_edges,
										'score'	=> '-',
										'seq'	=> $mini_cds_seq,
							});
						}
						$logger->debug("\tMini_trans_CDS_2:$mini_exon_edges\[$mini_cds_edges\]\n");
					}
					else {
						$logger->error("\tMini_trans_CDS_2:ERROR\n");
					}
				}
			}
			push(@{$cds_list},{
						'index'		=> $cds_index,
						'frame'		=> $cds_frame,
						'coord'		=> $exon_index,						
						'score'		=> '-',
						'mini_cds'	=> $mini_cds_list
			});					
		}
	}
	return $cds_list;	
}




main();



__END__

=head1 NAME

matador3d

=head1 DESCRIPTION

Run Matador3D program

=head1 SYNOPSIS

matador3d

=head2 Required arguments:

    --conf <Config file>

    --gff <GFF file that contains CDS information>

    --input <Fasta sequence file>

    --output <Annotation output file>

=head2 Optional arguments:

	--appris <Flag that enables the output for APPRIS (default: NONE)>

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
    

=head1 EXAMPLE

perl matador3d.pl

	--conf=../conf/pipeline.ini

	--gff=data/gencode.v3c.annotation.GRCh37.gtf
	
	--input=examples/ENSG00000016864.faa
	
	--output=examples/ENSG00000016864.output


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
