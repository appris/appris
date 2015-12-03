#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$DEFAULT_CONFIG_FILE
	
	$PROG_IN_SUFFIX
	$PROG_OUT_SUFFIX
	$WSPACE_BASE
	$WSPACE_CACHE
	$RUN_PROGRAM
	$PROG_DB
	$PROG_EVALUE
	$APPRIS_CUTOFF
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
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_BASE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('MATADOR3D_VARS', 'name').'/';
$WSPACE_CACHE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('CACHE_VARS', 'name').'/';
$RUN_PROGRAM		= $cfg->val( 'MATADOR3D_VARS', 'program');
$PROG_DB			= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$cfg->val('MATADOR3D_VARS', 'db');
$PROG_EVALUE		= $cfg->val('MATADOR3D_VARS', 'evalue');
$PROG_IN_SUFFIX		= 'faa';
$PROG_OUT_SUFFIX	= 'pdb';
$APPRIS_CUTOFF		= $cfg->val( 'MATADOR3D_VARS', 'cutoff');

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
sub _parse_blast($$);
sub _check_alignment($$$);
sub _run_blastpgp($$);
sub _get_best_cds_score($);
sub _get_biggest_cds($$$);
sub _get_biggest_mini_cds($$$);
sub _get_record_annotations($);
sub _get_appris_annotations($);
sub _get_cds_coordinates_from_gff($$);
sub _get_peptide_coordinate($$$$$$);
sub _get_mini_peptide_coordinate($$$$$$$);
sub _get_init_report($$$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Get the CDS coordinates info from gff file ---------------
	$logger->info("##Get cds coordinates from gff ---------------\n");
	my ($trans_cds_coords, $sort_total_cds_coords) = _get_cds_coordinates_from_gff($input_file, $gff_file);
    $logger->debug("##CDS coordenates ---------------\n".Dumper($trans_cds_coords));
    $logger->debug("##Sorted CDS coordenates ---------------\n".Dumper($sort_total_cds_coords));


	# For every sequence run the method ---------------
    my ($transcript_report);
    my $fasta_object = Bio::SeqIO->new(
                        -file => $input_file,
                        -format => 'Fasta'
    );
	while ( my $seq = $fasta_object->next_seq() )
	{
		if ( $seq->id=~/([^|]*)/ )
		{			
			my ($sequence_id) = $1;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
            my ($sequence) = $seq->seq;
            
            $logger->info("\n##$sequence_id ###############################\n");
            

		    # get the intersection of CDS coordinates creating "mini-cds"
    		my ($cds_reports) = _get_init_report($sequence_id, $trans_cds_coords, $sort_total_cds_coords);
			$logger->debug("##Init CDS report ---------------\n".Dumper($cds_reports));

			if(defined $cds_reports)
			{
                # Init transcript report from cds list
				$transcript_report->{$sequence_id}->{'score'} = 0;
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
	
	# Create temporal file for blast
	my ($fasta_sequence_file) = $WSPACE_BASE.'/'.$sequence_id.'.'.$PROG_IN_SUFFIX;
	unless(-e $fasta_sequence_file and (-s $fasta_sequence_file > 0) ) # Cached fasta
	{	
		my ($fasta_sequence_content_file) = ">$sequence_id\n$sequence";
		my ($print_fasta) = printStringIntoFile($fasta_sequence_content_file, $fasta_sequence_file);
		unless( defined $print_fasta ) {
			$logger->error("Can not create temporal file:$fasta_sequence_file: $!\n");
		}
	}
	
	# Run blast
	my ($blast_sequence_file) = $WSPACE_CACHE.'/'.$sequence_id.'.'.$PROG_OUT_SUFFIX;               
	unless(-e $blast_sequence_file and (-s $blast_sequence_file >0)) # Cached Blast
	{
		eval
		{
			# parameters to run blaspgp
			$logger->debug("$RUN_PROGRAM -d $PROG_DB -i $fasta_sequence_file -e$PROG_EVALUE -o $blast_sequence_file\n");                   
			system("$RUN_PROGRAM -d $PROG_DB -i $fasta_sequence_file -e$PROG_EVALUE -o $blast_sequence_file");			
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

    my ($loopstart) = 0;
	my ($loopstart2) = 0;
	my ($align_start) = 0;
	my ($align_end) = 0;

	foreach my $cds (@{$cds_list})
	{
    	my ($cds_index) = $cds->{'index'};
    	my ($cds_coord) = $cds->{'coord'};
    	my ($mini_cds_list) = $cds->{'mini_cds'};
    	my ($mini_cds_score_list);
		my ($mini_score) = 0;    	
		$logger->debug("#cds_trans:$cds_coord\[$cds_index\]\n");

		foreach my $mini_cds (@{$mini_cds_list})
		{
	    	my ($mini_cds_index) = $mini_cds->{'index'};
	    	my ($mini_cds_coord) = $mini_cds->{'coord'};
			my @boundaries = split ":", $mini_cds_index;
			my $mini_cds_start = $boundaries[0];
			my $mini_cds_end = $boundaries[1];
			$logger->debug("#mini_cds_trans:$mini_cds_coord\[$mini_cds_index\]\n");           

			# Finish when target-candidate don't cover exons or
			# Next until exon is within target-candidate or
			# Next if the mini range is less than 6 residues
			if (	($mini_cds_start > $targend) or
					($targstart > $mini_cds_end) or 
					($mini_cds_end - $mini_cds_start < 6)
			) {
				$mini_cds_score_list->{$mini_cds_index} = {
							'index'			=> $mini_cds_index,
							'coord'			=> $mini_cds_coord,
							'score'			=> '-',
				};
				next;
			}			
		        
			my $gapres = 0;
			my $identities = 0;
			my $res = 0;
			my $j = 0;
			my $n = 0;
	        
			# Get the index start depending the minus coordinate 
			if ($mini_cds_start < $targstart)
				{$loopstart2 = $targstart}
			else
				{$loopstart2 = $mini_cds_start}
	
			# Count and classify the residues      		
			for ($res=$loopstart,$j=$loopstart2;$j<=$mini_cds_end;$j++,$res++)
			{
				if(defined $target[$res] and defined $candidate[$res])
				{
					if ($target[$res] eq $candidate[$res])
						{$identities++}
					if ($target[$res] eq "-")
						{$gapres++;$j--}
					if ($candidate[$res]eq "-")
						{$gapres++;}
					$n++;           	
				}
			}
			
			$loopstart = $res;
			$align_start = $loopstart2;
			$align_end = $loopstart2 + $n - 1; 

			$logger->debug("#loopstart: $loopstart loopstart: $loopstart2 mini_cds_end: $mini_cds_end res: $res j: $j n: $n\n");

			$logger->debug("#mini_aligns_coord\[$align_start:$align_end\]\n");
	
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
			$logger->debug("\tIdentity: $identity = $identities/$mini_cds_residues * 100\n");
			$logger->debug("\tAlign_identity: $align_identity = $identities/$align_residues * 100\n");           
			$logger->debug("\tTotal identity: $totalidentity\n");
	
			# Value of gap score: the sort of conditions is important!!!
			my $totalgaps = 0;
			my $gaps = 0;
			$gaps = $gapres/$mini_cds_residues*100;
			if ($gaps <= 3)
				{$totalgaps = 1}
			elsif ($gaps <= 6)
				{$totalgaps = 0.80}
			elsif ($gaps <= 10)
				{$totalgaps = 0.50}
			elsif ($gaps <= 14)
				{$totalgaps = 0.33}
			elsif ($gaps <= 18)
				{$totalgaps = 0.20}
			elsif ($gaps > 40)
				{$totalgaps = -1}
			elsif ($gaps > 33)
				{$totalgaps = -0.5}
			elsif ($gaps > 25)
				{$totalgaps = 0}
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
			if ( ($totalgaps == -1 or $totalgaps == -0.5) and ($mini_cds_residues > 5 ) ) { $escore = $totalgaps }
			$logger->debug("\tScore: $escore\n");
									
			my ($mini_cds_report) = {
						'index'			=> $mini_cds_index,
						'coord'			=> $mini_cds_coord,
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
				if ( ($cds_pep_coord[1] - $cds_pep_coord[0]) >= 6 )
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
			if ( ($mini_cds_pep_coord[1] - $mini_cds_pep_coord[0]) >= 6 )
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
	my ($biggest_mini_pdb_cds_pdb) = $main_mini_pdb_cds_report->{'pdb'};
	$biggest_mini_pdb_cds_score = 0 if ( !defined $biggest_mini_pdb_cds_score or ($biggest_mini_pdb_cds_score eq '-') );
	$biggest_mini_pdb_cds_pdb = '' unless (defined $biggest_mini_pdb_cds_pdb );
	
	# If we found one CDS coming from the same exon and it had bigger score than 
	# the given cds, we get it
	while(my ($sequence_id,$trans_report) = each(%{$transcript_report}))
	{
		next if ($main_seq_id eq $sequence_id); # Jump for the same sequence
		
		if(exists $trans_report->{'cds'} and defined $trans_report->{'cds'})
		{
			foreach my $trans_pdb_cds_report (@{$trans_report->{'cds'}})
			{
				my ($trans_mini_cds_list) = $trans_pdb_cds_report->{'mini_cds'};				
				foreach my $trans_mini_cds (@{$trans_mini_cds_list})			
				{
					my (@trans_mini_cds) = split(':',$trans_mini_cds->{'coord'}); # Sticking paster for reverse strand
					my ($trans_mini_cds_coord_aux) = $trans_mini_cds[1].':'.$trans_mini_cds[0];
					#$logger->debug("(".$trans_mini_cds->{'coord'}." eq ".$biggest_mini_pdb_cds_coord.") or (".$trans_mini_cds_coord_aux." eq ".$biggest_mini_pdb_cds_coord.")) and (".$trans_mini_cds->{'score'}." > ".$biggest_mini_pdb_cds_score.") and (".$trans_mini_cds->{'pdb'}." ne ".$biggest_mini_pdb_cds_pdb.")\n");						
					if(
						(($trans_mini_cds->{'coord'} eq $biggest_mini_pdb_cds_coord) or
						($trans_mini_cds_coord_aux eq $biggest_mini_pdb_cds_coord))
						and 
						($trans_mini_cds->{'score'} ne '-') 
						and
						($trans_mini_cds->{'score'} > $biggest_mini_pdb_cds_score)
						and
						($trans_mini_cds->{'pdb'} ne $biggest_mini_pdb_cds_pdb)
					){
						$logger->debug("Into:_get_biggest_mini_cds\n");
						$external_id = $sequence_id;
						$biggest_mini_pdb_cds_score = $trans_mini_cds->{'score'};
						$mini_pdb_cds_report = $trans_mini_cds;
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
        if($seq->id=~/([^|]*)/)
        {
            my ($sequence_id) = $1;
            if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
			$transcript_codition .= ' $9 ~ /'.$sequence_id.'/ ||';
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
	my ($total_trans_cds_coords);
	my ($trans_cds_coords);
    my ($in2) = Bio::SeqIO->new(
						-file => $input_file,
						-format => 'Fasta'
	);
    while ( my $seq = $in2->next_seq() )
    {

        if($seq->id=~/([^|]*)/)
        {
            my ($sequence_id) = $1;
            if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }

			foreach my $trans_cds_info (@global_trans_cds_info)
			{
				#if (defined $trans_cds_info and ($trans_cds_info =~ /"$sequence_id"/))
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
	
						if(defined $strand) # Get the whole CDS coordinates by strand
						{
							$total_trans_cds_coords->{$strand}->{$start} = undef unless(exists $total_trans_cds_coords->{$strand}->{$start});
							$total_trans_cds_coords->{$strand}->{$end} = undef unless(exists $total_trans_cds_coords->{$strand}->{$end});						
						}
						
						my ($edges) = join ":", $start, $end; # Get the CDS coord for each transcript 
						push(@{$trans_cds_coords->{$sequence_id}->{'coord'}}, $edges);
					}					
				}
			}
		}
    }
    
    my ($sort_total_cds_coords);
    while(my ($strand,$total_cds_coords) = each(%{$total_trans_cds_coords}))
    {
    	if($strand eq '-')
    	{
    		$sort_total_cds_coords->{$strand}=join ":", sort {$b <=> $a} (keys %{$total_cds_coords});    		
    	} else {
    		$sort_total_cds_coords->{$strand}=join ":", sort {$a <=> $b} (keys %{$total_cds_coords});
    	}
    }

    return ($trans_cds_coords, $sort_total_cds_coords);	
}

# Get relative coordinates of peptide coming from transcriptome coordinates
sub _get_peptide_coordinate($$$$$$)
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
		
	return ($pep_cds_start,$pep_cds_end,$frame);
}

# Get relative coordinates of peptide coming from transcriptome coordinates
sub _get_mini_peptide_coordinate($$$$$$$)
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
		
	return ($pep_cds_start,$pep_cds_end,$frame);
}

# Get the intersection of CDS coordinates creating "mini-cds"
sub _get_init_report($$$)
{
	my ($sequence_id, $trans_cds_coords, $sort_total_cds_coords) = @_;

	# Get CDS coordinates relative to peptide
	my ($cds_list) = (); # CDS info
	my ($num_cds) = scalar(@{$trans_cds_coords->{$sequence_id}->{'coord'}});
	my ($pep_cds_start) = 0;
	my ($pep_cds_end) = 0;
	my ($frame) = 0;
	for ( my $cds_order_id = 0; $cds_order_id < $num_cds; $cds_order_id++ )
	{
		my ($trans_strand) = $trans_cds_coords->{$sequence_id}->{'strand'};
		my ($cds_coords) = $trans_cds_coords->{$sequence_id}->{'coord'}->[$cds_order_id];
		my (@boundaries) = split ":", $cds_coords;
		my ($trans_cds_start) = $boundaries[0];
		my ($trans_cds_end) = $boundaries[1];
		($pep_cds_start,$pep_cds_end,$frame) = _get_peptide_coordinate($trans_cds_start, $trans_cds_end, $pep_cds_end, $frame, $cds_order_id, $num_cds);
		my ($cds_index) = "$pep_cds_start:$pep_cds_end";
		my ($exon_index) = "$trans_cds_start:$trans_cds_end";
		
		$logger->debug("#Trans_CDS:$trans_cds_start-$trans_cds_end\[$pep_cds_start-$pep_cds_end:$frame\]\n");
				
		# Sort the list of CDS depending on strand
		my ($mini_cds_coord_range);
		if	(($trans_strand eq '-') and 
			($sort_total_cds_coords->{$trans_strand}=~/($trans_cds_end.*:$trans_cds_start)/)){
			$mini_cds_coord_range = $1;				
		} elsif (($trans_strand eq '+') and 
			($sort_total_cds_coords->{$trans_strand}=~/($trans_cds_start.*:$trans_cds_end)/)){
			$mini_cds_coord_range = $1;
		}
		if(defined $mini_cds_coord_range)
		{
			my ($mini_cds_list);
			my (@mini_cds_coord_list) = split ":", $mini_cds_coord_range;				
			if (scalar(@mini_cds_coord_list) == 2) {
				# Theres is not mini cds => the same CDS coordinates
				my ($mini_cds_edges) = join ":", $pep_cds_start,$pep_cds_end;
				my ($mini_exon_edges) = join ":", $trans_cds_start,$trans_cds_end;
				push(@{$mini_cds_list}, {
							'index'		=> $mini_cds_edges,
							'coord'		=> $mini_exon_edges,
							'score'		=> '-',
				});						
			}
			else
			{
				my ($mini_trans_cds_start) = 0;
				my ($mini_trans_cds_end) =  0;
				my ($mini_pep_cds_start) = $pep_cds_start;
				my ($mini_pep_cds_end) = $pep_cds_start;
				my ($mini_pep_frame) = 0;
				my ($mini_num_cds) = scalar(@mini_cds_coord_list)-1;
				for(my $k=1;$k<$mini_num_cds;$k++)
				{
					$mini_trans_cds_start = $mini_cds_coord_list[$k-1];
					$mini_trans_cds_end = $mini_cds_coord_list[$k];
					
					($mini_pep_cds_start,$mini_pep_cds_end,$mini_pep_frame) = _get_mini_peptide_coordinate($mini_trans_cds_start, $mini_trans_cds_end, $mini_pep_cds_start, $mini_pep_cds_end, $mini_pep_frame, $k, $mini_num_cds);
					
					my ($mini_cds_edges) = join ":", $mini_pep_cds_start,$mini_pep_cds_end;
					my ($mini_exon_edges);
					if ($trans_strand eq '-') {
						$mini_exon_edges = join ":", $mini_trans_cds_end,$mini_trans_cds_start;
					}
					elsif ($trans_strand eq '+') {
						$mini_exon_edges = join ":", $mini_trans_cds_start,$mini_trans_cds_end;
					}
					if (defined $mini_exon_edges) {
						push(@{$mini_cds_list}, {
									'index' => $mini_cds_edges,
									'coord' => $mini_exon_edges,
									'score'		=> '-',
						});
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
					my ($mini_exon_edges);
					if ($trans_strand eq '-') {
						$mini_exon_edges = join ":", $trans_cds_start,$mini_trans_cds_end;
					}
					elsif ($trans_strand eq '+') {
						$mini_exon_edges = join ":", $mini_trans_cds_end,$trans_cds_end;
					}
					if (defined $mini_exon_edges) {
						push(@{$mini_cds_list}, {
									'index' => $mini_cds_edges,
									'coord' => $mini_exon_edges,
									'score'		=> '-',
						});
						$logger->debug("\tMini_trans_CDS_2:$mini_exon_edges\[$mini_cds_edges\]\n");
					}
					else {
						$logger->error("\tMini_trans_CDS_2:ERROR\n");
					}
				}
			}
			push(@{$cds_list},{
						'index'		=> $cds_index,
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
