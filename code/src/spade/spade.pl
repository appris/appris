#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;
use List::Util qw(max minstr sum);

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile prepare_workspace );

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
my ($cfg) 			= new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_TMP			= $ENV{APPRIS_TMP_DIR};
$WSPACE_CACHE		= $ENV{APPRIS_PROGRAMS_CACHE_DIR};
$CACHE_FLAG			= $cfg->val('SPADE_VARS', 'cache');
$RUN_PROGRAM		= $cfg->val('SPADE_VARS', 'program');
$PROG_DB			= $cfg->val('SPADE_VARS', 'db');
$PROG_DB_DIR		= $ENV{APPRIS_PROGRAMS_DB_DIR}.'/'.$PROG_DB;
$PROG_EVALUE		= $cfg->val('SPADE_VARS', 'evalue');
$APPRIS_CUTOFF		= $cfg->val('SPADE_VARS', 'cutoff');


# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

my $EXP_CFG = new Config::IniFiles( -file => $ENV{APPRIS_EXP_CONF_FILE} );
my $recal_integrity_score = $EXP_CFG->val( 'spade_integrity', 'recal_integrity_score', 1 );

#####################
# Method prototypes #
#####################
sub _run_pfamscan($$);
sub _parse_pfamscan($$);
sub _get_best_score($$);
sub _get_best_domain($$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# For every sequence run the method ---------------
    my ($transcript_report);
    my ($sequence_report);
    my ($fasta_object) = Bio::SeqIO->new(
                        -file => $input_file,
                        -format => 'Fasta'
    );
	SEQ: while ( my $seq = $fasta_object->next_seq() )
	{
		if ( $seq->id=~/^([^|]*)\|([^|]*)/ )
		{			
			my ($sequence_id) = $2;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
            my ($sequence) = $seq->seq;
            $sequence_report->{$sequence_id} = $sequence;
            
            $logger->info("\n##$sequence_id ###############################\n");

			my ($stop_index) = index($sequence, '*');
			if ( $stop_index >= 0 && $stop_index < length($sequence) - 1 )
			{
				$logger->warning("skipping sequence $sequence_id due to internal translation stop site\n");
				next SEQ;
			}

			# Run blast
			$logger->info("\n##Running Pfamscan ---------------\n");
			my ($pfamscan_sequence_file) = _run_pfamscan($sequence_id, $sequence);
								  
			# Parse blast
			$logger->info("\n##Parsing Pfamscan ---------------\n");                
			my ($pfam_report) = _parse_pfamscan($pfamscan_sequence_file, $sequence_id);
			$transcript_report->{$sequence_id} = $pfam_report if (defined $pfam_report);
		}
    }

	if ( defined $transcript_report ) {
		$logger->debug("##CDS coordenates ---------------\n".Dumper($transcript_report));

		# If a transcript has not domains, we check if other transcript has the same sequence with domain (external domains)
		# But rememeber, we only accept the external domains that do not align with any region domain from the current transcript.
		_get_best_score($transcript_report, $sequence_report);
		$logger->debug("##The Best CDS scores ---------------\n".Dumper($transcript_report)."\n");

		# Get records by transcript ---------------
		my ($output_content) = _get_record_annotations($transcript_report);

		# Get the annotations for the main isoform /* APPRIS */ ----------------
		if ( defined $appris )
		{
			$output_content .= _get_appris_annotations_bitscore($transcript_report);
		}

		# Print records by transcript ---------------
		my ($print_out) = printStringIntoFile($output_content, $output_file);
		unless( defined $print_out ) {
			$logger->error("Can not create output file: $!\n");
		}
	}

	$logger->finish_log();
	
	exit 0;
}

# Get the matched domains from transcripts that has not domains
sub _get_best_score($$)
{
	my ($transcript_report, $sequence_report) = @_;
	
	my (%aux_transcript_report) = %{$transcript_report};
	while(my ($sequence_id, $trans_report) = each(%{$transcript_report}))
	{
		my ($sequence) = $sequence_report->{$sequence_id};
		$sequence = uc($sequence);
		my ($pfam_report) = _get_best_domain($sequence_id, $sequence, \%aux_transcript_report);
		if ( defined $pfam_report and defined $pfam_report->{'domains'} and (scalar(@{$pfam_report->{'domains'}}) > 0) ) {
			$transcript_report->{$sequence_id} = $pfam_report;
			$transcript_report->{$sequence_id}->{'result'} = $trans_report->{'result'};
		}
	}		
}

# calculate the domain integrity score
sub _calc_domain_integrity_score($$$$) {

	my ($num_domains, $num_possibly_damaged_domains, $num_damaged_domains, $num_wrong_domains) = @_;

	my $domain_integrity_score;
	if ($recal_integrity_score) {
		$domain_integrity_score = ($num_domains*0.75)
			+ ($num_possibly_damaged_domains*0.5)
			+ ($num_damaged_domains*0.25);
	} else {
		$domain_integrity_score = ($num_domains*1)
			+ ($num_possibly_damaged_domains*0.75)
			+ ($num_damaged_domains*0.5)
			+ ($num_wrong_domains*0.25);
	}

	return $domain_integrity_score
}

# Get the first domain that we have found
sub _get_best_domain($$$)
{
	my ($sequence_id, $sequence, $transcript_report) = @_;
	
	my ($pfam_report);
	my ($exist_domains);
	my ($best_domains);
	
	my ($alignment_list);
	my ($num_domains) = 0;
	my ($num_possibly_damaged_domains) = 0;
	my ($num_damaged_domains) = 0;
	my ($num_wrong_domains) = 0;
	my ($bitscore) = 0;
		
	# init the domains with the current domains
	if ( exists $transcript_report->{$sequence_id} ) {
		my ($pfam_report) = $transcript_report->{$sequence_id};
		if ( exists $pfam_report->{'domains'} and $pfam_report->{'domains'} )
		{
			foreach my $alignment_report (@{$pfam_report->{'domains'}})
			{
				my ($seq_domain) = $alignment_report->{'sequence'};
				my ($seq_domain_nogaps) = $seq_domain;
				$seq_domain_nogaps =~ s/\-*//g;
				my ($hmm_name) = $alignment_report->{'hmm_name'};
				my ($hmm_type) = $alignment_report->{'hmm_type'};
				my ($aln_start) = $alignment_report->{'alignment_start'};
				my ($aln_end) = $alignment_report->{'alignment_end'};
				my ($aln_length) = $alignment_report->{'alignment_length'};
				my ($bit_score) = $alignment_report->{'bit_score'};
				my ($aln_index) = $aln_start.":".$aln_end;
				$exist_domains->{$aln_index} = {					
					'hmm_name'		=> $hmm_name,
					'hmm_type'		=> $hmm_type,
					'seq_domain'	=> $seq_domain_nogaps,
					'aln_start'		=> $aln_start,
					'aln_end'		=> $aln_end,
					'aln_length'	=> $aln_length,
					'bit_score'		=> $bit_score,
					'report'		=> $alignment_report,
				};
			}
		}
	}
	
	# Discard original overlapping
	$best_domains = $exist_domains;
	if ( defined $exist_domains ) {
		my (@sorted_hmm_index) = sort { my ($a1,$b1); if($a=~/^([0-9]*)/){$a1=$1}; if($b=~/^([0-9]*)/){$b1=$1}; $a1<=>$b1} keys(%{$exist_domains});
		for (my $i=0; $i< scalar(@sorted_hmm_index); $i++) {
			my ($hmm_index_i) = $sorted_hmm_index[$i];
			if ( defined $exist_domains->{$hmm_index_i} ) {
				my ($domain_report_i) = $exist_domains->{$hmm_index_i};
				my ($hmm_name_i) = $domain_report_i->{'hmm_name'};
				my ($hmm_bitscore_i) = $domain_report_i->{'bit_score'};
				my ($aln_start_i) = $domain_report_i->{'aln_start'};
				my ($aln_end_i) = $domain_report_i->{'aln_end'};
				
				for (my $j=$i+1; $j < scalar(@sorted_hmm_index); $j++) {
					my ($hmm_index_j) = $sorted_hmm_index[$j];
					if ( defined $exist_domains->{$hmm_index_j} ) {
						my ($domain_report_j) = $exist_domains->{$hmm_index_j};
						my ($hmm_name_j) = $domain_report_j->{'hmm_name'};
						my ($hmm_bitscore_j) = $domain_report_j->{'bit_score'};
						my ($aln_start_j) = $domain_report_j->{'aln_start'};
						my ($aln_end_j) = $domain_report_j->{'aln_end'};

						# overlapping
						if ( ($aln_start_i <= $aln_end_j) and ($aln_start_j <= $aln_end_i) and ($hmm_name_i eq $hmm_name_j) ) {
							if ( $hmm_bitscore_j < $hmm_bitscore_i ) {
								$best_domains->{$hmm_index_j}->{'discarded'} = 1;
							}
							else {
								$best_domains->{$hmm_index_i}->{'discarded'} = 1;
							}
						}					
					}
				}	
				
			}
		}		
	}

	my (%candidates_left) = map { $_ => 1 } grep { $_ ne $sequence_id } keys(%{$transcript_report});

	my ($continuing_search) = 1;
	while ( $continuing_search && keys(%candidates_left) > 0 )
	{
		my (%cand_domain_sets);
		foreach my $sequence_id2 (keys %candidates_left)
		{
			my ($pfam_report) = $transcript_report->{$sequence_id2};
			if ( exists $pfam_report->{'domains'} and $pfam_report->{'domains'} )
			{
				my (@domains_by_dec_score) = sort {
					$b->{'bit_score'} <=> $a->{'bit_score'}
				} @{$pfam_report->{'domains'}};
				
				my ($cand_domains);
				foreach my $alignment_report (@domains_by_dec_score)
				{
					# check if domain exists in other sequence (without gaps)
					my ($seq_domain) = $alignment_report->{'sequence'};
					my ($seq_domain_nogaps) = $seq_domain;
					$seq_domain_nogaps =~ s/\-*//g;
					my ($seq_domain_num_gaps) = $seq_domain =~ tr/\-//;
					my ($index_domain) = rindex($sequence, $seq_domain_nogaps);

					# get the domain from the other sequence, saying from who
					# get the new range of alignment
					my ($hmm_name) = $alignment_report->{'hmm_name'};
					my ($hmm_type) = $alignment_report->{'hmm_type'};
					my ($aln_start) = $index_domain + 1;
					my ($aln_end) = $index_domain + length($seq_domain) - $seq_domain_num_gaps;
					my ($aln_length) = $alignment_report->{'alignment_length'};
					my ($hmm_bit_score) = $alignment_report->{'bit_score'};
					my ($hmm_index) = $aln_start.":".$aln_end;

					my (%alignment_report2) = %{$alignment_report};
					$alignment_report2{'external_id'} = $sequence_id2;
					$alignment_report2{'alignment_start'} = $aln_start;
					$alignment_report2{'alignment_end'} = $aln_end;

					if ( $index_domain != -1 )
					{
						my ($adding_domain) = 0;

						# the same alignment region
						if ( exists $best_domains->{$hmm_index} and !(exists $best_domains->{$hmm_index}->{'discarded'}) ) {
							# get the best bitscore for the same alignment
							my ($exists_hmm_name) = $best_domains->{$hmm_index}->{'hmm_name'};
							my ($exists_hmm_bitscore) = $best_domains->{$hmm_index}->{'bit_score'};
							if ( ($hmm_name eq $exists_hmm_name) and ($hmm_bit_score > $exists_hmm_bitscore) ) {
								# If there is no candidate domain in this region, take this one.
								# As we've sorted the domains by bitscore, we can be sure that
								# any existing candidate has a better bitscore than this one.
								if ( ! exists($cand_domains->{$hmm_index}) ) {
									$adding_domain = 1;
								}
							}
						}
						else {
							# different alignment region
							if ( ! _domain_overlaps_existing($hmm_index, $best_domains) &&
									! _domain_overlaps_existing($hmm_index, $cand_domains) ) {
								$adding_domain = 1;
							}
						}

						if ($adding_domain) {
							$cand_domains->{$hmm_index} = {
								'hmm_name'		=> $hmm_name,
								'hmm_type'		=> $hmm_type,
								'seq_domain'	=> $seq_domain_nogaps,
								'aln_start'		=> $aln_start,
								'aln_end'		=> $aln_end,
								'aln_length'	=> $aln_length,
								'bit_score'		=> $hmm_bit_score,
								'report'		=> \%alignment_report2,
							};
						}
					}
				}

				if ( keys(%{$cand_domains}) > 0 ) {
					my ($cand_set_bitscore) = sum(map { $_->{'bit_score'} } values(%{$cand_domains}) );
					$cand_domain_sets{$cand_set_bitscore}{$sequence_id2} = $cand_domains;
				}
			}
		}

		if ( keys(%cand_domain_sets) > 0 ) {
			my ($max_bitscore) = max(keys(%cand_domain_sets));

			# TODO: better tiebreaking between candidate domain sets
			my ($chosen_seq_id) = minstr(keys(%{$cand_domain_sets{$max_bitscore}}));
			delete $candidates_left{$chosen_seq_id};

			my $chosen_domains = $cand_domain_sets{$max_bitscore}{$chosen_seq_id};
			while (my ($hmm_index, $domain_report) = each(%{$chosen_domains}) ) {
				$best_domains->{$hmm_index} = $domain_report;
			}
		} else {
			$continuing_search = 0;
		}
	}

	# report the alignment list
	while (my ($hmm_index, $domain_report) = each(%{$best_domains}) ) {
		my ($alignment_report) = $domain_report->{'report'};
		if ( !(exists $best_domains->{$hmm_index}->{'discarded'}) ) {
			if( $alignment_report->{'type_domain'} eq 'domain' ) {
				$num_domains++;
			}
			elsif( $alignment_report->{'type_domain'} eq 'domain_possibly_damaged' ) {
				$num_possibly_damaged_domains++;
			}
			elsif ( $alignment_report->{'type_domain'} eq 'domain_damaged' ) {
				$num_damaged_domains++;
			}
			elsif ( $alignment_report->{'type_domain'} eq 'domain_wrong' ) {
				$num_wrong_domains++;
			}
			
			if ( exists $alignment_report->{'bit_score'} ) { $bitscore += $alignment_report->{'bit_score'} }				
		}
		else {
			$alignment_report->{'discarded'} = 1;
		}
		push(@{$alignment_list}, $alignment_report);
	}
	$pfam_report->{'bitscore'} = $bitscore;
	$pfam_report->{'domains'} = $alignment_list;
	$pfam_report->{'num_domains'} = $num_domains;
	$pfam_report->{'num_possibly_damaged_domains'} = $num_possibly_damaged_domains;
	$pfam_report->{'num_damaged_domains'} = $num_damaged_domains;
	$pfam_report->{'num_wrong_domains'} = $num_wrong_domains;
	$pfam_report->{'domain_integrity'} = _calc_domain_integrity_score($num_domains,
																	  $num_possibly_damaged_domains,
																	  $num_damaged_domains,
																	  $num_wrong_domains);
	
	return $pfam_report;	
}

# Check if domain overlaps any existing domain.
sub _domain_overlaps_existing($$)
{
	my ($hmm_index, $existing_domains) = @_;
	my ($overlap_found) = 0;

	my ($aln_start, $aln_end) = split(/:/, $hmm_index);

	foreach my $exists_hmm_rep (values %{$existing_domains}) {

		if ( !(exists $exists_hmm_rep->{'discarded'}) &&
			  ($exists_hmm_rep->{'aln_start'} <= $aln_end) &&
			  ($aln_start <= $exists_hmm_rep->{'aln_end'}) ) {
			$overlap_found = 1;
			last;
		}

	}

	return $overlap_found;
}

# Run pfamscan
sub _run_pfamscan($$)
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
	
	# Create temporal file
	my ($fasta_sequence_file) = $ws_cache.'/seq.faa';
	unless(-e $fasta_sequence_file and (-s $fasta_sequence_file > 0) ) # Cached fasta
	{
		my ($fasta_sequence_content_file) = ">Query\n$sequence";
		my ($print_fasta) = printStringIntoFile($fasta_sequence_content_file, $fasta_sequence_file);
		unless( defined $print_fasta ) {
			$logger->error("Can not create temporal file: $!\n");
		}
	}

	# Run pfamscan
	#my ($pfamscan_sequence_file) = $ws_cache.'/seq.pfam';
	my ($pfamscan_sequence_file) = $ws_cache.'/seq.'.$PROG_DB;
	unless(-e $pfamscan_sequence_file and (-s $pfamscan_sequence_file > 0) and ($CACHE_FLAG eq 'yes')) # Cached pfamscan
	{
		eval {
			my ($cmd) = "rm -rf $pfamscan_sequence_file && $RUN_PROGRAM -as -align -dir $PROG_DB_DIR -fasta $fasta_sequence_file -outfile $pfamscan_sequence_file &> /dev/null";
			$logger->debug("$cmd\n");                        
			my @out = `$cmd`;
			if (!(-e $pfamscan_sequence_file) or (-s $pfamscan_sequence_file <= 0) ) {
				$logger->error("Running pfamscan\n");
			}
		};
		$logger->error("Running pfamscan\n") if($@);
	}
	return $pfamscan_sequence_file;		
}

#E-value MUST be less than 0.00001
#Only Pfam A domains (this is already OK in pfamscan)
#The domains shouldnt overlap (this is already OK in pfamscan)
#
#The decision on whether a domain is damaged or not comes from:
#<hmm start> <hmm end> <hmm length> and the number of gaps
sub _parse_pfamscan($$)
{
	my ($pfamscan_sequence_file, $seq_id) = @_;
	
	my ($cutoffs);
	my ($alignment_list);
	my ($transcript_result) = '';
	my ($bitscore) = 0;
	my ($num_domains) = 0;
	my ($num_possibly_damaged_domains) = 0;
	my ($num_damaged_domains) = 0;
	my ($num_wrong_domains) = 0;

	# Get pfamscan output
	my ($result) = getStringFromFile($pfamscan_sequence_file);
	$result =~ s/^Query/$seq_id/mg; # changet the 'Query' identifier from firestar to the current identifier
	unless( defined $result ) {
		$logger->error("Can not open pfamscan result: $!\n");
	}
	my (@results) = split(/(\Q$seq_id\E[^\n]*\n#HMM[^\n]*\n#MATCH[^\n]*\n#PP[^\n]*\n#SEQ[^\n]*\n+)/,$result);
	
	foreach my $line (@results)
	{
		my ($alignment_report);

		if ( $line =~ /^(\Q$seq_id\E[^\n]*)/ )
		{
			$transcript_result .= '>'.$line;
				
			my ($value_line) = $1;
			# Check if value line matches PfamScan (Pfam-A) output, as
			# specified in Bio::Pfam::HMM::HMMResultsIO::write_ascii_out,
			# allowing for most non-numeric fields to be empty strings.
			if( $value_line=~/^
					(?<seq_id>\S+)
					\s
					\s{0,5}(?<alignment_start>[0-9]+)
					\s
					\s{0,5}(?<alignment_end>[0-9]+)
					\s
					\s{0,5}(?<envelope_start>[0-9]+)
					\s
					\s{0,5}(?<envelope_end>[0-9]+)
					\s
					(?<hmm_acc>\S*)\s{0,11}
					\s
					(?<hmm_name>\S*)\s{0,16}
					\s
					\s{0,7}(?<hmm_type>\S*)
					\s
					\s{0,4}(?<hmm_start>[0-9]+)
					\s
					\s{0,4}(?<hmm_end>[0-9]+)
					\s
					\s{0,4}(?<hmm_length>[0-9]+)
					\s
					\s{0,7}(?<bit_score>\S+)
					\s
					\s{0,8}(?<e_value>\S+)
					\s
					\s{0,2}(?<significance>[0-9]+)
					\s
					(?<clan>\S+)\s{0,7}
					\s
					\s(?<predicted_active_site_residues>predicted_active_site)?
					$/x )
			{
				my ($alignment_start) = $+{'alignment_start'};
				my ($alignment_end) = $+{'alignment_end'};
				my ($hmm_name) = $+{'hmm_name'};
				my ($hmm_start) = $+{'hmm_start'};
				my ($hmm_end) = $+{'hmm_end'};
				my ($hmm_length) = $+{'hmm_length'};
				my ($bit_score) = $+{'bit_score'};
				my ($e_value) = $+{'e_value'};

				# We can't use a PfamScan result that has an empty HMM
				# name, and it's a sign that something has gone wrong.
				if ( $hmm_name eq '' ) {
					$logger->error("Cannot process pfamscan result due to empty hmm name: '$value_line'\n");
				}

				# In a small proportion of PfamScan output files the HMM accession and HMM type
				# are parsed out as empty strings. We assign placeholder values to avoid confusion.
				my ($hmm_acc) = $+{'hmm_acc'} ne '' ? $+{'hmm_acc'} : 'Unknown_accession' ;
				my ($hmm_type) = $+{'hmm_type'} ne '' ? $+{'hmm_type'} : 'Unknown_type' ;

				if( $e_value < $PROG_EVALUE )
				{
					$alignment_report->{'alignment_start'} = $alignment_start;
					$alignment_report->{'alignment_end'} = $alignment_end;
					$alignment_report->{'alignment_length'} = ($alignment_end - $alignment_start)+1;
					$alignment_report->{'envelope_start'} = $+{'envelope_start'};
					$alignment_report->{'envelope_end'} = $+{'envelope_end'};
					$alignment_report->{'hmm_acc'} = $hmm_acc;
					$alignment_report->{'hmm_name'} = $hmm_name;
					$alignment_report->{'hmm_type'} = $hmm_type;
					$alignment_report->{'hmm_start'} = $hmm_start;
					$alignment_report->{'hmm_end'} = $hmm_end;
					$alignment_report->{'hmm_length'} = $hmm_length;
					$alignment_report->{'bit_score'} = $bit_score;
					$alignment_report->{'e_value'} = $e_value;
					$alignment_report->{'significance'} = $+{'significance'} if($+{'significance'});
					$alignment_report->{'clan'} = $+{'clan'} if($+{'clan'});

					if ( defined $+{'predicted_active_site_residues'} ) {
						$alignment_report->{'predicted_active_site_residues'} = $+{'predicted_active_site_residues'};
					}

					if ( $line =~ /(#HMM[^\n]*)\n*(#MATCH[^\n]*)\n*(#PP[^\n]*)\n*(#SEQ[^\n]*)\n*/ )			
					{
						my ($raw_hmm) = $1;
						my ($raw_seq) = $4;						
						
						# Get raw sequence
						my ($seq_hmm);
						if ( $raw_hmm =~ /#HMM\s*([^\n]*)/ ) {
							$seq_hmm = uc($1);
						}						my ($seq);
						if ( $raw_seq =~ /#SEQ\s*([^\n]*)/ ) {
							$seq = uc($1);
							$alignment_report->{'sequence'} = $seq; # uppercase for gap alignment
						}					
												
						# Get domain score
						my ($domain_score) = 0;
						my ($domain_score_c) = $hmm_length - $hmm_end;
						my ($domain_score_n) = $hmm_start - 1;
						if ( $domain_score_n > $domain_score_c ) {
							$domain_score = $domain_score_n;
						}
						else {
							$domain_score = $domain_score_c;
						}
						$alignment_report->{'score'} = $domain_score;
						my (@D_SCORE) = ( 4, 7, 12 );
						
						# Get the number of gaps that are together
						my ($get_gaps) = sub {
							my ($c) = shift;
							my ($str) = shift;
							my ($max) = 0;
							while ($str =~ /($c)(\1+)/g) {
								my ($num) = length($1.$2);
								if (length($1.$2) > $max) {
									$max = $num;
								}
							}
							return $max;
						};
						my ($num_gaps) = &$get_gaps('\-',$seq);
						my ($num_gaps_hmm) = &$get_gaps('\.',$seq_hmm);						
						my (@N_GAPS);
						if ( $hmm_length <= 100 )    { @N_GAPS = ( 3, 5, 7 ) }
						elsif ( $hmm_length <= 200 ) { @N_GAPS = ( 3, 6, 9 ) }
						else                         { @N_GAPS = ( 3, 7, 11 ) }
						
						# Get type of domain from hmm/seq length and gaps
						if ( ($domain_score <= $D_SCORE[0]) and ($num_gaps <= $N_GAPS[0]) and ($num_gaps_hmm <= $N_GAPS[0]) )
						{
							$alignment_report->{'type_domain'} = 'domain';
							$num_domains++;
						}
						elsif ( ($domain_score <= $D_SCORE[1]) and ($num_gaps <= $N_GAPS[1]) and ($num_gaps_hmm <= $N_GAPS[1]) )
						{
							$alignment_report->{'type_domain'} = 'domain_possibly_damaged';
							$num_possibly_damaged_domains++;
						}
						elsif ( ($domain_score <= $D_SCORE[2]) and ($num_gaps <= $N_GAPS[2]) and ($num_gaps_hmm <= $N_GAPS[2]) )
						{
							$alignment_report->{'type_domain'} = 'domain_damaged';
							$num_damaged_domains++;
						}
						else
						{
							$alignment_report->{'type_domain'} = 'domain_wrong';
							$num_wrong_domains++;								
						}
						
					}
					# Increase the total bitscore
					if (defined $bit_score) { $bitscore += $bit_score }
					
					# Save alignments
					push(@{$alignment_list}, $alignment_report);																
				}
				
			}
		}			
	}
	$cutoffs->{'bitscore'} = $bitscore;
	$cutoffs->{'domains'} = $alignment_list;
	$cutoffs->{'num_domains'} = $num_domains;
	$cutoffs->{'num_possibly_damaged_domains'} = $num_possibly_damaged_domains;
	$cutoffs->{'num_damaged_domains'} = $num_damaged_domains;
	$cutoffs->{'num_wrong_domains'} = $num_wrong_domains;
	$cutoffs->{'domain_integrity'} = _calc_domain_integrity_score($num_domains,
																  $num_possibly_damaged_domains,
																  $num_damaged_domains,
																  $num_wrong_domains);

	# Save result for each transcript
	$cutoffs->{'result'} = $transcript_result;		

	return $cutoffs;
}

# Get records by transcript
sub _get_record_annotations($)
{
	my ($transcript_report) = @_;
	my ($output_content) = '';

	$output_content .= "/*\n";
	$output_content .= " * SPADE\n";
	$output_content .= " *   predicts the presence of protein domains thanks to Pfamscan tool\n";
	$output_content .= " *   ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/README .\n";
	$output_content .= " * \n";
	$output_content .= " * The Pfam database (https://pfam.xfam.org/) is a large collection of protein families,\n";
	$output_content .= " * each represented by multiple sequence alignments and hidden Markov models (HMMs).\n";	
	$output_content .= " * \n";
	$output_content .= " */\n";

	$output_content .= "\n";
	$output_content .= "# ========================== #\n";
	$output_content .= "# Presence of protein domain #\n";
	$output_content .= "# ========================== #\n";
		
	while (my ($sequence_id, $trans_report) = each(%{$transcript_report}) )
	{
		if(	exists $trans_report->{'num_domains'} and defined $trans_report->{'num_domains'} and
			exists $trans_report->{'num_possibly_damaged_domains'} and defined $trans_report->{'num_possibly_damaged_domains'} and
			exists $trans_report->{'num_damaged_domains'} and defined $trans_report->{'num_damaged_domains'} and
			exists $trans_report->{'num_wrong_domains'} and defined $trans_report->{'num_wrong_domains'}			
		)
		{
			$output_content .= ">".$sequence_id."\t".
									$trans_report->{'domain_integrity'}."\t".
									$trans_report->{'bitscore'}."\t".
									$trans_report->{'num_domains'}."\t".
									$trans_report->{'num_possibly_damaged_domains'}."\t".
									$trans_report->{'num_damaged_domains'}."\t".
									$trans_report->{'num_wrong_domains'}."\n";
		
		}
		if ( exists $trans_report->{'domains'} and defined $trans_report->{'domains'} and scalar(@{$trans_report->{'domains'}}) > 0 ) {
			my (@sort_domains_trans_report) = sort {$a->{'alignment_start'} <=> $b->{'alignment_start'}} @{$trans_report->{'domains'}};
			foreach my $alignment_report (@sort_domains_trans_report)
			{
				# <type_domain> <domain score>
				# <alignment start> <alignment end> <envelope start> <envelope end>
				# <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value>
				# <significance> <clan> <predicted_active_site_residues> -optional values-
				# <external id, if applied> -optional values-
				$output_content .= 	$alignment_report->{'type_domain'}."\t".
									#$alignment_report->{'score'}."\t".
									$alignment_report->{'alignment_start'}."\t".
									$alignment_report->{'alignment_end'}."\t".
									$alignment_report->{'envelope_start'}."\t".
									$alignment_report->{'envelope_end'}."\t".
									$alignment_report->{'hmm_acc'}."\t".
									$alignment_report->{'hmm_name'}."\t".
									$alignment_report->{'hmm_type'}."\t".
									$alignment_report->{'hmm_start'}."\t".
									$alignment_report->{'hmm_end'}."\t".
									$alignment_report->{'hmm_length'}."\t".
									$alignment_report->{'bit_score'}."\t".
									$alignment_report->{'e_value'}."\t";
	
				if (exists $alignment_report->{'significance'} and defined $alignment_report->{'significance'} )
					{ $output_content .= $alignment_report->{'significance'}."\t"; }
				if (exists $alignment_report->{'clan'} and defined $alignment_report->{'clan'} )
					{ $output_content .= $alignment_report->{'clan'}."\t"; }
				if (exists $alignment_report->{'predicted_active_site_residues'} and defined $alignment_report->{'predicted_active_site_residues'} )
					{ $output_content .= $alignment_report->{'predicted_active_site_residues'}."\t"; }
	
				if(	exists $alignment_report->{'discarded'} and defined $alignment_report->{'discarded'} ){
					$output_content .= "[discarded]"."\t";
				}
				else {
					if (exists $alignment_report->{'external_id'} and defined $alignment_report->{'external_id'} )
						{ $output_content .= "[ext:".$alignment_report->{'external_id'}."]"."\t"; }
				}
	
				$output_content =~ s/\t$/\n/;								
			}
		}	
	}
	
	$output_content .= "\n";
	$output_content .= "# ============ #\n";
	$output_content .= "# Pfam reports #\n";
	$output_content .= "# ============ #\n";
		
	while (my ($sequence_id, $trans_report) = each(%{$transcript_report}) )
	{
		if ( exists $trans_report->{'result'} and defined $trans_report->{'result'} )
		{
			$output_content .= $trans_report->{'result'};	
		}
	}
	return $output_content;
}

# Get the annotations for the main isoform using the bitscore /* APPRIS */ ----------------
sub _get_appris_annotations_bitscore($)
{
	my ($transcript_report) = @_;

	my ($score_domains_transcripts);	
	my (@unknow_tag);
	my (@no_tag);	
	my ($output_content) = '';

	$output_content .= "\n";
	$output_content .= "# ============= #\n";
	$output_content .= "# Whole domains #\n";
	$output_content .= "# ============= #\n";

	# get the whole list of domains
	while (my ($sequence_id, $trans_report) = each(%{$transcript_report}) ) {
		if(	exists $trans_report->{'bitscore'} and defined $trans_report->{'bitscore'} ){	
			my ($bitscore) = $trans_report->{'bitscore'};
									
			push(@{$score_domains_transcripts->{$bitscore}}, $sequence_id);
		}
	}
	$logger->debug("DOMAIN:\n".Dumper($score_domains_transcripts)."\n");

	if ( defined $score_domains_transcripts )
	{
		# sort by descending order
		my (@sorted_domain_list) = sort { $b <=> $a } keys (%{$score_domains_transcripts}); 
		

		if ( scalar(@sorted_domain_list) > 0 )
		{			
			# Get the first and second big domain
			my ($biggest_num_domain) = $sorted_domain_list[0];
			my ($bigger_num_domain) = -1;
			if ( scalar(@sorted_domain_list) > 1 )
				{ $bigger_num_domain = $sorted_domain_list[1]; }

			while ( my ($num_domain, $domain_transcripts) = each(%{$score_domains_transcripts}) )
			{
				# We tag the transcript as UNKOWN whose num domains are biggest					
				if ( ($num_domain == $biggest_num_domain) or (($biggest_num_domain - $num_domain) <= $APPRIS_CUTOFF) )
				{
					foreach my $trans_id (@{$domain_transcripts})
					{
						push(@unknow_tag, $trans_id);
					}						
				}
				else
				{
					foreach my $trans_id (@{$domain_transcripts})
					{
						push(@no_tag, $trans_id);
					}
				}
			}
			
			# Add labels
			if ( scalar(@unknow_tag) == 1 )
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
	}	
	return $output_content;
}


main();



__END__

=head1 NAME

spade

=head1 DESCRIPTION

Run SPADE from PfamScan program

=head1 SYNOPSIS

spade

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

perl spade.pl

	--conf=../conf/pipeline.ini

	--input=examples/ENSG00000160796.faa
	
	--output=examples/ENSG00000160796.output

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
