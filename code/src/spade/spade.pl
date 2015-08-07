#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use Bio::SeqIO;
use Config::IniFiles;
use Data::Dumper;

use APPRIS::Utils::Logger;
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$WSPACE_BASE
	$WSPACE_CACHE
	$RUN_PROGRAM
	$PROG_DB_DIR
	$PROG_IN_SUFFIX
	$PROG_OUT_SUFFIX
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
my ($cfg) = new Config::IniFiles( -file =>  $config_file );
$LOCAL_PWD			= $FindBin::Bin;
$WSPACE_BASE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('SPADE_VARS', 'name').'/';
$WSPACE_CACHE		= $cfg->val('APPRIS_PIPELINE', 'workspace').'/'.$cfg->val('CACHE_VARS', 'name').'/';
$RUN_PROGRAM		= $cfg->val( 'SPADE_VARS', 'program');
$PROG_DB_DIR		= $ENV{APPRIS_PROGRAMS_DB_DIR};
$PROG_IN_SUFFIX		= 'faa';
$PROG_OUT_SUFFIX	= 'pfam';
$PROG_EVALUE		= $cfg->val( 'SPADE_VARS', 'evalue');
$APPRIS_CUTOFF		= $cfg->val( 'SPADE_VARS', 'cutoff');

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
	while ( my $seq = $fasta_object->next_seq() )
	{
		if ( $seq->id=~/([^|]*)/ )
		{			
			my ($sequence_id) = $1;
			if ( $sequence_id =~ /^ENS/ ) { $sequence_id =~ s/\.\d*$// }
            my ($sequence) = $seq->seq;
            $sequence_report->{$sequence_id} = $sequence;
            
            $logger->info("\n##$sequence_id ###############################\n");

			# Run blast
			$logger->info("\n##Running Pfamscan ---------------\n");
			my ($pfamscan_sequence_file) = _run_pfamscan($sequence_id, $sequence);
								                    			                
			# Parse blast
			$logger->info("\n##Parsing Pfamscan ---------------\n");                
			my ($pfam_report) = _parse_pfamscan($pfamscan_sequence_file, $sequence_id);
			$transcript_report->{$sequence_id} = $pfam_report if (defined $pfam_report);
		}
    }
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

# Get the first domain that we have found
sub _get_best_domain($$$)
{
	my ($sequence_id, $sequence, $transcript_report) = @_;
	
	my ($exits_domain);
	my ($exits_domain_no_redundancy);	
	my ($pfam_report);
	my ($alignment_list);
	my ($num_domains) = 0;
	my ($num_possibly_damaged_domains) = 0;
	my ($num_damaged_domains) = 0;
	my ($num_wrong_domains) = 0;
	
	# init the domains with the current domains
	if ( exists $transcript_report->{$sequence_id} ) {
		my ($pfam_report) = $transcript_report->{$sequence_id};
		if ( exists $pfam_report->{'domains'} and $pfam_report->{'domains'} )
		{
			foreach my $alignment_report (@{$pfam_report->{'domains'}})
			{
				my ($seq_domain) = $alignment_report->{'sequence'};
				$exits_domain->{$seq_domain}->{'e_value'} = $alignment_report->{'e_value'};
				$exits_domain->{$seq_domain}->{'report'} = $alignment_report;
			}
		}
	}
	
	# scan the rest of sequences
	while ( my ($sequence_id2, $pfam_report2) = each(%{$transcript_report}) )
	{
		next if ($sequence_id eq $sequence_id2); # Jump

		if ( exists $pfam_report2->{'domains'} and $pfam_report2->{'domains'} )
		{
			foreach my $alignment_report2 (@{$pfam_report2->{'domains'}})
			{
				my ($seq_domain2) = $alignment_report2->{'sequence'};
				my ($seq_domain2_nogaps) = $seq_domain2;
				$seq_domain2_nogaps =~ s/\-*//g;
				my ($seq_domain2_num_gaps) = $seq_domain2 =~ tr/\-//;
				
				# check if domain exists in other sequence (without gaps)
				my ($index_domain) = rindex($sequence, $seq_domain2_nogaps);
				if ( $index_domain != -1 )
				{
					if ( !(exists $exits_domain->{$seq_domain2}) )
					{ # first domain found that has a better e_value
						$exits_domain->{$seq_domain2}->{'e_value'} = $alignment_report2->{'e_value'};
						
						# get the domain from saying from who
						my (%alignment_report) = %{$alignment_report2};					
						$alignment_report{'external_id'} = $sequence_id2;
	
						# get the new range of alignment
						$alignment_report{'alignment_start'} = $index_domain + 1;
						$alignment_report{'alignment_end'} = $index_domain + length($seq_domain2) - $seq_domain2_num_gaps;
	
						$exits_domain->{$seq_domain2}->{'report'} = \%alignment_report;
					}
					elsif ( exists $exits_domain->{$seq_domain2} and 
							($alignment_report2->{'e_value'} < $exits_domain->{$seq_domain2}->{'e_value'})
					)
					{ # next domains with better evalue
						$exits_domain->{$seq_domain2}->{'e_value'} = $alignment_report2->{'e_value'};
						
						# get the domain from saying from who
						my (%alignment_report) = %{$alignment_report2};					
						$alignment_report{'external_id'} = $sequence_id2;
	
						# get the new range of alignment
						$alignment_report{'alignment_start'} = $index_domain + 1;
						$alignment_report{'alignment_end'} = $index_domain + length($seq_domain2) - $seq_domain2_num_gaps;
	
						$exits_domain->{$seq_domain2}->{'report'} = \%alignment_report;						
					}
				}
			}
		}
	}
	
	# delete the redundancy between alignments (we take the best domain)
	my (@aux_list_seq_domain) = sort { length($a) <=> length($b) } keys (%{$exits_domain});
	my (@aux_list_seq_domain2) = sort { length($b) <=> length($a) } keys (%{$exits_domain});
	foreach my $aux_seq_domain (@aux_list_seq_domain) {
		my ($aux_alignment_report) = $exits_domain->{$aux_seq_domain}->{'report'};
		my ($aux_type_domain) = 0;
		if ( $aux_alignment_report->{'type_domain'} eq 'domain_wrong') { $aux_type_domain = 1 }
		elsif ( $aux_alignment_report->{'type_domain'} eq 'domain_damaged') { $aux_type_domain = 2 }
		elsif ( $aux_alignment_report->{'type_domain'} eq 'domain_possibly_damaged') { $aux_type_domain = 3 }
		elsif ( $aux_alignment_report->{'type_domain'} eq 'domain') { $aux_type_domain = 4 }		
		my ($better_domain) = {
						'seq'	=> $aux_seq_domain,
						'type'	=> $aux_type_domain
		};
		foreach my $aux_seq_domain2 (@aux_list_seq_domain2) {
			my ($aux_alignment_report2) = $exits_domain->{$aux_seq_domain2}->{'report'};
			my ($aux_s_domain) = $better_domain->{'seq'};
			my ($aux_s_domain2) = $aux_seq_domain2;			
			# the seq domain has not to be the same
			if ( $aux_s_domain ne $aux_s_domain2 ) {

				# the seq domain has to be bigger
				if ( $aux_s_domain2 =~ /$aux_s_domain/ ) {
					# the seq domain has the same domain
					if ( $aux_alignment_report->{'hmm_name'} eq $aux_alignment_report2->{'hmm_name'} ) {
						# the seq domain has to be better type of aligment
						my ($aux_t_domain) = $better_domain->{'type'};
						my ($aux_t_domain2) = 0;	
						if ( $aux_alignment_report2->{'type_domain'} eq 'domain_wrong') { $aux_t_domain2 = 1 }
						elsif ( $aux_alignment_report2->{'type_domain'} eq 'domain_damaged') { $aux_t_domain2 = 2 }
						elsif ( $aux_alignment_report2->{'type_domain'} eq 'domain_possibly_damaged') { $aux_t_domain2 = 3 }
						elsif ( $aux_alignment_report2->{'type_domain'} eq 'domain') { $aux_t_domain2 = 4 }
						if ( $aux_t_domain2 >= $aux_t_domain ) {
							$better_domain->{'seq'} = $aux_s_domain2;
							$better_domain->{'type'} = $aux_t_domain2;
						}
					}
				}				
			}
		}
		my ($best_seq_domain) = $better_domain->{'seq'};
		$exits_domain_no_redundancy->{$best_seq_domain}->{'report'} = $exits_domain->{$best_seq_domain}->{'report'};
	}
		
	# report the alignment list
	while (my ($seq_domain2, $domain_report2) = each(%{$exits_domain_no_redundancy}) ) {
		my ($alignment_report2) = $domain_report2->{'report'};
		
		if( $alignment_report2->{'type_domain'} eq 'domain' ) {
			$num_domains++;
		}
		elsif( $alignment_report2->{'type_domain'} eq 'domain_possibly_damaged' ) {
			$num_possibly_damaged_domains++;
		}
		elsif ( $alignment_report2->{'type_domain'} eq 'domain_damaged' ) {
			$num_damaged_domains++;
		}
		elsif ( $alignment_report2->{'type_domain'} eq 'domain_wrong' ) {
			$num_wrong_domains++;
		}
		push(@{$alignment_list}, $alignment_report2);
	}
	$pfam_report->{'domains'} = $alignment_list;
	$pfam_report->{'num_domains'} = $num_domains;
	$pfam_report->{'num_possibly_damaged_domains'} = $num_possibly_damaged_domains;
	$pfam_report->{'num_damaged_domains'} = $num_damaged_domains;
	$pfam_report->{'num_wrong_domains'} = $num_wrong_domains;
	
	return $pfam_report;	
}

# Run pfamscan
sub _run_pfamscan($$)
{
	my ($sequence_id, $sequence) = @_;

	# Create temporal file
	my ($fasta_sequence_file) = $WSPACE_BASE.'/'.$sequence_id.'.'.$PROG_IN_SUFFIX;
	unless(-e $fasta_sequence_file and (-s $fasta_sequence_file > 0) ) # Cached fasta
	{
		my ($fasta_sequence_content_file) = ">$sequence_id\n$sequence";
		my ($print_fasta) = printStringIntoFile($fasta_sequence_content_file, $fasta_sequence_file);
		unless( defined $print_fasta ) {
			$logger->error("Can not create temporal file: $!\n");
		}
	}
	

	# Run pfamscan
	my ($pfamscan_sequence_file) = $WSPACE_CACHE.'/'.$sequence_id.'.'.$PROG_OUT_SUFFIX;
	my ($pfamscan_err_file) = $WSPACE_BASE.'/'.$sequence_id.'.'.$PROG_OUT_SUFFIX.'.err';              
	unless(-e $pfamscan_sequence_file and (-s $pfamscan_sequence_file > 0) ) # Cached pfamscan
	{
		eval {
			my ($cmd) = "$RUN_PROGRAM -as -align -dir $PROG_DB_DIR -fasta $fasta_sequence_file -outfile $pfamscan_sequence_file &> /dev/null";
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
	my ($num_domains) = 0;
	my ($num_possibly_damaged_domains) = 0;
	my ($num_damaged_domains) = 0;
	my ($num_wrong_domains) = 0;

	# Get pfamscan output
	my ($result) = getStringFromFile($pfamscan_sequence_file);
	unless( defined $result ) {
		$logger->error("Can not open pfamscan result: $!\n");
	}
	my (@results) = split(/($seq_id[^\n]*\n#HMM[^\n]*\n#MATCH[^\n]*\n#PP[^\n]*\n#SEQ[^\n]*\n+)/,$result);
	
	foreach my $line (@results)
	{
		my ($alignment_report);

		if ( $line =~ /^($seq_id[^\n]*)/ )
		{
			$transcript_result .= '>'.$line;
				
			my ($value_line) = $1;
			my (@value_list) = split(/\s+/,$value_line);
			if ( scalar(@value_list) > 12 )
			{
				my ($id) = $value_list[0];
				my ($alignment_start) = $value_list[1];
				my ($alignment_end) = $value_list[2];
				my ($envelope_start) = $value_list[3];
				my ($envelope_end) = $value_list[4];
				my ($hmm_acc) = $value_list[5];
				my ($hmm_name) = $value_list[6];
				my ($hmm_type) = $value_list[7];
				my ($hmm_start) = $value_list[8];
				my ($hmm_end) = $value_list[9];
				my ($hmm_length) = $value_list[10];
				my ($bit_score) = $value_list[11];
				my ($e_value) = $value_list[12];
					
				if($e_value < $PROG_EVALUE)
				{
					$alignment_report->{'alignment_start'} = $alignment_start;
					$alignment_report->{'alignment_end'} = $alignment_end;
					$alignment_report->{'envelope_start'} = $envelope_start;
					$alignment_report->{'envelope_end'} = $envelope_end;
					$alignment_report->{'hmm_acc'} = $hmm_acc;
					$alignment_report->{'hmm_name'} = $hmm_name;
					$alignment_report->{'hmm_type'} = $hmm_type;
					$alignment_report->{'hmm_start'} = $hmm_start;
					$alignment_report->{'hmm_end'} = $hmm_end;
					$alignment_report->{'hmm_length'} = $hmm_length;
					$alignment_report->{'bit_score'} = $bit_score;
					$alignment_report->{'e_value'} = $e_value;
					$alignment_report->{'significance'} = $value_list[13] if(defined $value_list[13]);
					$alignment_report->{'clan'} = $value_list[14] if(defined $value_list[14]);
					$alignment_report->{'predicted_active_site_residues'} = $value_list[15] if(defined $value_list[15]);

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
					push(@{$alignment_list}, $alignment_report);																
				}
			}
		}			
	}
	$cutoffs->{'domains'} = $alignment_list;
	$cutoffs->{'num_domains'} = $num_domains;
	$cutoffs->{'num_possibly_damaged_domains'} = $num_possibly_damaged_domains;
	$cutoffs->{'num_damaged_domains'} = $num_damaged_domains;
	$cutoffs->{'num_wrong_domains'} = $num_wrong_domains;

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
	$output_content .= " * The Pfam database (http://pfam.sanger.ac.uk) is a large collection of protein families,\n"; 
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
									$alignment_report->{'score'}."\t".
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
	
				if (exists $alignment_report->{'external_id'} and defined $alignment_report->{'external_id'} )
					{ $output_content .= "[".$alignment_report->{'external_id'}."]"."\t"; }
	
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

# Get the annotations for the main isoform /* APPRIS */ ----------------
sub _get_appris_annotations($)
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
	while (my ($sequence_id, $trans_report) = each(%{$transcript_report}) )
	{
		if(	exists $trans_report->{'num_domains'} and defined $trans_report->{'num_domains'} and
			exists $trans_report->{'num_possibly_damaged_domains'} and defined $trans_report->{'num_possibly_damaged_domains'} and
			exists $trans_report->{'num_damaged_domains'} and defined $trans_report->{'num_damaged_domains'} and
			exists $trans_report->{'num_wrong_domains'} and defined $trans_report->{'num_wrong_domains'}
		){	
			my ($domain_score) = 	($trans_report->{'num_domains'}*1) +
									($trans_report->{'num_possibly_damaged_domains'}*0.75) +
									($trans_report->{'num_damaged_domains'}*0.5) +
									($trans_report->{'num_wrong_domains'}*0.25);
									
			push(@{$score_domains_transcripts->{$domain_score}}, $sequence_id);
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
