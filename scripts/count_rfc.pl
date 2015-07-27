#!/usr/bin/perl -W

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Bio::SeqIO;

use APPRIS::Registry;
use APPRIS::Utils::File qw( printStringIntoFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info warning throw );

use Data::Dumper;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_APPRIS_DB_FILE
	$ALIGN_VARS
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_APPRIS_DB_FILE	= $LOCAL_PWD.'/conf/apprisdb.ini';
$ALIGN_VARS = {
	'ucsc' => {
		'main_specie'		=> 'hg19',
		'alter_species'		=> ['panTro2','rheMac2','mm9','canFam2'],
		'alter_sp_val'		=> [3,2,1,1],
		'alter_sp_list'		=> {'panTro2'=>'alter','rheMac2'=>'alter','mm9'=>'alter','canFam2'=>'alter'},
		'types'				=> ['filter','kalign','prank']
	},
	'compara' => {
		'main_specie'		=> 'homo_sapiens',
		'alter_species'		=> ['pan_troglodytes','macaca_mulatta','mus_musculus','canis_familiaris'],
		'alter_sp_val'		=> [3,2,1,1],
		'alter_sp_list'		=> {'pan_troglodytes'=>'alter','macaca_mulatta'=>'alter','mus_musculus'=>'alter','canis_familiaris'=>'alter'}
	}
};

# Input parameters
my ($position) = undef;
my ($inpath) = undef;
my ($output_file) = undef;
my ($t_align) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'position=s'		=> \$position,
	'inpath=s'			=> \$inpath,
	't-align=s'			=> \$t_align,
	'output=s'			=> \$output_file,
	'apprisdb-conf=s'	=> \$apprisdb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);
unless (
	defined $inpath and
	defined $output_file
){
	print `perldoc $0`;
	exit 1;
}

# Optional arguments
# get vars of appris db
unless ( defined $apprisdb_conf_file ) {
	$apprisdb_conf_file = $CONFIG_INI_APPRIS_DB_FILE;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();


#####################
# Method prototypes #
#####################
sub get_aligns_from_pos($);
sub get_aligns_from_dir($);
sub get_align_scores($$$\$);
sub get_rfc_scores($);
sub get_best_rfc_scores($);
sub get_specie_aligns($);
sub count_rfc_scores($);
sub count_rfc($$);
sub get_codons($$$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# get aligns from ...
	my ($gene_report);
	if ( defined $position ) {
		$position = undef if ( $position eq 'genome' );
		$logger->debug("-- get aligns from position -------\n");
		$gene_report = get_aligns_from_pos($position);		
	}
	else {
		$logger->debug("-- get aligns from inpath -------\n");
		$gene_report = get_aligns_from_dir($inpath);		
	}
#$logger->debug("CHR_REPORT:\n".Dumper($gene_report)."\n");
	
	# Get the RFC score for principal isoform per gene
	$logger->debug("-- get the RFC score for principal isoform per gene -------\n");
	my ($rfc_scores) = get_rfc_scores($gene_report);
#$logger->debug("RFC_SCORES:\n".Dumper($rfc_scores)."\n");

	# Get the RFC output
	$logger->debug("-- get the best RFC scores -------\n");
	my ($p_log) = get_best_rfc_scores($rfc_scores);
	unless ( defined $p_log) {
		die ("ERROR: printing ");		
	}

	$logger->finish_log();
	
	exit 0;	
	
}

sub get_aligns_from_pos($)
{
	my ($ipos) = @_;
	my ($report);

	# APPRIS Registry
	my ($registry) = APPRIS::Registry->new();
	$registry->load_registry_from_ini(-file	=> $apprisdb_conf_file);	
	my ($chromosome) = $registry->fetch_basic_by_region($ipos);
		
	foreach my $gene (@{$chromosome}) {
		my ($gene_id) = $gene->stable_id;
		my ($g_version) = $gene->version;
		my ($chr) = $gene->chromosome;
		my ($biotype) = $gene->biotype;		
		#$logger->debug("-- $gene_id: ");
		
		my ($g_report);
		my ($g_c_report);
		my ($pi_id);
		my ($pi_c_score);
		foreach my $transcript (@{$gene->transcripts}) {	
			my ($transcript_id) = $transcript->stable_id;
			#$logger->debug("$transcript_id ");
			
			if ($transcript->translate) {
				my ($translation_length) = length($transcript->translate->sequence);
				my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id,'appris');
				my ($c_analysis) = $registry->fetch_analysis_by_stable_id($transcript_id,'corsair');
				my ($corsair_score) = 0;
				if ( $c_analysis and $c_analysis->corsair and $c_analysis->corsair->score ) {
					$corsair_score = $c_analysis->corsair->score;
					if ( defined $corsair_score ) {
						push(@{$g_c_report->{$corsair_score}}, {
							'chr'						=> $chr,
							'gene_id'					=> $gene_id,
							'version'					=> $g_version,
							'gene_type'					=> $biotype,
							'transcript_id'				=> $transcript_id,
							'translation_length'		=> $translation_length,
							'cosair_score'				=> $corsair_score,
							'appris_annot'				=> 'BEST_CORSAIR',
						});
					}
				}
				# get the principal isoform or the longest btter isoform
				if ( $analysis and $analysis->appris and $analysis->appris->principal_isoform_signal ) {
					my ($appris_annot) = $analysis->appris->principal_isoform_signal;
					if ( $appris_annot eq 'YES' ) {
						$g_report->{'chr'}					= $chr;
						$g_report->{'gene_id'}				= $gene_id;
						$g_report->{'version'}				= $g_version;
						$g_report->{'gene_type'}			= $biotype;
						$g_report->{'transcript_id'}		= $transcript_id;
						$g_report->{'translation_length'}	= $translation_length;
						$g_report->{'appris_annot'}			= 'PRINCIPAL';
						$pi_id = $transcript_id;
						if ( defined $corsair_score ) {
							$g_report->{'cosair_score'}		= $corsair_score;
							$pi_c_score = $corsair_score;
						}
					}
					elsif ( $appris_annot eq 'UNKNOWN' ) {
						unless ( defined $g_report ) {
							$g_report->{'chr'}					= $chr;
							$g_report->{'gene_id'}				= $gene_id;
							$g_report->{'version'}				= $g_version;
							$g_report->{'gene_type'}			= $biotype;
							$g_report->{'transcript_id'}		= $transcript_id;
							$g_report->{'translation_length'}	= $translation_length;
							$g_report->{'appris_annot'}			= 'LONGEST';
							$pi_id = $transcript_id;
							if ( defined $corsair_score ) {
								$g_report->{'cosair_score'}		= $corsair_score;
								$pi_c_score = $corsair_score;
							}							
						}
						else {
							if ( $translation_length > $g_report->{'translation_length'} ) { # maximum length
								$g_report->{'chr'}					= $chr;
								$g_report->{'gene_id'}				= $gene_id;
								$g_report->{'version'}				= $g_version;
								$g_report->{'gene_type'}			= $biotype;
								$g_report->{'transcript_id'}		= $transcript_id;
								$g_report->{'translation_length'}	= $translation_length;
								$g_report->{'appris_annot'}			= 'LONGEST';
								$pi_id = $transcript_id;
								if ( defined $corsair_score ) {
									$g_report->{'cosair_score'}		= $corsair_score;
									$pi_c_score = $corsair_score;
								}
							}
						}
												
					}
				}
			}
		}
		if ( defined $g_report ) {
			# add the principal isoform or the longest sequence
			push(@{$report}, $g_report);
			
			# add the isoforms with better corsair score
			my (@sort_cosair_scores) = sort({$b <=> $a} keys(%{$g_c_report}));
			if ( scalar(@sort_cosair_scores) > 0 ) {
				my ($best_corsair) = $sort_cosair_scores[0];
				if ( defined $pi_c_score and defined $best_corsair ) {
					if ( $best_corsair >= $pi_c_score ) {
						foreach my $c_report (@{$g_c_report->{$best_corsair}}) {
							if ( $c_report->{'transcript_id'} ne $pi_id ) {
								push(@{$report}, $c_report);
							}
						}	
					}					
				}
			}
		}		
		#$logger->debug("\n");
	}
	return $report;
}

sub get_aligns_from_dir($)
{
	my ($dir) = @_;
	my ($report);
		
	opendir(INPUT_DIR, $dir) || die "can't opendir $dir: $!";
    my (@gene_list) = grep { !/^\./ && -d "$dir/$_" } readdir(INPUT_DIR);
    closedir INPUT_DIR;
     
	foreach my $g_id (@gene_list) {
		my ($g_report);
		my ($g_dir) = "$dir/$g_id";
		if ( $g_id =~ /^([^\.]*)\.([^\$])$/ ) {
			my ($gene_id) = $1;
			my ($g_version) = $2;
			
			my (@chr_list) = `cut -f 1 $g_dir/$g_id.annot.gtf | sort | uniq`;
			my ($chr) = $chr_list[0];
			$chr =~ s/^\s*//;$chr =~ s/\s*$//;$chr =~ s/^chr//;
			
			my (@biotype_list) = `grep "gene_type" $g_dir/$g_id.annot.gtf`;
			my ($biotype);
			if ( $biotype_list[0] =~ /gene_type "([^\"]*)"/ ) {
				$biotype = $1;
			}
						
			opendir(G_DIR, $g_dir) || die "can't opendir $g_dir: $!";
		    my (@transc_list) = grep { !/^\./ && /\.align\.faa$/ } readdir(G_DIR);
		    closedir G_DIR;

		    foreach my $t_id (@transc_list) {
				if ( $t_id =~ /^([^\.]*)\.align\.faa$/ ) {
					my ($transc_id) = $1;
					my ($g_report) = {
						'chr'				=> $chr,
						'gene_id'			=> $gene_id,
						'version'			=> $g_version,
						'gene_type'			=> $biotype,
						'transcript_id'		=> $transc_id,
					};
					push(@{$report}, $g_report);					
				}
		    }
		#$logger->debug("-- $gene_id: ");			
		}
		#$logger->debug("\n");
	}
	return $report;
}

sub get_rfc_scores($)
{
	my ($report) = @_;
	my ($rfc_scores);
		
	if ( defined $report ) {
		foreach my $g_rep (@{$report}) {
			my ($chr) = $g_rep->{'chr'};
			my ($gene_id) = $g_rep->{'gene_id'};
			my ($g_version) = $g_rep->{'version'};
			my ($gene_idv) = $gene_id.'.'.$g_version;
			my ($transc_id) = $g_rep->{'transcript_id'};
			my ($genetype) = $g_rep->{'gene_type'};
			
			my ($align_file);
			my ($type_align);
			my ($aln_report);
			if ( $t_align eq 'ucsc') {
				
				#$align_file = $inpath.'/'."chr$chr".'/inertia/'.$transc_id.'.filter.faa';
				#$rfc_scores->{$chr}->{$gene_id}->{$transc_id} = get_align_scores($g_rep, $t_align, $align_file);
				foreach my $t_a (@{$ALIGN_VARS->{$t_align}->{'types'}} ) {					
					if ( $t_a eq 'filter') {
						$align_file = $inpath.'/'."chr$chr".'/inertia/'.$transc_id.'.filter.faa';					
					}
					elsif ( $t_a eq 'kalign') {
						$align_file = $inpath.'/'."chr$chr".'/inertia/'.$transc_id.'.kalign.faa';
					}					
					elsif ( $t_a eq 'prank') {
						$align_file = $inpath.'/'."chr$chr".'/inertia/'.$transc_id.'.prank.output.2.fas';					
					}
					$logger->debug("$align_file\n");					
					get_align_scores($g_rep, $t_a, $align_file, $aln_report);
				}							
			}
			elsif ( $t_align eq 'compara') {
				$type_align = 'compara';
				$align_file = $inpath.'/'.$gene_idv.'/'.$transc_id.'.align.faa';
				$logger->debug("$align_file\n");
				get_align_scores($g_rep, $t_align, $align_file, $aln_report);
			}
			if ( defined $report ) {
				$rfc_scores->{$chr}->{$gene_id}->{$transc_id} = $aln_report;
			}			
		}
	}
	return $rfc_scores;
}

sub get_align_scores($$$\$)
{
	my ($g_rep, $type, $file, $report) = @_;
	
	my ($align_seqs) = get_specie_aligns($file);
	#$logger->debug("ALIGN_SEQS:\n".Dumper($align_seqs)."\n");
	
	my ($scores) = count_rfc_scores($align_seqs);
	#$logger->debug("RFC_SCORES:\n".Dumper($scores)."\n");
	
	foreach my $alter_id (@{$ALIGN_VARS->{$t_align}->{'alter_species'}} ) { # ucsc or compara
		my ($rfc_score) = 'NA';
		my ($nbp) = 'NA';					
		if ( defined $scores ) {
			if ( exists $scores->{$alter_id} ) {
				if ( exists $scores->{$alter_id}->{'score'}) {
					$rfc_score = $scores->{$alter_id}->{'score'};
				}												
				if ( exists $scores->{$alter_id}->{'nbp'}) {
					$nbp  = $scores->{$alter_id}->{'nbp'};
				}												
			}												
		}
		$$report->{'alter_species'}->{$alter_id}->{$type} = { 'rfc_score' => $rfc_score, 'nbp' => $nbp };
	}
	$$report->{'nbp'} = $align_seqs->{'main'}->{'len'};
	$$report->{'gene_type'} = $g_rep->{'gene_type'};
	return $report;
}

sub get_specie_aligns($)
{
	my ($file) = @_;
	my ($data) = {
		'main'		=> undef,
		'alter'		=> undef, 
	};

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			my ($seq_id) = $seq->id;
			if ( $seq_id eq $ALIGN_VARS->{$t_align}->{'main_specie'} ) {
				my ($s) = $seq->seq;
				(my $s2 = $s) =~ s/-*//g;
				my ($n) = length($s2);				
				$data->{'main'} = { 'id' => $seq_id, 'seq' => $s, 'len' => $n };
			}
			else {
				if ( exists $ALIGN_VARS->{$t_align}->{'alter_sp_list'}->{$seq_id} ) {
					$data->{'alter'}->{$seq_id} = $seq->seq;
				}				
			}			
		}		
	}
	return $data;
}

sub count_rfc_scores($)
{
	my ($data) = @_;
	my ($rfc_scores);
	
	if ( defined $data->{'main'} ) {
		my ($main_report) = $data->{'main'};		
		my ($main_id) = $main_report->{'id'};
		my ($main_seq) = $main_report->{'seq'};
		my ($alter_report) = $data->{'alter'};
		foreach my $alter_id (@{$ALIGN_VARS->{$t_align}->{'alter_species'}} ) {
#$logger->debug(">ALTER_ID: $alter_id\n");			
			my ($score) = 'NA';
			my ($ntlen) = '-';
			if ( exists $alter_report->{$alter_id} ) {
				my ($alter_seq) = $alter_report->{$alter_id};
#$logger->debug(">MAIN:$main_id".":".length($main_seq)."\n$main_seq\n");
#$logger->debug(">MAIN_NOGAPS:$main_id:$nbp\n$main_seq_nogaps\n");
#$logger->debug(">ALTER:$alter_id".":".length($alter_seq)."\n$alter_seq\n");

				# get codon for every start offset
				my ($count) = 0;
				my ($nbp) = 0;
				foreach my $a_frame (1,2,3) {
					my ($main_codons) = undef;
					my ($alter_codons) = undef;		
					($nbp, $main_codons, $alter_codons) = get_codons($main_seq, $alter_seq, $a_frame);
					#(my $main_seq_nogaps = $main_seq) =~ s/-*//g;
					#my ($nbp) = length($main_seq_nogaps);
					if ( ($main_codons ne '') and ($alter_codons ne '') ) {
#$logger->debug("ALTER:$alter_id:$a_frame\n");
						my ($c) = count_rfc($main_codons, $alter_codons);
#$logger->debug("COUNT_RFC: $c/$nbp\n");
#$logger->debug(">$main_id\n$main_codons\n\n");
#$logger->debug(">$alter_id\n$alter_codons\n\n");
						unless ( defined $count ) {
							$count = $c;
						}
						else {
							if ( $c > $count ) {
								$count = $c;
							}
						}
						$ntlen = 'EQUAL';
					}
					else {
						warning("NOT EQUAL THE SEQUENCES\n");
						$ntlen = 'NOT';
					}
				}
				if ( defined $nbp and ($nbp != 0) and defined $count and ($count != 0) ) {
					$score = sprintf "%.2f", ($count / $nbp);
					if ( $score == 0 ) { $score = 'NA'; }			
#$logger->debug("SCORE:$alter_id: $count / $nbp = $score\n");
				}
				else {
					$score = 'NA';
				}
			}
			$rfc_scores->{$alter_id}->{'score'} = $score;
			$rfc_scores->{$alter_id}->{'nbp'} = $ntlen;
		}
	}
	return $rfc_scores;
}

sub get_den_length($$)
{
	my ($main, $alter) = @_;
	my ($nbp) = 0;
	my (@main_chars) = split(//,$main);
	my (@alter_chars) = split(//,$alter);
	if ( scalar(@main_chars) == scalar(@alter_chars) ) {
		for (my $i=0; $i < scalar(@main_chars); $i++) {
			my ($main_nt) = $main_chars[$i];
			my ($alter_nt) = $alter_chars[$i];
			$nbp++ if ( ($main_nt ne '-') and ($alter_nt ne '-') );
		}
	}	
	return $nbp;	
}

sub get_codons($$$)
{
	my ($main, $alter, $a_frame) = @_;
	my ($m_codon_list) = '';
	my ($a_codon_list) = '';
	my ($t_nbp) = 0;
	
	
	my (@main_chars) = split(//,$main);
	my (@alter_chars) = split(//,$alter);
	if ( scalar(@main_chars) == scalar(@alter_chars) ) {
		my ($nb) = 1;
		my ($main_codon) = 1;
		my ($alter_codon) = $a_frame;
		for (my $i=0; $i < scalar(@main_chars); $i++) {
			# nt
			my ($main_nt) = $main_chars[$i];
			my ($alter_nt) = $alter_chars[$i];
			$t_nbp++ if ( ($main_nt ne '-') and ($alter_nt ne '-') );
			
			# init num. codon
			if ( $main_codon >= 4 ) { $main_codon = 1; }
			if ( $alter_codon >= 4 ) { $alter_codon = 1; }
			my ($m_codon) = $main_codon;
			my ($a_codon) = $alter_codon;
			if ( $main_nt eq '-' ) { $m_codon = '-'; }
			if ( $alter_nt eq '-' ) { $a_codon = '-'; }
			
			$m_codon_list .= "$m_codon";
			$a_codon_list .= "$a_codon";
						
			# increase counters
			if ( $main_nt ne '-' ) { $main_codon++; }
			if ( $alter_nt ne '-' ) { $alter_codon++; }
			$nb++;				
		}	
	}
	return ($t_nbp, $m_codon_list,$a_codon_list);
}

sub count_rfc($$)
{
	my ($main, $alter) = @_;
	my ($count) = 0;
		
	while ( $main =~ /([^\-]*)\-*/g ) {
		my ($main_init) = $-[1];
		my ($main_end) = $+[1];
		my ($main_length) = $main_end - $main_init;
		my ($gap_init) = $main_end;
		my ($gap_end) = $+[0];
		my ($gap_length) = $gap_end - $gap_init;
		
#$logger->debug("MAIN_INIT: ".$main_init."\n");
#$logger->debug("MAIN_END: ".$main_end."\n");
#$logger->debug("MAIN_LENGTH: ".$gap_init."\n");

		my ($main_seq) = substr($main,$main_init,$main_length);
		my ($alter_seq) = substr($alter,$main_init,$main_length);
		my ($gap_seq) = substr($main,$gap_init, $gap_length);
#$logger->debug("MAIN_SEQ: ".$main_seq."\n");
#$logger->debug("ALTE_SEQ: ".$alter_seq."\n");

		while ( $alter_seq =~ /([^\-]*)/g ) {
			my ($alter_init) = $-[1];
			my ($alter_end) = $+[1];
			my ($alter_length) = $alter_end - $alter_init;
#$logger->debug("\talter_init: ".$alter_init."\n");
#$logger->debug("\talter_end: ".$alter_end."\n");
#$logger->debug("\talter_len: ".$alter_length."\n");
			my ($main_read) = substr($main_seq,$alter_init,$alter_length);
			my ($alter_read) = substr($alter_seq,$alter_init,$alter_length);
			if ( $alter_read ne '' ) {
#$logger->debug("\tmain_read: ".$main_read."\n");
#$logger->debug("\talte_read: ".$alter_read."\n");
				if ( $alter_read =~ /^$main_read$/ ) {
					my ($num_rfc) = length($alter_read);
#$logger->debug("\tNUM_RFC: ".$num_rfc."\n");
					$count += $num_rfc;				
				}
			}
		}
#$logger->debug("\n");	
	}
	return $count;
}

sub get_best_rfc_scores($)
{
	my ($scores) = @_;
	my ($rfc_output) = '';
	my ($p_log);
	if ( defined $scores ) {
		my ($rfc_header) = '';
		foreach my $alter_id (@{$ALIGN_VARS->{$t_align}->{'alter_species'}} ) {
			$rfc_header .= $alter_id."\t";
		}
		$rfc_header =~ s/\t$//;
		$rfc_output .= '#'.'chr'."\t".'gene_id'."\t".'transc_id'."\t".$rfc_header."\t".'gene_type'."\t".'gene_len'."\n";		
		foreach my $chr (sort {$a <=> $b} keys(%{$scores}) ) {
			while (my ($gene_id, $g_rep) = each(%{$scores->{$chr}}) ) {
				my ($best_rfc_rep);
				my ($best_rfc_rep2);
				my ($best_transc_rfc);
				my ($gene_types);
				my ($nbps);
				while (my ($transc_id, $transc_rep) = each(%{$g_rep}) ) {
					$gene_types->{$gene_id}->{$transc_id} = $transc_rep->{'gene_type'};
					$nbps->{$gene_id}->{$transc_id} = $transc_rep->{'nbp'};
					$best_transc_rfc->{$transc_id} = 0;
					foreach my $alter_id (@{$ALIGN_VARS->{$t_align}->{'alter_species'}} ) {
						
						if ( $t_align eq 'ucsc' ) {
							my ($rfc_max) = '-';
							my ($nbp_max) = '-';
							while ( my ($t_a, $rfc_rep) = each(%{$transc_rep->{'alter_species'}->{$alter_id}}) ) {
								my ($rfc_score) = $rfc_rep->{'rfc_score'};
								my ($nbp) = $rfc_rep->{'nbp'};
								if ( defined $rfc_score and ($rfc_score ne 'NA') ) {
									if ( defined $rfc_max ) {
										if ( $rfc_score > $rfc_max ) {
											$rfc_max = $rfc_score;
											$nbp_max = $nbp;
										}									
									}
									else {
										$rfc_max = $rfc_score;
										$nbp_max = $nbp;
									}
								}								
							}			
							#unless ( defined $rfc_max ) {
							#	$rfc_max = 'NA';
							#}
							$best_rfc_rep->{$transc_id}->{$alter_id}->{'score'} = $rfc_max;
							$best_rfc_rep->{$transc_id}->{$alter_id}->{'nbp'} = $nbp_max;
							push(@{$best_rfc_rep2->{$alter_id}->{$rfc_max}}, $transc_id);
						}
						elsif ( $t_align eq 'compara' ) {
							my ($rfc_rep) = $transc_rep->{'alter_species'}->{$alter_id}->{$t_align};						
							my ($rfc_score) = $rfc_rep->{'rfc_score'};
							my ($nbp) = $rfc_rep->{'nbp'};
							$best_rfc_rep->{$transc_id}->{$alter_id}->{'score'} = $rfc_score;
							$best_rfc_rep->{$transc_id}->{$alter_id}->{'nbp'} = $nbp;
							push(@{$best_rfc_rep2->{$alter_id}->{$rfc_score}}, $transc_id);							
						}
					}
				}
#$logger->debug("BEST_RFC_SCORES:\n".Dumper($best_rfc_rep)."\n");
#$logger->debug("BEST_RFC_SCORES_2:\n".Dumper($best_rfc_rep2)."\n");
				
				# get the best transc with the best rfc scores
				my ($found) = 0;
				for (my $i=0; $i<scalar(@{$ALIGN_VARS->{$t_align}->{'alter_species'}});$i++) {
					my ($alter_id) = $ALIGN_VARS->{$t_align}->{'alter_species'}->[$i];
					my ($best_rfc);
					my (@sort_rfc_rep2) = sort({$b <=> $a} keys(%{$best_rfc_rep2->{$alter_id}}));
					if ( scalar(@sort_rfc_rep2) > 0) {
						my ($rfc_max) = $sort_rfc_rep2[0];
#$logger->debug("RFC_SCORES:$alter_id: ");						
						foreach my $t_id (@{$best_rfc_rep2->{$alter_id}->{$rfc_max}}) {
#$logger->debug("$t_id ");							
							$best_transc_rfc->{$t_id} += $ALIGN_VARS->{$t_align}->{'alter_sp_val'}->[$i];
						}
					}
				}				
#$logger->debug("BEST_RFC_SCORES:\n".Dumper($best_transc_rfc)."\n");
				if ( defined $best_transc_rfc ) {
					my (@best_val) = sort { $best_transc_rfc->{$b} <=> $best_transc_rfc->{$a} } keys(%{$best_transc_rfc});
					my ($best_t_id) = $best_val[0]; 
#$logger->debug("BEST_VALS:\n".Dumper(@best_val)."\n");
					$rfc_output .= $chr."\t".$gene_id."\t".$best_t_id."\t";
					foreach my $alter_id (@{$ALIGN_VARS->{$t_align}->{'alter_species'}} ) {
						my ($rfc_score) = $best_rfc_rep->{$best_t_id}->{$alter_id}->{'score'};
						if ( defined $rfc_score ) {
							#$rfc_output .= $rfc_score."\t";
							$rfc_output .= $rfc_score."|";
							$rfc_output .= $best_rfc_rep->{$best_t_id}->{$alter_id}->{'nbp'}."\t";
						}
						else {
							#$rfc_output .= 'NA'."\t";
							$rfc_output .= 'NA'."|";
							$rfc_output .= $best_rfc_rep->{$best_t_id}->{$alter_id}->{'nbp'}."\t";
						}
					}
					#$rfc_output =~ s/\t$/\n/;
					$rfc_output .= $gene_types->{$gene_id}->{$best_t_id}."\t";
					$rfc_output .= $nbps->{$gene_id}->{$best_t_id}."\n";					
				}
				else {
					warning("NO PRINCIPAL: $gene_id");
				}
			}
		}
	}
	if ($rfc_output ne '') {
		$p_log = printStringIntoFile($rfc_output,$output_file);
	}
	return $p_log;
}

main();


1;

__END__

=head1 NAME

count_rfc

=head1 DESCRIPTION

Count the rfc of pairwise alignment of principal isoforms

=head1 SYNOPSIS

count_rfc

=head2 Required arguments:

	--chr  <Genomic region>	
	or
	--inpath  <Alignment dir>

	--output <Output file that has the main isoforms>	

=head2 Optional arguments:

	--t-align= <Type of alignment: [ucsc,compara]>	

=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl count_rfc.pl
	
		--inpath=/home/jmrodriguez/projects/Encode/appris_rel12/features/rfc/rand_e_genes
		
		--output=/home/jmrodriguez/projects/Encode/appris_rel12/features/rfc/rand_e_genes.txt

=head1 EXAMPLE

	perl count_rfc.pl
	
		--position=21
		
		--inpath=/home/jmrodriguez/projects/Encode/appris_rel12/features/annotations
		
		--output=/home/jmrodriguez/projects/Encode/appris_rel12/features/rfc/rand_e_genes.txt
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
