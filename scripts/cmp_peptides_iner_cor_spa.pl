#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use List::Util qw(min max sum);
use APPRIS::Utils::File qw( getStringFromFile getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;
use Data::Dumper;

###################
# Global variable #
###################
use vars qw(
	$CONFIG_INI_APPRIS_DB_FILE
	$WRONG_GENES_FOR_INERTIA
	$CDS_GENCODE12
	$READTHROUGH_GENES_GENCODE12
	$GENE_MODEL_PROBLEMS_GENCODE12
	$BONAFIDE_AS_GENES_GENCODE12
	$NUM_PEPS_GENES_GENCODE12
	$NUM_SPECIES_MAF_ALIGN_GENCODE12
);

$CONFIG_INI_APPRIS_DB_FILE	= $ENV{APPRIS_SCRIPTS_CONF_DIR}.'/apprisdb.ini';
$WRONG_GENES_FOR_INERTIA = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/wrong_genes_for_inertia.txt');
$CDS_GENCODE12 = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/features/num_cds.rel12.txt');
$READTHROUGH_GENES_GENCODE12 = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/AllG12_RT.txt');
$GENE_MODEL_PROBLEMS_GENCODE12 = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/GeneModelProblems.txt');
$BONAFIDE_AS_GENES_GENCODE12 = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/BonaFideASGenes.txt');
my ($file_num_peps_genes_gencode12) = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/num_peps_genes.gencode12.txt');
while ( $file_num_peps_genes_gencode12 =~ /([^\t]*)\t+([^\t]*)\t+(\d*)\n+/mg ) {
	my ($gene_name) = $1;
	my ($gene_id) = $2;
	my ($num_peptides) = $3;
	$NUM_PEPS_GENES_GENCODE12->{$gene_id} = $num_peptides;
}
my ($file_num_aligned_species_gencode12) = getStringFromFile('/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/num_species_in_maf_aligns.gencode12.txt');
while ( $file_num_aligned_species_gencode12 =~ /([^\t]*)\t+([^\t]*)\t+(\d*)\n+/mg ) {
	my ($chr) = $1;
	my ($transc_id) = $2;
	my ($num_aligned_species) = $3;
	$NUM_SPECIES_MAF_ALIGN_GENCODE12->{$transc_id} = $num_aligned_species;
}

# Input parameters
my ($in_proteo) = undef;
my ($in_inertia) = undef;
my ($in_corsair) = undef;
my ($in_spade) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'in-proteo=s'		=> \$in_proteo,
	'in-inertia=s'		=> \$in_inertia,
	'in-corsair=s'		=> \$in_corsair,
	'in-spade=s'		=> \$in_spade,
	'outfile=s'			=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless( defined $in_proteo and defined $in_inertia and defined $in_corsair and defined $in_spade and defined $outfile )
{
	print `perldoc $0`;
	exit 1;
}

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	
	my ($in_proteo_txt) = getTotalStringFromFile($in_proteo);
	my ($in_inertia_txt) = getTotalStringFromFile($in_inertia);
	my ($in_corsair_txt) = getTotalStringFromFile($in_corsair);	
	my ($in_spade_txt) = getTotalStringFromFile($in_spade);	
	
	
	#--------#
	# PROTEO #
	#--------# 
	$logger->info("-- create proteomics report\n");
	# exon_chr   exon_start   exon_end   gene_id;transc_id     pep_chr  pep_start  pep_end  transc_id;pep_id   peplen-numexp    overlaping_pep_within_exon
	# CHR	STAR	END	GENE;TRANS	PEPITDE_LEN-NUM_EXPERIMENTS BP_OVERLAPING
	# 1       877939  878438  ENSG00000187634.6;ENST00000342066.3     1       878321  878371  ENST00000342066;XTAPEP070095346 17-1    50
	# 1       877939  878438  ENSG00000187634.6;ENST00000342066.3     1       877940  877966  ENST00000342066;XTAPEP070095347 9-1     26
	# Y       27184956        27185061        ENSG00000185894.3;ENST00000382287.1
	my ($pro_report);
	foreach my $in_line (@{$in_proteo_txt}) {
		if ( defined $in_line ) {
			my (@fields) = split ('\t', $in_line);
			my ($chr);
			my ($start);
			my ($end);
			my ($gene_transc);
			my ($transc_peps);
			my ($score) = '';
			my ($bp_overlap) = '';
			$chr = $fields[0] if ( defined $fields[0] );
			$start = $fields[1] if ( defined $fields[1] );
			$end = $fields[2] if ( defined $fields[2] );
			$gene_transc = $fields[3] if ( defined $fields[3] );
			$transc_peps = $fields[7] if ( defined $fields[7] );
			
			if ( defined $fields[8] and ($fields[8] ne '') ) {
				$score = $fields[8];
				$score =~ s/\s*//g;
			}
			if ( defined $fields[9] and ($fields[9] ne '') ) {
				$bp_overlap = $fields[9];
				$bp_overlap =~ s/\s*//g;
				$bp_overlap += 1 if ( $bp_overlap ne '' ); # fix the problem with bedtool
			}
			if ( $gene_transc =~ /([^\;]*)\;([^\$]*)$/ ) {
				my ($gene_id) = $1;
				my ($transc_id) = $2;
				my ($pep_id);
				$gene_id =~ s/\.\d*$//; $transc_id =~ s/\.\d*$//; $gene_id =~ s/\n*//g; $transc_id =~ s/\n*//g;
				if ( defined $transc_peps and $transc_peps =~ /([^\;]*)\;([^\$]*)$/ ) {
					$pep_id = $2;
				}
				
				# exons bigger or equal than 42bp. Discard chrM. Discard wrong genes
				if ( (abs($end-$start) < 42) || ($chr eq 'M') ) {
					next;
				}

				# Discard exons than mapped less than a peptide
				if ( ($bp_overlap ne '') and ($bp_overlap < 3) ) {
					next;
				}

				my ($exon_index) = $start.'-'.$end;
				unless ( exists $pro_report->{$exon_index} ) {
					$pro_report->{$exon_index}->{$gene_id} = { 'transc' => $transc_id };
					if ( defined $pep_id ) {						
						$pro_report->{$exon_index}->{$gene_id}->{'peps'} = $pep_id;
					}
					if ( $score ne '' ) {						
						$pro_report->{$exon_index}->{$gene_id}->{'score'} = $score;
					}
					if ( $bp_overlap ne '' ) {
						$pro_report->{$exon_index}->{$gene_id}->{'bp_overlap'} = $bp_overlap;
					}
				}
				else {
					if ( exists $pro_report->{$exon_index}->{$gene_id} ) {
						$pro_report->{$exon_index}->{$gene_id}->{'transc'} .= ";".$transc_id;
						if ( exists $pro_report->{$exon_index}->{$gene_id}->{'peps'} ) {
							$pro_report->{$exon_index}->{$gene_id}->{'peps'} .= ";".$pep_id;
						}
						if ( exists $pro_report->{$exon_index}->{$gene_id}->{'score'} ) {
							$pro_report->{$exon_index}->{$gene_id}->{'score'} .= ";".$score;
						}				
						if ( exists $pro_report->{$exon_index}->{$gene_id}->{'bp_overlap'} ) {
							$pro_report->{$exon_index}->{$gene_id}->{'bp_overlap'} .= ";".$bp_overlap;
						}						
					}
					else {
						#$logger->warning("\tPROTEO: Different genes with common exon...read-through?:\n".$in_line."\n");
						$pro_report->{$exon_index}->{$gene_id} = { 'transc' => $transc_id };
						if ( defined $pep_id ) {						
							$pro_report->{$exon_index}->{$gene_id}->{'peps'} = $pep_id;
						}
						if ( $score ne '' ) {						
							$pro_report->{$exon_index}->{$gene_id}->{'score'} = $score;
						}
						if ( $bp_overlap ne '' ) {
							$pro_report->{$exon_index}->{$gene_id}->{'bp_overlap'} = $bp_overlap;
						}
					}
				}
			}			
		}
	}
	$logger->debug(Dumper($pro_report)."\n");


	#---------#
	# INERTIA #
	#---------# 
	$logger->info("-- create appris/inertia report\n");
	# CHR	EXON_START:EXON_END		GENE_ID		TRANSC_ID;TRANSC_ID2;TRANSC_ID2		(OVERLAP_DIFF_FRAME|NO_PRINC)		MIN_OMEGA_REG_SCORE	
	#1       151584678:151584857     ENSG00000143376 ENST00000368841 NO_PRINC      0.230423333333333
	#1       151611332:151611595     ENSG00000143376 ENST00000368838 NO_PRINC      0.0620738636363636
	#Y       59342487:59342520       ENSGR0000124334 ENSTR0000540897;ENSTR0000369423 OVERLAP_DIFF_FRAME      0.941127272727273
	#Y       59342487:59343077       ENSGR0000124334 ENSTR0000244174;ENSTR0000424344 OVERLAP_DIFF_FRAME      0.979612755102041
	my ($iner_report);
	foreach my $in_line (@{$in_inertia_txt}) {
		if ( $in_line =~ /([^\t]*)\t+([^\t]*)\t+([^\t]*)\t+([^\t]*)\t+([^\t]*)\t+([^\n]*)\n/ ) {
			my ($chr) = $1;
			my ($exon_index) = $2;
			my ($gene_id) = $3;
			my ($transc_list) = $4;
			my ($label_annot) = $5;
			my ($cor_annot) = undef; #$6;
			my ($iner_annot) = undef; #$7;
			my ($omega_score) = $6; #$8;
			$gene_id =~ s/\.\d*$//; $gene_id =~ s/\n*//g; $transc_list =~ s/\n*//g;
			
			# exons bigger or equal than 42bp. Discard chrM. Discard wrong genes
			if ( $WRONG_GENES_FOR_INERTIA =~ /$gene_id/ ) {
				next;
			}

			# We discard exons/genes when at least one their transcripts have less than 4 aligned species (counting human)
			my ($discard) = 0;
			foreach my $transc_id (split(';',$transc_list)) {
				if ( exists $NUM_SPECIES_MAF_ALIGN_GENCODE12->{$transc_id} and 
					defined $NUM_SPECIES_MAF_ALIGN_GENCODE12->{$transc_id} and 
					($NUM_SPECIES_MAF_ALIGN_GENCODE12->{$transc_id} <= 4 ) ) {
						$discard = 1;
						last;						
				}
			}
			next if ( $discard == 1 );
		
			unless ( exists $iner_report->{$exon_index} ) {
				$iner_report->{$exon_index}->{$gene_id} = {
					'transcs'		=> $transc_list,
					'label'			=> $label_annot,
					'cor_annot'		=> $cor_annot,
					'iner_annot'	=> $iner_annot,
					'omega_sc'		=> $omega_score,
				}				
			}
			else {
				#$logger->warning("\tINERTIA: Different genes with common exons...read-through?:\n".$in_line."\n");
				$iner_report->{$exon_index}->{$gene_id} = {
					'transcs'		=> $transc_list,
					'label'			=> $label_annot,
					'cor_annot'		=> $cor_annot,
					'iner_annot'	=> $iner_annot,
					'omega_sc'		=> $omega_score,
				}
			}			
		}
	}
	$logger->debug(Dumper($iner_report)."\n");
	
			
	
	#---------#
	# CORSAIR #
	#---------# 
	$logger->info("-- create corsair report\n");
	my ($cor_report);
#	foreach my $in_line (@{$in_corsair_txt}) {
#		#45940956:45941132       14
#		#45869728:45869877       0.5 {14-ENST00000374391.2}
#		if ( defined $in_line ) {
#			my (@fields) = split ('\s', $in_line);
#			my ($exon_index);
#			my ($score) = '';
#			my ($score2) = '';
#			$exon_index = $fields[0] if ( defined $fields[0] );
#			$score = $fields[1] if ( defined $fields[1] );
#			$score2 = $fields[2] if ( defined $fields[2] );
#			$exon_index =~ s/\s*//g;
#			unless ( exists $cor_report->{$exon_index} ) {
#				if ( $score ne '' ) {						
#					$cor_report->{$exon_index}->{'score'} = $score;
#				}
#				if ( ($score2 ne '') and ($score2=~/\{([^\-]*)/) ) {
#					my ($sc) = $1;
#					$cor_report->{$exon_index}->{'max_sc'} = $sc;
#				}
#			}
#			else {
#				if ( exists $cor_report->{$exon_index}->{'score'} ) {
#					$cor_report->{$exon_index}->{'score'} .= ";".$score;
#				}				
#				if ( exists $cor_report->{$exon_index}->{'max_sc'} ) {
#					if ( $cor_report->{$exon_index}->{'max_sc'}=~/\{([^\-]*)/ ) {
#						my ($sc) = $1;
#						$cor_report->{$exon_index}->{'max_sc'} .= ";".$sc;
#					}
#				}
#			}
#		}
#	}
#	$logger->debug(Dumper($cor_report)."\n");	


	#-------#
	# SPADE #
	#-------# 
	$logger->info("-- create spade report\n");
	# HEAD: dom_chr dom_start dom_end gene_id;transc_id     pep_chr  pep_start  pep_end  transc_id;pep_id   peplen-numexp    overlaping_pep_within_exon
	#1       151584810       151611448       ENSG00000143376;ENST00000368843 1       151584822       151584851       ENST00000368843;XTAPEP070127136 10-2    29
	#1       151584810       151611448       ENSG00000143376;ENST00000368843 1       151584882       151584953       ENST00000368843;XTAPEP070193024 24-5    71
	#1       151584810       151611448       ENSG00000143376;ENST00000368843 1       151584969       151611397       ENST00000368843;XTAPEP070193022 18-1    26428 # DOMAIN WITHIN SEVERAL EXONS
	
	
	
	
	
	
	
	#-------#	

	$logger->info("-- sort all cds of gencode12\n");
	my (@all_cds) = split('\n', $CDS_GENCODE12);
	my (@sorted_coord_exons) = sort { # sort and uniq exon coords
					my ($a2,$b2);
					$a =~ /^([^\:]*)\:/; $a2=$1;
					$b=~ /^([^\:]*)\:/; $b2=$1;
					$a2 <=> $b2
				} @all_cds;
	my (@pro_exons) = keys(%{$pro_report});
	my (@iner_exons) = keys(%{$iner_report});
	my (@cor_exons) = keys(%{$cor_report});				
	#$logger->debug(Dumper(@sorted_coord_exons)."\n");
	
	
	
	$logger->info("-- create joined report\n");
	my ($joined_report);
	foreach my $exon_index (@sorted_coord_exons) {
		my (@exon_coord) = split( '-', $exon_index );
		my ($e_start) = $exon_coord[0];
		my ($e_end) = $exon_coord[1];
		my ($exon_len) = abs($e_end-$e_start);
		my ($peplen_numexp) = '-';
		my ($p_g_list, $p_tr_list, $p_pe_list) = ('-','-','-');
		
		# Create report with a subset of exons-genes-omega_sc
		if ( exists $iner_report->{$exon_index} and defined $iner_report->{$exon_index} ) {
			
			while (my ($gene_id, $rep) = each(%{$iner_report->{$exon_index}}) )
			{
				# We discard genes that have been labeled as readthrough
				if ( $READTHROUGH_GENES_GENCODE12 =~ /$gene_id/ ) {
					next;
				}				
				# We discard genes that have problems with the model
				if ( $GENE_MODEL_PROBLEMS_GENCODE12 =~ /$gene_id/ ) {
					next;
				}				
				# We discard genes that have been labeled as alternative splincing
				if ( $BONAFIDE_AS_GENES_GENCODE12 =~ /$gene_id/ ) {
					next;
				}

				# We discard gene that has less than 2 founded peptides
				if ( exists $NUM_PEPS_GENES_GENCODE12->{$gene_id} and 
					defined $NUM_PEPS_GENES_GENCODE12->{$gene_id} and 
					($NUM_PEPS_GENES_GENCODE12->{$gene_id} < 2) ) {
						next;
				}				

				# We use alternative transcripts whose exons dont overlap
				my ($alter_label) = '-';
				if ( exists $rep->{'label'} ) {
					$alter_label = $rep->{'label'};
					# DEPRECATED 
					#if ( $rep->{'label'} eq 'PRINCIPAL-OVERLAP' ) {
					#	$alter_label = $rep->{'label'};
					#} else {
					#	next;
					#}
				}
				my ($trid_list) = '-';
				if ( exists $rep->{'transcs'} ) {
					$trid_list = $rep->{'transcs'};
				}
				my ($inertia_sc) = '-';
				if ( exists $rep->{'omega_sc'} ) {
					$inertia_sc = $rep->{'omega_sc'};
				}				
				$joined_report->{$exon_index}->{$gene_id} = {
						'alter_label'		=> $alter_label,
						'alter_nooverlap'	=> $trid_list,
						'inertia_sc'		=> $inertia_sc,
						'alter_pepevents'	=> '-',
						'proteo_sc'			=> '-',
						'pep_bp_over'		=> '-',
						'pepid_list'		=> '-',
				}
			}
		}
		
		# Add into report the subset of exons-peptides
		if ( exists $pro_report->{$exon_index} and defined $pro_report->{$exon_index} ) {
						
			while (my ($gene_id, $rep) = each(%{$pro_report->{$exon_index}}) )
			{
				my ($trid_list) = '-';
				if ( exists $rep->{'transc'} ) {
					$trid_list = $rep->{'transc'};
				}
				my ($proteo_sc) = '-';
				if ( exists $rep->{'score'} ) {
					$proteo_sc = $rep->{'score'};
				}
				my ($pep_bp_over) = '-';
				if ( exists $rep->{'bp_overlap'} ) {
					$pep_bp_over = $rep->{'bp_overlap'};
				}
				my ($pepid_list) = '-';
				if ( exists $rep->{'peps'} ) {
					$pepid_list = $rep->{'peps'};
				}
				
				$joined_report->{$exon_index}->{$gene_id}->{'alter_pepevents'} = $trid_list;
				$joined_report->{$exon_index}->{$gene_id}->{'proteo_sc'} = $proteo_sc;
				$joined_report->{$exon_index}->{$gene_id}->{'pep_bp_over'} = $pep_bp_over;
				$joined_report->{$exon_index}->{$gene_id}->{'pepid_list'} = $pepid_list;
				
			}
		}
	}
	$logger->debug(Dumper($joined_report)."\n");
	
		
	$logger->info("-- print joined report\n");
	my ($output) = '';
	foreach my $exon_index (@sorted_coord_exons) {
		if ( exists $joined_report->{$exon_index}) {
			while (my ($gene_id,$g_report) = each(%{$joined_report->{$exon_index}}) ) {
				if ( exists $g_report->{'alter_label'} and exists $g_report->{'alter_nooverlap'} and exists $g_report->{'inertia_sc'} ) {
					$output .= 	$exon_index."\t".
								$gene_id."\t".
								$g_report->{'alter_label'}."\t".
								$g_report->{'alter_nooverlap'}."\t".
								$g_report->{'inertia_sc'}."\t".
								$g_report->{'alter_pepevents'}."\t".
								$g_report->{'pepid_list'}."\t".
								$g_report->{'proteo_sc'}."\t".
								$g_report->{'pep_bp_over'}."\n";
								
				}
			}			
		}
	}
	if ($output ne '') {
		my ($printing_file_log) = printStringIntoFile($output,$outfile);
		$logger->error("printing") unless(defined $printing_file_log);		
	}

	$logger->finish_log();
	
	exit 0;
	
}

main();

1;


__END__

=head1 NAME

cmp_peptides_iner_cor_spa.pl

=head1 DESCRIPTION

=head2 Required arguments:

	--infile= <Input file>
	
	--outfile= <Output file>

=head2 Optional arguments (logs):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl convert_output_scoreNR.pl
			
		--in-proteo=/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/proteo.exons_VS_peplen-numexp.rel12.bed
		
		--in-inertia=/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/appris_data.princ_noover_difframe_exons.gen12.v5.txt
		
		--in-corsair=/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/corsair_exon.rel12.txt
		
		--in-spade=/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/proteo.domains_VS_peplen-numexp.rel12.bed
		
		--outfile=/Users/jmrodriguez/projects/APPRIS/workspaces/gencode12/exon_cmp_inertia_corsair_proteo/cmp_exons.princ_VS_peplen-numexp_VS_inertia.tsv
		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut