#! /usr/bin/perl -W

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;

use APPRIS::Utils::File qw( prepare_workspace getStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	
	$USED_ORGANISMS
	$ALL_ORGANISMS
	
	
	$UCSC_46WAY_TREE
	$ALIGN_IN_SUFFIX
	$ALIGN_SUFFIX
	$TREE_SUFFIX
);

$LOCAL_PWD		= $FindBin::Bin;

$USED_ORGANISMS	= {
	'hg19'		=> 'homo_sapiens',
	'panTro2'	=> 'pan_troglodytes',
	'ponAbe2'	=> 'pongo_abelii',
	'rheMac2'	=> 'macaca_mulatta',
	'mm9'		=> 'mus_musculus',
	'rn4'		=> 'rattus_norvegicus',
	'cavPor3'	=> 'cavia_porcellus',
	'oryCun2'	=> 'oryctolagus_cuniculus',
	'bosTau4'	=> 'bos_taurus',
	'equCab2'	=> 'equus_caballus',
	'canFam2'	=> 'canis_familiaris',
	'loxAfr3'	=> 'loxodonta_africana',
	'monDom5'	=> 'monodelphis_domestica', # differ from ECompara
	'ornAna1'	=> 'anolis_carolinensis', # differ from ECompara
	'galGal3'	=> 'gallus_gallus', # differ from ECompara
	'taeGut1'	=> 'taeniopygia_guttata',  # differ from ECompara
	'xenTro2'	=> 'xenopus_tropicalis', # differ from ECompara
	'danRer6'	=> 'danio_rerio', # differ from ECompara
};

$ALL_ORGANISMS	= {
	'hg19'		=> 'homo_sapiens',
	'panTro2'	=> 'pan_troglodytes',
	'gorGor1'	=> 'gorilla_gorilla',
	'ponAbe2'	=> 'pongo_abelii',
	'rheMac2'	=> 'macaca_mulatta',
	'papHam1'	=> 'papio_hamadryas',
	'calJac1'	=> 'callithrix_jacchus',
	'tarSyr1'	=> 'tarsius_syrichta',
	'micMur1'	=> 'microcebus_murinus',
	'otoGar1'	=> 'otolemur_garnettii',
	'tupBel1'	=> 'tupaia_belangeri',
	'mm9'		=> 'mus_musculus',
	'rn4'		=> 'rattus_norvegicus',
	'dipOrd1'	=> 'dipodomys_ordii',
	'cavPor3'	=> 'cavia_porcellus',
	'speTri1'	=> 'spermophilus_tridecemlineatus',	
	'oryCun2'	=> 'oryctolagus_cuniculus',
	'ochPri2'	=> 'ochotona_princeps',
	'vicPac1'	=> 'vicugna_pacos',
	'turTru1'	=> 'tursiops_truncatus',
	'bosTau4'	=> 'bos_taurus',
	'equCab2'	=> 'equus_caballus',
	'felCat3'	=> 'felis_catus',
	'canFam2'	=> 'canis_familiaris',
	'myoLuc1'	=> 'myotis_lucifugus',
	'pteVam1'	=> 'pteropus_vampyrus',
	'eriEur1'	=> 'erinaceus_europaeus',
	'sorAra1'	=> 'sorex_araneus',
	'loxAfr3'	=> 'loxodonta_africana',
	'proCap1'	=> 'procavia_capensis',
	'echTel1'	=> 'echinops_telfairi',
	'dasNov2'	=> 'dasypus_novemcinctus',
	'choHof1'	=> 'choloepus_hoffmanni',
	'macEug1'	=> 'macropus_eugenii',
	'monDom5'	=> 'monodelphis_domestica',
	'ornAna1'	=> 'ornithorhynchus_anatinus',
	'galGal3'	=> 'gallus_gallus',
	'taeGut1'	=> 'taeniopygia_guttata',
	'anoCar1'	=> 'anolis_carolinensis',
	'xenTro2'	=> 'xenopus_tropicalis',
	'tetNig2'	=> 'tetraodon_nigroviridis',
	'fr2'		=> 'takifugu_rubripes',
	'gasAcu1'	=> 'gasterosteus_aculeatus',
	'oryLat2'	=> 'oryzias_latipes',
	'danRer6'	=> 'danio_rerio',
	'petMar1'	=> 'petromyzon_marinus',	
};
  

$UCSC_46WAY_TREE		= $LOCAL_PWD.'/46way.nh';
$ALIGN_IN_SUFFIX		= '.ucsc.init';
$ALIGN_SUFFIX			= '.ucsc.faa';
$TREE_SUFFIX			= '.ucsc.nh';


# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($species) = undef;
my ($gff_file) = undef;
my ($transl_file) = undef;
my ($outpath) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'species=s'			=> \$species,
	'maf=s'				=> \$maf_file,
	'gff=s'				=> \$gff_file,
	'transl_file=s'		=> \$transl_file,	
	'outpath=s'		=> \$outpath,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined $target and 
	defined $gff_file and
	defined $transl_file and
	defined $maf_file and defined $outpath
)
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
$logger->init_log($str_params);

##############
# Prototypes #
##############
sub get_cds_list_for_transcript($$$$);
sub get_cds_sequences($$);
sub extract_cds_sequence($$$);
sub reverse_dna_complement($);
sub create_cds_alignments($$$$$);
sub create_transcript_alignments($$$);
sub create_filter_transcript_alignments($$$);
sub get_cds_alignment_from_file($);
sub filter_cds_alignment($);
sub print_cds_filter_aligment($$);


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	# Prepare workspace
	$logger->info("-- prepare workspace\n");
	$outpath = prepare_workspace($outpath);
	$logger->error("ERROR: creating workspace ") unless (defined $outpath);
	
	# gff parsing
	$logger->info("-- get cds coords for transc\n");
	my ($transcript_cds_list, $cds_list, $strand_cds_list) = get_cds_list_for_transcript($chr, $gff_file, $transl_file, $outpath);
	$logger->debug("TRANS:\n".Dumper($transcript_cds_list)."\n");
	$logger->debug("CDS_LIST:\n".Dumper($cds_list)."\n");
	$logger->debug("STRAND_CDS:\n".Dumper($strand_cds_list)."\n");

	# Maf parsing
	$logger->info("-- get cds sequences for transc\n");
	my ($cds_2_sequences, $cds_2_organims) = get_cds_sequences($maf_file, $cds_list);
	$logger->debug("CDS_2_SEQ:\n".Dumper($cds_2_sequences)."\n");
	$logger->debug("CDS_2_ORG:\n".Dumper($cds_2_organims)."\n");

	# Reverse complement	
	$logger->info("-- get cds rev sequences for transc\n");
	get_cds_reverse_sequences($cds_2_sequences,$strand_cds_list);
	$logger->debug("CDS_2_SEQ_REV:\n".Dumper($cds_2_sequences)."\n");
		
	# Create transcript alignments from every cds
	$logger->info("-- create cds aligns for transc\n");
	my ($transcrip_multi_alignment) = create_cds_alignments($chr, $transcript_cds_list, $cds_2_sequences, $cds_2_organims, $outpath);
		
	$logger->finish_log();
	
	exit 0;
}

#################
# .gff2 parsing #
#################

# Get sequence
sub _parse_sequence($)
{
	my ($file) = @_;
	my ($data);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() )
		{
			if ( $seq->id=~/([^|]*)/ )
			{
				my ($sequence_id) = $1;
				if(exists $data->{$sequence_id}) {
					throw("Duplicated sequence: $sequence_id");
				}
				else {
					$data->{$sequence_id} = $seq->seq;				
				}
			}
		}		
	}
	return $data;
}

# Get id and version from attributes list
sub _get_ids_version($)
{
	my ($attributes) = @_;
	my ($t_id, $g_id) = (undef,undef);
	
	my (@add_attributes) = split(";", $attributes);				
	for ( my $i=0; $i<scalar @add_attributes; $i++ )
	{
		$add_attributes[$i] =~ /^(.+)\s(.+)$/;
		my ($c_type) = $1;
		my ($c_value) = $2;
		if(	defined $c_type and !($c_type=~/^\s*$/) and
			defined $c_value and !($c_value=~/^\s*$/))
		{
			$c_type =~ s/^\s//;
			$c_value =~ s/"//g;
			if ($c_type eq 'transcript_id')
			{
				if ( $c_value =~ /^([^\.]+\.\d+)$/ ) {
					$t_id = $1;
				}				
			}
			if ($c_type eq 'gene_id')
			{
				if ( $c_value =~ /^([^\.]+\.\d+)$/ ) {
					$g_id = $1;
				}				
			}			
		}
	}
	return ($t_id, $g_id);
}

## De las lineas que contienen informacion sobre los exones sacamos coordenadas, strand y nombre del transcrito.
## Al final obtenemos un hash de arrays (%transcripts) que por cada transcrito nos dice por cuales exones 
## (el "nombre" del exon son sus coordenadas con un _ en el medio) esta constituido, 
## un hash %strand que guarda con clave {nombre del exon} el sentido (+/-), y otro hash %start_coord que con 
## clave {nombre del exon} guarda la coordenada de inicio
sub get_cds_list_for_transcript($$$$)
{
	my ($chr, $gff_file, $transl_file, $work_dir) = @_;
	my ($transcripts);
	my (@cds_list);
	my ($strand_cds_list);
	my ($aux_cds_list);
	my ($aux_trans_strand);
	
	# Get protein sequences
	my ($pep_sequences) = _parse_sequence($transl_file);

	# Get the cds info
	my($exons);
	my(@chr_global_info)=`awk '{if ( \$1=="chr$chr" ) {print \$0} }' $gff_file`;
	foreach my $chr_info_line (@chr_global_info)
	{
		if ($chr_info_line =~ /^[^\t]*\t+[^\t]*\t+transcript\t+[^\t]*\t+[^\t]*\t+[^\t]*\t+([^\t]*)\t+[^\t]*\t+([^\n]*)/)
		{
			#store ids and additional information in second hash
			my ($trans_strand) = $1;			
			my ($attributes) = $2;
			my ($transc_id, $gene_id) = _get_ids_version($attributes);
			if (exists $pep_sequences->{$transc_id} and defined $pep_sequences->{$transc_id}) { # check if protein sequence already exists 
				$transcripts->{$transc_id}->{'strand'} = $trans_strand;
				$transcripts->{$transc_id}->{'gene_id'} = $gene_id;
			}
		}
		elsif ($chr_info_line =~ /^[^\t]*\t+[^\t]*\t+CDS\t+([^\t]*)\t+([^\t]*)\t+[^\t]*\t+([^\t]*)\t+[^\t]*\t+([^\n]*)/)
		{
			my ($cds_start) = $1;
			my ($cds_end) = $2;
			my ($trans_strand) = $3;
			my ($attributes) = $4;
			my ($transc_id, $gene_id) = _get_ids_version($attributes);
			my ($cds_range) = $cds_start."_".$cds_end;
			if (exists $pep_sequences->{$transc_id} and defined $pep_sequences->{$transc_id}) { # check if protein sequence exists
				if (exists $strand_cds_list->{$cds_range}) {
					push(@{$transcripts->{$transc_id}->{'cds_list'}},$cds_range);
				}
				else {
					push(@{$transcripts->{$transc_id}->{'cds_list'}}, $cds_range);
					push(@{$aux_cds_list}, $cds_range);
					$strand_cds_list->{$cds_range} = $trans_strand;
				}
			}
		}
	}
	
	# Numerical sort CDS list (using the start coordinate -the first value-)
	if (defined $aux_cds_list) {
		@cds_list = sort { ($a=~/([^\_]*)/)[0] <=> ($b=~/([^\_]*)/)[0] } @{$aux_cds_list};		
	}

	return ($transcripts,\@cds_list,$strand_cds_list);
}

################
# .maf parsing #
################

## empezamos a parsear el fichero linea por linea porque por su dimensión sería imposible cargarlo en un array.
## La idea es utilizando variables de flag 
sub get_cds_sequences($$)
{
	my ($file, $cds_list) = @_;
	
	# Init variables
	my ($mafline_2_sequences);
	my ($cds_2_sequences);
	my ($cds_2_organims);
	my ($control_cds_list);
	foreach my $cds_coord (@{$cds_list}) {
		$control_cds_list->{$cds_coord} = undef;
		$cds_2_sequences->{$cds_coord} = undef;
	}

	open (FILE, $file);
	while (<FILE>) {

		if ($_ =~ /^s\s(\w+)\.\w+\s+(\d+)\s+(\d+)\s+.\s+\d+\s+(.+)/)
		{
			$logger->debug("ENTRA_1 ---------------\n");
	
			# Patron de las lineas con las secuencias en el formato Maf. 
			# Para mas informacion en http://genome.ucsc.edu/FAQ/FAQformat#format5
			my ($organism) = $1;
			my ($maf_start) = $2;
			my ($maf_length) = $3;
			my ($maf_sequence) = $4;
			my ($maf_gap) = 0;

			$mafline_2_sequences->{$organism}->{'sequence'} .= $maf_sequence;
			$mafline_2_sequences->{$organism}->{'length'} = $maf_length;
			$mafline_2_sequences->{$organism}->{'start'} = $maf_start + 1;
			$mafline_2_sequences->{$organism}->{'end'} = $maf_start + $maf_length;
		}
		elsif ( $_ =~ /^\n$/)
		{
			$logger->debug("ENTRA_2 ---------------\n");

			# Get coordinates of target alignment
			my ($target_maf_sequence) = $mafline_2_sequences->{$target}->{'sequence'};
			my ($target_maf_length) = $mafline_2_sequences->{$target}->{'length'};
			my ($target_maf_start_coord) = $mafline_2_sequences->{$target}->{'start'};
			my ($target_maf_end_coord) = $mafline_2_sequences->{$target}->{'end'};			
			
			# Scan every CDS taking into account the several cases that one CDS fall down
			# into the Maf Alignment
			foreach my $cds_coord (@{$cds_list}) {
				my (@split_coord) = split('_',$cds_coord);
				my ($cds_start_coord) = $split_coord[0];
				my ($cds_end_coord) = $split_coord[1];
				my ($cds_organisms);
				my ($init_sequence);
				my ($length_sequence);
				
				# Open CDS alignment 
				if ( ($target_maf_start_coord <= $cds_start_coord) and ($target_maf_end_coord >= $cds_start_coord) and
					 ($target_maf_start_coord <= $cds_end_coord) and ($target_maf_end_coord <= $cds_end_coord))
				{
					$logger->debug("ENTRA_2_1 ---------------\n");

					# Get the CDS coord (forgetting the gaps)
					$control_cds_list->{$cds_coord} = 'open';
					$init_sequence = $cds_start_coord - $target_maf_start_coord;
					$length_sequence = ($target_maf_end_coord - $cds_start_coord) + 1;
				}
				# Continue CDS alignment
				elsif ( ($target_maf_start_coord >= $cds_start_coord) and ($target_maf_end_coord >= $cds_start_coord) and
						($target_maf_start_coord <= $cds_end_coord) and ($target_maf_end_coord <= $cds_end_coord) )
				{
					$logger->debug("ENTRA_2_2 ---------------\n");

					# Get the CDS coord (forgetting the gaps)
					$control_cds_list->{$cds_coord} = 'open';
					$init_sequence = 0;
					$length_sequence = ($target_maf_end_coord - $target_maf_start_coord) + 1;
				}
				# Open and Close CDS alignment
				elsif ( ($target_maf_start_coord <= $cds_start_coord) and ($target_maf_end_coord >= $cds_end_coord) )
				{
					$logger->debug("ENTRA_2_3 ---------------\n");
					
					# Get the CDS coord (forgetting the gaps)
					$control_cds_list->{$cds_coord} = 'close';					
					$init_sequence = $cds_start_coord - $target_maf_start_coord;
					$length_sequence = ($cds_end_coord - $cds_start_coord) + 1;
				}
				# Close CDS alignment
				elsif ( ($target_maf_start_coord >= $cds_start_coord) and ($target_maf_end_coord >= $cds_start_coord) and
					 ($target_maf_start_coord <= $cds_end_coord) and ($target_maf_end_coord >= $cds_end_coord) )
				{
					$logger->debug("ENTRA_2_4 ---------------\n");

					# Get the CDS coord (forgetting the gaps)
					$control_cds_list->{$cds_coord} = 'close';					
					$init_sequence = 0;
					$length_sequence = ($cds_end_coord - $target_maf_start_coord) + 1;
				}
				
				# Extract the alignments of every organism from CDS coord (forgetting the gaps)
				if (defined $control_cds_list->{$cds_coord} and defined $init_sequence and defined $length_sequence) {
					my ($target_cds_sequence,$target_init_cds_sequence,$target_length_cds_sequence) = extract_target_cds_sequence($target_maf_sequence, $init_sequence, $length_sequence); # ... first for target (human)
					if (defined $target_cds_sequence and ($target_cds_sequence ne '')) {
						$cds_2_sequences->{$cds_coord}->{$target} .= $target_cds_sequence;
						unless (exists $cds_2_organims->{$cds_coord}->{$target}) {
							$cds_2_organims->{$cds_coord}->{$target} = undef;
						}
						
						$logger->debug("Maf_start: $target_maf_start_coord Maf_end: $target_maf_end_coord\n");
						$logger->debug("Cds_start: $cds_start_coord Cds_end: $cds_end_coord\n");
						$logger->debug("Init_seq:  $init_sequence\n");
						$logger->debug("Length_seq: ($target_maf_end_coord - $cds_start_coord) + 1 = $length_sequence\n");
						$logger->debug("Seq:$target: $target_cds_sequence\n");
					}
					foreach my $organism ( keys(%{$USED_ORGANISMS}) ) { # ... then for the rest of organism
						if ($organism ne $target) {
							my ($cds_sequence);
							if (exists $mafline_2_sequences->{$organism}) {								
								my ($maf_sequence) = $mafline_2_sequences->{$organism}->{'sequence'};
								$cds_sequence = extract_cds_sequence($maf_sequence,$target_init_cds_sequence,$target_length_cds_sequence);
							}
							else {
								my ($num_residues) = length($target_cds_sequence);
								$cds_sequence = sprintf('-' x $num_residues);
							}
							if (defined $cds_sequence and ($cds_sequence ne '')) {
								$cds_2_sequences->{$cds_coord}->{$organism} .= $cds_sequence;
								unless (exists $cds_2_organims->{$cds_coord}->{$organism}) {
									$cds_2_organims->{$cds_coord}->{$organism} = undef;
								}
							}
							
							$logger->debug("Target_Init_seq:  $target_init_cds_sequence  \n");
							$logger->debug("Target_Length_seq:  $target_length_cds_sequence\n");
							$logger->debug("Seq:$organism: $cds_sequence\n");
						}
					}
				}				
			}
			# Init MAF line 
			$mafline_2_sequences = undef; 
		}
	}
	close(FILE);

	return ($cds_2_sequences, $cds_2_organims);
}
sub extract_target_cds_sequence($$$)
{
	my ($seq, $init_seq, $length_seq) = @_;
	
	my ($num_included_residues) = 0;
	my ($num_excluded_residues) = 0;
	my ($num_total_init_residues) = 0;
	my ($target_seq) = '';
	my ($target_init_seq) = $init_seq;
	my ($target_length_seq) = 0;
	my (@seq_array) = split(//,$seq);
	for (my $index = 0; $index < scalar(@seq_array); $index++) {
		my ($residue) = $seq_array[$index];
	
		if ($num_excluded_residues >= $init_seq) {
			if ($num_excluded_residues == $init_seq) {
				$target_init_seq = $num_total_init_residues;
			}
			$target_seq .= $residue;
			if ($residue ne '-') {
				$num_included_residues++;
			}
			if ($num_included_residues >= $length_seq) {
				last;
			}
		}
		elsif ($num_excluded_residues < $init_seq) {
			$num_total_init_residues++;
			if ($residue ne '-') {
				$num_excluded_residues++;
			}
		}
	}
	return ($target_seq, $target_init_seq, length($target_seq));
}
sub extract_cds_sequence($$$)
{
	my ($seq, $init_seq, $length_seq) = @_;
	
	my ($num_residues) = 0;
	my ($target_sequence) = '';
	my (@seq_array) = split(//,$seq);
	my ($index) = $init_seq;
	while ( ($index < scalar(@seq_array)) and ($num_residues < $length_seq) ) {
		my ($residue) = $seq_array[$index];
		$target_sequence .= $residue;
		$index++;
		$num_residues++;
	}
	return $target_sequence;
}
sub get_cds_reverse_sequences($$)
{
	my ($cds_2_sequences,$strand_cds_list) = @_;
	
	while (my($cds_coord,$cds_strand)=each(%{$strand_cds_list})) {
		if($cds_strand eq '-') {
			while (my($maf_organism,$maf_sequence)=each(%{$cds_2_sequences->{$cds_coord}})) {
				$cds_2_sequences->{$cds_coord}->{$maf_organism} = reverse_dna_complement($maf_sequence);
			}
		}		
	}
}
sub reverse_dna_complement($) {
  my ($dna) = shift;
  my ($revcomp) = scalar reverse($dna);

  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}


#############################################
# Create aligment file for every transcript #
#############################################

sub create_cds_alignments($$$$$)
{
	my ($chr, $transcript_cds_list, $cds_2_sequences, $cds_2_organims, $work_dir) = @_;	
	my ($all_trans_alignment);
	my ($global_output_content) = '';
	my ($global_aligment_file) = $work_dir.'/'."chr$chr.trans.cdsNuc".'.txt';

	while ( my ($transcript_id, $trans_report) = each(%{$transcript_cds_list}))
	{
		if (exists $trans_report->{'cds_list'} and defined $trans_report->{'cds_list'})
		{
			my ($trans_alignment);
			my ($gene_id) = $trans_report->{'gene_id'};
			
			# Numerical sort CDS list (using the start coordinate -the first value-)
			my (@cds_list);
			if ($trans_report->{'strand'} eq '+') {
				@cds_list = sort { ($a=~/([^\_]*)/)[0] <=> ($b=~/([^\_]*)/)[0] } @{$trans_report->{'cds_list'}};
			}
			else {
				@cds_list = sort { ($b=~/([^\_]*)/)[0] <=> ($a=~/([^\_]*)/)[0] } @{$trans_report->{'cds_list'}};
			}
			# Get the list of the whole organisms (except human) for each CDS
			my ($trans_org_list);
			foreach my $cds_coord (@cds_list) {
				my ($cds_org_list) = $cds_2_organims->{$cds_coord};
				foreach my $cds_org (keys(%{$cds_org_list})) {
					unless (exists $trans_org_list->{$cds_org} or ($cds_org eq $target)) {
						$trans_org_list->{$cds_org} = undef;				
					}				
				}
			}
			
			# Create alignment from human alignment
			$global_output_content .= '>'.$transcript_id."\t"."cds:".scalar(@cds_list)."\n";
			$all_trans_alignment->{$transcript_id}->{$target} = '';
			$trans_alignment->{$transcript_id}->{$target} = '';

			foreach my $cds_coord (@cds_list) {
				$global_output_content .= '#'.$cds_coord."\n";
						
				# for human
				my ($human_cds_alignment) = $cds_2_sequences->{$cds_coord}->{$target};
				$global_output_content .= $target."\t".$human_cds_alignment."\n"; # for print
				$all_trans_alignment->{$transcript_id}->{$target} .= $human_cds_alignment; # keep cds align
				$trans_alignment->{$transcript_id}->{$target} .= $human_cds_alignment; # keep cds align
				
				# for the rest of organism
				foreach my $org (keys(%{$trans_org_list})) {
					my ($org_cds_alignment) = '';
					if (exists $cds_2_sequences->{$cds_coord}->{$org}) {
						$org_cds_alignment .= $cds_2_sequences->{$cds_coord}->{$org};
					}
					else { # fill with gaps
						my ($num_residues) = length($human_cds_alignment); 
						$org_cds_alignment .= sprintf('-' x $num_residues);
					}
					$global_output_content .= $org."\t".$org_cds_alignment."\n"; # for print
					$all_trans_alignment->{$transcript_id}->{$org} .= $org_cds_alignment; # keep cds align
					$trans_alignment->{$transcript_id}->{$org} .= $org_cds_alignment; # keep cds align
				}
			}
			# Print transcript CDS alignment
			create_transcript_alignments($gene_id, $trans_alignment, $work_dir);

			# Print transcript with filtered CDS alignment
			create_filter_transcript_alignments($gene_id, $trans_alignment, $work_dir);

			# Print transcript Phylogenetic tree
			create_transcript_tree($gene_id, $transcript_id, $work_dir);
		}
	}
	if ($global_output_content ne '') {
		my ($printing_file_log) = printStringIntoFile($global_output_content, $global_aligment_file);
		$logger->error("Printing output") unless ( defined $printing_file_log );		
	}
	
	return $all_trans_alignment;
}

sub create_transcript_tree($$$)
{
	my ($gene_id, $transcript_id, $work_dir) = @_;
	
	# get the 46Way phylogenetic tree of UCSC
	my ($ucsc_46way_cont) = getStringFromFile($UCSC_46WAY_TREE);

	# change the alias name for each specie
	while ( my ($id,$name) = each(%{$ALL_ORGANISMS}) ) {
		$ucsc_46way_cont =~ s/$id/$name/mg;
	}
	my ($g_work_dir) = $work_dir.'/'.$gene_id.'/data';
	$g_work_dir = prepare_workspace($g_work_dir);
	$logger->error("ERROR: creating workspace ") unless (defined $g_work_dir);
	
	my ($trans_output_file) = $g_work_dir.'/'.$transcript_id.$TREE_SUFFIX;	
	my ($printing_file_log) = printStringIntoFile($ucsc_46way_cont, $trans_output_file);
	$logger->error("Printing output") unless ( defined $printing_file_log );
}

sub create_transcript_alignments($$$)
{
	my ($gene_id, $trans_alignment, $work_dir) = @_;

	# for each transcript
	while ( my ($transcript_id, $trans_report) = each(%{$trans_alignment}))
	{		
		if ( exists $trans_report->{$target} and defined $trans_report->{$target} and $trans_report->{$target} ne '' )
		{
			my ($trans_target_align) = '';
			my ($trans_org_align) = '';
			
			# Prepare workspace
			my ($g_work_dir) = $work_dir.'/'.$gene_id.'/data';
			$g_work_dir = prepare_workspace($g_work_dir);
			$logger->error("ERROR: creating workspace ") unless (defined $g_work_dir);			
			
			# print the multiple alignment
			while ( my ($org, $org_aligment) = each(%{$trans_report}))
			{
				if ($org eq $target) {
					$trans_target_align .= '>'.$org."\n".$org_aligment."\n";
				}
				else {
					$trans_org_align .= '>'.$org."\n".$org_aligment."\n";
				}				
			}
			my ($trans_output_content) = $trans_target_align.$trans_org_align;
			my ($trans_output_file) = $g_work_dir.'/'.$transcript_id.$ALIGN_IN_SUFFIX;

			my ($printing_file_log) = printStringIntoFile($trans_output_content, $trans_output_file);
			$logger->error("Printing output") unless ( defined $printing_file_log );		
		}		
	}
}

sub create_filter_transcript_alignments($$$)
{
	my ($gene_id, $trans_alignment, $work_dir) = @_;

	# for each transcript
	while ( my ($transcript_id, $trans_report) = each(%{$trans_alignment}))
	{		
		if ( exists $trans_report->{$target} and defined $trans_report->{$target} and $trans_report->{$target} ne '' )
		{
			my ($g_work_dir) = $work_dir.'/'.$gene_id.'/data';
			my ($alignmet_file) = $g_work_dir.'/'.$transcript_id.$ALIGN_IN_SUFFIX;
			my ($cds_alignment_report) = get_cds_alignment_from_file($alignmet_file);			
			$logger->error("Getting the report of cds alignment\n") unless (defined $cds_alignment_report);
		
			# Get CDS alignment that its number of gaps is less than the half of nucleotides
			# And delete both "anc_aa" and "ancestor" alignment
			my ($filtered_alignment_report) = filter_cds_alignment($cds_alignment_report);
			$logger->error("Getting the report of filtered cds alignment\n") unless (defined $filtered_alignment_report);
		
			# Print CDS aligments sorting by exons as FASTA format
			my ($filtered_alignment_file) = $g_work_dir.'/'.$transcript_id.$ALIGN_SUFFIX;
			my ($filtered_alignment_fasta) = print_cds_filter_aligment($filtered_alignment_report,'fasta');
			$logger->error("Creating filtered cds alignment\n") unless (defined $filtered_alignment_fasta);
		
			my ($printing_file_log) = printStringIntoFile($filtered_alignment_fasta, $filtered_alignment_file);
			$logger->error("Printing output") unless ( defined $printing_file_log );
		}
	}
}

# Get CDS alignment from file
sub get_cds_alignment_from_file($)
{
	my ($alignment_file) = @_;
	my ($cds_alignment);

	$/=undef;
	local(*ALIGNMENT_FILE);
	open(ALIGNMENT_FILE,$alignment_file) or $logger->error("Can not open alignment directory: $!\n");
	my ($alignmet_content) = <ALIGNMENT_FILE>;
	close(ALIGNMENT_FILE);
	$/='\n';
	
	if ( defined $alignmet_content )
	{
		while ( $alignmet_content =~ />([^\n]+)\n+([a-zA-Z\-]+)/g )
		{
			my ($specie_id) = $1;
			my ($cds_sequence) = $2;
			my ($specie_name) = $specie_id;
			if ( exists $USED_ORGANISMS->{$specie_id} ) {
				$specie_name = $USED_ORGANISMS->{$specie_id}; # For the case we want the name
				my ($test_gaps) = $cds_sequence;
				my ($num_gaps) = ($test_gaps=~tr/\-/\-/);
				$cds_alignment->{$specie_name}->{'sequence'} = $cds_sequence;
				$cds_alignment->{$specie_name}->{'num_gaps'} = $num_gaps;
				$cds_alignment->{$specie_name}->{'length'} = length($cds_sequence);
			}
		}
	}
	return $cds_alignment;	
}

# Filter CDS alignment. We get the sequence alignment bigger than 50% exons
sub filter_cds_alignment($)
{
	my ($cds_alignment) = @_;
	my ($cds_filter_alignment);
	my ($target_name) = $USED_ORGANISMS->{$target}; # For the case we want the name
	
	my ($human_length) = $cds_alignment->{$target_name}->{'length'};
	while ( my ($specie_name,$specie_cds) = each(%{$cds_alignment}) )
	{
		if ( $specie_name eq $target_name )
		{
			$cds_filter_alignment->{$specie_name} = $specie_cds;
		}
		elsif ( $specie_cds->{'num_gaps'} < ($human_length/2) )
		{
			$cds_filter_alignment->{$specie_name} = $specie_cds; 
		}
	}
	return $cds_filter_alignment;
}

# Print CDS aligments as phylip or fasta format
sub print_cds_filter_aligment($$)
{
	my ($cds_filter_alignment,$format) = @_;
	my ($cds_alignment) = '';
	my ($target_name) = $USED_ORGANISMS->{$target}; # For the case we want the name

	if ( $format eq 'phylip' )
	{
		my ($num_species) = 0;
		my ($length_sequences) = length($cds_filter_alignment->{$target_name}->{'sequence'});
		my ($sequences)='';
		while ( my ($specie_name,$specie_cds) = each(%{$cds_filter_alignment}) )
		{
			$sequences .= $specie_name."\n".$specie_cds->{'sequence'}."\n";
			$num_species++;
		}
		$cds_alignment = $num_species.' '.$length_sequences."\n".$sequences;
	}
	elsif ($format eq 'fasta')
	{
		my ($cds_target_alignment_fasta) = '';
		my ($cds_alignment_fasta) = '';
		while( my ($specie_name,$specie_cds) = each(%{$cds_filter_alignment}) )
		{
			if ($specie_name eq $target_name) {
				$cds_target_alignment_fasta .= '>'.$specie_name."\n".$specie_cds->{'sequence'}."\n";
			}
			else {
				$cds_alignment_fasta .= '>'.$specie_name."\n".$specie_cds->{'sequence'}."\n";
			}
		}
		$cds_alignment = $cds_target_alignment_fasta.$cds_alignment_fasta;
	}
	return $cds_alignment;
}

main();

1;

__END__

=head1 NAME

getUCSCAlign

=head1 DESCRIPTION

Acquire Multiple Alignments from UCSC

=head1 SYNOPSIS

getUCSCAlign

=head2 Required arguments:

	--species= <Name of species -mammals->
	
	--data <Gencode data file>
	
	--translations=  <Gencode translations file>
	
	--outpath <Dir where will save the cds alignments>

=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
	
	
=head1 EXAMPLE

	perl getUCSCAlign.pl
	
		--species='Homo sapiens'
		
		--data=../../examples/ENSG00000140416/ENSG00000140416.annot.gtf
		
		--translations=../../examples/ENSG00000140416/ENSG00000140416.transl.fa
		
		--outpath=../../examples/ENSG00000140416/ENSG00000140416/		  

=head1 AUTHOR

Created by:
	Paolo Maietta -pmaietta@ext.cnio.es- 15-June-2009

Developed by:
	Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut