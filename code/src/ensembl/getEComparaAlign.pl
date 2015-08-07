#! /usr/bin/perl -W

use strict;
use warnings;
use FindBin;
use Config::IniFiles;
use Getopt::Long;
use Data::Dumper;
use Bio::SimpleAlign;
use Bio::AlignIO;
use IO::String;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use Bio::EnsEMBL::Registry;

use APPRIS::Parser qw( parse_gencode );
use APPRIS::Utils::Exception qw( throw warning deprecate );
use APPRIS::Utils::File qw( prepare_workspace printStringIntoFile );
use APPRIS::Utils::Logger;


###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD	
	$ENSEMBL_CONFIG_FILE
	$ORGANISMS
	$ALIGN_SPECIES_SEY_NAME
	$PROG_ALIGN_IN_SUFFIX
	$PROG_TREE_IN_SUFFIX
);

$LOCAL_PWD				= $FindBin::Bin;
#$ENSEMBL_CONFIG_FILE	= $LOCAL_PWD.'/../../conf/ensembldb.ini';
$ENSEMBL_CONFIG_FILE	= $ENV{APPRIS_CODE_CONF_DIR}.'/ensembl.ini';

# http://www.ensembl.org/info/genome/compara/analyses.html
$ORGANISMS = {
# 39 eutherian mammals EPO_LOW_COVERAGE
'homo_sapiens' => undef,
'gorilla_gorilla' => undef,
'pan_troglodytes' => undef,
'pongo_abelii' => undef,
'nomascus_leucogenys' => undef,
'macaca_mulatta' => undef,
'papio_anubis' => undef,
'chlorocebus_sabaeus' => undef,
'callithrix_jacchus' => undef,
'tarsius_syrichta' => undef,
'microcebus_murinus' => undef,
'otolemur_garnettii' => undef,
'mus_musculus' => undef,
'rattus_norvegicus' => undef,
'dipodomys_ordii' => undef,
'ictidomys_tridecemlineatus' => undef,
'cavia_porcellus' => undef,
'ochotona_princeps' => undef,
'oryctolagus_cuniculus' => undef,
'tupaia_belangeri' => undef,
'choloepus_hoffmanni' => undef,
'dasypus_novemcinctus' => undef,
'echinops_telfairi' => undef,
'loxodonta_africana' => undef,
'procavia_capensis' => undef,
'erinaceus_europaeus' => undef,
'sorex_araneus' => undef,
'myotis_lucifugus' => undef,
'pteropus_vampyrus' => undef,
'equus_caballus' => undef,
'felis_catus' => undef,
'canis_familiaris' => undef,
'ailuropoda_melanoleuca' => undef,
'mustela_putorius_furo' => undef,
'tursiops_truncatus' => undef,
'sus_scrofa' => undef,
'bos_taurus' => undef,
'ovis_aries' => undef,
'vicugna_pacos' => undef,
# 11 fish EPO_LOW_COVERAGE 
'danio_rerio' => undef,
'astyanax_mexicanus' => undef,
'gadus_morhua' => undef,
'oreochromis_niloticus' => undef,
'oryzias_latipes' => undef,
'xiphophorus_maculatus' => undef,
'poecilia_formosa' => undef,
'takifugu_rubripes' => undef,
'tetraodon_nigroviridis' => undef,
'gasterosteus_aculeatus' => undef,
'lepisosteus_oculatus' => undef,
};

$ALIGN_SPECIES_SEY_NAME = defined($ENV{APPRIS_TYPE_ALIGN_SPECIES_NAME}) ? $ENV{APPRIS_TYPE_ALIGN_SPECIES_NAME} : "mammals";


$PROG_ALIGN_IN_SUFFIX		= 'compara.faa';
$PROG_TREE_IN_SUFFIX		= 'compara.nh';

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($atype) = undef;
my ($e_version) = undef;
my ($species) = undef;
my ($data_file) = undef;
my ($transcripts_file) = undef;
my ($translations_file) = undef;
my ($outpath) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'atype=s'			=> \$atype,
	'species=s'			=> \$species,
	'e-version=s'		=> \$e_version,
	'data=s'			=> \$data_file,
	'transcripts=s'		=> \$transcripts_file,
	'translations=s'	=> \$translations_file,	
	'outpath=s'			=> \$outpath,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless(
	defined $species and
	defined $e_version and 
	defined $data_file and
	defined $outpath
){
	print `perldoc $0`;
	exit 1;
}

# Optional arguments
# get type of alignment: exon, cds
if ( defined $atype and ( (lc($atype) eq 'exon') or (lc($atype) eq 'cds') ) ) {
	$atype = lc($atype);
}
else {
	$atype = 'cds';
}

# Get conf file
my ($cfg) = new Config::IniFiles( -file =>  $ENSEMBL_CONFIG_FILE );

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
sub get_align_by_trans($$$);
sub create_phylo_clustalw($);
sub get_num_tree_nodes($);
sub remove_ghost_node_tree($);
sub remove_excess_node_tree($);
sub get_coords_by_trans($$$);
sub get_ensembl_specie_align_from_coord($$$$$);
sub get_ensembl_specie_aligntree_from_coord($$$$$);
sub get_align_strings($$$\$);

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
	
	# Get conf file
	$logger->info("-- -- load config file\n");
	my ($e_core_param) = {
			'-host'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host'),
			'-user'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user'),
			'-verbose'    => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose'),
			'-species'    => $species,
			'-db_version' => $e_version,			
	};
	my ($e_compara_param) = {
			'-host'		  => $cfg->val( 'ENSEMBL_COMPARA_REGISTRY', 'host'),
			'-dbname'     => $cfg->val( 'ENSEMBL_COMPARA_REGISTRY', 'dbname'),
			'-user'       => $cfg->val( 'ENSEMBL_COMPARA_REGISTRY', 'user'),
			'-verbose'    => $cfg->val( 'ENSEMBL_COMPARA_REGISTRY', 'verbose'),
			'-db_version' => $e_version,
	};
			
	# Ensembl connection details
	$logger->info("-- -- load ensembl compara registry\n");
	my ($e_compara_registry) = 'Bio::EnsEMBL::Registry';
	eval {
		$e_compara_registry->load_registry_from_db(%{$e_compara_param});
	};
	return undef if $@;	

	# Get cds data from input file
	$logger->info("-- get gtf data\n");
	my ($genes_data) = parse_gencode($data_file, $transcripts_file, $translations_file);
	#$logger->debug("-- gencode_data:\n".Dumper($genes_data)."\n");	
	unless ( defined $genes_data ) {
		$logger->error("can not parse data: $!\n");
	}
		
	# get the transcripts and their own cds's from the list of genes 
	$logger->info("-- get $atype coords by trans\n");
	my ($transc_align_list, $transc_aligntree_list) = get_coords_by_trans($e_compara_registry, $genes_data, $atype);
	$logger->debug("-- transc_align_list:\n".Dumper($transc_align_list)."\n");
	$logger->debug("-- transc_aligntree_list:\n".Dumper($transc_aligntree_list)."\n");

	# Get the alignment for cds regions
	$logger->info("-- get $atype aligns\n");
	if ( defined $transc_align_list  ) {
		get_align_by_trans($e_compara_registry, $transc_align_list, $transc_aligntree_list);
	}
		
	$logger->finish_log();
	
	exit 0;
}

sub get_align_by_trans($$$)
{
	my ($e_compara_registry, $transc_align_list, $transc_aligntree_list) = @_;		
		
	while ( my ($transcript_id, $cds_list) = each(%{$transc_align_list}) ) {
		$logger->info("-- -- id: $transcript_id\n");
		my ($trans_align);
		my ($trans_align_fasta_file) = $outpath.'/'.$transcript_id.'.'.$PROG_ALIGN_IN_SUFFIX;
		my ($trans_aligntree_nh_file) = $outpath.'/'.$transcript_id.'.'.$PROG_TREE_IN_SUFFIX;			
		
		foreach my $cds (@{$cds_list}) {
			my ($chr) = $cds->{'chr'};
			my ($start) = $cds->{'start'};
			my ($end) = $cds->{'end'};
			my ($strand) = 1;
			$strand = -1 if ( $cds->{'strand'} eq '-' );
			my ($seq_len) = abs($end - $start) +1;				
			$logger->info("-- -- -- get_ensembl_specie_align_from_coord ".$cds->{'id'}.">$chr:$start-$end:$strand);\n");
			my ($align_seq, $align_blocks) = get_ensembl_specie_align_from_coord($e_compara_registry, $chr, $start, $end, $strand);			
			#$logger->debug(Dumper($align_blocks)."\n");
			
			$logger->info("-- -- -- get_align_strings\n");
			get_align_strings($e_compara_registry, $align_seq, $align_blocks, $trans_align);
			#$logger->debug(Dumper($trans_align)."\n");
		}
		
		$logger->info("-- -- -- convert alignment and phylogenetic tree\n");
		my ($main_specie) = lc($species); $main_specie =~ s/^\s*//; $main_specie =~ s/\s*$//; $main_specie =~ s/\s/\_/;
		my ($trans_align_fasta_cont) = '';
		my ($trans_aligntree_nh_cont) = $transc_aligntree_list->{$transcript_id};
		while ( my ($name, $seq) = each(%{$trans_align}) ) {
			if ( $name =~ /^([\w{1}])(.)*\_(\w{3})[^\$]*$/ ) {
				if ( $seq !~ /^([\.|\-]*)$/ ) {
					my ($ref) = uc($1).$3;
					$trans_aligntree_nh_cont =~ s/$ref[^\:]*/$name/;
					$seq =~ s/N/-/mg; $seq =~ s/\./-/mg;
					$trans_align_fasta_cont .= ">".$name."\n".$seq."\n";
				}
				else {
					$logger->info("empty alignment: $name\n");
				}
			}
			else {
				$logger->error("don't match the name\n");
			}
		}
		
		# There is a bug in Bio::TreeIO (version of 2013-05-30). Works in the version 2013-02-20 of Bio::TreeIO
		$logger->info("-- -- -- remove nodes that exceed from the tree\n");
		$trans_aligntree_nh_cont = remove_excess_node_tree($trans_aligntree_nh_cont);
		
		# THIS METHOD IS A SHIT!!!!!
		#$logger->info("-- -- -- remove ghost node from the tree\n");
		#$trans_aligntree_nh_cont = remove_ghost_node_tree($trans_aligntree_nh_cont);

		if ( ($trans_align_fasta_cont ne '') and ($trans_aligntree_nh_cont ne '') ) {
			my ($printing_file_log2) = printStringIntoFile($trans_align_fasta_cont, $trans_align_fasta_file);
			$logger->error("Printing output") unless ( defined $printing_file_log2 );

			my ($printing_file_log3) = printStringIntoFile($trans_aligntree_nh_cont, $trans_aligntree_nh_file);
			$logger->error("Printing tree") unless ( defined $printing_file_log3 );	
			
			#$logger->info("-- -- -- phylogenetic tree generation using the clustalw\n");
			#create_phylo_clustalw($trans_align_fasta_file);
		}		
	}	
}

sub create_phylo_clustalw($)
{
	my ($align_fasta_file) = @_;
	eval {
		my ($cmd) = "clustalw2 -tree -clustering=Neighbour-joining -infile=$align_fasta_file 1>&2 ";				
		$logger->info("\t** script: $cmd\n");
		my ($cmd_stderr_txt) = '';
		local *STDERR;
		open(STDERR, "+<", \$cmd_stderr_txt);
		system($cmd);
		$logger->info("$cmd_stderr_txt\n");
	};
	$logger->error("creating phylogenetic tree") if $@;	
}

sub get_num_tree_nodes($)
{
	my ($treein_txt) = @_;
	my ($num) = 0;
	my ($io) = IO::String->new($treein_txt);
	my ($treein) = Bio::TreeIO->new(
								-format => 'newick',
								-fh => $io,
	);
	while ( my ($tin) = $treein->next_tree ) {
		for my $node ( $tin->get_nodes ) {
			if ( defined $node ) {
				if ( defined $node->id ) {
					$num++ if ( ($node->id ne '') and ($node->id ne '0') );
				}
			}
		}
	}
	return $num;
}

sub remove_ghost_node_tree($)
{
	my ($treein_txt) = @_;	
	my ($treeout_txt) = $treein_txt;
	my ($tmp_treeout_txt) = '';
		
	eval {
		my ($io) = IO::String->new($treein_txt);
		my ($treein) = Bio::TreeIO->new(
									-format => 'newick',
									-fh => $io,
		);
		open(my $fake_fh, "+<", \$tmp_treeout_txt);
		my ($treeout) = Bio::TreeIO->new(
									-format => 'newick',
									-fh => $fake_fh,
		);
		while ( my ($tin) = $treein->next_tree ) {
			for my $node ( $tin->get_nodes ) {
				if ( defined $node ) {
					if ( !$node->is_Leaf and defined $node->each_Descendent ) {
						my (@descendants) = $node->each_Descendent;
						if ( scalar(@descendants) == 1 ) {
							$tin->remove_Node($node);
						}					
					}
				}	
			}	
			$treeout->write_tree($tin);
		}
		close($fake_fh);		
	};
	return $treeout_txt if ($@);
	
	if ( $tmp_treeout_txt ne '' ) {
		$treeout_txt = $tmp_treeout_txt;
	}
	return $treeout_txt;
}

sub remove_excess_node_tree($)
{
	my ($treein_txt) = @_;
	my ($treeout_txt) = $treein_txt;
	my ($tmp_treeout_txt) = '';
	
	eval {
		my ($io) = IO::String->new($treein_txt);
		my ($treein) = Bio::TreeIO->new(
									-format => 'newick',
									-fh => $io,
		);
		open(my $fake_fh, "+<", \$tmp_treeout_txt);
		my ($treeout) = Bio::TreeIO->new(
									-format => 'newick',
									-fh => $fake_fh,
		);
		while ( my ($tin) = $treein->next_tree ) {
			for my $node ( $tin->get_nodes ) {
				if ( defined $node ) {
					if ( defined $node->id ) {
						my ($n_id) = $node->id;
						if ( ($n_id ne '') and ($n_id ne '0') ) {
							unless ( exists $ORGANISMS->{$n_id} ) { # remove a node of "excess" specie
								$tin->remove_Node($node);				
							}
						}					
					}
				}
			}	
			$treeout->write_tree($tin);
		}
		close($fake_fh);
	};
	return $treeout_txt if ($@);
	
	if ( $tmp_treeout_txt ne '' ) {
		$treeout_txt = $tmp_treeout_txt;
	}
	return $treeout_txt;
}

sub get_coords_by_trans($$$)
{
	my ($e_compara_registry, $genes_data, $atype) = @_;
	
	# get the transcripts and their coord from the list of genes 
	my ($align_list);
	my ($aligntree_list);
	foreach my $gene (@{$genes_data}) {
		my ($gene_id) = $gene->stable_id;
		my ($chr) = $gene->chromosome;
		$chr =~ s/^M$/MT/; # for ensembl api
		my ($g_start) = $gene->start;
		my ($g_end) = $gene->end;
		my ($g_strand) = 1;
		$g_strand = -1 if ( $gene->strand eq '-' );		
		$logger->info("-- id: $gene_id\n");
		
		if ( defined $gene->transcripts ) {

			# get gene tree
			$logger->info("-- -- -- get_ensembl_specie_aligntree_from_coord>$chr:$g_start-$g_end:$g_strand\n");
			my ($trans_aligntree_nh_cont) = get_ensembl_specie_aligntree_from_coord($e_compara_registry, $chr, $g_start, $g_end, $g_strand);
			foreach my $transcript (@{$gene->transcripts})
			{
				my ($transcript_id) = $transcript->stable_id;
				my ($t_start) = $transcript->start;
				my ($t_end) = $transcript->end;
				my ($t_strand) = 1;
				$t_strand = -1 if ( $transcript->strand eq '-' );				
				$logger->info("-- -- id: $transcript_id\n");
				
				if ( defined $atype ) {
					$logger->info("-- -- -- collect $atype list\n");
					if ( $atype eq 'exon' ) {
						if ( defined $transcript->exons ) {
							for (my $i = 0; $i < scalar(@{$transcript->exons}); $i++) {
								my ($exon) = $transcript->exons->[$i];
								my ($parameters) = {
												#id			=> $i+1,
												id			=> $exon->stable_id,
												chr             => $chr,
												start			=> $exon->start,
												end				=> $exon->end,
												strand			=> $exon->strand
								};
								push(@{$align_list->{$transcript_id}}, $parameters);
							}			
						}					
					}
					elsif ( $atype eq 'cds' ) {
						if ( defined $transcript->translate ) {
							my ($translation) = $transcript->translate;
							if ( defined $translation->cds ) {
								for (my $i = 0; $i < scalar(@{$translation->cds}); $i++) {
									my ($cds) = $translation->cds->[$i];
									my ($parameters) = {
													id			=> $i+1,
													chr             => $chr,
													start			=> $cds->start,
													end				=> $cds->end,
													strand			=> $cds->strand,
													phase			=> $cds->phase
									};
									push(@{$align_list->{$transcript_id}}, $parameters);
								}			
							}							
						}					
					}
					# print align tree
					#$logger->info("-- -- -- get_ensembl_specie_aligntree_from_coord>$chr:$t_start-$t_end:$t_strand\n");
					#my ($trans_aligntree_nh_cont) = get_ensembl_specie_aligntree_from_coord($e_compara_registry, $chr, $t_start, $t_end, $t_strand);
					#$aligntree_list->{$transcript_id} = $trans_aligntree_nh_cont;
					# save gene tree
					$logger->info("-- -- -- save gene tree\n");
					$aligntree_list->{$transcript_id} = $trans_aligntree_nh_cont;
				}		
			}			
		}	
	}
	
	return ($align_list, $aligntree_list);
}

sub get_ensembl_specie_align_from_coord($$$$$)
{
	my ($registry, $chr, $start, $end, $strand) = @_;
	my ($genomic_align_blocks);
	
	$logger->info("-- -- -- get compara adaptor\n");
	my ($genome_db_adaptor) = $registry->get_adaptor('Multi', 'compara', 'GenomeDB');
	throw("Cannot connect to Compara") if (!$genome_db_adaptor);
	
	$logger->info("-- -- -- get method link adaptor\n");
	my ($method_link_species_set_adaptor) = $registry->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');
	throw("Registry configuration file has no data for connecting ") if (!$method_link_species_set_adaptor);
	
	$logger->info("-- -- -- fetch method link by specie\n");	
	my ($method_link_species_set) = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name("EPO_LOW_COVERAGE", $ALIGN_SPECIES_SEY_NAME);
	
	$logger->info("-- -- -- get slice adapator\n");
	my ($slice_adaptor) = $registry->get_adaptor($species, 'core', 'Slice');
	
	$logger->info("-- -- -- get genomic align adaptor\n");
	my ($genomic_align_block_adaptor) = $registry->get_adaptor('Multi', 'compara', 'GenomicAlignBlock');

	$logger->info("-- -- -- fetch by region\n");
	my ($query_slice) = $slice_adaptor->fetch_by_region('toplevel', $chr, $start, $end, $strand);
	my ($query_length) = $query_slice->length(); # not used
	my ($query_seq) = $query_slice->seq();
	
	$logger->info("-- -- -- fetch all method link from specie\n");
	my ($restrict) = 1;	
	$genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set, $query_slice, undef, undef, $restrict);
	
	return ($query_seq, $genomic_align_blocks);
}

sub get_ensembl_specie_aligntree_from_coord($$$$$)
{
	my ($registry, $chr, $start, $end, $strand) = @_;
	my ($align_tree_nh);
	
	$logger->info("-- -- -- get compara adaptor\n");
	my ($genome_db_adaptor) = $registry->get_adaptor('Multi', 'compara', 'GenomeDB');
	throw("Cannot connect to Compara") if (!$genome_db_adaptor);
	
	$logger->info("-- -- -- get method link adaptor\n");
	my ($method_link_species_set_adaptor) = $registry->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');
	throw("Registry configuration file has no data for connecting ") if (!$method_link_species_set_adaptor);
	
	$logger->info("-- -- -- fetch method link by specie\n");	
	my ($method_link_species_set) = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name("EPO_LOW_COVERAGE", $ALIGN_SPECIES_SEY_NAME);
	
	$logger->info("-- -- -- get slice adapator\n");
	my ($slice_adaptor) = $registry->get_adaptor($species, 'core', 'Slice');
	
	$logger->info("-- -- -- fetch by region\n");
	my ($query_slice) = $slice_adaptor->fetch_by_region('toplevel', $chr, $start, $end, $strand);	

	$logger->info("-- -- -- get genomic tree adaptor\n");
	my ($genomic_align_tree_adaptor) = $registry->get_adaptor('Multi', 'compara', 'GenomicAlignTree');
	
	$logger->info("-- -- -- fetch all method link from specie\n");
	my ($restrict) = 1;
	my ($genomic_align_trees) = $genomic_align_tree_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set, $query_slice, undef, undef, $restrict);
	
	# We obtain the tree with max depth... If there are trees with the same depth, we obtain the last one... THIS WORKS???
	my ($max_depth);
	foreach my $genomic_align_tree (@{$genomic_align_trees}) {
		if (UNIVERSAL::isa($genomic_align_tree, "Bio::EnsEMBL::Compara::GenomicAlignTree")) {		
			my ($aln_tree_nh) = $genomic_align_tree->newick_format('simple');
			$aln_tree_nh =~ s/\[[+|-]+\]//mg; # remove strange characters (strand label)
			
			#my ($m_depth) = $genomic_align_tree->max_depth;
			my ($num_species) = scalar( @{[ $aln_tree_nh =~ /([A-Z]{1}[a-z]{3}\_)/g ]});
			my ($m_depth) = $num_species;
			unless ( defined $max_depth ) {
				$max_depth = $m_depth;
				$align_tree_nh = $aln_tree_nh;
			}
			else {
				if ( $m_depth >= $max_depth ) {
					$max_depth = $m_depth;
					$align_tree_nh = $aln_tree_nh;
				}				
			}
			
			#$align_tree_nh = $aln_tree_nh;
			#last; 	# Note: we obtain the first "align-tree"... Why? Because I THINK is the one that is the multiple alignment
		}
	}
	return $align_tree_nh;
}

sub get_align_strings($$$\$)
{
	my ($registry, $align_seq, $align_blocks, $ref_report) = @_;
	my ($main_specie) = lc($species); $main_specie =~ s/^\s*//; $main_specie =~ s/\s*$//; $main_specie =~ s/\s/\_/;	
	my ($report);
	my ($loc_rep);

	if ( defined $align_blocks and scalar(@{$align_blocks}) > 0 ) {
		# we have align blocks
		foreach my $align_block (@{$align_blocks}) {
			if (UNIVERSAL::isa($align_block, "Bio::EnsEMBL::Compara::GenomicAlignBlock")) {			
				my ($simple_align) = $align_block->get_SimpleAlign('uc');
				my ($join_rep);
				
				foreach my $aln ($simple_align->each_seq) {
					if ( $aln->id() =~ /^([^\/]*)\// ) {
						my ($spe) = $1;
						my ($id) = $aln->id();
						my ($seq) = $aln->seq();
						
						# join alignments with the same specie
						if ( exists $join_rep->{$spe} ) {				
							my $s1 = $join_rep->{$spe};
							my $s2 = $seq;
							my @s1 = split(//,$s1);
							my @s2 = split(//,$s2);
							my $s = '';
							for (my $i=0; $i < @s1; $i++) {
								if ( ($s1[$i] eq '.') and ($s2[$i] ne '.') )
									{ $s .= $s2[$i]; }
								elsif ( ($s1[$i] ne '.') and ($s2[$i] eq '.') )
									{ $s .= $s1[$i]; }
								elsif ( ($s1[$i] eq '.') and ($s2[$i] eq '.') )
									{ $s .= $s1[$i]; }
							}
							if ( $s ne '' ) {
								$join_rep->{$spe} = $s;
								$loc_rep->{$spe} = $s;
							}
						}
						else {
							$join_rep->{$spe} = $seq;
							$loc_rep->{$spe} = $seq;
						}
					}
				}	
				# if one specie has not alignment, then we print gaps
				foreach my $spe ( keys(%{$ORGANISMS}) ) {
					if ( exists $loc_rep->{$spe} ) {
						my ($s) = '';
						$s = $report->{$spe} if ( exists $report->{$spe} );
						$report->{$spe} = $loc_rep->{$spe}.$s;
					}
					else {
						my ($s) = '-' x length($loc_rep->{$main_specie});
						$report->{$spe} = $s;
					}
				}
			}	
		}
	}
	else {
		# in the case where the sequence of main specie has not alignment,
		# then we insert the seq and for the rest, we print gaps
		foreach my $spe ( keys(%{$ORGANISMS}) ) {
			if ( defined $align_seq and ($spe eq $main_specie) ) {
				my ($s) = '';
				$s = $report->{$spe} if ( exists $report->{$spe} );
				$report->{$spe} = $align_seq.$s;
			}
			else {
				my ($s) = '-' x length($align_seq);
				$report->{$spe} = $s;
			}			
		}
	}
				
	# save block aligns
	foreach my $spe ( keys(%{$ORGANISMS}) ) {
		if ( exists $report->{$spe} ) {
			${$ref_report}->{$spe} .= $report->{$spe};
		}
	}
}


main();

1;

__END__

=head1 NAME

getEComparaAlign

=head1 DESCRIPTION

Retrieves compara alignments of Ensembl from chromosome.

Acquire Phylogenetic tree using ClustalW program.

=head1 SYNOPSIS

getEComparaAlign

=head2 Required arguments:

	--species= <Name of species>
	
	--e-version= <Number of Ensembl version of identifier>
	
	--data <Gencode data file>
	
	--transcripts=  <Gencode transcript file>
	
	--translations=  <Gencode translations file>
	
	--outpath <Dir where will save the cds alignments>

=head2 Optional arguments:

	--atype=  <Type of alignment: [exon,cds] (default: cds)>

=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>
	
	
=head1 EXAMPLE

	perl getEComparaAlign.pl
				
		--species='Homo sapiens'
		
		--e-version=69
	
		--data=../../examples/ENSG00000140416/ENSG00000140416.annot.gtf
		
		--transcripts=../../examples/ENSG00000140416/ENSG00000140416.transc.fa
		
		--translations=../../examples/ENSG00000140416/ENSG00000140416.transl.fa
		
		--outpath=../../examples/ENSG00000140416/ENSG00000140416/
		  

=head1 AUTHOR

Created by:
	Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut