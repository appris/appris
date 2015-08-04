#!/usr/bin/perl -W

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Registry;
use APPRIS::Utils::File qw( prepare_workspace printStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
	$CONFIG_INI_APPRIS_DB_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;
$CONFIG_INI_APPRIS_DB_FILE	= $LOCAL_PWD.'/conf/apprisdb.ini';

# Input parameters
my ($species) = undef;
my ($method) = undef;
my ($position) = undef;
my ($type) = undef;
my ($output_file) = undef;
my ($apprisdb_conf_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'species=s'			=> \$species,
	'method=s'  		=> \$method,
	'position=s'		=> \$position,
	'output=s'			=> \$output_file,
	'apprisdb-conf=s'	=> \$apprisdb_conf_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless( defined $species and defined $method and defined $position and defined $output_file )
{
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
sub get_firestar_annot($$$);
sub get_matador3d_annot($$$);
sub get_spade_annot($$$);
sub get_inertia_annot($$$);


#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	# Output vars
	my ($methods);
	my ($out_reports);
	
	my ($output_content) = '';
	my ($output_content2) = '';
	
	# Get the list of methods
	foreach my $m ( split(',',$method) ) {
		$methods->{$m} = 1;
	}

	# APPRIS Registry from given specie
	my ($cfg) = new Config::IniFiles( -file => $apprisdb_conf_file );
	my ($spe) = $species; $spe =~ s/^\s*//; $spe =~ s/\s*$//; $spe =~ s/\s/\_/;	
	my ($specie_db) = uc($spe.'_db');
	my ($param) = {
			'-dbhost'       => $cfg->val($specie_db, 'host'),
			'-dbname'       => $cfg->val($specie_db, 'db'),
			'-dbuser'       => $cfg->val($specie_db, 'user'),
			'-dbpass'       => $cfg->val($specie_db, 'pass'),
			'-dbport'       => $cfg->val($specie_db, 'port'),
	};
	$logger->debug(Dumper($param)."\n");
	my ($registry) = APPRIS::Registry->new();
	$registry->load_registry_from_db(%{$param});
	
	my ($region);
	if ( $position =~ /^([^\:]*):([^\-]*)-([^\$]*)$/ ) {
		my ($position_chr) = $1;
		my ($position_start) = $2;
		my ($position_end) = $3;
		$region = $registry->fetch_basic_by_region($position_chr, $position_start, $position_end);			
	}
	elsif ( $position =~ /^chr([^\$]*)$/ ) {
		my ($position_chr) = $1;
		$region = $registry->fetch_basic_by_region($position_chr);
	}
	else {
		$region = $registry->fetch_basic_by_region($position);
	}
	$logger->info("OK\n");
	
	foreach my $gene (@{$region}) {
		my ($gene_id) = $gene->stable_id;
		$logger->info("$gene_id ---------------\n");

		foreach my $transcript (@{$gene->transcripts}) {			
			my ($transcript_id) = $transcript->stable_id;

			if ($transcript->translate and $transcript->exons) {
				my ($exons) = $transcript->exons;

				if ($transcript->translate) {
					my ($translate) = $transcript->translate;
					my ($translate_seq) = $translate->sequence;
										
					if ( defined $methods and (exists $methods->{'firestar'}) ) {
						my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id, 'firestar');
						if ( defined $analysis ) {
							$out_reports->{'firestar'} .= get_firestar_annot($gene_id,$transcript_id,$analysis);	
						}														
					}
					if ( defined $methods and (exists $methods->{'matador3d'}) ) {
						my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id, 'matador3d');
						if ( defined $analysis ) {
							$out_reports->{'matador3d'} .= get_matador3d_annot($gene_id,$transcript_id,$analysis);	
						}														
					}
					if ( defined $methods and (exists $methods->{'spade'}) ) {
						my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id, 'spade');
						if ( defined $analysis ) {
							$out_reports->{'spade'} .= get_spade_annot($gene_id,$transcript_id,$analysis);	
						}														
					}
					if ( defined $methods and (exists $methods->{'inertia'}) ) {
						my ($analysis) = $registry->fetch_analysis_by_stable_id($transcript_id, 'inertia');
						if ( defined $analysis ) {
							$out_reports->{'inertia'} .= get_inertia_annot($gene_id,$transcript_id,$analysis);	
						}														
					}
				}
			}
		}
	}
	
	# print outputs
	while ( my ($method,$output_content) = each(%{$out_reports}) ) {
		if ( $output_content ne '' ) {
			my ($printing_file_log) = printStringIntoFile($output_content, $output_file);
			$logger->error("Printing output") unless ( defined $printing_file_log );
		}		
	}
	
	$logger->info("retrieve_method_annots:$position:finished ---------------\n");
		
	$logger->finish_log();
	
	exit 0;
}

sub get_firestar_annot($$$)
{
	my ($gene_id,$transcript_id,$analysis) = @_;
	my ($report) = '';
	
	my ($method) = $analysis->firestar;
	
	if ( defined $method->residues ) {
		foreach my $region (@{$method->residues}) {
			if ( defined $region->ligands and defined $region->residue and defined $region->domain and 
				 defined $region->start and defined $region->end and defined $region->strand )
			{
			 	my ($ligand_list) = $region->ligands;
			 	my ($residue) = $region->residue;
			 	my ($domain) = $region->domain;
			 	#my ($start) = $region->start;
			 	#my ($end) = $region->end;
			 	#my ($strand) = $region->strand;
				my ($ligands) = '';				 	
			 	my (@ligand_list2) = split(/\|/,$ligand_list);
			 	foreach my $ligand (@ligand_list2) {
				 	$ligand =~ s/\[[^\]]*\]$//g;
				 	$ligands .= $ligand.'|';
			 	}
			 	if ( $ligands ne '' ) {
			 		$ligands =~ s/\|$//g;
				 	$report .= $ligands."\t".
							 	$residue."\t".
							 	$gene_id."\t".
							 	$transcript_id."\n";			 		
			 	}
			}
		}
	}
	return $report;
} # end get_firestar_annot

sub get_matador3d_annot($$$)
{
	my ($gene_id,$transcript_id,$analysis) = @_;
	my ($report) = '';
	
	my ($method) = $analysis->matador3d;
	
	if ( defined $method->alignments ) {
		foreach my $region (@{$method->alignments}) {
			if ( defined $region->pdb_id and defined $region->pstart and defined $region->pend and defined $region->score and 
				 defined $region->start and defined $region->end and defined $region->strand ) {
				 	my ($pstart) = $region->pstart;
				 	my ($pend) = $region->pend;
				 	my ($score) = $region->score;
				 	my ($align_start) = $region->alignment_start;
				 	my ($align_end) = $region->alignment_end;
				 	#my ($start) = $region->start;
				 	#my ($end) = $region->end;
				 	#my ($strand) = $region->strand;
				 	my ($pdb_id) = $region->pdb_id;
				 	my ($identity) = $region->identity;				 	
				 	
				 	$report .= $pdb_id."\t".
							 	$pstart."\t".
							 	$pend."\t".
							 	$align_start."\t".
							 	$align_end."\t".
							 	$identity."\t".
							 	$gene_id."\t".
							 	$transcript_id."\n";							 	
			}
		}
	}
	return $report;
} # end get_matador3d_annot

sub get_spade_annot($$$)
{
	my ($gene_id,$transcript_id,$analysis) = @_;
	my ($report) = '';
	
	my ($method) = $analysis->spade;
	
	if ( defined $method->regions ) {
		foreach my $region (@{$method->regions}) {
			if ( defined $region->hmm_name and defined $region->hmm_acc and
				 defined $region->alignment_start and defined $region->alignment_end and defined $region->evalue and 
				 defined $region->start and defined $region->end and defined $region->strand ) {
				 	my ($hmm_name) = $region->hmm_name;
				 	my ($hmm_acc) = $region->hmm_acc;
				 	my ($pstart) = $region->alignment_start;
				 	my ($pend) = $region->alignment_end;
				 	my ($evalue) = $region->evalue;
				 	my ($start) = $region->start;
				 	my ($end) = $region->end;
				 	my ($strand) = $region->strand;
				 	
					 	$report .= $hmm_acc."\t".
							 	$pstart."\t".
							 	$pend."\t".
							 	$start."\t".
							 	$end."\t".
							 	$hmm_name."\t".
							 	$evalue."\t".
							 	$gene_id."\t".
							 	$transcript_id."\n";							 	
			}
		}
	}
	return $report;
} # end get_spade_annot

sub get_inertia_annot($$$)
{
	my ($gene_id,$transcript_id,$analysis) = @_;
	my ($report) = '';
	my ($rep_regions);
	
	my ($method) = $analysis->inertia;
	
	if ( defined $method ) {
		foreach my $type ('filter','prank','kalign') {				
	
			# get the type of alignment
			my ($alignment);
			if ( ($type eq 'filter') and $method->mafft_alignment) {
				$alignment = $method->mafft_alignment;
			}
			elsif ( ($type eq 'prank') and $method->prank_alignment) {
				$alignment = $method->prank_alignment;
			}
			elsif ( ($type eq 'kalign') and $method->kalign_alignment) {
				$alignment = $method->kalign_alignment;
			}
			else {
				next; # we don't have aligment
			}
			
			if ( defined $alignment->result and defined $alignment->unusual_evolution ) {
				
				# get regions
				if ( defined $alignment->regions ) {
					for (my $icds = 0; $icds < scalar(@{$alignment->regions}); $icds++) {
						my ($region2) = $alignment->regions->[$icds];
						if ( defined $region2->start and defined $region2->end and defined $region2->unusual_evolution and defined $region2->omega_mean ) {
							my ($reg_index) = $region2->start.':'.$region2->end;
							# get the first omega mean
							unless ( exists $rep_regions->{$reg_index} ) {
								$rep_regions->{$reg_index} = $region2->omega_mean;
							}
							else {
								# if omega mean already exists, we get the smaller one
								my ($old_omega_mean) = $rep_regions->{$reg_index};
								if ( $region2->omega_mean < $old_omega_mean ) { 
									$rep_regions->{$reg_index} = $region2->omega_mean;
								}
							}
						}
						else {
							$logger->warning("No inertia-omega $type residues");
						}
					}
				}			
			}
			else {
				$logger->warning("No inertia-omega $type analysis");
			}
		}
	}
	
	# print report
	foreach my $r_index (sort { $rep_regions->{$b} <=> $rep_regions->{$a} } keys %{$rep_regions}) {
		my ($r_omega) = $rep_regions->{$r_index};
		my (@index) = split(':',$r_index);
		my ($r_start) = $index[0];
		my ($r_end) = $index[1];
		$report .=	$r_omega."\t".
					$r_start."\t".
					$r_end."\t".
				 	$gene_id."\t".
				 	$transcript_id."\n";		
	}		
	
	return $report;
} # end get_inertia_annot

main();

1;


__END__

=head1 NAME

retrieve_method_annots

=head1 DESCRIPTION

Exports the annotations for each method:

	* firestar: list of ligands/compounds.
	
	* Matador3D: list of pdb's, and coordinate alignments (at the moment, it's reported cds coordinates).
	
	* SPADE: list of domains, and alignment coordinates.
	
	* INERTIA: list of exons and the best score of slr (omega), the smaller score.
	

=head1 SYNOPSIS

retrieve_method_annots

=head2 Required arguments:

	--species= <Name of species: Homo sapiens, Mus musculus, etc>
	
	--position= <Genome position> (chr21 or chr1:109102711-109187522)
	
	--method= <Name of Method(s)> (firestar,matador3d,spade,inertia)	
	
	--output= <Output file>
	
=head2 Optional arguments:

	--apprisdb-conf <Config file of APPRIS database (default: 'conf/apprisdb.ini' file)>
		
=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

	perl retrieve_method_annots.pl
		
		--species='Homo sapiens'
		
		--position=chr21
		
		--method=firestar,matador3d,spade,inertia
				
		--output=data/appris_data.annots.txt
	

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut