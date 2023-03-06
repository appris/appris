=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::GTF - Utility functions for info exporting

=head1 SYNOPSIS

  use APPRIS::Export::GTF
    qw(
       get_trans_annotations
     );

  or to get all methods just

  use APPRIS::Export::GTF;

  eval { get_trans_annotations($feature,$version) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Retrieves information of transcript as GTF format.

=head1 METHODS

=cut

package APPRIS::Export::GTF;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::Constant qw(
	$OK_LABEL
	$NO_LABEL
	$UNKNOWN_LABEL
);

###################
# Global variable #
###################
use vars qw(
	$GTF_CONSTANTS
);

$GTF_CONSTANTS = {
	'appris'=>{
		'source'=>'APPRIS',
		'type'=>'principal_isoform',
		'annotation'=>['Principal Isoform','Possible Principal Isoform','No Principal Isoform']		
	},
	'firestar'=>{
		'source'=>'FIRESTAR',
		'type'=>'functional_residue',
	},
	'matador3d'=>{
		'source'=>'MATADOR3D',
		'type'=>'homologous_structure',
	},
	'matador3d2'=>{
		'source'=>'MATADOR3D2',
		'type'=>'homologous_structure',
	},
	'spade'=>{
		'source'=>'SPADE',
		'type'=>'functional_domain',
	},
	'corsair'=>{
		'source'=>'CORSAIR',
		'type'=>'vertebrate_conservation',
		'annot'=>['conservation', 'doubtful_conservation', 'no_conservation']
	},	
	'corsair_alt'=>{
		'source'=>'CORSAIR_ALT',
		'type'=>'vertebrate_conservation',
		'annot'=>['conservation', 'doubtful_conservation', 'no_conservation']
	},	
	'inertia'=>{
		'source'=>'INERTIA',
		'type'=>'neutral_evolution',
		'annot'=>['neutral_evolution','unusual_evolution'],
	},
	'crash'=>{
		'source'=>'CRASH',
		'type'=>['signal_peptide','mitochondrial_signal'],
		'annot'=>{
			'signal_peptide' => ['pep_signal','doubtful_pep_signal','no_pep_signal'],
			'mitochondrial_signal' => ['mit_signal','doubtful_mit_signal','no_mit_signal']
		}
	},
	'thump'=>{
		'source'=>'THUMP',
		'type'=>'transmembrane_signal',
		'annot'=>['tmh_signal','damaged_tmh_signal']
	},
	'proteo'=>{
		'source'=>'PROTEO',
		'type'=>'proteomic_evidence',
	},
	'trifid'=>{
		'source'=>'TRIFID',
		'type'=>'functional_importance',
	},
};

=head2 get_trans_annotations

  Arg [1]    : APPRIS::Transcript or undef
  Arg [2]    : String - $soure 
               List of sources
  Arg [3]    : String - $version
  Example    : $annot = get_trans_annotations($feature,$version);
  Description: Retrieves tabular information of transcript.
  Returntype : String or undef

=cut

sub get_trans_annotations {
    my ($feature, $source) = @_;
    my ($output) = '';

    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {
   	    
		if ($feature->stable_id) {
			my ($gene_id);
			my ($transcript_id) = $feature->stable_id;
			my ($external_id) = $feature->external_name;
			if ($feature->xref_identify) {
				foreach my $xref_identify (@{$feature->xref_identify}) {
					if ( $xref_identify->dbname eq 'Gene_Id' ) {
						$gene_id = $xref_identify->id;							
					}
				}		
			}
			if ( ($source =~ /appris/) or ($source eq 'all') ) { # print appris annotation for all transcripts
				$output .= get_appris_annotations(	$transcript_id,
   		        										$gene_id,
   		        										$external_id,
           												$feature);
			}
			if ($feature->translate and $feature->translate->sequence) { # proteing coding methods
				if ( ($source =~ /firestar/) or ($source eq 'all') ) {
					$output .= get_firestar_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /^matador3d$/) or ($source eq 'all') ) {
					$output .= get_matador3d_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /^matador3d2$/) or ($source eq 'all') ) {
					$output .= get_matador3d2_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /corsair/) or ($source eq 'all') ) {
					$output .= get_corsair_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /corsair_alt/) or ($source eq 'all') ) {
					$output .= get_corsair_alt_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /spade/) or ($source eq 'all') ) {
					$output .= get_spade_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /inertia/) or ($source eq 'all') ) {
					$output .= get_inertia_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /thump/) or ($source eq 'all') ) {
					$output .= get_thump_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /crash/) or ($source eq 'all') ) {
					$output .= get_crash_annotations(	$transcript_id,
	   		        										$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /proteo/) or ($source eq 'all') ) {
					$output .= get_proteo_annotations(	$transcript_id,
						   		        					$gene_id,
	   		        										$external_id,
	           												$feature);
				}
				if ( ($source =~ /trifid/) or ($source eq 'all') ) {
					$output .= get_trifid_annotations(	$transcript_id,
														$gene_id,
														$external_id,
														$feature);
				}
			}
		}
    }
    else {
		throw('Argument must be an APPRIS::Transcript');
   	}
	return $output;
}

=head2 print_annotations

  Arg [1]    : String - common attributes
  Arg [2]    : String - optional attributes
  Example    : $annot = print_annotations($feature,$version);
  Description: Print the (common and optional) annotations.
  Returntype : String or undef

=cut

sub print_annotations {
	my ($common, $optional) = @_;
	
	my ($output)='';

	# <seqname> <source> <feature> <start> <end>  <score> <strand> <frame> [attributes] [comments]
	#  Print feature: Common feature
	$output .= 	$common->{'seqname'}."\t".
				$common->{'source'}."\t".
				$common->{'type'}."\t".
				$common->{'start'}."\t".
				$common->{'end'}."\t".
				$common->{'score'}."\t".
				$common->{'strand'}."\t".
				$common->{'phase'}."\t";

	#  Print feature
	while ( my ($key,$value) = each(%{$optional}) )
	{
		$output .= $key.' "'.$value.'"; ' if ( defined $value );
	}
	$output=~s/\; $/\n/;
	
	return $output;
}

=head2 get_appris_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_appris_annotations($trans_id, $gen_id, $ext_id, $feat, $v);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_appris_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;
	
    my ($output) = '';
 	my ($biotype_annot) = $feature->biotype;
 	my ($status_annot) = $feature->status;
 	my ($ccds_id);
 	my ($xref_ids);
	if ($feature->xref_identify) {
		foreach my $xref_identify (@{$feature->xref_identify}) {
			if ( $xref_identify->dbname eq 'CCDS' ) {
				$ccds_id = $xref_identify->id;							
			}
			elsif ( $xref_identify->dbname =~ /_Gene_Id$/ or $xref_identify->dbname =~ /_Transcript_Id$/ ) {
				my ($db) = lc($xref_identify->dbname);
				$xref_ids->{$db} = $xref_identify->id;
			}
		}		
	}
	my ($tsl_annot) = $feature->tsl;
	my ($tag_annot) = $feature->tag;	
	my ($length_na);
	if ($feature->sequence) {
		$length_na = length($feature->sequence);
	}
	my ($length_aa);
	if ($feature->translate and $feature->translate->sequence) {
		$length_aa = length($feature->translate->sequence);
	}
		
	my ($no_codons) = 'start/stop';
	if ( $feature->translate and $feature->translate->codons ) {
		my ($start_found) = 0;
		my ($stop_found) = 0;
		foreach my $codon (@{$feature->translate->codons}) {
			if ( $codon->type eq 'start' ) {
				$start_found = 1;
			} elsif ( $codon->type eq 'stop' ) {
				$stop_found = 1;
			}
		}
		if ( $start_found && $stop_found ) {
			undef $no_codons;
		} elsif ( ! $start_found ) {
			$no_codons = 'start';
		} elsif ( ! $stop_found ) {
			$no_codons = 'stop';
		}
	}

	# Get annotations
	my ($method_annot);
	my ($reliability);
	my ($all_rejected);
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->appris ) {
	 		my ($method) = $analysis->appris;
	 		
	 		# print principal isoform feature
			if ( $method->principal_isoform_signal ) {
				my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
				my ($method_source) = $GTF_CONSTANTS->{'appris'}->{'source'};
				my ($method_type) = $GTF_CONSTANTS->{'appris'}->{'type'};
				my ($method_start,$method_end,$method_strand);			
				if ( $feature->start and $feature->end and $feature->strand ) {
					$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
				}
				else {
					$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
				}
				my ($method_score) = '.';
				my ($method_phase) = '.';
				my ($appris_annot) = $method->principal_isoform_signal;		
				if ($appris_annot eq $APPRIS::Utils::Constant::NO_LABEL) {		
					$method_annot = $GTF_CONSTANTS->{'appris'}->{'annotation'}->[2];
				}
				elsif ($appris_annot eq $APPRIS::Utils::Constant::UNKNOWN_LABEL) {
					$method_annot = $GTF_CONSTANTS->{'appris'}->{'annotation'}->[1];		
				}
				elsif ($appris_annot eq $APPRIS::Utils::Constant::OK_LABEL) {
					$method_annot = $GTF_CONSTANTS->{'appris'}->{'annotation'}->[0];
				}
				$reliability = $method->reliability if ( $method->reliability );
				my ($common) = {
						'seqname'	=> $method_seqname,
						'source'	=> $method_source,
						'type'		=> $method_type,
						'start'		=> $method_start,
						'end'		=> $method_end,
						'score'		=> $method_score,
						'strand'	=> $method_strand,
						'phase'		=> $method_phase
				};
				my ($optional) = {
						'gene_id'			=> $gene_id,
						'transcript_id'		=> $transcript_id,
						'transcript_name'	=> $external_id,
						'biotype'			=> $biotype_annot,
						'status'			=> $status_annot,
						'tsl'				=> $tsl_annot,
						'tag'				=> $tag_annot,				
				};
				$optional->{'ccds_id'}			= $ccds_id if (defined $ccds_id);
				$optional->{'length_na'}		= $length_na if (defined $length_na);
				$optional->{'length_aa'}		= $length_aa if (defined $length_aa);
				$optional->{'no_codons'}		= $no_codons if (defined $no_codons);
				$optional->{'annotation'}		= $method_annot if (defined $method_annot);
				$optional->{'reliability'}		= $reliability if (defined $reliability);
				$optional->{'all_rejected'}		= $all_rejected if (defined $all_rejected);
				while ( my ($x_id, $x_val) = each(%{$xref_ids}) ) { $optional->{$x_id} = $x_val }
				
				if (defined $common and defined $optional) { # print output
					$output .= print_annotations($common,$optional);			
				}
			}
			# functional residues feature
			if ( defined $method->functional_residues_signal ) {				
				my ($method_score) = 0;
				my ($method_type) = $GTF_CONSTANTS->{'firestar'}->{'type'};
				my ($annot) = $method->functional_residues_signal;
				if ( defined $method->functional_residues_score ) {
					$method_score = $method->functional_residues_score;
				}
				# print output
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			}
			# homolog structure feature
			if ( defined $method->homologous_structure_signal ) {				
				my ($method_score) = 0;
				my ($method_type) = $GTF_CONSTANTS->{'matador3d'}->{'type'};
				my ($annot) = $method->homologous_structure_signal;				
				if ( defined $method->homologous_structure_score ) {
					$method_score = $method->homologous_structure_score;
				}
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			}	
			# vertebrate conservation features
			if ( defined $method->vertebrate_conservation_signal ) {
				my ($method_score) = 0;
				my ($method_type) = $GTF_CONSTANTS->{'corsair'}->{'type'};							
				my ($annot) = $method->vertebrate_conservation_signal;
				if ( defined $method->vertebrate_conservation_score ) {
					$method_score = $method->vertebrate_conservation_score;
				}
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			} 		
			# pfam whole domains features	
			if ( defined $method->domain_signal ) {
			    my ($method_score) = '.';
				my ($method_type) = $GTF_CONSTANTS->{'spade'}->{'type'};
				my ($annot) = $method->domain_signal;
				if ( defined $method->domain_score ) {
					$method_score = $method->domain_score;	
				}				
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			}	
			# unsual exon features
			if ( defined $method->unusual_evolution_signal ) {
				my ($method_score) = '.';
				my ($method_type) = $GTF_CONSTANTS->{'inertia'}->{'type'};
				my ($annot) = $method->unusual_evolution_signal;
				if ( defined $method->unusual_evolution_score ) {
					$method_score = $method->unusual_evolution_score;	
				}				
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			}	
 		
			# helix transmembrane features
			if ( defined $method->transmembrane_helices_signal ) {
				my ($method_score) = '.';
				my ($method_type) = $GTF_CONSTANTS->{'thump'}->{'type'};
				my ($annot) = $method->transmembrane_helices_signal;
				if ( defined $method->transmembrane_helices_score ) {
					$method_score = $method->transmembrane_helices_score;	
				}				
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			}	
 		
			# signal peptide features
			if ( defined $method->peptide_signal and defined $method->mitochondrial_signal ) {
				my ($method_score) = '.';
				if ( defined $method->peptide_signal ) {	
					my ($annot) = $method->peptide_signal;
					my ($method_type) = $GTF_CONSTANTS->{'crash'}->{'type'}->[0];
					$output .= get_appris_annotations2(	$transcript_id,
														$gene_id,
														$external_id,
														$feature,
														$method_type,
														$method_score,
														$annot);				
				}	
				if ( defined $method->mitochondrial_signal ) {				
					my ($annot) = $method->mitochondrial_signal;
					my ($method_type) = $GTF_CONSTANTS->{'crash'}->{'type'}->[1];
					$output .= get_appris_annotations2(	$transcript_id,
														$gene_id,
														$external_id,
														$feature,
														$method_type,
														$method_score,
														$annot);				
				}				
			}
			# peptide evidence feature
			if ( defined $method->peptide_evidence_score ) {
				my ($method_score) = 0;
				my ($method_type) = $GTF_CONSTANTS->{'proteo'}->{'type'};
				my ($annot) = $method->peptide_evidence_signal;
				if ( defined $method->peptide_evidence_score ) {
					$method_score = $method->peptide_evidence_score;
				}
				# print output
				$output .= get_appris_annotations2(	$transcript_id,
													$gene_id,
													$external_id,
													$feature,
													$method_type,
													$method_score,
													$annot);				
			} 		
 		
 		}
	 	else
	 	{
	 		# There is not APPRIS annotations => Retrieves common features
	 		my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
			my ($method_source) = $GTF_CONSTANTS->{'appris'}->{'source'};
			my ($method_type) = $GTF_CONSTANTS->{'appris'}->{'type'};
			my ($method_start,$method_end,$method_strand);			
			if ( $feature->start and $feature->end and $feature->strand ) {
				$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
			}
			else {
				$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
			}
			my ($method_score) = '.';
			my ($method_phase) = '.';			
			my ($common) = {
					'seqname'	=> $method_seqname,
					'source'	=> $method_source,
					'type'		=> $method_type,
					'start'		=> $method_start,
					'end'		=> $method_end,
					'score'		=> $method_score,
					'strand'	=> $method_strand,
					'phase'		=> $method_phase
			};
			my ($optional) = {
					'gene_id'			=> $gene_id,
					'transcript_id'		=> $transcript_id,
					'transcript_name'	=> $external_id,
					'biotype'			=> $biotype_annot,
					'status'			=> $status_annot,				
			};
			if (defined $common and defined $optional) { # print output
				$output .= print_annotations($common,$optional);			
			}
	 	}
 	}
 	
	return $output;
	
} # End get_appris_annotations

sub get_appris_annotations2 {
	my ($transcript_id, $gene_id, $external_id, $feature, $method_type, $method_score, $appris_annot) = @_;
	
	my ($output) = '';
	my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_start,$method_end,$method_strand);			
	if ( $feature->start and $feature->end and $feature->strand ) {
		$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
	}
	else {
		my ($length_aa);
		if ($feature->translate and $feature->translate->sequence) {
			$length_aa = length($feature->translate->sequence);
		}
		$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
	}
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'appris'}->{'source'};

	# Common attributes
	my ($common) = {
			'seqname'	=> $method_seqname,
			'source'	=> $method_source,
			'type'		=> $method_type,
			'start'		=> $method_start,
			'end'		=> $method_end,
			'score'		=> $method_score,
			'strand'	=> $method_strand,
			'phase'		=> $method_phase
	};
	
	# Optinal attributes
	my($optional);
	$optional->{'gene_id'}			= $gene_id;
	$optional->{'transcript_id'}	= $transcript_id;
	$optional->{'transcript_name'}	= $external_id;
	$optional->{'annotation'}		= $appris_annot if ( defined $appris_annot );
	if (defined $common and defined $optional) {
		$output .= print_annotations($common,$optional);			
	}
	
	return $output;
	
} # End get_appris_annotations2

=head2 get_firestar_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version  
  Example    : $annot = get_firestar_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_firestar_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = 0;
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'firestar'}->{'source'};
	my ($method_type) = $GTF_CONSTANTS->{'firestar'}->{'type'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->firestar ) {
	 		my ($method) = $analysis->firestar;
	 		# get residue annotations
			if ( defined $method->residues ) {
				foreach my $region (@{$method->residues}) {
					if ( defined $region->residue ) {
						# get coords
						my ($method_start)  = ( defined $region->start )  ? $region->start  : $region->residue;
						my ($method_end)    = ( defined $region->end )    ? $region->end    : $region->residue;
						my ($method_strand) = ( defined $region->strand ) ? $region->strand : '.';
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						$optional->{'note'}				= "pep_position:".$region->residue if ($region->residue);

						my (@ligands);
						my (@ligand_fields) = split(/\|/, $region->ligands);
						foreach my $ligand_field (@ligand_fields) {
							if ( $ligand_field =~ /^([^[]+)\[/ ) {
								my ($ligand) = $1;
								if ( ! grep { $_ eq $ligand } @ligands ) {
									push(@ligands, $ligand);
								}
							}
						}
						if (@ligands) {
							$optional->{'note'}			.= ",ligands:".join('|', @ligands);
						}

						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_firestar_annotations

=head2 get_matador3d_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version  
  Example    : $annot = get_matador3d_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_matador3d_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'matador3d'}->{'source'};
	my ($method_type) = $GTF_CONSTANTS->{'matador3d'}->{'type'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->matador3d ) {
	 		my ($method) = $analysis->matador3d;
	 		
	 		# get method annotations
	 		if ( defined $method->score ) {
	 			# get coords
				my ($method_start,$method_end,$method_strand);			
				if ( $feature->start and $feature->end and $feature->strand ) {
					$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
				}
				else {
					my ($length_aa);
					if ($feature->translate and $feature->translate->sequence) {
						$length_aa = length($feature->translate->sequence);
					}
					$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
				}
				my ($method_score)  = ( defined $method->score )  ? $method->score  : 0;
				# common attributes
				my ($common) = {
						'seqname'	=> $method_seqname,
						'source'	=> $method_source,
						'type'		=> $method_type,
						'type'		=> $method_type,
						'start'		=> $method_start,
						'end'		=> $method_end,
						'score'		=> $method_score,
						'strand'	=> $method_strand,
						'phase'		=> $method_phase
				};
				# optinal attributes
				my($optional);
				$optional->{'gene_id'}			= $gene_id;
				$optional->{'transcript_id'}	= $transcript_id;
				if (defined $common and defined $optional) {
#					$output .= print_annotations($common,$optional);			
				}
	 		}
	 		
	 		# get residue annotations
			if ( defined $method->alignments ) {
				foreach my $region (@{$method->alignments}) {
					if ( defined $region->score and $region->pdb_id ) {
						# get coords
						my ($method_start)  = ( defined $region->start )  ? $region->start  : $region->score;
						my ($method_end)    = ( defined $region->end )    ? $region->end    : $region->score;
						my ($method_score)  = ( defined $region->score )  ? $region->score  : 0;
						my ($method_strand) = ( defined $region->strand ) ? $region->strand : '.';
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'type'		=> $method_type,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						$optional->{'note'}				= "pdb_id:".$region->pdb_id if ($region->pdb_id);
						$optional->{'note'}				.= ",identity:".$region->identity if ($region->identity);
						if ($region->pstart and $region->pend) {
							$optional->{'note'}			.= ",pep_start:".$region->pstart.",pep_end:".$region->pend;
						}						
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_matador3d_annotations
sub get_matador3d2_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'matador3d2'}->{'source'};
	my ($method_type) = $GTF_CONSTANTS->{'matador3d2'}->{'type'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->matador3d2 ) {
	 		my ($method) = $analysis->matador3d2;
	 		
	 		# get method annotations
	 		if ( defined $method->score and !defined $method->alignments ) {
	 			# get coords
				my ($method_start,$method_end,$method_strand);			
				if ( $feature->start and $feature->end and $feature->strand ) {
					$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
				}
				else {
					my ($length_aa);
					if ($feature->translate and $feature->translate->sequence) {
						$length_aa = length($feature->translate->sequence);
					}
					$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
				}
				my ($method_score)  = ( defined $method->score )  ? $method->score  : 0;
				# common attributes
				my ($common) = {
						'seqname'	=> $method_seqname,
						'source'	=> $method_source,
						'type'		=> $method_type,
						'start'		=> $method_start,
						'end'		=> $method_end,
						'score'		=> $method_score,
						'strand'	=> $method_strand,
						'phase'		=> $method_phase
				};
				# optinal attributes
				my($optional);
				$optional->{'gene_id'}			= $gene_id;
				$optional->{'transcript_id'}	= $transcript_id;
				if (defined $common and defined $optional) {
#					$output .= print_annotations($common,$optional);			
				}
	 		}
	 		
	 		# get residue annotations
			if ( defined $method->alignments ) {
				foreach my $region (@{$method->alignments}) {
					if ( defined $region->score and $region->pdb_id ) {
						# get coords
						my ($method_start)  = ( defined $region->start )  ? $region->start  : ( ( defined $region->pstart )  ? $region->pstart  : '.' );
						my ($method_end)    = ( defined $region->end   )  ? $region->end    : ( ( defined $region->pend   )  ? $region->pend    : '.' );
						my ($method_score)  = ( defined $region->score )  ? $region->score  : 0;
						my ($method_strand) = ( defined $region->strand ) ? $region->strand : '.';
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'type'		=> $method_type,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						$optional->{'note'}				= "pdb_id:".$region->pdb_id if ($region->pdb_id);
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_matador3d2_annotations

=head2 get_corsair_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_corsair_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_corsair_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = 0;
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'corsair'}->{'source'};

	# Get annotations
	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->corsair ) {
	 		my ($method) = $analysis->corsair;	 		
			if ( defined $method->score and $method->score != 0 ) {
				# get coords
				my ($method_start,$method_end,$method_strand);			
				if ( $feature->start and $feature->end and $feature->strand ) {
					$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
				}
				else {
					my ($length_aa);
					if ($feature->translate and $feature->translate->sequence) {
						$length_aa = length($feature->translate->sequence);
					}
					$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
				}			
				# get type
				my ($method_type) = $GTF_CONSTANTS->{'corsair'}->{'type'};				
				# common attributes
				my ($common) = {
						'seqname'	=> $method_seqname,
						'source'	=> $method_source,
						'type'		=> $method_type,
						'start'		=> $method_start,
						'end'		=> $method_end,
						'score'		=> $method->score,
						'strand'	=> $method_strand,
						'phase'		=> $method_phase
				};
				# optinal attributes
				my($optional);
				$optional->{'gene_id'}			= $gene_id;
				$optional->{'transcript_id'}	= $transcript_id;
				if (defined $common and defined $optional) {
					$output .= print_annotations($common,$optional);			
				}													 	
			}		
 		}
 	}

	return $output;
	
} # End get_corsair_annotations

=head2 get_corsair_alt_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_corsair_alt_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_corsair_alt_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = 0;
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'corsair_alt'}->{'source'};

	# Get annotations
	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->corsair_alt ) {
	 		my ($method) = $analysis->corsair_alt;	 		
			if ( defined $method->score and $method->score != 0 ) {
				# get coords
				my ($method_start,$method_end,$method_strand);			
				if ( $feature->start and $feature->end and $feature->strand ) {
					$method_start = $feature->start;$method_end = $feature->end;$method_strand = $feature->strand;			
				}
				else {
					my ($length_aa);
					if ($feature->translate and $feature->translate->sequence) {
						$length_aa = length($feature->translate->sequence);
					}
					$method_start = 1;$method_end = ( defined $length_aa ) ? $length_aa : 1; $method_strand = '.';
				}			
				# get type
				my ($method_type) = $GTF_CONSTANTS->{'corsair_alt'}->{'type'};				
				# common attributes
				my ($common) = {
						'seqname'	=> $method_seqname,
						'source'	=> $method_source,
						'type'		=> $method_type,
						'start'		=> $method_start,
						'end'		=> $method_end,
						'score'		=> $method->score,
						'strand'	=> $method_strand,
						'phase'		=> $method_phase
				};
				# optinal attributes
				my($optional);
				$optional->{'gene_id'}			= $gene_id;
				$optional->{'transcript_id'}	= $transcript_id;
				if (defined $common and defined $optional) {
					$output .= print_annotations($common,$optional);			
				}													 	
			}		
 		}
 	}

	return $output;
	
} # End get_corsair_alt_annotations

=head2 get_spade_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_spade_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_spade_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'spade'}->{'source'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->spade ) {
	 		my ($method) = $analysis->spade;	 		
			# get residue annotations
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {
					if ( defined $region->type_domain and defined $region->bit_score) {
						# get coords
						my ($method_start)  = ( defined $region->start )      ? $region->start      : $region->alignment_start;
						my ($method_end)    = ( defined $region->end )        ? $region->end        : $region->alignment_end;
						my ($method_score)  = ( defined $region->bit_score )  ? $region->bit_score  : '.';
						my ($method_strand) = ( defined $region->strand )     ? $region->strand     : '.';					
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $region->type_domain,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						$optional->{'note'}				= "hmm_name:".$region->hmm_name if ($region->hmm_name);
						$optional->{'note'}				.= ",evalue:".$region->evalue if ($region->evalue);
						if ($region->alignment_start and $region->alignment_end) {
							$optional->{'note'}			.= ",pep_start:".$region->alignment_start.",pep_end:".$region->alignment_end;
						}						
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_spade_annotations

=head2 get_inertia_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_inertia_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_inertia_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = '.';
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'inertia'}->{'source'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->inertia ) {
	 		my ($method) = $analysis->inertia;
	 		# get residue annotations
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {
					if ( defined $region->unusual_evolution and 
						 defined $region->start and defined $region->end and defined $region->strand ) {
						# get type
						my ($method_type) = '.';
						if ( $region->unusual_evolution eq $APPRIS::Utils::Constant::UNKNOWN_LABEL ) {
							$method_type = $GTF_CONSTANTS->{'inertia'}->{'annot'}->[0];
						}
						elsif ( $region->unusual_evolution eq $APPRIS::Utils::Constant::NO_LABEL ) {
							$method_type = $GTF_CONSTANTS->{'inertia'}->{'annot'}->[1];	
						}
						else {
							$method_type = $GTF_CONSTANTS->{'inertia'}->{'annot'}->[0];
						}							 	
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'start'		=> $region->start,
								'end'		=> $region->end,
								'score'		=> $method_score,
								'strand'	=> $region->strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_inertia_annotations

=head2 get_thump_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_thump_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_thump_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = '.';
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'thump'}->{'source'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->thump ) {
	 		my ($method) = $analysis->thump;
	 		# get residue annotations
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {
					#if ( defined $region->start and defined $region->end and defined $region->strand ) {
						# get coords
						my ($method_start)  = ( defined $region->start )      ? $region->start      : $region->pstart;
						my ($method_end)    = ( defined $region->end )        ? $region->end        : $region->pend;
						my ($method_strand) = ( defined $region->strand )     ? $region->strand     : '.';					
						# get type
						my ($method_type) = '.';
						if ( $region->damaged and ($region->damaged eq '0') ) {
							$method_type = $GTF_CONSTANTS->{'thump'}->{'annot'}->[0];
						}
						elsif ( $region->damaged and ($region->damaged eq '1') ) {
							$method_type = $GTF_CONSTANTS->{'thump'}->{'annot'}->[1];	
						}
						else {
							$method_type = $GTF_CONSTANTS->{'thump'}->{'annot'}->[0];
						}
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						if ($region->pstart and $region->pend) {
							$optional->{'note'}			= "pep_start:".$region->pstart.",pep_end:".$region->pend;
						}						
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					#}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_thump_annotations

=head2 get_crash_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_crash_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_crash_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = '.';
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'crash'}->{'source'};

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->crash ) {
	 		my ($method) = $analysis->crash;
	 		# get residue annotations
			if ( defined $method->regions ) {
				foreach my $region (@{$method->regions}) {
					if ( defined $region->s_mean and defined $region->s_prob and 
						 defined $region->d_score and defined $region->c_max and
						 defined $region->reliability and defined $region->localization and
						 defined $region->pstart and defined $region->pend ) {
						# get type
						my ($method_type) = '.';
						if ( $method->peptide_signal eq $APPRIS::Utils::Constant::UNKNOWN_LABEL ) {
							$method_type = $GTF_CONSTANTS->{'crash'}->{'annot'}->{'signal_peptide'}->[1];
						}
						elsif ( $method->peptide_signal eq $APPRIS::Utils::Constant::NO_LABEL ) {
							$method_type = $GTF_CONSTANTS->{'crash'}->{'annot'}->{'signal_peptide'}->[2];
						}
						elsif ( $method->peptide_signal eq $APPRIS::Utils::Constant::OK_LABEL ) {
							$method_type = $GTF_CONSTANTS->{'crash'}->{'annot'}->{'signal_peptide'}->[0];
						}
						# get coords
						my ($method_start)  = ( defined $region->start )      ? $region->start      : $region->pstart;
						my ($method_end)    = ( defined $region->end )        ? $region->end        : $region->pend;
						my ($method_strand) = ( defined $region->strand )     ? $region->strand     : '.';					
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						$optional->{'note'}				= "s_mean:".$region->s_mean.",".
										  				  "s_prob:".$region->s_prob.",".
														  "d_score:".$region->d_score.",".
														  "c_max:".$region->c_max."";
						if ($region->pstart and $region->pend) {
							$optional->{'note'}			.= ",pep_start:".$region->pstart.",pep_end:".$region->pend;
						}														  
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}
						# get type
						my ($method_type2) = '.';
						if ( $method->mitochondrial_signal eq $APPRIS::Utils::Constant::UNKNOWN_LABEL ) {
							$method_type2 = $GTF_CONSTANTS->{'crash'}->{'annot'}->{'mitochondrial_signal'}->[1];
						}
						elsif ( $method->mitochondrial_signal eq $APPRIS::Utils::Constant::NO_LABEL ) {
							$method_type2 = $GTF_CONSTANTS->{'crash'}->{'annot'}->{'mitochondrial_signal'}->[2];
						}
						elsif ( $method->mitochondrial_signal eq $APPRIS::Utils::Constant::OK_LABEL ) {
							$method_type2 = $GTF_CONSTANTS->{'crash'}->{'annot'}->{'mitochondrial_signal'}->[0];
						}
						# common attributes
						my ($common2) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type2,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional2);
						$optional2->{'gene_id'}			= $gene_id;
						$optional2->{'transcript_id'}	= $transcript_id;
						$optional2->{'note'}			= "reliability:".$region->reliability.",".
														  "localization:".$region->localization."";
						if (defined $common2 and defined $optional2) {
							$output .= print_annotations($common2,$optional2);			
						}					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_crash_annotations


=head2 get_proteo_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : String - $version
  Example    : $annot = get_proteo_annotations($trans_id, $gen_id, $ext_id, $feat, $v);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut
sub get_proteo_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;

    my ($output) = '';
    my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
	my ($method_score) = '.';
	my ($method_phase) = '.';
	my ($method_source) = $GTF_CONSTANTS->{'proteo'}->{'source'};
	my ($method_type) = $GTF_CONSTANTS->{'proteo'}->{'type'};
	
	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->proteo ) {
	 		my ($method) = $analysis->proteo;
	 		# get residue annotations
			if ( defined $method->peptides ) {
				foreach my $region (@{$method->peptides}) {
					if ( defined $region->sequence and defined $region->num_experiments ) {
						# get coords
						my ($method_start)  = ( defined $region->start )      ? $region->start      : $region->pstart;
						my ($method_end)    = ( defined $region->end )        ? $region->end        : $region->pend;
						my ($method_strand) = ( defined $region->strand )     ? $region->strand     : '.';
						# common attributes
						my ($common) = {
								'seqname'	=> $method_seqname,
								'source'	=> $method_source,
								'type'		=> $method_type,
								'start'		=> $method_start,
								'end'		=> $method_end,
								'score'		=> $method_score,
								'strand'	=> $method_strand,
								'phase'		=> $method_phase
						};
						# optinal attributes
						my($optional);
						$optional->{'gene_id'}			= $gene_id;
						$optional->{'transcript_id'}	= $transcript_id;
						if ($region->pstart and $region->pend and $region->sequence) {
							$optional->{'note'}			= "pep_start:".$region->pstart.",pep_end:".$region->pend.",pep_seq:".$region->sequence;
						}
						if (defined $common and defined $optional) {
							$output .= print_annotations($common,$optional);			
						}													 	
					}
				}
			} 		
 		}
 	}

	return $output;
	
} # End get_proteo_annotations

=head2 get_trifid_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [2]    : String - the stable identifier of gene
  Arg [3]    : String - the external database name associated with transcript
  Arg [4]    : APPRIS::Transcript
  Example    : $annot = get_trifid_annotations($trans_id, $gen_id, $ext_id, $feat);
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut
sub get_trifid_annotations {
	my ($transcript_id, $gene_id, $external_id, $feature) = @_;
	my ($output) = '';

	if ( $feature->analysis ) {
		my ($analysis) = $feature->analysis;
		if ( $analysis->trifid ) {
			if ( defined($analysis->trifid->norm_trifid_score) ) {
				my ($norm_trifid_score) = $analysis->trifid->norm_trifid_score;
				my ($trifid_score) = $analysis->trifid->trifid_score;

				my ($method_seqname) = ( $feature->chromosome ) ? $feature->chromosome : $gene_id;
				my ($method_source) = $GTF_CONSTANTS->{'trifid'}->{'source'};
				my ($method_type) = $GTF_CONSTANTS->{'trifid'}->{'type'};
				my ($method_start) = $feature->start;
				my ($method_end) = $feature->end;
				my ($method_score) = $norm_trifid_score;
				my ($method_strand) = $feature->strand;
				my ($method_phase) = '.';

				# Common attributes
				my ($common) = {
						'seqname'	=> $method_seqname,
						'source'	=> $method_source,
						'type'		=> $method_type,
						'start'		=> $method_start,
						'end'		=> $method_end,
						'score'		=> $method_score,
						'strand'	=> $method_strand,
						'phase'		=> $method_phase
				};

				# Optional attributes
				my($optional);
				$optional->{'gene_id'}			= $gene_id;
				$optional->{'transcript_id'}	= $transcript_id;
				$optional->{'transcript_name'}	= $external_id;
				$optional->{'note'}				= "trifid_score:$trifid_score";
				if (defined $common and defined $optional) {
					$output .= print_annotations($common,$optional);
				}
			}
		}
	}

	return $output;

} # End get_trifid_annotations

1;
