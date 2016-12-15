=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::RES - Utility functions for info exporting

=head1 SYNOPSIS

  use APPRIS::Export::RES
    qw(
       get_trans_annotations
     );

  or to get all methods just

  use APPRIS::Export::RES;

  eval { get_trans_annotations($feature,$params) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Retrieves sequences of transcript as fasta format.

=head1 METHODS

=cut

package APPRIS::Export::RES;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

use APPRIS::Utils::Exception qw(throw warning deprecate);

use APPRIS::Utils::Constant qw(
	$OK_LABEL
	$METHOD_DESC
);

###################
# Global variable #
###################
use vars qw(
	$OK_LABEL
	$METHOD_DESC
	$METHOD_HEADS
	$METHOD_LABEL_DESC
);

$OK_LABEL = $APPRIS::Utils::Constant::OK_LABEL;
$METHOD_DESC = $APPRIS::Utils::Constant::METHOD_DESC;
$METHOD_LABEL_DESC = $APPRIS::Utils::Constant::METHOD_LABEL_DESC;
$METHOD_HEADS = {
	'appris'	=> "principal isoform (reliability, 1-->5)",
	'firestar'	=> "PDB ligand (reliability, 1-->6)",
	'matador3d'	=> "best PDB template (%ID)",
	'spade'		=> "best Pfam domain name (e-value)",
	'corsair'	=> "nearest homologue (%ID)",
	'crash'		=> "type signal",
	'thump'		=> "type signal",
	'inertia'	=> "slr_omega_score",
	'proteo'	=> "peptide (no. experiments found)",
};

=head2 get_trans_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : String - $soure List of sources  
  Arg [3]    : Int    - $res Residue position
  Example    : $annot = $exporter->get_trans_annotations($feature,'aa');
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_trans_annotations {
    my ($feature, $source, $inres) = @_;
    my ($report);
	
    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {    	
		my ($id);
		if ($feature->stable_id) {
			$id = $feature->stable_id;
			my ($len) = 0;
			if ( $feature->translate and $feature->translate->sequence ) {
				$len = length($feature->translate->sequence);
			}
			
		 	if ( $feature->analysis ) {		 		
		 		my ($methods);
		 		my ($analysis) = $feature->analysis;
		 		if ( ($source =~ /firestar/) or ($source eq 'all') ) {
			 		if ( $analysis->firestar and $analysis->firestar->result ) {		 			
						my ($method) = $analysis->firestar;	 			
				 		my ($residues) = parser_firestar_residues($method, $inres);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'firestar',
								'id'		=> $METHOD_DESC->{'firestar'},
								'label'		=> $METHOD_LABEL_DESC->{'firestar'},							
								'title'		=> $METHOD_HEADS->{'firestar'},
								'residues'	=> $residues
							});								
						}
			 		}
		 		}
		 		if ( ($source =~ /matador3d/) or ($source eq 'all') ) {
			 		if ( $analysis->matador3d and $analysis->matador3d->result ) {		 			
						my ($method) = $analysis->matador3d;	 			
				 		my ($residues) = parser_matador3d_residues($method, $inres);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'matador3d',
								'id'		=> $METHOD_DESC->{'matador3d'},
								'label'		=> $METHOD_LABEL_DESC->{'matador3d'},
								'title'		=> $METHOD_HEADS->{'matador3d'},
								'residues'	=> $residues
							});
						}
			 		}		 			
		 		}
		 		if ( ($source =~ /spade/) or ($source eq 'all') ) {	 		
			 		if ( $analysis->spade and $analysis->spade->result ) {		 			
						my ($method) = $analysis->spade;
				 		my ($residues) = parser_spade_residues($method, $inres);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'spade',
								'id'		=> $METHOD_DESC->{'spade'},
								'label'		=> $METHOD_LABEL_DESC->{'spade'}->[2],
								'title'		=> $METHOD_HEADS->{'spade'},
								'residues'	=> $residues
							});
						}
			 		}
		 		}
		 		if ( ($source =~ /corsair/) or ($source eq 'all') ) {
			 		if ( $analysis->corsair and $analysis->corsair->result ) {		 			
						my ($method) = $analysis->corsair;
				 		my ($residues) = parser_corsair_residues($method, $len);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'corsair',
								'id'		=> $METHOD_DESC->{'corsair'},
								'label'		=> $METHOD_LABEL_DESC->{'corsair'},
								'title'		=> $METHOD_HEADS->{'corsair'},
								'residues'	=> $residues
							});
						}
			 		}
		 		}
		 		if ( ($source =~ /thump/) or ($source eq 'all') ) {
			 		if ( $analysis->thump and $analysis->thump->result ) {
			 			my ($method) = $analysis->thump;		 			
				 		my ($residues) = parser_thump_residues($method, $inres);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'thump',
								'id'		=> $METHOD_DESC->{'thump'},
								'label'		=> $METHOD_LABEL_DESC->{'thump'}->[0],
								'title'		=> $METHOD_HEADS->{'thump'},
								'residues'	=> $residues
							});
						}
			 		}
		 		}
		 		if ( ($source =~ /crash/) or ($source eq 'all') ) {		 		
			 		if ( $analysis->crash and $analysis->crash->result ) {		 			
						my ($method) = $analysis->crash;	 			
				 		my ($label,$residues) = parser_crash_residues($method, $inres);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'crash',
								'id'		=> $METHOD_DESC->{'crash'},
								'label'		=> $METHOD_LABEL_DESC->{'crash'},
								'title'		=> $METHOD_HEADS->{'crash'},
								'residues'	=> $residues
							});
						}
			 		}
		 		}
		 		if ( ($source =~ /proteo/) ) {
			 		if ( $analysis->proteo and $analysis->proteo->result ) {
			 			my ($method) = $analysis->proteo;		 			
				 		my ($residues) = parser_proteo_residues($method, $len);
						if ( defined $residues and scalar($residues) > 0 ) {
							push(@{$methods}, {
								'name'		=> 'proteo',
								'id'		=> $METHOD_DESC->{'proteo'},
								'label'		=> $METHOD_LABEL_DESC->{'proteo'},
								'title'		=> $METHOD_HEADS->{'proteo'},
								'residues'	=> $residues
							});
						}
			 		}
		 		}
		 		if ( defined $methods ) {
					$report = {
						'id'		=> $id,
						'methods'	=> $methods
					};		 			
		 		}
		 	}
		}		
    }
    else {
		throw('Argument must be an APPRIS::Transcript');
   	}
	return $report;
}

sub parser_firestar_residues {
	my ($method, $inres) = @_;
	my ($residues);
	if ( defined $method->residues and defined $method->result ) {

		foreach my $region (@{$method->residues}) {
			if ( defined $region->residue and defined $region->ligands and $region->domain ) {
				my ($dom) = '';
				my ($lig) = '';
				my ($sc) = '';
				my (@ligands) = split(/\|/, $region->ligands);
				foreach my $ligs (@ligands) {
					if ( $ligs ne '' ) {
						if ( $ligs =~ /^([^\[]*)\[/ ) {
							$lig .= $1.',';
						}
						$ligs =~ s/^[^\[]*\[[^\,]*\,//;
						$ligs =~ s/\,[^\]]*\]$//;					
						$sc .= $ligs.',';						
					}
				}
				$lig =~ s/\,$//;
				$sc =~ s/\,$//;
				if ( $region->domain ) {
					$dom = $region->domain;
					$dom =~ s/^[\w|\-]{6}//;
					if (  length($dom) > 6 ) {
						$dom =~ s/[\w|\-]{6}$//;
					}
					else {
						my ($s) = $dom;
						$dom = substr $s, 0, 1;
					}
				}				
				if ( !defined $inres or ( defined $inres and ($region->residue == $inres) ) ) {
					my ($res) = {
						'start'		=> $region->residue,
#						'end'		=> $region->residue,
					};
					if ( ($dom ne '') and ($lig ne '') and ($sc ne '') ) {
						#$res->{'annot'} = $dom . "\t" . $lig." (".$sc.")";
						$res->{'annot'} = $lig." (".$sc.")";
					}
					push(@{$residues},$res);					
				}
			}
		}
	}
	return $residues;
		
} # end parser_firestar_residues

#sub parser_matador3d_residues {
#	my ($method, $inres) = @_;
#	my ($residues);
#	
#	if ( defined $method->result ) {
#		if ( defined $method->alignments ) {
#			foreach my $region (@{$method->alignments}) {	
#				if ( ($region->type eq 'mini-exon') and 
#					defined $region->score and 
#					defined $region->pstart and defined $region->pend and 
#					defined $region->pdb_id and defined $region->identity ) {
#						if ( !defined $inres or ( defined $inres and ($region->pstart <= $inres) and ($region->pend >= $inres) ) ) {
#							my ($res) = {
#								'start'		=> $region->pstart,
#								'end'		=> $region->pend,
#								'annot'		=> $region->pdb_id ." (" . $region->identity . ")",
#							};
#							push(@{$residues},$res);
#						}
#				}
#			}
#		}
#	}
#	
#	return $residues;
#	
#} # end parser_matador3d_residues
sub parser_matador3d_residues {
	my ($method, $inres) = @_;
	my ($residues);
	
	if ( defined $method->result ) {
		if ( defined $method->alignments ) {
			foreach my $region (@{$method->alignments}) {	
				if ( 
					defined $region->pstart and defined $region->pend and 
					defined $region->pdb_id and defined $region->score
				){
					if ( !defined $inres or ( defined $inres and ($region->pstart <= $inres) and ($region->pend >= $inres) ) ) {
						my ($res) = {
							'start'		=> $region->pstart,
							'end'		=> $region->pend,
							'annot'		=> $region->pdb_id ." (" . $region->score . ")",
						};
						push(@{$residues},$res);
					}
				}
			}
		}
	}
	
	return $residues;
	
} # end parser_matador3d_residues

sub parser_spade_residues {
	my ($method, $inres) = @_;
	my ($residues);
	
	if ( defined $method->result ) {
		if ( defined $method->regions ) {
			foreach my $region (@{$method->regions}) {	
				if ( defined $region->alignment_start and defined $region->alignment_end and defined $region->hmm_name and defined $region->evalue and defined $region->type_domain ) {
					if ( !defined $inres or ( defined $inres and ($region->alignment_start <= $inres) and ($region->alignment_end >= $inres) ) ) {
						my ($res) = {
							'type'		=> $region->type_domain,
							'start'		=> $region->alignment_start,
							'end'		=> $region->alignment_end,
							'annot'		=> $region->hmm_name . " (" . $region->evalue . ")",
						};						
					 	if ( ($region->type_domain eq 'domain_damaged') or ($region->type_domain eq 'domain_wrong') ) {
					 		$res->{'damaged'} = 1;
					 	}
					 	push(@{$residues},$res);
					}
				}
			}
		}
	}
	
	return $residues;
	
} # end parser_spade_residues

sub parser_corsair_residues {
	my ($method, $len) = @_;
	my ($residues);
	
	if ( defined $method->result ) {
		my ($s_report) = $method->result;
		#>ENST00000511833        0.5
		#Homo sapiens    100.00  0.5
		#        - 472575:472608[1:12]   0.5
		#                Homo sapiens    100.00  0.5
		#        - 441698:444432[13:923] 0.5
		#                Homo sapiens    100.00  0.5		
		$s_report =~ s/^\>[^\n]*\n//;
		my ($annot) = '';
		while ( $s_report =~ /([^\n]*)\n/g ) {
			my ($line) = $1;
			my (@cols) = split(/\t/,$line);
			if ( ($cols[0] ne '') and ($cols[1] ne '') and ($cols[2] ne '0') ) {
				my ($sp) = $cols[0];
				my ($id) = $cols[1];
				$annot .=	$sp . " (" . $id . ")" . "\n";
			}
			elsif ( ($cols[0] eq '') and ($cols[1] ne '') and ($cols[2] ne '0') ) {
				#unless ( $cols[1] =~ /^\s*\-/ ) {
				#	my ($id) = sprintf ("%.2f",$cols[1]);
				#	$annot .=	$sp . " (" . $id . ")";
				#}				
			}
			elsif ( ($cols[0] eq '') and ($cols[1] eq '') and ($cols[2] ne '0') ) {
				
			}
		}
		if ( $annot ne '' ) {
			$annot =~ s/\n$//;
			my ($res) = {
				'start'		=> "1",
				'end'		=> "$len",							
				'annot'		=> $annot,
			};						
		 	push(@{$residues},$res);			
		}
	}
	
	return $residues;
	
} # end parser_corsair_residues

sub parser_crash_residues {
	my ($method, $inres) = @_;
	my ($residues);
	my ($type);
	
	if ( defined $method->result and defined $method->peptide_signal and defined $method->mitochondrial_signal ) {
		if ( ($method->peptide_signal eq $OK_LABEL) or ($method->mitochondrial_signal eq $OK_LABEL) ) {
			if ( defined $method->regions ) {
				if ( ($method->peptide_signal eq $OK_LABEL) and ($method->mitochondrial_signal eq $OK_LABEL) ) {
					$type = 'peptide_mitochondrial_signal';
				}
				elsif ( $method->peptide_signal eq $OK_LABEL ) {
					$type = 'peptide_signal';
				}					
				elsif ( $method->mitochondrial_signal eq $OK_LABEL ) {
					$type = 'mitochondrial_signal';
				}					
				
				foreach my $region (@{$method->regions}) {
					if ( defined $region->pstart and defined $region->pend ) {
						if ( !defined $inres or ( defined $inres and ($region->pstart <= $inres) and ($region->pend >= $inres) ) ) {
							my ($res) = {
								'type'	=> $type,
								'start'	=> $region->pstart,
								'end'	=> $region->pend,
								'annot'	=> $type,
							};
							push(@{$residues},$res);
						}
					}
				}
			}						
		}
	}
	
	return $residues;
	
} # end parser_crash_residues

sub parser_thump_residues {
	my ($method, $inres) = @_;
	my ($residues);
	
	if ( defined $method->result ) {
		if ( defined $method->regions ) {
			foreach my $region (@{$method->regions}) {
				if ( defined $region->pstart and defined $region->pend ) {
					if ( !defined $inres or ( defined $inres and ($region->pstart <= $inres) and ($region->pend >= $inres) ) ) {			
						my ($res) = {
							'start'	=> $region->pstart,
							'end'	=> $region->pend,
							#'start'		=> "1",
							#'end'		=> "$len",
							'annot'		=> "transmembrane signal",			
						};
					 	if ( defined $region->damaged ) {
					 		$res->{'damaged'} = 1;
					 		$res->{'type'} = 'transmembrane_wrong';
					 	}
					 	else {
					 		$res->{'type'} = 'transmembrane';				 		
					 	}
					 	push(@{$residues},$res);
					}
				}
			}
		}
	}
		
	return $residues;
	
} # end parser_thump_residues

sub parser_appris_residues {
	my ($method, $len) = @_;
	my ($residues);
	
	if ( defined $method->principal_isoform_signal and defined $method->result ) {
		if ( defined $method->reliability ) {
			my ($annot) = $method->reliability;
			if ( $annot ne '' ) {
				my ($res) = {
					'start'		=> "1",
					'end'		=> "$len",							
					'annot'		=> $annot,
				};						
			 	push(@{$residues},$res);			
			}			
		}
	}
		
	return $residues;
	
} # end parser_appris_residues

sub parser_proteo_residues {
	my ($method, $inres) = @_;
	my ($residues);
	
	if ( defined $method->result ) {
		if ( defined $method->peptides ) {
			foreach my $region (@{$method->peptides}) {	
				if ( defined $region->sequence and defined $region->num_experiments and defined $region->pstart and defined $region->pend ) {
					my ($res) = {
						'start'		=> $region->pstart,
						'end'		=> $region->pend,
						'annot'		=> $region->sequence . " (".$region->num_experiments.")",
					};
					push(@{$residues},$res);
				}
			}
		}
	}
	
	return $residues;
	
} # end parser_proteo_residues

1;
