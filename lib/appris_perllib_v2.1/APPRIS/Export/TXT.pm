=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::TXT - Utility functions for info exporting

=head1 SYNOPSIS

  use APPRIS::Export::TXT
    qw(
       get_tsv_annotations
     );

  or to get all methods just

  use APPRIS::Export::TXT;

  eval { get_tsv_annotations($feature) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Retrieves information of transcript as TXT format.

=head1 METHODS

=cut

package APPRIS::Export::TXT;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::Constant qw(
        $OK_LABEL
);

###################
# Global variable #
###################
use vars qw(
	$METHOD_HEADS
	$OK_LABEL
	$UNKNOWN_LABEL
	$NO_LABEL
	$APPRIS_PRINC_LABEL
	$APPRIS_CANDI_LABEL
	$APPRIS_ALTER_LABEL
);

$METHOD_HEADS = {
	'appris'=>{
		'tsv'=>
			"== principal_isoform: appris_isoform\tappris_reliability\n",
		'raw'=>
			"== principal_isoform: text_result\n",
	},
	'firestar'=>{
		'tsv'=>
			"== functional_residue: residue\tamino_acid\tligand\treliability_score (1-->6)\n",
		'raw'=>
			"== functional_residue: text_result\n",
	},
	'matador3d'=>{
		'tsv'=>
			"== homologous_structure: residues\tbest_template\t%ID\n",
		'raw'=>
			"== homologous_structure: text_result\n",
	},
	'matador3d2'=>{
		'tsv'=>
			"== homologous_structure: residues\tbest_template\tbitscore\n",
		'raw'=>
			"== homologous_structure: text_result\n",
	},
	'spade'=>{
		'tsv'=>
			"== functional_domain: start\tend\tdomain_acc\tdomain_name\tbest_e-value\n",
		'raw'=>
			"== functional_domain: text_result\n",
	},
	'corsair'=>{ 
		'tsv'=>
			"== vertebrate_conservation: nearest_homologue\t%ID\n",
		'raw'=>
			"== vertebrate_conservation: text_result\n",
	},
	'corsair_alt'=>{ 
		'tsv'=>
			"== vertebrate_conservation: nearest_homologue\t%ID\n",
		'raw'=>
			"== vertebrate_conservation: text_result\n",
	},
	'crash'=>{
		'tsv'=>
			"== signal_peptide_mitochondrial: type_signal\tstart\tend\n",
		'raw'=>
			"== signal_peptide_mitochondrial: text_result\n",
	},
	'thump'=>{
		'tsv'=>
			"== transmembrane_signal: helix_start\thelix_end\tdamaged\n",	
		'raw'=>
			"== transmembrane_signal: text_result\n",
	},
	'proteo'=>{
		'tsv'=>
			"== proteomic_evidence: peptides\tsequence\tno.experiments_found\n",
		'raw'=>
			"== proteomic_evidence: text_result\n",
	},
	'inertia'=>{
		'tsv'=>
			"== INERTIA: slr_omega_score\texon_start\texon_end\n",
		'raw'=>
			"== INERTIA: text_result\n",
	}
};
$OK_LABEL					= $APPRIS::Utils::Constant::OK_LABEL;
$UNKNOWN_LABEL				= $APPRIS::Utils::Constant::UNKNOWN_LABEL;
$NO_LABEL					= $APPRIS::Utils::Constant::NO_LABEL;
$APPRIS_PRINC_LABEL			= $APPRIS::Utils::Constant::APPRIS_PRINC_LABEL;
$APPRIS_CANDI_LABEL			= $APPRIS::Utils::Constant::APPRIS_CANDI_LABEL;
$APPRIS_ALTER_LABEL			= $APPRIS::Utils::Constant::APPRIS_ALTER_LABEL;


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
    my ($feature, $source, $res) = @_;
    my ($output) = '';

    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {
   	    
		if ($feature->stable_id) {
			if ($feature->translate and $feature->translate->sequence) { # proteing coding methods
				$output .= ">".$feature->stable_id."\n";
		 		my (%sc) = map { $_ => 1 } split(',', $source);				
				if ( (exists $sc{firestar}) or ($source eq 'all') ) {
					$output .= get_trans_firestar_tsv_annot($feature, $res);
				}
				if ( (exists $sc{matador3d}) or ($source eq 'all') ) {
					$output .= get_trans_matador3d_tsv_annot($feature, $res);
				}
				if ( (exists $sc{matador3d2}) or ($source eq 'all') ) {
					$output .= get_trans_matador3d2_tsv_annot($feature, $res);
				}
				if ( (exists $sc{corsair}) or ($source eq 'all') ) {
					$output .= get_trans_corsair_tsv_annot($feature, $res);
				}
				if ( (exists $sc{corsair_alt}) or ($source eq 'all') ) {
					$output .= get_trans_corsair_alt_tsv_annot($feature, $res);
				}
				if ( (exists $sc{spade}) or ($source eq 'all') ) {
					$output .= get_trans_spade_tsv_annot($feature, $res);
				}
				if ( (exists $sc{inertia}) or ($source eq 'all') ) {
					$output .= get_trans_inertia_tsv_annot($feature, $res);
				}
				if ( (exists $sc{thump}) or ($source eq 'all') ) {
					$output .= get_trans_thump_tsv_annot($feature, $res);
				}
				if ( (exists $sc{crash}) or ($source eq 'all') ) {
					$output .= get_trans_crash_tsv_annot($feature, $res);
				}
				if ( (exists $sc{proteo}) ) {
					$output .= get_trans_proteo_tsv_annot($feature,$res);
				}
#				if ( ($source =~ /appris/) or ($source eq 'all') ) {
#					$output .= get_trans_appris_tsv_annot($feature,$res);
#					#if ( $type eq 'tsv' ) {
#					#	$output .= get_trans_appris_tsv_annot($feature);
#					#}
#					#elsif ( $type eq 'raw' ) {
#					#	$output .= get_trans_appris_raw_annot($feature);
#					#}
#				}
			}
		}
    }
    else {
		throw('Argument must be an APPRIS::Transcript');
   	}
	return $output;
}

=head2 get_tsv_annotations

  Arg [1]    : APPRIS::Transcript or undef
  Arg [2]    : String - $methods 
               List of sources
  Example    : $annot = get_tsv_annotations($feature,$method);
  Description: Retrieves tabular information of transcript.
  Returntype : String or undef

=cut

sub get_tsv_annotations {
    my ($feature, $methods) = @_;
    my ($type) = 'tsv';
    my ($output) = '';

    if ( defined $feature ) {
    	foreach my $method ( split(',', $methods) ) {
			if ( ($method eq 'firestar') or ($method eq 'all') ) {
				$output .= get_firestar_annot($method, $type, $feature);
			}
			if ( ($method eq 'matador3d') or ($method eq 'all') ) {
				$output .= get_matador3d_annot($method, $type, $feature);
			}
			if ( ($method eq 'matador3d2') or ($method eq 'all') ) {
				$output .= get_matador3d2_annot($method, $type, $feature);
			}
			if ( ($method eq 'corsair') or ($method eq 'all') ) {
				$output .= get_corsair_annot($method, $type, $feature);
			}
			if ( ($method eq 'spade') or ($method eq 'all') ) {
				$output .= get_spade_annot($method, $type, $feature);
			}
			if ( ($method eq 'thump') or ($method eq 'all') ) {
				$output .= get_thump_annot($method, $type, $feature);
			}
			if ( ($method eq 'crash') or ($method eq 'all') ) {
				$output .= get_crash_annot($method, $type, $feature);
			}
			if ( ($method eq 'inertia') or ($method eq 'all') ) {
				$output .= get_inertia_annot($method, $type, $feature);
			}
			if ( ($method eq 'appris') or ($method eq 'all') ) {
				$output .= get_appris_annot($method, $type, $feature);
			}
    	}
    }
    else {
		throw('Argument must be defined');
   	}
	return $output;	
}

=head2 get_raw_annotations

  Arg [1]    : APPRIS::Transcript or undef
  Arg [2]    : String - $methods 
               List of sources
  Example    : $annot = get_raw_annotations($feature,$method);
  Description: Retrieves tabular information of transcript.
  Returntype : String or undef

=cut

sub get_raw_annotations {
    my ($feature, $methods) = @_;
	my ($type) = 'raw';
    my ($output) = '';
    
    if ( defined $feature ) {
    	foreach my $method ( split(',', $methods) ) {
			if ( ($method eq 'firestar') or ($method eq 'all') ) {
				$output .= get_firestar_annot($method, $type, $feature);
			}
			if ( ($method eq 'matador3d') or ($method eq 'all') ) {
				$output .= get_matador3d_annot($method, $type, $feature);
			}
			if ( ($method eq 'matador3d2') or ($method eq 'all') ) {
				$output .= get_matador3d2_annot($method, $type, $feature);
			}
			if ( ($method eq 'corsair') or ($method eq 'all') ) {
				$output .= get_corsair_annot($method, $type, $feature);
			}
			if ( ($method eq 'corsair_alt') or ($method eq 'all') ) {
				$output .= get_corsair_alt_annot($method, $type, $feature);
			}
			if ( ($method eq 'spade') or ($method eq 'all') ) {
				$output .= get_spade_annot($method, $type, $feature);
			}
			if ( ($method eq 'thump') or ($method eq 'all') ) {
				$output .= get_thump_annot($method, $type, $feature);
			}
			if ( ($method eq 'crash') or ($method eq 'all') ) {
				$output .= get_crash_annot($method, $type, $feature);
			}
			if ( ($method eq 'inertia') or ($method eq 'all') ) {
				$output .= get_inertia_annot($method, $type, $feature);
			}
			if ( ($method eq 'appris') or ($method eq 'all') ) {
				$output .= get_appris_annot($method, $type, $feature);
			}
			if ( ($method eq 'proteo') ) {
				$output .= get_proteo_annot($method, $type, $feature);
			}
    	}
    }
    else {
		throw('Argument must be defined');
   	}
	return $output;	
}

=head2 get_tsv_results

  Arg [1]    : Hash report or undef
  Arg [2]    : String - $methods 
               List of sources  
  Example    : $annot = get_results($report);
  Description: Retrieves tabular information of transcript.
  Returntype : String or undef

=cut

sub get_tsv_results {
    my ($report,$methods) = @_;
    my ($output) = '';
    if ( defined $report ) {
    	foreach my $method ( split(',', $methods) ) {
    		if ( exists $report->{$method} ) {
    			my ($m_report) = $report->{$method};
				if ( ($method eq 'firestar') or ($method eq 'all') ) {
					$output .= get_firestar_tsv_result($method, $m_report);
				}
				if ( ($method eq 'matador3d') or ($method eq 'all') ) {
					$output .= get_matador3d_tsv_result($method, $m_report);
				}
				if ( ($method eq 'corsair') or ($method eq 'all') ) {
					$output .= get_corsair_tsv_result($method, $m_report);
				}
				if ( ($method eq 'corsair_alt') or ($method eq 'all') ) {
					$output .= get_corsair_alt_tsv_result($method, $m_report);
				}
				if ( ($method eq 'spade') or ($method eq 'all') ) {
					$output .= get_spade_tsv_result($method, $m_report);
				}
				if ( ($method eq 'thump') or ($method eq 'all') ) {
					$output .= get_thump_tsv_result($method, $m_report);
				}
				if ( ($method eq 'crash') or ($method eq 'all') ) {
					$output .= get_crash_tsv_result($method, $m_report);
				}
				if ( ($method eq 'appris') or ($method eq 'all') ) {
					$output .= get_appris_tsv_result($method, $m_report);
				}				
    		}
    	}    	
    }
    else {
		throw('Argument must be defined');
   	}
	return $output;
}

=head2 get_raw_results

  Arg [1]    : Hash report or undef
  Arg [2]    : String - $methods 
               List of sources  
  Example    : $annot = get_results($report);
  Description: Retrieves tabular information of transcript.
  Returntype : String or undef

=cut

sub get_raw_results {
    my ($report,$methods) = @_;
    my ($output) = '';
    if ( defined $report ) {
    	foreach my $method ( split(',', $methods) ) {
    		if ( exists $report->{$method} ) {
				$output .= "== APPRIS > ".uc($method)."\n";    			
    			my ($m_report) = $report->{$method};
				while (my ($seq_id, $seq_report) = each(%{$m_report}) ) {
					if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
						$output .= $seq_report->{'result'}."\n"; 
					}
				}
    		}
    	}
    }
    else {
		throw('Argument must be defined');
   	}
	return $output;	
}

sub res_is_inside {
	my ($res,$start,$end) = @_;
	my ($inside);
	
	if ( defined $res ) {
		my ($ires); map { $ires->{$_} = 1 } split(',',$res);
		if ( defined $start and defined $end ) {
			foreach my $i ( keys(%{$ires}) ) {
				if ( ($start <= $i) and ($end >= $i) ) {
					$inside = 1;
					last;
				}
			}
		}
		elsif  ( defined $start ) {
			$inside = 1 if ( exists $ires->{$start} ); 
		}		
	}
	else {
		$inside = 1;
	}
	return $inside;
}


# FIRESTAR METHODS -------
sub get_firestar_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_firestar_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_firestar_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_firestar_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_firestar_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_firestar_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_firestar_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_firestar_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_firestar_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;		
}
# firestar: annotation input. tabular output
sub get_trans_firestar_tsv_annot {
	my ($transcript, $res) = @_;
	my ($imet) = 'firestar';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->firestar and $analysis->firestar->result ) {
					my ($method) = $analysis->firestar;
					if ( defined $method->residues and defined $method->result ) {
						foreach my $region (@{$method->residues}) {
							if ( defined $region->residue ) {
								next unless ( res_is_inside($res, $region->residue) );
								my (%parameters) = (
												peptide_position	=> $region->residue,
								);
								if ( $region->domain ) {
									$parameters{domain} = $region->domain;
									$parameters{domain} =~ s/^[\w|\-]{6}//;
									if (  length($parameters{domain}) > 6 ) {
										$parameters{domain} =~ s/[\w|\-]{6}$//;
									}
									else {
										my ($s) = $parameters{domain};
										$parameters{domain} = substr $s, 0, 1;
									}
								}				
								if ( $region->ligands ) {
									my ($lig) = '';
									my ($sc) = '';
									my (@ligands) = split(/\|/, $region->ligands);
									foreach my $ligs (@ligands) {
										if ( $ligs =~ /^([^\[]*)\[/ ) {
											$lig .= $1.',';
										}
										$ligs =~ s/^[^\[]*\[[^\,]*\,//;
										$ligs =~ s/\,[^\]]*\]$//;					
										$sc .= $ligs.',';						
									}
									$lig =~ s/\,$//;
									$sc =~ s/\,$//;
									$parameters{ligands} = $lig if ( $lig ne '' );
									$parameters{rel_score} = $sc if ( $sc ne '' );
									
								}				
								$output .=	$parameters{peptide_position}."\t".
											$parameters{domain}."\t".
											$parameters{ligands}."\t".
											$parameters{rel_score}."\n";
							}
						}
					}
								
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
	return $output;	
}
# firestar: annotation input. raw output
sub get_trans_firestar_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'firestar';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->firestar and $analysis->firestar->result ) {
					my ($rst) = $analysis->firestar->result;
					$rst =~ s/ACCEPT:([^\n]+)\n*//mg;
					$rst =~ s/REJECT:([^\n]+)\n*//mg;
					$output .= $rst."\n";
				}
		 	}
		}
    }
	if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
	return $output;	
}
# firestar: result input. tabular output
sub get_firestar_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'residues'}) ) {
			foreach my $region (@{$seq_report->{'residues'}}) {
				if (defined $region->{'residue'} ) {				
					my ($residue) = $region->{'residue'};
					my (%params) = ( residue => $residue );						
					if ( $region->{'domain'} ) {
						$params{domain} = $region->{'domain'};
						$params{domain} =~ s/^[\w|\-]{6}//;
						if (  length($params{domain}) > 6 ) {
							$params{domain} =~ s/[\w|\-]{6}$//;
						}
						else {
							my ($s) = $params{domain};
							$params{domain} = substr $s, 0, 1;
						}
					}
					if ( $region->{'ligands'} ) {
						my ($lig) = '';
						my ($sc) = '';
						my (@ligands) = split(/\|/, $region->{'ligands'});
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
						$params{ligands} = $lig if ( $lig ne '' );
						$params{rel_score} = $sc if ( $sc ne '' );
						
					}
					$output .=	$params{residue}."\t".
								$params{domain}."\t".
								$params{ligands}."\t".
								$params{rel_score}."\n";
				}
			}
						
		}		
	}
	return $output;
}





# MATADOR3D METHODS -------
sub get_matador3d_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_matador3d_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_matador3d_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_matador3d_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_matador3d_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_matador3d_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_matador3d_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_matador3d_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_matador3d_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# matador3d: annotation input. tabular output
#sub get_trans_matador3d_tsv_annot {
#	my ($transcript,$res) = @_;
#	my ($imet) = 'matador3d';
#	my ($output) = '';
#	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
#		if ($transcript->stable_id) {
#			my ($transcript_id) = $transcript->stable_id;
#			if ( $transcript->analysis ) {
#				my ($analysis) = $transcript->analysis;				
#				if ( $analysis->matador3d and $analysis->matador3d->result ) {
#					my ($method) = $analysis->matador3d;					
#					if ( defined $method->result ) {
#						if ( defined $method->alignments ) {
#							foreach my $region (@{$method->alignments}) {
#								if ( ($region->type eq 'mini-exon') and 
#									defined $region->pstart and defined $region->pend and 
#									defined $region->pdb_id and defined $region->identity ) {
#										next unless ( res_is_inside($res, $region->pstart, $region->pend) );								
#										$output .=	$region->pstart.':'.$region->pend."\t".
#													$region->pdb_id."\t".
#													$region->identity."\n";
#								}
#							}
#						}
#					}
#										
#				}
#		 	}
#		}
#    }
#    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
#	return $output;	
#}
sub get_trans_matador3d_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'matador3d';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->matador3d and $analysis->matador3d->result ) {
					my ($method) = $analysis->matador3d;					
					if ( defined $method->result ) {
						if ( defined $method->alignments ) {
							foreach my $region (@{$method->alignments}) {
								if ( 
									defined $region->pstart and defined $region->pend and 
									defined $region->pdb_id and defined $region->score
								) {
									next unless ( res_is_inside($res, $region->pstart, $region->pend) );								
									$output .=	$region->pstart.':'.$region->pend."\t".
												$region->pdb_id."\t".
												$region->score."\n";
								}
							}
						}
					}
										
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
	return $output;	
}
# matador3d: annotation input. raw output
sub get_trans_matador3d_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'matador3d';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->matador3d and $analysis->matador3d->result ) {
					$output .= $analysis->matador3d->result."\n";
				}
		 	}
		}
    }
	if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }    
	return $output;	
}
# matador3d: result input. tabular output
#sub get_matador3d_tsv_result
#{
#	my ($method, $report) = @_;
#	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
#	
#	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
#		if ( ($seq_id ne 'result') and (exists $seq_report->{'alignments'}) ) {
#			foreach my $region (@{$seq_report->{'alignments'}}) {
#				if ( ($region->{'type'} eq 'mini-exon') and 
#					defined $region->{'pstart'} and defined $region->{'pend'} and 
#					defined $region->{'pdb_id'} and defined $region->{'identity'} ) {
#						$output .=	$region->{'pstart'}.':'.$region->{'pend'}."\t".
#									$region->{'pdb_id'}."\t".
#									$region->{'identity'}."\n";
#				}
#			}
#						
#		}
#	}
#	return $output;
#}
sub get_matador3d_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'alignments'}) ) {
			foreach my $region (@{$seq_report->{'alignments'}}) {
				if ( 
					defined $region->{'pstart'} and defined $region->{'pend'} and 
					defined $region->{'pdb_id'} and defined $region->{'score'}
				) {
					$output .=	$region->{'pstart'}.':'.$region->{'pend'}."\t".
								$region->{'pdb_id'}."\t".
								$region->{'score'}."\n";
				}
			}
		}
	}
	return $output;
}



# MATADOR3D2 METHODS -------
sub get_matador3d2_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_matador3d2_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_matador3d2_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_matador3d2_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_matador3d2_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_matador3d2_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_matador3d2_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_matador3d2_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_matador3d2_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# matador3d2: annotation input. tabular output
sub get_trans_matador3d2_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'matador3d2';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->matador3d2 and $analysis->matador3d2->result ) {
					my ($method) = $analysis->matador3d2;					
					if ( defined $method->result ) {
						if ( defined $method->alignments ) {
							foreach my $region (@{$method->alignments}) {
								if ( 
									defined $region->pstart and defined $region->pend and 
									defined $region->pdb_id and defined $region->score
								) {
									next unless ( res_is_inside($res, $region->pstart, $region->pend) );								
									$output .=	$region->pstart.':'.$region->pend."\t".
												$region->pdb_id."\t".
												$region->score."\n";
								}
							}
						}
					}
										
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
	return $output;	
}
# matador3d2: annotation input. raw output
sub get_trans_matador3d2_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'matador3d2';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->matador3d2 and $analysis->matador3d2->result ) {
					$output .= $analysis->matador3d2->result."\n";
				}
		 	}
		}
    }
	if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }    
	return $output;	
}
# matador3d2: result input. tabular output
sub get_matador3d2_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'alignments'}) ) {
			foreach my $region (@{$seq_report->{'alignments'}}) {
				if ( 
					defined $region->{'pstart'} and defined $region->{'pend'} and 
					defined $region->{'pdb_id'} and defined $region->{'score'}
				) {
					$output .=	$region->{'pstart'}.':'.$region->{'pend'}."\t".
								$region->{'pdb_id'}."\t".
								$region->{'score'}."\n";
				}
			}
		}
	}
	return $output;
}



# SPADE METHODS -------
sub get_spade_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_spade_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_spade_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_spade_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_spade_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_spade_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_spade_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_spade_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_spade_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# spade: annotation input. tabular output
sub get_trans_spade_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'spade';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->spade and $analysis->spade->result ) {
					my ($method) = $analysis->spade;
					if ( defined $method->result ) {
						if ( defined $method->regions ) {
							foreach my $region (@{$method->regions}) {					
								if ( defined $region->alignment_start and defined $region->alignment_end and
									 defined $region->hmm_acc and defined $region->hmm_name and defined $region->evalue ) {
									next unless ( res_is_inside($res, $region->alignment_start, $region->alignment_end) );
									my (%parameters) = (
														alignment_start		=> $region->alignment_start,
														alignment_end		=> $region->alignment_end,
														hmm_acc			=> $region->hmm_acc,
														hmm_name			=> $region->hmm_name,
														evalue				=> $region->evalue,
									);
									$output .=	$parameters{alignment_start}."\t".
												$parameters{alignment_end}."\t".
												$parameters{hmm_acc}."\t".							
												$parameters{hmm_name}."\t".							
												$parameters{evalue}."\n";				
								}
							}
						}
						#if ( $output eq '' ) {
						#	$output .= "No domains\n";
						#}
						#else {
						#	my ($rst) = $method->result;
						#	$rst =~ s/^\>[^\>]*//;
						#	$output .= 	"\n";
						#	$output .= "### PfamScan alignments ###\n";
						#	$output .= "----------------------------------------------------------------------\n";
						#	$output .= 	$rst;
						#}
					}

										
				}
		 	}
		}
    }
	if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
	return $output;	
}
# spade: annotation input. raw output
sub get_trans_spade_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'spade';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->spade and $analysis->spade->result ) {
					$output .= $analysis->spade->result."\n";
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }
	return $output;	
}
# spade: result input. tabular output
sub get_spade_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
			foreach my $region (@{$seq_report->{'domains'}}) {				
				if ( defined $region->{'alignment_start'} and defined $region->{'alignment_end'} and
					 defined $region->{'hmm_acc'} and defined $region->{'hmm_name'} and defined $region->{'e_value'} ) {
					my (%parameters) = (
										alignment_start		=> $region->{'alignment_start'},
										alignment_end		=> $region->{'alignment_end'},
										hmm_acc			=> $region->{'hmm_acc'},
										hmm_name			=> $region->{'hmm_name'},
										evalue				=> $region->{'e_value'},
					);
					$output .=	$parameters{alignment_start}."\t".
								$parameters{alignment_end}."\t".
								$parameters{hmm_acc}."\t".							
								$parameters{hmm_name}."\t".							
								$parameters{evalue}."\n";				
				}
			}
			#if ( $output eq '' ) {
			#	$output .= "No domains\n";
			#}
			#else {
			#	my ($rst) = $method->{'result'};
			#	$rst =~ s/^\>[^\>]*//;
			#	$output .= 	"\n";
			#	$output .= "### PfamScan alignments ###\n";
			#	$output .= "----------------------------------------------------------------------\n";
			#	$output .= 	$rst;
			#}
					
		}
	}
	return $output;
}



# CORSAIR METHODS -------
sub get_corsair_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_corsair_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_corsair_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_corsair_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_corsair_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_corsair_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_corsair_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_corsair_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_corsair_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# corsair: annotation input. tabular output
sub get_trans_corsair_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'corsair';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->corsair and $analysis->corsair->result ) {
					#my ($method) = $analysis->corsair;
					#if ( defined $method->vertebrate_signal and defined $method->result ) {
					#	if ( defined $method->alignments ) {
					#		foreach my $region (@{$method->alignments}) {
					#
					#			if ( defined $region->pstart and defined $region->pend and
					#				 defined $region->score ) {
					#				my (%parameters) = (
					#									pstart		=> $region->pstart,
					#									pend		=> $region->pend,
					#									score		=> $region->score,
					#				);
					#				$output .=	$parameters{pstart}."\t".
					#							$parameters{pend}."\t".
					#							$parameters{score}."\n";				
					#			}
					#		}
					#	}
					#}
					my ($s_report) = $analysis->corsair->result;
					$s_report =~ s/^\>[^\n]*\n//;
					while ( $s_report =~ /([^\n]*)\n/g ) {
						my ($line) = $1;
						my (@cols) = split(/\t/,$line);
						if ( ($cols[2] ne '0') and ($cols[1] ne '') ) {
							my ($sp) = $cols[0];
							unless ( $cols[1] =~ /^\s*\-/ ) {
								my ($id) = sprintf ("%.2f",$cols[1]);
								$output .=	$sp."\t".
											$id."\n";

							}							
						}
					}
			
										
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }
	return $output;	
}
# corsair: annotation input. raw output
sub get_trans_corsair_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'corsair';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->corsair and $analysis->corsair->result ) {
					$output .= $analysis->corsair->result."\n";
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }
	return $output;	
}
# corsair: result input. tabular output
sub get_corsair_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
			my ($s_report) = $seq_report->{'result'};
			$s_report =~ s/^\>[^\n]*\n//;
			while ( $s_report =~ /([^\n]*)\n/g ) {
				my ($line) = $1;
				my (@cols) = split(/\t/,$line);
				if ( ($cols[2] ne '0') and ($cols[1] ne '') ) {
					my ($sp) = $cols[0];
					unless ( $cols[1] =~ /^\s*\-/ ) {
						my ($id) = sprintf ("%.2f",$cols[1]);
						$output .=	$sp."\t".
									$id."\n";
					}
				}
			}
			
		}
	}
	return $output;
}




# corsair_alt METHODS -------
sub get_corsair_alt_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_corsair_alt_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_corsair_alt_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_corsair_alt_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_corsair_alt_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_corsair_alt_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_corsair_alt_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_corsair_alt_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_corsair_alt_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# corsair_alt: annotation input. tabular output
sub get_trans_corsair_alt_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'corsair_alt';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->corsair_alt and $analysis->corsair_alt->result ) {
					my ($s_report) = $analysis->corsair_alt->result;
					$s_report =~ s/^\>[^\n]*\n//;
					while ( $s_report =~ /([^\n]*)\n/g ) {
						my ($line) = $1;
						my (@cols) = split(/\t/,$line);
						if ( ($cols[2] ne '0') and ($cols[1] ne '') ) {
							my ($sp) = $cols[0];
							unless ( $cols[1] =~ /^\s*\-/ ) {
								my ($id) = sprintf ("%.2f",$cols[1]);
								$output .=	$sp."\t".
											$id."\n";

							}							
						}
					}
			
										
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }
	return $output;	
}
# corsair_alt: annotation input. raw output
sub get_trans_corsair_alt_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'corsair_alt';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->corsair_alt and $analysis->corsair_alt->result ) {
					$output .= $analysis->corsair_alt->result."\n";
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }
	return $output;	
}
# corsair_alt: result input. tabular output
sub get_corsair_alt_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
			my ($s_report) = $seq_report->{'result'};
			$s_report =~ s/^\>[^\n]*\n//;
			while ( $s_report =~ /([^\n]*)\n/g ) {
				my ($line) = $1;
				my (@cols) = split(/\t/,$line);
				if ( ($cols[2] ne '0') and ($cols[1] ne '') ) {
					my ($sp) = $cols[0];
					unless ( $cols[1] =~ /^\s*\-/ ) {
						my ($id) = sprintf ("%.2f",$cols[1]);
						$output .=	$sp."\t".
									$id."\n";
					}
				}
			}
			
		}
	}
	return $output;
}




# THUMP METHODS -------
sub get_thump_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_thump_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_thump_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_thump_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_thump_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_thump_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_thump_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_thump_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_thump_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# thump: annotation input. tabular output
sub get_trans_thump_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'thump';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->thump and $analysis->thump->result ) {					
					my ($method) = $analysis->thump;					
					if ( defined $method->result ) {
						if ( defined $method->regions ) {
							foreach my $region (@{$method->regions}) {
								if ( defined $region->pstart and defined $region->pend ) {
									next unless ( res_is_inside($res, $region->pstart, $region->pend) );
									my (%parameters) = (
													start				=> $region->pstart,
													end					=> $region->pend,
									);
									my ($dam) = '-';
									$dam = 'damaged' if ( defined $region->damaged );
									$output .=	$parameters{start}."\t".
												$parameters{end}."\t".
												$dam."\n";
								}
							}
						}
					}
	
			
										
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }
	return $output;	
}
# thump: annotation input. raw output
sub get_trans_thump_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'thump';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->thump and $analysis->thump->result ) {
					$output .= $analysis->thump->result."\n";
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }
	return $output;	
}
# thump: result input. tabular output
sub get_thump_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
			foreach my $region (@{$seq_report->{'tmhs'}}) {
				if ( defined $region->{'start'} and defined $region->{'end'} ) {
					my (%parameters) = (
									start				=> $region->{'start'},
									end					=> $region->{'end'},
					);
					my ($dam) = '-';
					$dam = 'damaged' if ( defined $region->{'damaged'} );
					$output .=	$parameters{start}."\t".
								$parameters{end}."\t".
								$dam."\n";
				}
			}
	
					
		}
	}
	return $output;
}




# CRASH METHODS -------
sub get_crash_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_crash_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_crash_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_crash_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_crash_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_crash_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_crash_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_crash_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_crash_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# crash: annotation input. tabular output
sub get_trans_crash_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'crash';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->crash and $analysis->crash->result ) {					
					my ($method) = $analysis->crash;					
					if ( defined $method->result and defined $method->peptide_signal and defined $method->mitochondrial_signal ) {
						if ( ($method->peptide_signal eq $OK_LABEL) or ($method->mitochondrial_signal eq $OK_LABEL) ) {
							if ( defined $method->regions ) {
								if ( ($method->peptide_signal eq $OK_LABEL) and ($method->mitochondrial_signal eq $OK_LABEL) ) {
									$output .=	"Peptide-Mitochondrial signal\t";					
								}
								elsif ( $method->peptide_signal eq $OK_LABEL ) {
									$output .=	"Peptide signal\t";
								}					
								elsif ( $method->mitochondrial_signal eq $OK_LABEL ) {
									$output .=	"Mitochondrial signal\t";
								}
								foreach my $region (@{$method->regions}) {
									if ( defined $region->pstart and defined $region->pend ) {
										next unless ( res_is_inside($res, $region->pstart, $region->pend) );
										my (%parameters) = (
															start				=> $region->pstart,
															end					=> $region->pend,
										);
										$output .=	$parameters{start}."\t".
													$parameters{end}."\n";
									}					
								}
							}						
						}
					}
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }
	return $output;	
}
# crash: annotation input. raw output
sub get_trans_crash_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'crash';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->crash and $analysis->crash->result ) {
					$output .= $analysis->crash->result."\n";
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }
	return $output;	
}
# crash: result input. tabular output
sub get_crash_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
			my ($method) = $seq_report;
			if ( defined $method->{'peptide_signal'} and defined $method->{'mitochondrial_signal'} ) {
				if ( ($method->{'peptide_signal'} eq $OK_LABEL) or ($method->{'mitochondrial_signal'} eq $OK_LABEL) ) {
					if ( defined $method->{'regions'} ) {
						if ( ($method->{'peptide_signal'} eq $OK_LABEL) and ($method->{'mitochondrial_signal'} eq $OK_LABEL) ) {
							$output .=	"Peptide-Mitochondrial signal\t";					
						}
						elsif ( $method->{'peptide_signal'} eq $OK_LABEL ) {
							$output .=	"Peptide signal\t";
						}					
						elsif ( $method->{'mitochondrial_signal'} eq $OK_LABEL ) {
							$output .=	"Mitochondrial signal\t";
						}
						foreach my $region (@{$method->{'regions'}}) {
							if ( defined $region->{'pstart'} and defined $region->{'pend'} ) {
								my (%parameters) = (
													start				=> $region->{'pstart'},
													end					=> $region->{'pend'},
								);
							$output .=	$parameters{start}."\t".
										$parameters{end}."\n";
							}					
						}
					}						
				}
			}
		}
	}
	return $output;
}







# APPRIS METHODS -------
sub get_appris_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_appris_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_appris_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_appris_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_appris_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_appris_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_appris_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_appris_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_appris_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# appris: annotation input. tabular output
sub get_trans_appris_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'appris';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->appris and $analysis->appris->result ) {
					my ($method) = $analysis->appris;
					if ( defined $method->principal_isoform_signal and defined $method->reliability and defined $method->result ) {
						my ($principal_isoform_signal) = '-';
						if ( $method->principal_isoform_signal eq $OK_LABEL ) {
							$principal_isoform_signal = $APPRIS_PRINC_LABEL;
						}
						elsif ( $method->principal_isoform_signal eq $UNKNOWN_LABEL ) {
							$principal_isoform_signal = $APPRIS_CANDI_LABEL;
						}
						elsif ( $method->principal_isoform_signal eq $NO_LABEL ) {
							$principal_isoform_signal = $APPRIS_ALTER_LABEL;
						}						
						my (%parameters) = (
											principal_isoform_signal	=> $principal_isoform_signal,
											reliability					=> $method->reliability,
						);
						$output .=	$parameters{principal_isoform_signal}."\t".
									$parameters{reliability}."\n";		
					}
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }
	return $output;	
}
# appris: annotation input. raw output
sub get_trans_appris_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'appris';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->appris and $analysis->appris->result ) {
					$output .= $analysis->appris->result."\n";
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }
	return $output;	
}
# appris: result input. tabular output
sub get_appris_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'result'}) ) {
			my ($method) = $seq_report;
			if ( defined $method->{'principal_isoform_signal'} and defined $method->{'reliability'} and defined $method->{'result'} ) {
				my ($principal_isoform_signal) = '-';
				if ( $method->{'principal_isoform_signal'} eq $OK_LABEL ) {
					$principal_isoform_signal = $APPRIS_PRINC_LABEL;
				}
				elsif ( $method->{'principal_isoform_signal'} eq $UNKNOWN_LABEL ) {
					$principal_isoform_signal = $APPRIS_CANDI_LABEL;
				}
				elsif ( $method->{'principal_isoform_signal'} eq $NO_LABEL ) {
					$principal_isoform_signal = $APPRIS_ALTER_LABEL;
				}
				my (%parameters) = (
									principal_isoform_signal	=> $principal_isoform_signal,
									reliability					=> $method->{'reliability'},
				);
				$output .=	$parameters{principal_isoform_signal}."\t".
							$parameters{reliability}."\n";		
			}
					
			
		}
	}
	return $output;
}



# PROTEO METHODS -------
sub get_proteo_annot {
	my ($imethod, $itype, $feature) = @_;
	my ($output) = '';
	if ( $itype eq 'tsv' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'tsv'};	
	}
	elsif ( $itype eq 'raw' ) {
		$output .= $METHOD_HEADS->{$imethod}->{'raw'};
	}
	if ($feature and (ref($feature) ne 'ARRAY')) { # one object
	   	if ($feature->isa("APPRIS::Gene")) {
			foreach my $transcript (@{$feature->transcripts}) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_proteo_tsv_annot($transcript);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_proteo_raw_annot($transcript);
				}
			}
 		}
    	elsif ($feature->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_proteo_tsv_annot($feature);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_proteo_raw_annot($feature);
				}
    	}
    	else {
			throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
    	}
	}
	elsif ($feature and (ref($feature) eq 'ARRAY') ) { # listref of objects
	   	foreach my $feat (@{$feature}) {
	    	if ($feat->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feat->transcripts}) {
					if ( $itype eq 'tsv' ) {
						$output .= get_trans_proteo_tsv_annot($transcript);
					}
					elsif ( $itype eq 'raw' ) {
						$output .= get_trans_proteo_raw_annot($transcript);
					}
				}
	    	}
	   		elsif ($feat->isa("APPRIS::Transcript")) {
				if ( $itype eq 'tsv' ) {
					$output .= get_trans_proteo_tsv_annot($feat);
				}
				elsif ( $itype eq 'raw' ) {
					$output .= get_trans_proteo_raw_annot($feat);
				}
			}
			else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
			}    		
		}
	}
	return $output;
}
# proteo: annotation input. tabular output
sub get_trans_proteo_tsv_annot {
	my ($transcript,$res) = @_;
	my ($imet) = 'proteo';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->proteo and $analysis->proteo->result ) {
					my ($method) = $analysis->proteo;					
					if ( defined $method->result ) {
						if ( defined $method->peptides ) {
							foreach my $region (@{$method->peptides}) {
								if ( defined $region->sequence and defined $region->num_experiments and defined $region->pstart and defined $region->pend ) {
									$output .=	$region->pstart.':'.$region->pend."\t".
												$region->sequence."\t".
												$region->num_experiments."\n";
								}
							}
						}
					}				
				}
		 	}
		}
    }
    if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'tsv'} . $output }    
	return $output;	
}
# proteo: annotation input. raw output
sub get_trans_proteo_raw_annot {
	my ($transcript) = @_;
	my ($imet) = 'proteo';
	my ($output) = '';
	if (ref($transcript) and $transcript->isa("APPRIS::Transcript")) {   	    
		if ($transcript->stable_id) {
			my ($transcript_id) = $transcript->stable_id;
			if ( $transcript->analysis ) {
				my ($analysis) = $transcript->analysis;				
				if ( $analysis->proteo and $analysis->proteo->result ) {
					$output .= $analysis->proteo->result."\n";
				}
		 	}
		}
    }
	if ( $output ne '' ) { $output = $METHOD_HEADS->{$imet}->{'raw'} . $output }    
	return $output;	
}
# proteo: result input. tabular output
sub get_proteo_tsv_result
{
	my ($method, $report) = @_;
	my ($output) = $METHOD_HEADS->{$method}->{'tsv'};
	
	while (my ($seq_id, $seq_report) = each(%{$report}) ) {
		if ( ($seq_id ne 'result') and (exists $seq_report->{'alignments'}) ) {
			foreach my $region (@{$seq_report->{'alignments'}}) {
				if ( ($region->{'type'} eq 'mini-exon') and 
					defined $region->{'pstart'} and defined $region->{'pend'} and 
					defined $region->{'pdb_id'} and defined $region->{'identity'} ) {
						$output .=	$region->{'pstart'}.':'.$region->{'pend'}."\t".
									$region->{'pdb_id'}."\t".
									$region->{'identity'}."\n";
				}
			}
						
		}
	}
	return $output;
}
1;
