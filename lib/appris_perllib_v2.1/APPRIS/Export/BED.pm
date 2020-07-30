=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::BED - Utility functions for error handling

=head1 SYNOPSIS

  use APPRIS::Export::BED qw(get_annotations);
  
  or to get all methods just

  use APPRIS::Export::BED;

  eval { get_annotations("text to file",file_path) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to export database values as BED format.

=head1 METHODS

=cut

package APPRIS::Export::BED;

use strict;
use warnings;
use Data::Dumper;

use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::Constant qw(
	$METHOD_LABEL_DESC
);

###################
# Global variable #
###################
use vars qw(
	$LABEL_TRACKS
	$METHOD_LABEL_DESC
);

$LABEL_TRACKS = {
	'unknown'	=> $APPRIS::Utils::Constant::UNKNOWN_LABEL,
	'ok'		=> $APPRIS::Utils::Constant::OK_LABEL,
	'no'		=> $APPRIS::Utils::Constant::NO_LABEL,
};

$METHOD_LABEL_DESC = $APPRIS::Utils::Constant::METHOD_LABEL_DESC;

=head2 get_annotations

  Arg [1]    : Listref of APPRIS::Gene or undef
  Arg [2]    : String - $position
               genome position (chr22:20116979-20137016)
  Arg [3]    : String - $methods
               List of sources ('all', ... )
  Example    : get_annotations($feature,'chr22:20116979-20137016','no','appris');  
  Description: Retrieves text as BED format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_annotations {
    my ($source, $feature, $position, $methods, $typebed) = @_;
	my ($chromosome,$start,$end) = (undef,undef,undef);
	my ($pos);
    my ($output) = '';
    my ($track) = [
		# APPRIS
		[
			'appris',
			{
				'title' => "track name=APPRIS_principal_isoform description='".$METHOD_LABEL_DESC->{'appris'}."' visibility=2 color='180,95,4' ",
				'body' => '' 
			},
		],
		# Firestar
		[
			'firestar',
			{
				'title' => "track name=Known_functional_residues description='".$METHOD_LABEL_DESC->{'firestar'}."' visibility=2 color='244,169,39' ",			
				'body' => '' 
			},
		],		
		# Matador3D
		[
			'matador3d',
			{
				'title' => "track name=Known_3D_structure description='".$METHOD_LABEL_DESC->{'matador3d'}."' visibility=2 color='128,0,0' ",
				'body' => '' 
			},
		],		
		# SPADE
		[
			'spade',
			{
				'title' => "track name=Functional_domains description='".$METHOD_LABEL_DESC->{'spade'}->[0]."' visibility=2 color='164,189,83' ",
				'body' => '' 
			},
			{
				'title' => "track name=Damaged_functional_domains description='".$METHOD_LABEL_DESC->{'spade'}->[1]."' visibility=2 color='86,191,85' ",
				'body' => '' 
			},
		],
		# CORSAIR
		[	
			'corsair',
			{
				'title' => "track name=Cross_species_evidence description='".$METHOD_LABEL_DESC->{'corsair'}."'  visibility=2 color='4,95,180' ",
				'body' => '' 
			},	
		],	
		# INERTIA
		[
			'inertia',
			{
				'title' => "track name=Neutrally_evolving_exons description='".$METHOD_LABEL_DESC->{'inertia'}->[0]."'  visibility=2 color='190,129,247' ",
				'body' => '' 
			},
			{
				'title' => "track name=Unusually_evolving_exons description='".$METHOD_LABEL_DESC->{'inertia'}->[1]."'  visibility=2 color='190,129,247' ",
				'body' => '' 
			}
		],
		# CRASH
		[
			'crash',
			{
				'title' => "track name=Signal_peptide_sequence description='".$METHOD_LABEL_DESC->{'crash'}->[0]."'  visibility=2 color='153,153,102' ",
				'body' => '' 
			},
			{
				'title' => "track name=Mitochondrial_signal_sequence description='".$METHOD_LABEL_DESC->{'crash'}->[1]."'  visibility=2 color='153,153,102' ",
				'body' => '' 
			}
		],
		# THUMP
		[
			'thump',
			{
				'title' => "track name=Transmembrane_helices description='".$METHOD_LABEL_DESC->{'thump'}->[0]."'  visibility=2 color='245,169,208' ",
				'body' => '' 
			},
			{
				'title' => "track name=Damaged_transmembrane_helices description='".$METHOD_LABEL_DESC->{'thump'}->[1]."'  visibility=2 color='245,169,208' ",
				'body' => '' 
			}
		],	
		# PROTEO
		[
			'proteo',
			{
				'title' => "track name=CNIO_Proteomic_Evidence description='".$METHOD_LABEL_DESC->{'proteo'}."' visibility=2 color='153,102,0' ",
				'body' => '' 
			}
		],	
	];
    
    # Convert position value for BED format
    if ( defined $position and ($position =~ /^([^\:]*):([^\-]*)-([^\$]*)$/) ) {
		($chromosome, $start, $end) = ($1,$2,$3);
		if ( $chromosome =~ /^NC_/ ) {
			$chromosome =~ s/NC_[0]*//g;
			$chromosome =~ s/\.[0-9]*$//g;
			if ( $chromosome eq '23' ) { $chromosome = 'X' }
			if ( $chromosome eq '24' ) { $chromosome = 'Y' }
			if ( $chromosome eq '12920' ) { $chromosome = 'M' }				
		}
		$pos = $chromosome;
		
		if ( !($pos =~ /^chr/) ) {
			$pos = 'chr'.$pos;
		}
		else {
			$pos =~ s/\.([0-9]*)/\-$1/g;
		}
		$position = "$pos:$start\-$end";
    }        

	# Get the bed annotations
    if ( defined $methods and ($methods ne '') ) {
		if ($feature and (ref($feature) ne 'ARRAY')) {
	    	if ($feature->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feature->transcripts}) {
					get_trans_annotations($typebed, $track, $transcript, $methods);
				}
				#if ( $methods =~ /proteo/ ) {
				#	get_gen_annotations($track, $feature, $methods);
				#}				
	    	}
	    	elsif ($feature->isa("APPRIS::Transcript")) {
	    		get_trans_annotations($typebed, $track, $feature, $methods);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}
	    }
		elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
	    	foreach my $feat (@{$feature}) {
		    	if ($feat->isa("APPRIS::Gene")) {
					foreach my $transcript (@{$feat->transcripts}) {
			    		get_trans_annotations($typebed, $track, $transcript, $methods);
					}
					#if ( $methods =~ /proteo/ ) {
					#	get_gen_annotations($track, $feat, $methods);
					#}					
		    	}
		    	elsif ($feat->isa("APPRIS::Transcript")) {
		    		get_trans_annotations($typebed, $track, $feat, $methods);
		    	}
		    	else {
					throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
		    	}    		
	    	}
		}
	    else {
			throw('Argument must be define');
	   	}
		$output .= print_annotations($typebed, $track);
    }
    else {
		throw('Arguments must be define');
    }
    
    # print head or not
    # Print output
	if ( $output ne '' ) {
		if ( $typebed eq 'bed' ) {
			my ($ucscGeneTrack) = '';
			if ( defined $source and $source eq 'ensembl' ) {
				$ucscGeneTrack =
					"browser full knownGene"."\n".
					"browser full wgEncodeGencodeVM4"."\n". # HARD-CORE!!! due UCSC does not update correctly					
					"browser full ensGene"."\n".
					"browser full ccdsGene"."\n";				
			}
			elsif ( defined $source and $source eq 'refseq' ) {
				$ucscGeneTrack =
					"browser full refSeqComposite"."\n".
					"browser full refGene"."\n". # deprecated
					"browser full ccdsGene"."\n";				
			}
			else {
				$ucscGeneTrack =
					"browser full knownGene"."\n".
					"browser full ensGene"."\n".
					"browser full refSeqComposite"."\n".
					"browser full refGene"."\n". # deprecated
					"browser full ccdsGene"."\n";	
			}
			$output =
				"browser position $position"."\n".
				"browser hide all"."\n".
				$ucscGeneTrack."\n".
				$output;			
		}
	}
	return $output;
}

=head2 get_trans_annotations

  Arg [1]    : APPRIS::Transcript or undef
  Arg [2]    : String - $position
               genome position (chr22:20116979-20137016)
  Arg [3]    : String - $methods
               List of sources ('all', ... )
  Example    : get_annotations($feature,'chr22:20116979-20137016','appris');  
  Description: Retrieves bed information of transcript.
  Returntype : Nothing (reference output)

=cut

sub get_gen_annotations {
    my ($track, $gene, $methods) = @_;
    
    if (ref($gene) and $gene->isa("APPRIS::Gene")) {
		if ( ($methods =~ /proteo/) or ($methods eq 'all') ) {
			get_g_proteo_annotations( $gene,
									 \$track->[8]
			);
		}
    	
    }
    else {
		throw('Argument must be correct');
   	}    
}

=head2 get_trans_annotations

  Arg [1]    : APPRIS::Transcript or undef
  Arg [2]    : String - $position
               genome position (chr22:20116979-20137016)
  Arg [3]    : String - $methods
               List of sources ('all', ... )
  Example    : get_annotations($feature,'chr22:20116979-20137016','appris');  
  Description: Retrieves bed information of transcript.
  Returntype : Nothing (reference output)

=cut

sub get_trans_annotations {
    my ($typebed, $track, $feature, $methods) = @_;

    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {
   	    
		if ($feature->stable_id) {
			if ($feature->translate and $feature->translate->sequence) {
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
				my (%sc) = map { $_ => 1 } split(',', $methods);				
				if ( (exists $sc{appris}) or ($methods eq 'all') ) {				
					get_appris_annotations(	$typebed,
											$transcript_id,
	           								$feature,
	           								\$track->[0]
					);
				}
				if ( (exists $sc{firestar}) or ($methods eq 'all') ) {				
					get_firestar_annotations( 	$typebed,
												$transcript_id,
	           									$feature,
	           									\$track->[1]
					);
				}
				if ( (exists $sc{matador3d}) or ($methods eq 'all') ) {				
					get_matador3d_annotations(	$typebed,
												$transcript_id,
	           									$feature,
	           									\$track->[2]
					);
				}
				if ( (exists $sc{matador3d2}) or ($methods eq 'all') ) {				
					get_matador3d2_annotations(	$typebed,
												$transcript_id,
	           									$feature,
	           									\$track->[2]
					);
				}
				if ( (exists $sc{spade}) or ($methods eq 'all') ) {
					if ( $methods =~ /spade\-([^\,\$]*)/ ) {
						while ( $methods =~ /spade\-([^\,\$]*)/mg ) {
							my ($s_name) = $1;
							if ( $s_name eq 'domain' ) {		
								get_spade_annotations(	$typebed,
														$transcript_id,
				           								$feature,
				           								\$track->[3],
				           								$s_name
								);
							}
							elsif ( $s_name eq 'damaged_domain' ) {		
								get_spade_annotations(	$typebed,
														$transcript_id,
				           								$feature,
				           								\$track->[3],
				           								$s_name
								);
							}
						}
					} else {
						get_spade_annotations(	$typebed,
												$transcript_id,
		           								$feature,
		           								\$track->[3],
		           								$methods
						);
						get_spade_annotations(	$typebed,
												$transcript_id,
		           								$feature,
		           								\$track->[3],
		           								'domain'
						);
						get_spade_annotations(	$typebed,
												$transcript_id,
		           								$feature,
		           								\$track->[3],
		           								'damaged_domain'
						);
					}	
				}
				if ( (exists $sc{corsair}) or ($methods eq 'all') ) {				
					get_corsair_annotations(	$typebed,
												$transcript_id,
	           									$feature,
	           									\$track->[4]
					);
				}
				if ( (exists $sc{corsair_alt}) or ($methods eq 'all') ) {				
					get_corsair_alt_annotations(	$typebed,
												$transcript_id,
	           									$feature,
	           									\$track->[4]
					);
				}
				if ( (exists $sc{inertia}) or ($methods eq 'all') ) {				
					get_inertia_annotations(	$typebed,
												$transcript_id,
	           									$feature,
	           									\$track->[5]
					);
				}
				if ( (exists $sc{crash}) or ($methods eq 'all') ) {				
					get_crash_annotations(	$typebed,
											$transcript_id,
	           								$feature,
	           								\$track->[6]
					);
				}
				if ( (exists $sc{thump}) or ($methods eq 'all') ) {				
					get_thump_annotations(	$typebed,
											$transcript_id,
	           								$feature,
	           								\$track->[7]
					);
				}
				# Now, the tracks of PROTEO are printed per peptide (gene)
				if ( (exists $sc{proteo}) or ($methods eq 'all') ) {
					get_proteo_annotations(	$typebed,
											$transcript_id,
	           								$feature,
	           								\$track->[8]
					);
				}
			}
		}
    }
    else {
		throw('Argument must be correct');
   	}
}

=head2 print_track

  Arg [1]    : Hast - data attributes of bed
               $data = {
                   'chr'          => ''
                   'start'        => ''
                   'end'          => ''
                   'name'         => ''
                   'score'        => ''
                   'strand'       => ''
                   'thick_start'  => ''
                   'thick_end'    => ''
                   'color'        => ''
                   'blocks'       => ''
                   'block_sizes'  => ''
                   'block_starts' => ''
               };
  Example    : $annot = print_data($data);
  Description: Print the bed annotations.
               Chromosome Start End Name Score Strand 
                 ThickStart ThickEnd Color Blocks BlockSizes BlockStarts
                 
                 The start values are -1 due to UCSC Genome browser prints bad
  Returntype : String or ''

=cut

sub print_track {
	my ($typebed, $data) = @_;
	
	my ($output) = '';
    my ($pos) = $data->{'chr'};
	my ($chromStart) = ($data->{'start'}-1);
	my ($chromEnd) = $data->{'end'};
	my ($name) = $data->{'name'}; $name =~ s/\s/\_/g;
	my ($score) = (exists $data->{'score'} and defined $data->{'score'}) ? $data->{'score'} : undef;
	my ($strand) = (exists $data->{'strand'} and defined $data->{'strand'}) ? $data->{'strand'} : undef;
	my ($thick_start) = (exists $data->{'thick_start'} and defined $data->{'thick_start'}) ? ($data->{'thick_start'}-1) : undef;
	my ($thick_end) = (exists $data->{'thick_end'} and defined $data->{'thick_end'}) ? $data->{'thick_end'} : undef;
	my ($color) = (exists $data->{'color'} and defined $data->{'color'}) ? $data->{'color'} : undef;
	my ($blocks) = (exists $data->{'blocks'} and defined $data->{'blocks'}) ? $data->{'blocks'} : undef;
	my ($block_sizes) = (exists $data->{'block_sizes'} and defined $data->{'block_sizes'}) ? $data->{'block_sizes'} : undef;
	my ($block_starts) = (exists $data->{'block_starts'} and defined $data->{'block_starts'}) ? $data->{'block_starts'} : undef;
	my ($note) = (exists $data->{'note'} and defined $data->{'note'}) ? $data->{'note'} : undef; $note =~ s/\s/\_/g;
	
	my ($chromStarts, $blockSizes_last, $chromStarts_first, $chromStarts_last) = (0,0,0,0);	
	my ($block_ascending_notoverlap) = 1;

    # Convert position value for BED format
	if ( $pos =~ /^NC_/ ) {
		$pos =~ s/NC_[0]*//g;
		$pos =~ s/\.[0-9]*$//g;
		if ( $pos eq '23' ) { $pos = 'X' }
		if ( $pos eq '24' ) { $pos = 'Y' }
		if ( $pos eq '12920' ) { $pos = 'M' }				
	}    
    if ( !($pos =~ /^chr/) ) {
    	$pos = 'chr'.$pos;
	}
	else {
		$pos =~ s/\.([0-9]*)/\-$1/g;
	}
    
    # for BED12: Change the Notes for Names
    if ( $typebed eq 'bed12' ) {
		if ( defined $note and defined $name ) {
			my ($na) = $name;
			my ($no) = $note;
			$name = $no;
			$note = $na;
		}    	
    }
    
	$output .= 	$pos."\t".
				$chromStart."\t".
				$chromEnd."\t".
				$name."\t";

	if ( defined $score and defined $strand ){
		$output .= 	$score."\t".
					$strand."\t";
	}

	if (defined $thick_start and defined $thick_end and defined $color and defined $blocks ){
		$output .= 	$thick_start."\t".
					$thick_end."\t".
					$color."\t".
					$blocks."\t";
	}
	my (@block_sizes);
	if ( defined $block_sizes ) {
		if ( $strand eq '-') { @block_sizes=reverse @{$block_sizes} }
		else { @block_sizes = @{$block_sizes} }		
		foreach my $size (@block_sizes) {
			$output .= $size.',';
			$chromStarts += $size;
		}
		$blockSizes_last = $block_sizes[scalar(@block_sizes)-1];		
		$output =~ s/,$/\t/;
	}
	my (@block_starts);
	if ( defined $block_starts ) {
		if ( $strand eq '-') { @block_starts = reverse @{$block_starts} }
		else { @block_starts = @{$block_starts} }		
		for (my $i=0; $i < scalar(@block_starts); $i++ ) {
			$output .= $block_starts[$i].',';
			if ( $i > 0 ) {
				if ( ($block_starts[$i-1] + $block_sizes[$i-1]) > $block_starts[$i] ) {
					$block_ascending_notoverlap = 0
				}
			} 
		}
		$chromStarts_first = $block_starts[0];
		$chromStarts_last = $block_starts[scalar(@block_starts)-1];		
		$output =~ s/,$/\t/;
	}
	
	if ( $typebed eq 'bed12' ) {
		if ( defined $note ) {
			$output .= $note."\t";
		}
	}
	
	# BED blocks must be in ascending order without overlap.
	if ( $block_ascending_notoverlap == 0 ) {
		$output = '#1# '.$output;
	}	
	# BED chromStarts[0] = 1, must be 0 so that (chromStart + chromStarts[0]) equals chromStart.
	unless ( ($chromStart + $chromStarts_first) == $chromStart ) {
		$output = '#2# '.$output;
	}
	# BED chromStarts[i]+chromStart must be less or equal than chromEnd.
	unless ( ($chromStart + $chromStarts) <= $chromEnd ) {
		$output = '#3# '.$output;
	}
	# BED (chromStart + blockSizes[last] + chromStarts[last] ) must equal chromEnd.
	unless ( ($chromStart + $blockSizes_last + $chromStarts_last ) == $chromEnd ) {
		$output = '#4# '.$output;
	}

	$output =~ s/\t$/\n/;	
	return $output;
}

=head2 merge_annotations

  Arg [1]    : String - annotations by method
  Example    : merge_annotations($output);
  Description: Get the BED block from exon.
  Returntype : (Int- init, Int- size) or (undef,undef)
  Exceptions : none

=cut

sub merge_annotations {
	my ($inannot) = @_;
	my ($output) = '';
    my ($report);
	foreach my $line ( split("\n",$inannot) ) {
		my (@cols) = split("\t", $line);
		next if ( scalar(@cols) < 12);					
		my (
			$chrom,
			$chromStart,
			$chromEnd,
			$name,
			$score,
			$strand,
			$thickStart,
			$thickEnd,
			$reserved,
			$blockCount,
			$blockSizes,
			$chromStarts,
			$transcripts
		) = @cols;
		my ($idx) = $chrom.":".$chromStart."-".$chromEnd;
		unless ( $report->{$idx} ) {
			$report->{$idx} = {
					'chrom' 		=> $chrom,
					'chromStart'	=> $chromStart,
					'chromEnd'		=> $chromEnd,
					'name'			=> $name, 
					'score'			=> $score, 
					'strand'		=> $strand,
					'thickStart'	=> $thickStart,
					'thickEnd'		=> $thickEnd,
					'reserved'		=> $reserved,
					'blockCount'	=> $blockCount,
					'blockSizes'	=> $blockSizes,
					'chromStarts'	=> $chromStarts,
					'transcripts'	=> $transcripts
			}
		}
		else {
			unless ( ($report->{$idx}->{'name'} =~ /^($name)/) or ($report->{$idx}->{'name'} =~ /\;($name)\;/) ) {
				$report->{$idx}->{'name'}		.= ';'.$name;
			}			
			$report->{$idx}->{'transcripts'}	.= ';'.$transcripts if ( defined $transcripts );
		}
	}
	if ( defined $report ) {
		foreach my $idx ( sort(keys(%{$report})) ) {
			my ($rep) = $report->{$idx};
			$output .= $rep->{'chrom'}."\t".
						$rep->{'chromStart'}."\t".
						$rep->{'chromEnd'}."\t".
						$rep->{'name'}."\t".
						$rep->{'score'}."\t".
						$rep->{'strand'}."\t".
						$rep->{'thickStart'}."\t".
						$rep->{'thickEnd'}."\t".
						$rep->{'reserved'}."\t".
						$rep->{'blockCount'}."\t".
						$rep->{'blockSizes'}."\t".
						$rep->{'chromStarts'};
			$output .= "\t".$rep->{'transcripts'} if ( defined $rep->{'transcripts'} );
			$output .= "\n";
		}
	}
	return $output;
}

=head2 print_annotations

  Arg [1]    : Hash - internal hash of bed annotations
  Example    : $annot = print_annotations($report);
  Description: Print the bed annotations.
  Returntype : String or ''

=cut
sub print_annotations {
	my ($typebed, $report) = @_;
	my ($output) = '';
	
	foreach my $method_report (@{$report}) {
		if (ref($method_report) eq 'ARRAY') {
			my ($methods);
			for ( my $i=0; $i<scalar(@{$method_report}); $i++) {
				if ( $i == 0 ) {
					$methods = $method_report->[$i];
				}
				else {
					my ($track_report) = $method_report->[$i];
					if ($track_report->{'body'} ne '') {
						if ( $typebed eq 'bed' ) {
							$output .= $track_report->{'title'}."\n".$track_report->{'body'};	
						}
						elsif ( $typebed eq 'bed12' ) {
							$output .= merge_annotations($track_report->{'body'});
						}
					}					
				}
			}
		}
	}
	return $output;
}

=head2 get_block_from_exon

  Arg [1]    : Int - start position
  Arg [2]    : Int - end position
  Arg [3]    : Int - start position of exon
  Arg [4]    : Int - end position of exon
  Example    : get_block_from_exon($e_start, $e_end, $b_start, $b_end);
  Description: Get the BED block from exon.
  Returntype : (Int- init, Int- size) or (undef,undef)
  Exceptions : none

=cut

sub get_block_from_exon {
	my ($pos_start, $pos_end, $exon_strand, $block_start, $block_end) = @_;
	my ($init,$length);
		
	my ($residue_start);
	my ($residue_size);
	my ($residue_end);
	if ($exon_strand eq '-') {
		$residue_start = abs($pos_end - $block_start);
		$residue_size = abs($pos_end - $pos_start) +1;		
	} else {
		$residue_start = abs($pos_start - $block_start);
		$residue_size = abs($pos_start - $pos_end) +1;
	}		
	
	if ( defined $residue_start and defined $residue_size ) {
		$init = $residue_start;
		$length = $residue_size;
	}
	return ($init, $length);
}

sub extract_track_cds {
	my ($transcript_id, $feature, $aux_res_list) = @_;
	my ($data);
	
	if ( (ref($aux_res_list) eq 'ARRAY') and (scalar(@{$aux_res_list}) > 0) ) {	
		my ($num_exons) = scalar(@{$aux_res_list});
		my ($trans_chr) = $feature->chromosome;
		my ($trans_start) = $feature->start;
		my ($trans_end) = $feature->end;
		my ($trans_strand) = $feature->strand;
		my ($score) = 0;
		my ($thick_start) = $feature->start;
		my ($thick_end) = $feature->end;
		my ($color) = 0;
		my ($blocks) = $num_exons;					
		if ( $feature->translate and $feature->translate->cds ) {
			my ($cds_list) = $feature->translate->cds;
			my ($num_cds) = scalar(@{$cds_list});						
			if ($trans_strand eq '-') {
				$thick_start = $cds_list->[$num_cds-1]->start;
				$thick_end = $cds_list->[0]->end;
			}
			else {
				$thick_start = $cds_list->[0]->start; 
				$thick_end = $cds_list->[$num_cds-1]->end;
			}			
		}
		$data = {
				'chr'			=> $trans_chr,
				'name'			=> $transcript_id,
				'start'			=> $trans_start,
				'end'			=> $trans_end,
				'strand'		=> $trans_strand,
				'score'			=> $score,
				'thick_start'	=> $thick_start,
				'thick_end'		=> $thick_end,			
				'color'			=> $color,
				'blocks'		=> $blocks
		};
		# get block annotations
		foreach my $exon (@{$aux_res_list}) {
			my ($pos_start) = $exon->start;
			my ($pos_end) = $exon->end;
			my ($pos_strand) = $exon->strand;
			if ($trans_strand eq '-') {
				$pos_start = $exon->end;
				$pos_end = $exon->start;
			}
			else {
				$pos_start = $exon->start;
				$pos_end = $exon->end;
			}
			my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $data->{'start'}, $data->{'end'});
			push(@{$data->{'block_starts'}}, $init);
			push(@{$data->{'block_sizes'}}, $length);
		}
 	}
 	return $data;
}
sub extract_track_region {
	my ($transcript_id, $feature, $aux_res_list, $attribues) = @_;
	my ($data);

	if ( (ref($aux_res_list) eq 'ARRAY') and (scalar(@{$aux_res_list}) > 0) ) {
		
		# sort the list of alignments
		my ($res_list);
		if ($feature->strand eq '-') {
			@{$res_list} = sort { $b->start <=> $a->start } @{$aux_res_list};
		}
		else {
			@{$res_list} = sort { $a->start <=> $b->start } @{$aux_res_list};
		}
				
		# get initial data
		my ($num_res) = scalar(@{$res_list});		
		my ($trans_chr) = $feature->chromosome;
		my ($trans_start) = $feature->start;
		my ($trans_end) = $feature->end;
		my ($trans_strand) = $feature->strand;
		my ($score) = 0;
		my ($thick_start) = $feature->start;
		my ($thick_end) = $feature->end;
		my ($color) = 0;
		my ($blocks) = 0;
		if ( $trans_strand eq '-' ) {
			$trans_start = $res_list->[$num_res-1]->start;
			$trans_end = $res_list->[0]->end;
			$thick_start = $trans_start;
			$thick_end = $trans_end;
		}
		else {
			$trans_start = $res_list->[0]->start; 
			$trans_end = $res_list->[$num_res-1]->end;
			$thick_start = $trans_start;
			$thick_end = $trans_end;
		}
		$data = {
				'chr'			=> $trans_chr,
				'name'			=> $transcript_id,
				'start'			=> $trans_start,
				'end'			=> $trans_end,
				'strand'		=> $trans_strand,
				'score'			=> $score,
				'thick_start'	=> $thick_start,
				'thick_end'		=> $thick_end,			
				'color'			=> $color,
				'blocks'		=> $blocks
		};
		# get block annotations
		if ( $feature->translate and $feature->translate->cds ) {
			my ($translation) = $feature->translate;
			my ($cds_list) = $translation->cds;
			my ($num_cds) = scalar(@{$cds_list});						
			foreach my $res (@{$res_list}) {
				my ($contained_cds) = $translation->get_overlapping_cds($res->start, $res->end);
				my (@sorted_contained_cds) = @{$contained_cds};
				for (my $i = 0; $i < scalar(@sorted_contained_cds); $i++) {
						
					if ( scalar(@sorted_contained_cds) == 1 ) { # Within one CDS
						my ($pos_start) = $res->start;
						my ($pos_end) = $res->end;
						my ($pos_strand) = $res->strand;
						if ($trans_strand eq '-') {
							$pos_start = $res->end;
							$pos_end = $res->start;
						}
						my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $data->{'thick_start'}, $data->{'thick_end'});
						push(@{$data->{'block_starts'}}, $init);
						push(@{$data->{'block_sizes'}}, $length);
						$data->{'blocks'}++;	
						last;
					}								
					else { # Within several CDS
						my ($cds_out) = $sorted_contained_cds[$i];
						my ($pos_start) = $cds_out->start;
						my ($pos_end) = $cds_out->end;
						my ($pos_strand) = $cds_out->strand;
						if ( $trans_strand eq '-' ) {
							$pos_start = $cds_out->end;
							$pos_end = $cds_out->start;
						}

						if ( $i==0 ) {
							if ($trans_strand eq '-') {
								$pos_start = $res->end;
								$pos_end = $cds_out->start;
							}
							else {
								$pos_start = $res->start;
								$pos_end = $cds_out->end;
							}
						}
						elsif ( $i == scalar(@sorted_contained_cds)-1 ) {
							if ( $trans_strand eq '-' ) {
								$pos_start = $cds_out->end;
								$pos_end = $res->start;
							}
							else {
								$pos_start = $cds_out->start;
								$pos_end = $res->end;
							}
						}
						my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $data->{'thick_start'}, $data->{'thick_end'});
						push(@{$data->{'block_starts'}}, $init);
						push(@{$data->{'block_sizes'}}, $length);
						$data->{'blocks'}++;	
					}					
				}
				foreach my $attr (@{$attribues}) {
					my ($a_name) = $attr->{'name'};
					my ($v_name) = $attr->{'value'};
					if ( $res->$v_name ) {
						$data->{$a_name} = $res->$v_name;
					}
				}
			}
		}
	}
	return $data;
}

=head2 get_appris_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_appris_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_appris_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->appris ) {
	 		my ($method) = $analysis->appris;
			if(	$method->principal_isoform_signal and (
					($method->principal_isoform_signal eq $APPRIS::Utils::Constant::OK_LABEL) or
					($method->principal_isoform_signal eq $APPRIS::Utils::Constant::UNKNOWN_LABEL) 
				) and 
				$method->reliability and ($method->reliability =~ /^PRINCIPAL/)
			) {
				if ( $feature->exons ) {
					my ($res_list) = $feature->exons;
					my ($data) = extract_track_cds($transcript_id,
													$feature,
													$res_list);
					if (defined $data ) {
						if ( $method->reliability ) {
							$data->{'name'} = $method->reliability;
							$data->{'note'} = $transcript_id;	
						}
						if ( $typebed eq 'bed' ) {
							my ($na) = $data->{'name'}; my ($no) = $data->{'note'};
							$data->{'name'} = $no;
							$data->{'note'} = $na;						
						}
						${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
					}
				}
 			}
 		}
 	}
}

=head2 get_firestar_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_firestar_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_firestar_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->firestar ) { 			
	 		my ($method) = $analysis->firestar;
	 		# get residue annotations
			if ( defined $method->residues ) {
								
				# get the residues
				my ($res_list) = $method->residues;
				if ( $typebed eq 'bed' ) {
					my ($data) = _aux_get_firestar_annotations($transcript_id,
																$feature,
																$res_list);
					if (defined $data ) {
						${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
					}
				}
				elsif ( $typebed eq 'bed12' ) {
					foreach my $res (@{$res_list}) {
						my ($data) = _aux_get_firestar_annotations($transcript_id,
																	$feature,
																	[$res]);
						if (defined $data ) {
							${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
						}
					}
				}
			}
 		}
 	}
}
sub _aux_get_firestar_annotations {
	my ($transcript_id, $feature, $aux_res_list) = @_;
	my ($data);
	if ( (ref($aux_res_list) eq 'ARRAY') and (scalar(@{$aux_res_list}) > 0) ) {
		
		# sort the list of alignments
		my ($res_list);
		if ($feature->strand eq '-') {
			@{$res_list} = sort { $b->start <=> $a->start } @{$aux_res_list};
		}
		else {
			@{$res_list} = sort { $a->start <=> $b->start } @{$aux_res_list};
		}
				
		# get initial data
		my ($num_res) = scalar(@{$res_list});		
		my ($trans_chr) = $feature->chromosome;
		my ($trans_start) = $feature->start;
		my ($trans_end) = $feature->end;
		my ($trans_strand) = $feature->strand;
		my ($score) = 0;
		my ($thick_start) = $feature->start;
		my ($thick_end) = $feature->end;
		my ($color) = 0;
		my ($blocks) = 0;
		if ( $trans_strand eq '-' ) {
			$trans_start = $res_list->[$num_res-1]->end;
			$trans_end = $res_list->[0]->start;
		}
		else {
			$trans_start = $res_list->[0]->start; 
			$trans_end = $res_list->[$num_res-1]->end;
		}
		if ( $trans_start > $trans_end ) {
			my ($d) = $trans_end;
			$trans_end = $trans_start;
			$trans_start = $d;
		}
		$thick_start = $trans_start;
		$thick_end = $trans_end;
		$data = {
				'chr'			=> $trans_chr,
				'name'			=> $transcript_id,
				'start'			=> $trans_start,
				'end'			=> $trans_end,
				'strand'		=> $trans_strand,
				'score'			=> $score,
				'thick_start'	=> $thick_start,
				'thick_end'		=> $thick_end,			
				'color'			=> $color,
				'blocks'		=> $blocks
		};
		# get block annotations
		if ( $feature->translate and $feature->translate->cds ) {
			my ($translation) = $feature->translate;
			my ($cds_list) = $translation->cds;
			my ($num_cds) = scalar(@{$cds_list});
			foreach my $res (@{$res_list}) {
				my ($res_start) = $res->start;
				my ($res_end) = $res->end;
				my ($res_strand) = $res->strand;
				if ($res_strand eq '-') {
					$res_start = $res->end;
					$res_end = $res->start;
				}
				my ($contained_cds) = $translation->get_overlapping_cds($res_start, $res_end);						
				my (@sorted_contained_cds) = @{$contained_cds};
				for (my $i = 0; $i < scalar(@sorted_contained_cds); $i++) {
					my ($cds_out) = $sorted_contained_cds[$i];
					my ($cds_strand) = $cds_out->strand;
					my ($cds_phase) = $cds_out->phase;
					if ( scalar(@sorted_contained_cds) == 1 ) { # Residue fulldown in one CDS
						my ($pos_start) = $res->start;
						my ($pos_end) = $res->end;
						my ($pos_strand) = $res->strand;
						my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $thick_start, $thick_end);
						$init = 0; # HARD-CORE
						push(@{$data->{'block_starts'}}, $init);
						push(@{$data->{'block_sizes'}}, $length);
						$data->{'blocks'}++;	
						last;
					}												
					else { # Residue fulldown in multiple CDS
						my ($pos_start) = $cds_out->start;
						my ($pos_end) = $cds_out->end;
						my ($thick_start) = $data->{'thick_start'};
						my ($thick_end) = $data->{'thick_end'};
						if ( $i==0 ) {
							if ($trans_strand eq '-') {
								$pos_start = $cds_out->start;
								$pos_end = $res->start;
								# the residue falls down between two CDS
								if ( $cds_phase eq '0') {
									$thick_start = $thick_start;
								}
								elsif ( $cds_phase eq '2') {
									$thick_start = $thick_start;
								}
							}
							else {
								$pos_start = $res->start;
								$pos_end = $cds_out->end;
							}
							
						}
						elsif ( $i == scalar(@sorted_contained_cds)-1 ) {
							if ( $trans_strand eq '-' ) {
								$pos_start = $res->end;
								$pos_end = $cds_out->end;
								# the residue falls down between two CDS
								if ( $cds_phase eq '0') {
									$thick_start = $thick_start +1;
								}
								elsif ( $cds_phase eq '2') {
									$thick_start = $thick_start +1;
								}
							}
							else {
								$pos_start = $cds_out->start;
								$pos_end = $res->end;
							}
						}
						my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $cds_strand, $thick_start, $thick_end);
						push(@{$data->{'block_starts'}}, $init);
						push(@{$data->{'block_sizes'}}, $length);
						$data->{'blocks'}++;
					}
				}
				if ( $res->ligands ) {
					my ($lig) = '';
					my (@ligands) = split(/\|/, $res->ligands);
					foreach my $ligs (@ligands) {
						if ( $ligs ne '' ) {
							if ( $ligs =~ /^([^\[]*)\[/ ) {
								$lig .= $1.',';
							}
							$ligs =~ s/^[^\[]*\[[^\,]*\,//;
							$ligs =~ s/\,[^\]]*\]$//;
						}
					}
					$lig =~ s/\,$//;
					if ( $lig ne '' ) {
						$data->{'note'} = $lig;
					}
				}
			}
		}
	}
	return $data;	
}

=head2 get_matador3d_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_matador3d_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

#sub get_matador3d_annotations {
#	my ($typebed, $transcript_id, $feature, $ref_output) = @_;
#
#	# Get annotations
# 	if ( $feature->analysis ) {
# 		my ($analysis) = $feature->analysis;
# 		if ( $analysis->matador3d ) { 			
#	 		my ($method) = $analysis->matador3d;
#	 		# get residue annotations
#			if ( defined $method->alignments ) {
#				
#				# get the residues with 'mini-exon' info
#				my ($res_list) = $method->alignments;				
#				my ($num_res) = scalar(@{$res_list});
#				my ($res_exon);
#				my ($res_mini_exon);
#				foreach my $res (@{$res_list}) {
#					if ( $res->type eq 'exon' ) {
#						push(@{$res_exon}, $res);
#					}
#					elsif ( $res->type eq 'mini-exon' ) {
#						push(@{$res_mini_exon}, $res);
#					}
#				}
#				my ($data) = extract_track_region(	$transcript_id,
#													$feature,
#													$res_mini_exon,
#													[{
#														'name' => 'note',
#														'value' => 'pdb_id'
#													}]);
#				if (defined $data ) {
#					${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
#				}
#			}
# 		}
# 	}
#}
sub get_matador3d_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->matador3d ) { 			
	 		my ($method) = $analysis->matador3d;
	 		# get residue annotations
			if ( defined $method->alignments ) {
				
				# get the residues
				my ($res_list) = $method->alignments;				
				my ($num_res) = scalar(@{$res_list});
				my ($res_align);
				foreach my $res (@{$res_list}) {
					push(@{$res_align}, $res);
				}
				my ($data) = extract_track_region(	$transcript_id,
													$feature,
													$res_align,
													[{
														'name' => 'note',
														'value' => 'pdb_id'
													}]);
				if (defined $data ) {
					${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
				}
			}
 		}
 	}
}

=head2 get_matador3d2_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_matador3d2_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_matador3d2_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->matador3d2 ) { 			
	 		my ($method) = $analysis->matador3d2;
	 		# get residue annotations
			if ( defined $method->alignments ) {
				
				# get the residues
				my ($res_list) = $method->alignments;				
				my ($num_res) = scalar(@{$res_list});
				my ($res_align);
				foreach my $res (@{$res_list}) {
					push(@{$res_align}, $res);
				}
				my ($data) = extract_track_region(	$transcript_id,
													$feature,
													$res_align,
													[{
														'name' => 'note',
														'value' => 'pdb_id'
													}]);
				if (defined $data ) {
					${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
				}
			}
 		}
 	}
}

=head2 get_spade_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_spade_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

# PRINT the alignments together from Pfam 
#sub get_spade_annotations {
#	my ($transcript_id, $feature, $ref_output, $methods) = @_;
#
#	# Get annotations
# 	if ( $feature->analysis ) {
# 		my ($analysis) = $feature->analysis;
# 		if ( $analysis->spade ) { 			
#	 		my ($method) = $analysis->spade;
#	 		# get residue annotations
#			if ( defined $method->regions ) {
#								
#				# get the residues with 'domain', 'domain_possibly_damaged', 'domain_damaged', and 'domain_wrong' separetly
#				my ($res_list) = $method->regions;
#				my ($num_res) = scalar(@{$res_list});
#				my ($res_domains);
#				my ($res_damaged_domains);
#				foreach my $res (@{$res_list}) {
#					if ( $res->type_domain eq 'domain' ) {
#						#push(@{$res_domains}, $res);
#					}
#					elsif ( $res->type_domain eq 'domain_possibly_damaged' ) {
#						push(@{$res_domains}, $res);
#					}
#					elsif ( $res->type_domain eq 'domain_damaged' ) {
#						push(@{$res_damaged_domains}, $res);
#					}
#					elsif ( $res->type_domain eq 'domain_wrong' ) {
#						push(@{$res_damaged_domains}, $res);
#					}
#				}
#				if ( $methods eq 'domain' ) {
#					_aux_get_spade_annotations('domain',
#												$transcript_id,
#												$feature,
#												$res_domains,
#												$ref_output);
#				}
#				elsif ( $methods eq 'damaged_domain' ) {
#					_aux_get_spade_annotations('damaged_domain',
#												$transcript_id,
#												$feature,
#												$res_damaged_domains,
#												$ref_output);
#				}
#				else {
#					_aux_get_spade_annotations('domain',
#												$transcript_id,
#												$feature,
#												$res_domains,
#												$ref_output);
#					_aux_get_spade_annotations('damaged_domain',
#												$transcript_id,
#												$feature,
#												$res_damaged_domains,
#												$ref_output);
#				}
#			}
# 		}
# 	}
#}
# PRINT the alignments from Pfam separetly
sub get_spade_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output, $methods) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->spade ) { 			
	 		my ($method) = $analysis->spade;
	 		# get residue annotations
			if ( defined $method->regions ) {								
				# get the residues with 'domain', 'domain_possibly_damaged', 'domain_damaged', and 'domain_wrong' separetly
				my ($res_list) = $method->regions;
				my ($num_res) = scalar(@{$res_list});
				foreach my $res (@{$res_list}) {
					my ($res_domains);
					my ($res_damaged_domains);
					if ( $res->type_domain eq 'domain' ) {
						push(@{$res_domains}, $res);
					}
					elsif ( $res->type_domain eq 'domain_possibly_damaged' ) {
						push(@{$res_domains}, $res);
					}
					elsif ( $res->type_domain eq 'domain_damaged' ) {
						push(@{$res_damaged_domains}, $res);
					}
					elsif ( $res->type_domain eq 'domain_wrong' ) {
						push(@{$res_damaged_domains}, $res);
					}
					my ($data_domain);
					if ( ($methods eq 'domain') or ($methods eq 'spade') ) {
						$data_domain = extract_track_region( $transcript_id,
															$feature,
															$res_domains,
															[{
																'name' => 'note',
																'value' => 'hmm_name'
															}]);						
					}
					my ($data_damg_domain);
					if ( ($methods eq 'damaged_domain') or ($methods eq 'spade') ) {
						$data_damg_domain = extract_track_region( $transcript_id,
															$feature,
															$res_damaged_domains,
															[{
																'name' => 'note',
																'value' => 'hmm_name'
															}]);
					}
					if (defined $data_domain ) {
						${$ref_output}->[1]->{'body'} .= print_track($typebed, $data_domain);
					}
					if (defined $data_damg_domain ) {
						${$ref_output}->[2]->{'body'} .= print_track($typebed, $data_damg_domain);
					}
				}
			}
 		}
 	}
}
#sub _aux_get_spade_annotations {
#	my ($type, $transcript_id, $feature, $aux_res_list, $attribues) = @_;
#	my ($data);
#
#	if ( (ref($aux_res_list) eq 'ARRAY') and (scalar(@{$aux_res_list}) > 0) ) {
#		
#		# sort the list of alignments
#		my ($res_list);
#		if ($feature->strand eq '-') {
#			@{$res_list} = sort { $b->start <=> $a->start } @{$aux_res_list};
#		}
#		else {
#			@{$res_list} = sort { $a->start <=> $b->start } @{$aux_res_list};
#		}
#				
#		# get initial data
#		my ($num_res) = scalar(@{$res_list});		
#		my ($trans_chr) = $feature->chromosome;
#		my ($trans_start) = $feature->start;
#		my ($trans_end) = $feature->end;
#		my ($trans_strand) = $feature->strand;
#		my ($score) = 0;
#		my ($thick_start) = $feature->start;
#		my ($thick_end) = $feature->end;
#		my ($color) = 0;
#		my ($blocks) = 0;
#		if ( $trans_strand eq '-' ) {
#			$trans_start = $res_list->[$num_res-1]->start;
#			$trans_end = $res_list->[0]->end;
#			$thick_start = $trans_start;
#			$thick_end = $trans_end;
#		}
#		else {
#			$trans_start = $res_list->[0]->start; 
#			$trans_end = $res_list->[$num_res-1]->end;
#			$thick_start = $trans_start;
#			$thick_end = $trans_end;
#		}
#		$data = {
#				'chr'			=> $trans_chr,
#				'name'			=> $transcript_id,
#				'start'			=> $trans_start,
#				'end'			=> $trans_end,
#				'strand'		=> $trans_strand,
#				'score'			=> $score,
#				'thick_start'	=> $thick_start,
#				'thick_end'		=> $thick_end,			
#				'color'			=> $color,
#				'blocks'		=> $blocks
#		};
#		# get block annotations
#		if ( $feature->translate and $feature->translate->cds ) {
#			my ($translation) = $feature->translate;
#			my ($cds_list) = $translation->cds;
#			my ($num_cds) = scalar(@{$cds_list});						
#			foreach my $res (@{$res_list}) {
#				my ($contained_cds) = $translation->get_overlapping_cds($res->start, $res->end);
#				my (@sorted_contained_cds) = @{$contained_cds};
#				for (my $i = 0; $i < scalar(@sorted_contained_cds); $i++) {
#						
#					if ( scalar(@sorted_contained_cds) == 1 ) { # Within one CDS
#						my ($pos_start) = $res->start;
#						my ($pos_end) = $res->end;
#						my ($pos_strand) = $res->strand;
#						if ($trans_strand eq '-') {
#							$pos_start = $res->end;
#							$pos_end = $res->start;
#						}
#						my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $data->{'thick_start'}, $data->{'thick_end'});
#						push(@{$data->{'block_starts'}}, $init);
#						push(@{$data->{'block_sizes'}}, $length);
#						$data->{'blocks'}++;	
#						last;
#					}								
#					else { # Within several CDS
#						my ($cds_out) = $sorted_contained_cds[$i];
#						my ($pos_start) = $cds_out->start;
#						my ($pos_end) = $cds_out->end;
#						my ($pos_strand) = $cds_out->strand;
#						if ( $trans_strand eq '-' ) {
#							$pos_start = $cds_out->end;
#							$pos_end = $cds_out->start;
#						}
#
#						if ( $i==0 ) {
#							if ($trans_strand eq '-') {
#								$pos_start = $res->end;
#								$pos_end = $cds_out->start;
#							}
#							else {
#								$pos_start = $res->start;
#								$pos_end = $cds_out->end;
#							}
#						}
#						elsif ( $i == scalar(@sorted_contained_cds)-1 ) {
#							if ( $trans_strand eq '-' ) {
#								$pos_start = $cds_out->end;
#								$pos_end = $res->start;
#							}
#							else {
#								$pos_start = $cds_out->start;
#								$pos_end = $res->end;
#							}
#						}
#						my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $data->{'thick_start'}, $data->{'thick_end'});
#						push(@{$data->{'block_starts'}}, $init);
#						push(@{$data->{'block_sizes'}}, $length);
#						$data->{'blocks'}++;	
#					}					
#				}
#				if ( $res->hmm_name ) {
#					$data->{'note'} = $res->hmm_name;
#				}				
#			}
#		}
#	}
#	return $data;
#}

=head2 get_corsair_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_corsair_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_corsair_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->corsair ) {
	 		my ($method) = $analysis->corsair;
	 		if ( defined $method->score and $method->score != 0 ) {  	 		
				if ( $feature->exons ) {
					my ($res_list) = $feature->exons;					
					my ($data) = extract_track_cds($transcript_id,
													$feature,
													$res_list);
					if (defined $data ) {
						$data->{'note'} = 'Species Conserv';
						${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
					}
				}
	 		}
 		}
 	}
}

=head2 get_corsair_alt_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_corsair_alt_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_corsair_alt_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->corsair_alt ) {
	 		my ($method) = $analysis->corsair_alt;
	 		if ( defined $method->score and $method->score != 0 ) {  	 		
				if ( $feature->exons ) {
					my ($res_list) = $feature->exons;					
					my ($data) = extract_track_cds($transcript_id,
													$feature,
													$res_list);
					if (defined $data ) {
						$data->{'note'} = 'Species Conserv';
						${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
					}
				}
	 		}
 		}
 	}
}

=head2 get_inertia_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_inertia_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_inertia_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->inertia ) {
	 		my ($method) = $analysis->inertia;
	 		# get residue annotations
			if ( defined $method->regions ) {
				# get the residues with 'neutral_evolution', and 'unusual_evolution' separetly
				my ($res_list) = $method->regions;
				my ($num_res) = scalar(@{$res_list});
				my ($res_evol);
				my ($res_u_evol);
				foreach my $res (@{$res_list}) {
					if ( $res->unusual_evolution eq $APPRIS::Utils::Constant::UNKNOWN_LABEL ) {
						push(@{$res_evol}, $res);
					}
					elsif ( $res->unusual_evolution eq $APPRIS::Utils::Constant::NO_LABEL ) {
						push(@{$res_u_evol}, $res);
					}
				}
				my ($data_evol) = extract_track_cds(	$transcript_id,
														$feature,
														$res_evol);
				my ($data_u_evol) = extract_track_cds(	$transcript_id,
														$feature,
														$res_u_evol);
				if (defined $data_evol ) {
					${$ref_output}->[1]->{'body'} .= print_track($typebed, $data_evol);
				}
				if (defined $data_u_evol ) {
					${$ref_output}->[2]->{'body'} .= print_track($typebed, $data_u_evol);
				}
			}
 		}
 	}
}

=head2 get_crash_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_crash_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_crash_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->crash ) { 			
	 		my ($method) = $analysis->crash;
	 		# get residue annotations
			if ( defined $method->regions ) {
								
				# get the residues
				my ($res_list) = $method->regions;
				my ($num_res) = scalar(@{$res_list});
				foreach my $res (@{$res_list}) {
					my ($data) = extract_track_region(	$transcript_id,
														$feature,
														[$res]);
					if (defined $data ) {
						$data->{'note'} = 'Signal Peptide';
						if ( $method->peptide_signal and ($method->peptide_signal eq $APPRIS::Utils::Constant::OK_LABEL) ){
							${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
						}
						if ( $method->mitochondrial_signal and ($method->mitochondrial_signal eq $APPRIS::Utils::Constant::OK_LABEL) ){
							${$ref_output}->[2]->{'body'} .= print_track($typebed, $data);			
						}						
					}					
				}
			}
 		}
 	}
}

=head2 get_thump_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_thump_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_thump_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->thump ) {
	 		my ($method) = $analysis->thump;
	 		# get residue annotations
			if ( defined $method->regions ) {
				# get the residues with 'transmembrane_helixes', and 'damaged_transmembrane_helixes' separetly
				my ($res_list) = $method->regions;
				my ($num_res) = scalar(@{$res_list});
				my ($res_helix);
				my ($res_helix_damaged);
				foreach my $res (@{$res_list}) {
					if (defined $res->damaged and ($res->damaged eq '1') ) {
						push(@{$res_helix_damaged}, $res);
					}
					else {
						push(@{$res_helix}, $res);
					}
				}
				my ($data_helix) = extract_track_region(	$transcript_id,
															$feature,
															$res_helix);
				my ($data_damg_helix) = extract_track_region(	$transcript_id,
																$feature,
																$res_helix_damaged);
				if (defined $data_helix ) {
					$data_helix->{'note'} = 'TMHelices';					
					${$ref_output}->[1]->{'body'} .= print_track($typebed, $data_helix);
				}					
				#if (defined $data_damg_helix ) {
				#	$data_damg_helix->{'note'} = 'Damaged TMHelices';					
				#	${$ref_output}->[2]->{'body'} .= print_track($typebed, $data_damg_helix);
				#}
			}
 		}
 	}
}

=head2 get_proteo_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_proteo_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_proteo_annotations {
	my ($typebed, $transcript_id, $feature, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->proteo ) { 			
	 		my ($method) = $analysis->proteo;
	 		# get residue annotations
			if ( defined $method->peptides ) {
								
				# get the residues
				my ($res_list) = $method->peptides;
				my ($num_res) = scalar(@{$res_list});
				foreach my $res (@{$res_list}) {
					my ($data) = extract_track_region($transcript_id,
													$feature,
													[$res],
													[{
														'name' => 'score',
														'value' => 'num_experiments'
													},
													{
														'name' => 'note',
														'value' => 'sequence'
													}]);
					if (defined $data ) {
						# BED format does not support scores > 1000.
						if ( exists($data->{'score'}) &&
								defined($data->{'score'}) &&
								$data->{'score'} > 1000 ) {
							$data->{'score'} = 1000;
						}
						${$ref_output}->[1]->{'body'} .= print_track($typebed, $data);
					}													
				}
			}
 		}
 	}
}

1;