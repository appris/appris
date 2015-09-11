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
  Arg [3]    : String - $source
               List of sources ('all', ... )
  Arg [4]    : String - $typebed
               flag of head title ('yes','no','only')
  Arg [5]    : String - $version
  Arg [6]    : String - $date
  Example    : get_annotations($feature,'chr22:20116979-20137016','no','appris');  
  Description: Retrieves text as BED format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_annotations {
    my ($feature, $position, $source, $typebed) = @_;
	my ($chromosome,$start,$end) = (undef,undef,undef);
	my ($pos);
    my ($output) = '';
    my ($title_typebed) = '';
    if ( defined $typebed and ($typebed eq 'bedDetail') ) {
    	$title_typebed = 'type='.$typebed;    	
    }    
    my ($track) = [
		# APPRIS
		[
			'appris',
			{
				'title' => "track name=APPRIS_principal_isoform description='".$METHOD_LABEL_DESC->{'appris'}."' $title_typebed visibility=2 color='180,95,4' ",
				'body' => '' 
			},
		],
		# Firestar
		[
			'firestar',
			{
				'title' => "track name=Known_functional_residues description='".$METHOD_LABEL_DESC->{'firestar'}."' $title_typebed visibility=2 color='244,169,39' ",			
				'body' => '' 
			},
		],		
		# Matador3D
		[
			'matador3d',
			{
				'title' => "track name=Known_3D_structure description='".$METHOD_LABEL_DESC->{'matador3d'}."' $title_typebed visibility=2 color='128,0,0' ",
				'body' => '' 
			},
		],		
		# SPADE
		[
			'spade',
			{
				'title' => "track name=Functional_domains description='".$METHOD_LABEL_DESC->{'spade'}->[0]."' $title_typebed visibility=2 color='164,189,83' ",
				'body' => '' 
			},
			{
				'title' => "track name=Damaged_functional_domains description='".$METHOD_LABEL_DESC->{'spade'}->[1]."' $title_typebed visibility=2 color='86,191,85' ",
				'body' => '' 
			},
		],
		# CORSAIR
		[	
			'corsair',
			{
				'title' => "track name=Cross_species_evidence description='".$METHOD_LABEL_DESC->{'corsair'}."' $title_typebed  visibility=2 color='4,95,180' ",
				'body' => '' 
			},	
		],	
		# INERTIA
		[
			'inertia',
			{
				'title' => "track name=Neutrally_evolving_exons description='".$METHOD_LABEL_DESC->{'inertia'}->[0]."' $title_typebed  visibility=2 color='190,129,247' ",
				'body' => '' 
			},
			{
				'title' => "track name=Unusually_evolving_exons description='".$METHOD_LABEL_DESC->{'inertia'}->[1]."' $title_typebed  visibility=2 color='190,129,247' ",
				'body' => '' 
			}
		],
		# CRASH
		[
			'crash',
			{
				'title' => "track name=Signal_peptide_sequence description='".$METHOD_LABEL_DESC->{'crash'}->[0]."' $title_typebed  visibility=2 color='153,153,102' ",
				'body' => '' 
			},
			{
				'title' => "track name=Mitochondrial_signal_sequence description='".$METHOD_LABEL_DESC->{'crash'}->[1]."' $title_typebed  visibility=2 color='153,153,102' ",
				'body' => '' 
			}
		],
		# THUMP
		[
			'thump',
			{
				'title' => "track name=Transmembrane_helices description='".$METHOD_LABEL_DESC->{'thump'}->[0]."' $title_typebed  visibility=2 color='245,169,208' ",
				'body' => '' 
			},
			{
				'title' => "track name=Damaged_transmembrane_helices description='".$METHOD_LABEL_DESC->{'thump'}->[1]."' $title_typebed  visibility=2 color='245,169,208' ",
				'body' => '' 
			}
		],	
		# PROTEO
		[
			'proteo',
			{
				#'title' => "track type=bedGraph description='".$METHOD_LABEL_DESC->{'proteo'}."' visibility=full color='213,181,117' altColor='153,102,0' maxHeightPixels='32'",
				'title' => "track name=CNIO_Proteomic_Evidence description='".$METHOD_LABEL_DESC->{'proteo'}."' $title_typebed visibility=2 color='153,102,0' ",
				'body' => '' 
			}
		],	
	];
    
    # Convert position value for BED format
    if ( $position =~ /^([^\:]*):([^\-]*)-([^\$]*)$/ ) {
		($chromosome, $start, $end) = ($1,$2,$3);
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
    if ( defined $source and ($source ne '') ) {
		if ($feature and (ref($feature) ne 'ARRAY')) {
	    	if ($feature->isa("APPRIS::Gene")) {
				foreach my $transcript (@{$feature->transcripts}) {
					get_trans_annotations($typebed, $track, $transcript, $position, $source);
				}
				if ( $source =~ /proteo/ ) {
					get_gen_annotations($typebed, $track, $feature, $position, $source);
				}				
	    	}
	    	elsif ($feature->isa("APPRIS::Transcript")) {
	    		get_trans_annotations($typebed, $track, $feature, $position, $source);
	    	}
	    	else {
				throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
	    	}
	    }
		elsif ($feature and (ref($feature) eq 'ARRAY') ) { # in the case that we have a list of objects
	    	foreach my $feat (@{$feature}) {
		    	if ($feat->isa("APPRIS::Gene")) {
					foreach my $transcript (@{$feat->transcripts}) {
			    		get_trans_annotations($typebed, $track, $transcript, $position, $source);
					}
					if ( $source =~ /proteo/ ) {
						get_gen_annotations($typebed, $track, $feat, $position, $source);
					}					
		    	}
		    	elsif ($feat->isa("APPRIS::Transcript")) {
		    		get_trans_annotations($typebed, $track, $feat, $position, $source);
		    	}
		    	else {
					throw('Argument must be an APPRIS::Gene or APPRIS::Transcript');
		    	}    		
	    	}
		}
	    else {
			throw('Argument must be define');
	   	}
		$output .= print_annotations($track, $typebed);
    }
    else {
		throw('Arguments must be define');
    }
    
    # print head or not
    # Print output
	if ( $output ne '' ) {
	    if ( !defined $typebed or (defined $typebed and ($typebed !~ /bed12\_/) ) ) {
			$output =
				"browser position $position"."\n".
				"browser hide all"."\n".
				"browser full knownGene"."\n".
				"browser full wgEncodeGencodeVM4"."\n". # HARD-CORE!!! due UCSC does not update correctly
				"browser full ensGene"."\n".
				"browser full ccdsGene"."\n".
				$output;    	
	    }
	}
	return $output;
}

=head2 get_trans_annotations

  Arg [1]    : APPRIS::Transcript or undef
  Arg [4]    : String - $typebed
               flag of head title ('bed4', 'bed12', 'bedDetail', 'bed12+1', 'yes','no','only')
  Arg [2]    : String - $position
               genome position (chr22:20116979-20137016)
  Arg [3]    : String - $source
               List of sources ('all', ... )
  Example    : get_annotations($feature,'chr22:20116979-20137016','appris');  
  Description: Retrieves bed information of transcript.
  Returntype : Nothing (reference output)

=cut

sub get_gen_annotations {
    my ($typebed, $track, $gene, $position, $source) = @_;
    
    if (ref($gene) and $gene->isa("APPRIS::Gene")) {
		if ( ($source =~ /proteo/) or ($source eq 'all') ) {
			get_g_proteo_annotations( $gene,
          								$typebed,
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
  Arg [4]    : String - $typebed
               flag of head title ('bed4', 'bed12', 'bedDetail', 'bed12+1', 'yes','no','only')
  Arg [2]    : String - $position
               genome position (chr22:20116979-20137016)
  Arg [3]    : String - $source
               List of sources ('all', ... )
  Example    : get_annotations($feature,'chr22:20116979-20137016','appris');  
  Description: Retrieves bed information of transcript.
  Returntype : Nothing (reference output)

=cut

sub get_trans_annotations {
    my ($typebed, $track, $feature, $position, $source) = @_;

    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {
   	    
		if ($feature->stable_id) {
			if ($feature->translate and $feature->translate->sequence) {
				my ($gene_id);
				my ($transcript_id) = $feature->stable_id;
				my ($external_id) = $feature->external_name;
				if ($feature->xref_identify) {
					foreach my $xref_identify (@{$feature->xref_identify}) {
						if ( $xref_identify->dbname eq 'Ensembl_Gene_Id' ) {
							$gene_id = $xref_identify->id;							
						}
					}		
				}
				if ( ($source =~ /appris/) or ($source eq 'all') ) {				
					get_appris_annotations(	$transcript_id,
	           									$feature,
	           									$typebed,
	           									\$track->[0]
					);
				}
				if ( ($source =~ /firestar/) or ($source eq 'all') ) {				
					get_firestar_annotations(	$transcript_id,
	           									$feature,
	           									$typebed,
	           									\$track->[1]
					);
				}
				if ( ($source =~ /matador3d/) or ($source eq 'all') ) {				
					get_matador3d_annotations(	$transcript_id,
	           									$feature,
	           									$typebed,
	           									\$track->[2]
					);
				}
				if ( ($source =~ /spade/) or ($source eq 'all') ) {
					if ( $source =~ /spade\-([^\,\$]*)/ ) {
						foreach my $s_name (split('\|', $1)) {
							if ( $s_name eq 'domain' ) {		
								get_spade_annotations(	$transcript_id,
				           								$feature,
				           								$typebed,
				           								\$track->[3],
				           								$s_name
								);
							}
							elsif ( $s_name eq 'damaged_domain' ) {		
								get_spade_annotations(	$transcript_id,
				           								$feature,
				           								$typebed,
				           								\$track->[3],
				           								$s_name
								);
							}
						}
					} else {
						get_spade_annotations(	$transcript_id,
		           								$feature,
		           								$typebed,
		           								\$track->[3],
		           								$source
						);						
					}	
				}
				if ( ($source =~ /corsair/) or ($source eq 'all') ) {				
					get_corsair_annotations(	$transcript_id,
	           									$feature,
	           									$typebed,
	           									\$track->[4]
					);
				}
				if ( ($source =~ /inertia/) or ($source eq 'all') ) {				
					get_inertia_annotations(	$transcript_id,
	           									$feature,
	           									$typebed,
	           									\$track->[5]
					);
				}
				if ( ($source =~ /crash/) or ($source eq 'all') ) {				
					get_crash_annotations(	$transcript_id,
	           									$feature,
	           									\$track->[6]
					);
				}
				if ( ($source =~ /thump/) or ($source eq 'all') ) {				
					get_thump_annotations(	$transcript_id,
	           									$feature,
	           									$typebed,
	           									\$track->[7]
					);
				}
				# Now, the tracks of PROTEO are printed per peptide (gene)
				#if ( ($source =~ /proteo/) or ($source eq 'all') ) {
				#	get_proteo_annotations(	$transcript_id,
	           	#								$feature,
	           	#								$typebed,
	           	#								\$track->[8]
				#	);
				#}
			}
		}
    }
    else {
		throw('Argument must be correct');
   	}
}

=head2 print_data

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

sub print_data {
	my ($data) = @_;
	
	my ($output) = '';
	my ($chromStart) = ($data->{'start'}-1);
	my ($chromEnd) = $data->{'end'};
	my ($chromStarts) = 0;

    # Convert position value for BED format
    my ($pos) = $data->{'chr'};
    if ( !($pos =~ /^chr/) ) {
    	$pos = 'chr'.$pos;
	}
	else {
		$pos =~ s/\.([0-9]*)/\-$1/g;
	}
    
	$output .= $pos."\t".
				($data->{'start'}-1)."\t".
				$data->{'end'}."\t".
				$data->{'name'}."\t";

	if (exists $data->{'score'} and defined $data->{'score'} and
		exists $data->{'strand'} and defined $data->{'strand'}
	){
		$output .= $data->{'score'}."\t".
							$data->{'strand'}."\t";
	}

	if (exists $data->{'thick_start'} and defined $data->{'thick_start'} and
		exists $data->{'thick_end'} and defined $data->{'thick_end'} and
		exists $data->{'color'} and defined $data->{'color'} and
		exists $data->{'blocks'} and defined $data->{'blocks'}
	){
		$output .= ($data->{'thick_start'}-1)."\t".
							$data->{'thick_end'}."\t".
							$data->{'color'}."\t".
							$data->{'blocks'}."\t";
	}
	if (exists $data->{'block_sizes'} and defined $data->{'block_sizes'}) {
		my (@block_sizes) = @{$data->{'block_sizes'}};
		if ($data->{'strand'} eq '-') {
			@block_sizes=reverse @{$data->{'block_sizes'}};
		}
		
		foreach my $size (@block_sizes)
		{
			$output .= $size.',';
			$chromStarts += $size;
		}
		$output =~ s/,$/\t/;
	}

	if (exists $data->{'block_starts'} and defined $data->{'block_starts'}) {
		my (@block_starts) = @{$data->{'block_starts'}};
		if($data->{'strand'} eq '-') {
			@block_starts = reverse @{$data->{'block_starts'}};
		}
		
		foreach my $start (@block_starts) {
			$output .= $start.',';
		}
		$output =~ s/,$/\t/;
	}
	
	if (exists $data->{'note'} and defined $data->{'note'}) {
		$output .= $data->{'note'}."\t";
	}
	
	# BED chromStarts[i]+chromStart must be less or equal than chromEnd:
	unless ( ($chromStart + $chromStarts) <= $chromEnd ) {
		$output = '#'.$output;
	}

	$output =~ s/\t$/\n/;	
	return $output;
}

=head2 print_annotations

  Arg [1]    : Hash - internal hash of bed annotations
  Example    : $annot = print_annotations($report);
  Description: Print the bed annotations.
  Returntype : String or ''

=cut
sub print_annotations {
	my ($report, $typebed) = @_;
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
						if ( !defined $typebed or (defined $typebed and ($typebed !~ /bed12\_/) ) ) {
							$output .= $track_report->{'title'}."\n".$track_report->{'body'};
						}
						else {
							$output .= $track_report->{'body'};
						}
						#if ( defined $typebed and $typebed =~ /no:name/ ) {
						#	$output .= $track_report->{'body'};
						#}
						#else {							
						#	$output .= $track_report->{'title'}."\n".$track_report->{'body'};
						#}
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

=head2 get_appris_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_appris_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_appris_annotations {
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

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
					# get initial data					
					my ($exon_list) = $feature->exons;
					my ($num_exons) = scalar(@{$exon_list});
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
					my ($data) = {
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
					foreach my $exon (@{$exon_list}) {
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
					if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {
						if ( $method->reliability ) {
							$data->{'note'} = $method->reliability;
						}						
					}
					${$ref_output}->[1]->{'body'} .= print_data($data);
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
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->firestar ) {
	 		my ($method) = $analysis->firestar;
	 		# get residue annotations
			if ( defined $method->residues ) {

				# get initial data
				my ($res_list) = $method->residues;				
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
					$thick_start = $trans_start;
					$thick_end = $trans_end;
				}
				else {
					$trans_start = $res_list->[0]->start; 
					$trans_end = $res_list->[$num_res-1]->end;
					$thick_start = $trans_start;
					$thick_end = $trans_end;
				}
				my ($data) = {
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
										#elsif ( $cds_phase eq '2') {
										#	$thick_start = $thick_start +1;
										#}
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
											$thick_start = $thick_start;	
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
						if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {
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
								$data->{'note'} = $lig if ( $lig ne '' );
							}							
						}						
					}
				}
				${$ref_output}->[1]->{'body'} .= print_data($data);
			}
 		}
 	}
}

=head2 get_matador3d_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_matador3d_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_matador3d_annotations {
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->spade ) { 			
	 		my ($method) = $analysis->matador3d;
	 		# get residue annotations
			if ( defined $method->alignments ) {
				
				# get the residues with 'mini-exon' info
				my ($res_list) = $method->alignments;				
				my ($num_res) = scalar(@{$res_list});
				my ($res_exon);
				my ($res_mini_exon);
				foreach my $res (@{$res_list}) {
					if ( $res->type eq 'exon' ) {
						push(@{$res_exon}, $res);
					}
					elsif ( $res->type eq 'mini-exon' ) {
						push(@{$res_mini_exon}, $res);
					}
				}
				
				_aux_get_matador3d_annotations('mini-exon',
											$transcript_id,
											$feature,
											$res_mini_exon,
											$typebed,
											$ref_output);
			}
 		}
 	}
}

sub _aux_get_matador3d_annotations {
	my ($type, $transcript_id, $feature, $aux_res_list, $typebed, $ref_output) = @_;
	
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
		my ($data) = {
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
					if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {
						if ( $res->pdb_id ) {
							$data->{'note'} = $res->pdb_id;
						}
					}										
				}
			}
		}
		if ( $type eq 'mini-exon' ) {
			${$ref_output}->[1]->{'body'} .= print_data($data);							
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
#	my ($transcript_id, $feature, $typebed, $ref_output, $source) = @_;
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
#				if ( $source eq 'domain' ) {
#					_aux_get_spade_annotations('domain',
#												$transcript_id,
#												$feature,
#												$res_domains,
#												$typebed,
#												$ref_output);
#				}
#				elsif ( $source eq 'damaged_domain' ) {
#					_aux_get_spade_annotations('damaged_domain',
#												$transcript_id,
#												$feature,
#												$res_damaged_domains,
#												$typebed,
#												$ref_output);
#				}
#				else {
#					_aux_get_spade_annotations('domain',
#												$transcript_id,
#												$feature,
#												$res_domains,
#												$typebed,
#												$ref_output);
#					_aux_get_spade_annotations('damaged_domain',
#												$transcript_id,
#												$feature,
#												$res_damaged_domains,
#												$typebed,
#												$ref_output);
#				}
#			}
# 		}
# 	}
#}
# PRINT the alignments from Pfam separetly
sub get_spade_annotations {
	my ($transcript_id, $feature, $typebed, $ref_output, $source) = @_;

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
					if ( $source eq 'domain' ) {
						_aux_get_spade_annotations('domain',
													$transcript_id,
													$feature,
													$res_domains,
													$typebed,
													$ref_output);
					}
					elsif ( $source eq 'damaged_domain' ) {
						_aux_get_spade_annotations('damaged_domain',
													$transcript_id,
													$feature,
													$res_damaged_domains,
													$typebed,
													$ref_output);
					}
					else {
						_aux_get_spade_annotations('domain',
													$transcript_id,
													$feature,
													$res_domains,
													$typebed,
													$ref_output);
						_aux_get_spade_annotations('damaged_domain',
													$transcript_id,
													$feature,
													$res_damaged_domains,
													$typebed,
													$ref_output);
					}
				}
			}
 		}
 	}
}

sub _aux_get_spade_annotations {
	my ($type, $transcript_id, $feature, $aux_res_list, $typebed, $ref_output) = @_;

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
		my ($data) = {
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
				if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {
					if ( $res->hmm_name ) {
						#$data->{'note'} = $res->hmm_name;
						$data->{'name'} = $res->hmm_name;
						$data->{'note'} = $transcript_id;
					}
				}				
			}
		}
		if ( $type eq 'domain' ) {
			${$ref_output}->[1]->{'body'} .= print_data($data);								
		}
		elsif ( $type eq 'damaged_domain' ) { # join the whole damaged domains
			${$ref_output}->[2]->{'body'} .= print_data($data);								
		}
	}	
}

=head2 get_corsair_annotations

  Arg [1]    : String - the stable identifier of transcript
  Arg [4]    : APPRIS::Transcript
  Arg [5]    : Object - internal BED variable 
  Example    : $annot = get_corsair_annotations($trans_id, $feat, $ref_out);  
  Description: Retrieves specific annotation.
  Returntype : String or undef

=cut

sub get_corsair_annotations {
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->corsair ) {
	 		my ($method) = $analysis->corsair;
			if ( $feature->exons ) {
				
				# get initial data					
				my ($exon_list) = $feature->exons;
				my ($num_exons) = scalar(@{$exon_list});
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
				my ($data) = {
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
				foreach my $exon (@{$exon_list}) {					
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
				if ( $method->vertebrate_signal and 
						( ($method->vertebrate_signal eq $APPRIS::Utils::Constant::OK_LABEL) or 
						($method->vertebrate_signal eq $APPRIS::Utils::Constant::UNKNOWN_LABEL) ) 
				) {
					${$ref_output}->[1]->{'body'} .= print_data($data);								
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
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

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
				_aux_get_inertia_annotations('neutral_evolution',
											$transcript_id,
											$feature,
											$res_evol,
											$typebed,
											$ref_output);
				_aux_get_inertia_annotations('unusual_evolution',
											$transcript_id,
											$feature,
											$res_u_evol,
											$typebed,
											$ref_output);
			}
 		}
 	}
}

sub _aux_get_inertia_annotations {
	my ($type, $transcript_id, $feature, $res_list2, $typebed, $ref_output) = @_;

	my ($res_list);
	if (defined $res_list2) {
		if($feature->strand eq '-') {
			@{$res_list} = sort { $b->{'start'} <=> $a->{'start'} } @{$res_list2};
		}
		else {
			@{$res_list} = sort { $a->{'start'} <=> $b->{'start'} } @{$res_list2};
		}		
	}

	if ( (ref($res_list) eq 'ARRAY') and (scalar(@{$res_list}) > 0) ) {
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
				
		my ($data) = {
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
		foreach my $res (@{$res_list}) {
			my ($pos_start) = $res->start;
			my ($pos_end) = $res->end;
			my ($pos_strand) = $res->strand;
			if ($trans_strand eq '-') {
				$pos_start = $res->end;
				$pos_end = $res->start;
			}
			else {
				$pos_start = $res->start;
				$pos_end = $res->end;
			}
			my ($init, $length) = get_block_from_exon($pos_start, $pos_end, $pos_strand, $data->{'thick_start'}, $data->{'thick_end'});
			push(@{$data->{'block_starts'}}, $init);
			push(@{$data->{'block_sizes'}}, $length);
			$data->{'blocks'}++;	

		}
		if ( $type eq 'neutral_evolution' ) {
			${$ref_output}->[1]->{'body'} .= print_data($data);								
		}
		elsif ( $type eq 'unusual_evolution' ) {
			${$ref_output}->[2]->{'body'} .= print_data($data);								
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
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

	# Get annotations
 	if ( $feature->analysis ) {
 		my ($analysis) = $feature->analysis;
 		if ( $analysis->crash ) {
	 		my ($method) = $analysis->crash;
	 		# get residue annotations
			if ( defined $method->regions ) {

				# get initial data
				my ($res_list) = $method->regions;				
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
				my ($data) = {
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
								if ( $i==0 ) {
									if ($trans_strand eq '-') {
										#$pos_start = $cds_out->start;
										#$pos_end = $res->end;
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
										#$pos_start = $res->start;
										#$pos_end = $cds_out->end;
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
					}
				}
				if ( $method->peptide_signal and ($method->peptide_signal eq $APPRIS::Utils::Constant::OK_LABEL) ){
					${$ref_output}->[1]->{'body'} .= print_data($data);								
				}
				if ( $method->mitochondrial_signal and ($method->mitochondrial_signal eq $APPRIS::Utils::Constant::OK_LABEL) ){
					${$ref_output}->[2]->{'body'} .= print_data($data);								
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
	my ($transcript_id, $feature, $typebed, $ref_output) = @_;

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
				_aux_get_thump_annotations('transmembrane_helixes',
											$transcript_id,
											$feature,
											$res_helix,
											$typebed,
											$ref_output);
				_aux_get_thump_annotations('damaged_transmembrane_helixes',
											$transcript_id,
											$feature,
											$res_helix_damaged,
											$typebed,
											$ref_output);
			}
 		}
 	}
}

sub _aux_get_thump_annotations {
	my ($type, $transcript_id, $feature, $res_list, $typebed, $ref_output) = @_;

	if ( (ref($res_list) eq 'ARRAY') and (scalar(@{$res_list}) > 0) ) {
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
						
		my ($data) = {
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
						if ( $i==0 ) {
							if ($trans_strand eq '-') {
								#$pos_start = $cds_out->start;
								#$pos_end = $res->end;
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
								#$pos_start = $res->start;
								#$pos_end = $cds_out->end;
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
			}
		}
		if ( $type eq 'transmembrane_helixes' ) {
			${$ref_output}->[1]->{'body'} .= print_data($data);								
		}
		elsif ( $type eq 'damaged_transmembrane_helixes' ) {
			${$ref_output}->[2]->{'body'} .= print_data($data);								
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

sub get_g_proteo_annotations {
	my ($gene, $typebed, $ref_output) = @_;

	# Get all peptides per gene
	my ($res_peptides);
	foreach my $transcript (@{$gene->transcripts}) {
		my ($transcript_id) = $transcript->stable_id;
	 	if ( $transcript->analysis ) {
	 		my ($analysis) = $transcript->analysis;
	 		if ( $analysis->proteo ) { 			
		 		my ($method) = $analysis->proteo;
		 		# get residue annotations
				if ( defined $method->peptides ) {
					# get the peptides
					foreach my $res (@{$method->peptides}) {
						my ($pep_seq) = $res->sequence;
						my ($pep_idx) = $res->start.'-'.$res->end.':'.$res->strand;
						my ($aux_output) = _aux_get_g_proteo_annotations('peptides',
																		$transcript_id,
																		$transcript,
																		[ $res ],
																		$typebed);		
						if ( !exists $res_peptides->{$pep_seq} ) {
							$res_peptides->{$pep_seq}->{$aux_output} = $transcript_id;
						}
						elsif ( exists $res_peptides->{$pep_seq} and !exists $res_peptides->{$pep_seq}->{$aux_output} ) {
							$res_peptides->{$pep_seq}->{$aux_output} = $transcript_id;
						}
						else {
							$res_peptides->{$pep_seq}->{$aux_output} .= ','.$transcript_id;
						}
					}					
				}
	 		}
	 	}
	}
	
	# Get annotations
	while ( my ($pep_seq, $pep_rep) = each(%{$res_peptides}) ) {
		while ( my ($track, $track_rep) = each(%{$pep_rep}) ) {
			if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {
				${$ref_output}->[1]->{'body'} .= $track."\t".$track_rep."\n";
			}
			else {
				${$ref_output}->[1]->{'body'} .= $track."\n";
			}
		}
	}
}

sub _aux_get_g_proteo_annotations {
	my ($type, $transcript_id, $feature, $aux_res_list, $typebed) = @_;
	my ($output) = '';

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
		my ($data) = {
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
				if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {				
					if ( $res->num_experiments ) {
						$data->{'score'} = $res->num_experiments;
					}				
					if ( $res->sequence ) {
						#$data->{'note'} = $res->sequence;
						$data->{'name'} = $res->sequence;
						#$data->{'note'} = $transcript_id;
					}
				}		
			}
		}
		if ( $type eq 'peptides' ) {
			$output .= print_data($data);
			$output =~ s/\n*$//;
		}
	}	
	return $output;
}

## DEPRECATED but STILL USEFULL
#sub get_proteo_annotations {
#	my ($transcript_id, $feature, $ref_output) = @_;
#
#	# Get annotations
# 	if ( $feature->analysis ) {
# 		my ($analysis) = $feature->analysis;
# 		if ( $analysis->proteo ) {
#	 		my ($method) = $analysis->proteo;	 		
#			if ( defined $method->peptides ) {
#				# get the peptides
#				${$ref_output}->[1]->{'body'} .= ${$ref_output}->[1]->{'title'}."  name=$transcript_id "."\n";
#				foreach my $res (@{$method->peptides}) {
#					${$ref_output}->[1]->{'body'} .= "chr".$feature->chromosome."\t".($res->start-1)."\t".$res->end."\t".$res->num_experiments."\n";
#				}
#			}
# 		}
# 	}
#}

## THIS VERSION PRINT ALL PEPTIDES TOGETHER FOR ONE TRANSCRIPT
#sub get_proteo_annotations {
#	my ($transcript_id, $feature, $typebed, $ref_output) = @_;
#
#	# Get annotations
# 	if ( $feature->analysis ) {
# 		my ($analysis) = $feature->analysis;
# 		if ( $analysis->proteo ) { 			
#	 		my ($method) = $analysis->proteo;
#	 		# get residue annotations
#			if ( defined $method->peptides ) {
#				# get the peptides
#				my ($res_list) = $method->peptides;
#				my ($num_res) = scalar(@{$res_list});
#				my ($res_peptides);
#				foreach my $res (@{$res_list}) {
#					push(@{$res_peptides}, $res);
#				}
#				_aux_get_proteo_annotations('peptides',
#											$transcript_id,
#											$feature,
#											$res_peptides,
#											$ref_output);				
#			}
# 		}
# 	}
#}

# THIS VERSION PRINTS PEPTIDES PER TRANSCRIPT
#sub get_proteo_annotations {
#	my ($transcript_id, $feature, $typebed, $ref_output) = @_;
#
#	# Get annotations
# 	if ( $feature->analysis ) {
# 		my ($analysis) = $feature->analysis;
# 		if ( $analysis->proteo ) { 			
#	 		my ($method) = $analysis->proteo;
#	 		# get residue annotations
#			if ( defined $method->peptides ) {
#				# get the peptides
#				my ($res_list) = $method->peptides;
#				my ($num_res) = scalar(@{$res_list});
#				foreach my $res (@{$res_list}) {
#					my ($res_peptides);					
#					push(@{$res_peptides}, $res);
#					_aux_get_proteo_annotations('peptides',
#												$transcript_id,
#												$feature,
#												$res_peptides,
#												$typebed,
#												$ref_output);
#				}
#			}
# 		}
# 	}
#}
#
#sub _aux_get_proteo_annotations {
#	my ($type, $transcript_id, $feature, $aux_res_list, $typebed, $ref_output) = @_;
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
#		my ($data) = {
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
#		
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
#				if ( defined $typebed and ( ($typebed eq 'bedDetail') or ($typebed =~ /bed12\_/)) ) {				
#					if ( $res->num_experiments ) {
#						$data->{'score'} = $res->num_experiments;
#					}				
#					if ( $res->sequence ) {
#						#$data->{'note'} = $res->sequence;
#						$data->{'name'} = $res->sequence;
#						$data->{'note'} = $transcript_id;
#					}
#				}		
#			}
#		}
#		if ( $type eq 'peptides' ) {
#			${$ref_output}->[1]->{'body'} .= print_data($data);
#		}
#	}	
#}

1;