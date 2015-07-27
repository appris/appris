=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Utils::ProCDS - Utility functions for protein handling

=head1 SYNOPSIS

  use APPRIS::Utils::ProCDS
    qw(
       get_protein_cds_sequence
     );

  or to get all methods just

  use APPRIS::Utils::ProCDS;

  eval { get_protein_cds_sequence(cds_list) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to handle files.

=head1 METHODS

=cut

package APPRIS::Utils::ProCDS;

use strict;
use warnings;

use APPRIS::CDS;
use APPRIS::ProCDS;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	sort_cds
	get_protein_cds_sequence
	get_contained_cds
	get_cds_init_length
	get_coords_from_residue
	get_cds_from_residue
);

=head2 sort_cds

  Arg [1]    : Listref of APPRIS::CDS $cds_list
  Arg [2]    : (optional) String - the aminoacid sequence 
               that for this peptide
  Example    : use APPRIS::Utils::ProCDS qw(get_protein_cds_sequence);
               get_protein_cds_sequence($cds_list);
  Description: Sort list of CDS depending orientation.
  Returntype : Listref of APPRIS::ProCDS or undef
  Exceptions : none

=cut

sub sort_cds($$)
{
	my ($cds_list, $strand) = @_;

	my ($sort_cds_list);

	# Sort the CDS depending orientation from transcript
	if ($cds_list) {
		if ($strand eq '-') {
			@{$sort_cds_list}= sort { $b->start <=> $a->start } @{$cds_list};			
		}
		else {
			@{$sort_cds_list}= sort { $a->start <=> $b->start } @{$cds_list};				
		}		
	}
	return $sort_cds_list;
}

=head2 get_protein_cds_sequence

  Arg [1]    : Listref of APPRIS::CDS $cds_list
  Arg [2]    : (optional) String - the aminoacid sequence 
               that for this peptide
  Example    : use APPRIS::Utils::ProCDS qw(get_protein_cds_sequence);
               get_protein_cds_sequence($cds_list);
  Description: Get the peptide coordinates and sequence from CDS list.
  Returntype : Listref of APPRIS::ProCDS or undef
  Exceptions : none

=cut

sub get_protein_cds_sequence($;$)
{
	my ($cds_list, $sequence) = @_;

	my ($protein_cds_list);
	
	my ($pro_cds_start) = 0;
	my ($pro_cds_end) = 0;
	my ($start_phase) = 0;
	my ($end_phase) = 0;

	foreach my $cds (@{$cds_list}) {

		my ($accumulate) = 0;
		$start_phase = $end_phase;
		if($start_phase == 1) {
			$accumulate = 2;
		}
		elsif($start_phase == 2) {
			$accumulate = 1;
		}
		else {
			$accumulate = 0;
		}
		my ($cds_start) = $cds->start;
		my ($cds_end) = $cds->end;
				
		$pro_cds_start = $pro_cds_end + 1;
		my ($pro_cds_end_div) = (abs($cds_end - $cds_start) + 1 - $accumulate) / 3;
		if ($pro_cds_end_div == 0) {
			$pro_cds_end = $pro_cds_start + $pro_cds_end_div;
			$end_phase = 0;
		}
		elsif ($pro_cds_end_div =~ /(\d+)\.(\d{2})/) {
			my ($pro_cds_end_int) = $1;
			$pro_cds_end = $pro_cds_start + $pro_cds_end_int;			
			my ($pro_cds_end_mod) = '0.'.$2;
			if ($pro_cds_end_mod == '0.33') {
				$end_phase = 1;
			}
			elsif ($pro_cds_end_mod == '0.66') {
				$end_phase = 2;
			}
		}
		else {
			$pro_cds_end = $pro_cds_start + $pro_cds_end_div - 1;
			$end_phase = 0;
		}
		my ($protein_cds) = APPRIS::ProCDS->new(
											-start				=> $pro_cds_start,
											-end				=> $pro_cds_end,
											-start_phase		=> $start_phase,
											-end_phase			=> $end_phase
										);
		if (defined $sequence) {
			my ($pro_cds_length) = ($pro_cds_end - $pro_cds_start) + 1;
			my ($pro_seq) = substr($sequence, ($pro_cds_start - 1), $pro_cds_length);
			$protein_cds->sequence($pro_seq) if (defined $pro_seq);
		}
		push(@{$protein_cds_list}, $protein_cds);
	}
	
	return $protein_cds_list;
}	

=head2 get_contained_cds

  Arg [1]    : Listref of APPRIS::CDS $cds_list
  Arg [2]    : Int - the start location of looking region
  Arg [3]    : Int - the end location of looking region  
  Example    : use APPRIS::Utils::ProCDS qw(get_contained_cds);
               get_contained_cds($cds_list,$start,$end);
  Description: Get the CDS list that are within given region
  Returntype : Listref of APPRIS::CDS or undef
  Exceptions : none

=cut

sub get_contained_cds($$$)
{
	my ($cds_list, $residue_start, $residue_end) = @_;
	
	my ($contained_cds);

#use Data::Dumper;
#open(DATA_7, ">>/tmp/data7.ccds.log");		

	# Sort the exons depending orientation from transcript
	# For that we take the orientation of the first CDS
	my ($trans_strand) = $cds_list->[0]->strand;
	my ($sort_cds_list) = sort_cds($cds_list, $trans_strand);

#print DATA_7 "SORT_CDS_LIST:$trans_strand:\n".Dumper($sort_cds_list)."\n";

	for (my $i = 0; $i < scalar(@{$sort_cds_list}); $i++) {
		my ($cds_start) = $sort_cds_list->[$i]->start;
		my ($cds_end) = $sort_cds_list->[$i]->end;
		my ($cds_strand) = $sort_cds_list->[$i]->strand;
		my ($cds_phase) = $sort_cds_list->[$i]->phase;

#print DATA_7 "RES_COORD: start: $residue_start end: $residue_end\n";
#print DATA_7 "CDS_COORD: cds_start: $cds_start cds_end: $cds_end cds_strand: $cds_strand cds_phase: $cds_phase\n";					


		# within one CDS (strand +)
		if ($cds_strand eq '+' and $residue_start >= $cds_start and $residue_end <= $cds_end) {

#print DATA_7 "ONE_CDS:$cds_strand: $residue_start <= $cds_end and $residue_end >= $cds_start\n";
			
			push(@{$contained_cds},
				APPRIS::CDS->new
				(
					-start		=> $cds_start,
					-end		=> $cds_end,
					-strand		=> $cds_strand,
					-phase		=> $cds_phase
				)
			);
			last;			
		}
		# within one CDS (strand -)
		elsif ($cds_strand eq '-' and $residue_start >= $cds_start and $residue_end <= $cds_end) {

#print DATA_7 "ONE_CDS:$cds_strand: $residue_end >= $cds_start and $residue_start <= $cds_end\n";
			
			push(@{$contained_cds},
				APPRIS::CDS->new
				(
					-start		=> $cds_start,
					-end		=> $cds_end,
					-strand		=> $cds_strand,
					-phase		=> $cds_phase
				)
			);
			last;			
		}
		# within several CDS's (strand +)
		elsif ($cds_strand eq '+' and $residue_start <= $cds_end and $residue_end >= $cds_end) {
			
#print DATA_7 "MORE_CCDS:$cds_strand: $residue_start <= $cds_end and $residue_end > $cds_end\n";

			push(@{$contained_cds},
				APPRIS::CDS->new
				(
					-start		=> $cds_start,
					-end		=> $cds_end,
					-strand		=> $cds_strand,
					-phase		=> $cds_phase
				)
			);
			
			for (my $j = $i+1; $j < scalar(@{$sort_cds_list}); $j++) {
				my ($next_cds_start) = $sort_cds_list->[$j]->start;
				my ($next_cds_end) = $sort_cds_list->[$j]->end;
				my ($next_cds_strand) = $sort_cds_list->[$j]->strand;
				my ($next_cds_phase) = $sort_cds_list->[$j]->phase;

				if ($residue_end >= $next_cds_start and $residue_end >= $next_cds_end) {
					push(@{$contained_cds},
						APPRIS::CDS->new
						(
							-start		=> $next_cds_start,
							-end		=> $next_cds_end,
							-strand		=> $next_cds_strand,
							-phase		=> $next_cds_phase
						)
					);
				}
				elsif($residue_end >= $next_cds_start and $residue_end <= $next_cds_end) {
					push(@{$contained_cds},
						APPRIS::CDS->new
						(
							-start		=> $next_cds_start,
							-end		=> $next_cds_end,
							-strand		=> $next_cds_strand,
							-phase		=> $next_cds_phase
						)
					);					
					last;
				}
			}
			last;			
		}
		# within several CDS's (strand -)
		elsif ($cds_strand eq '-' and $residue_start <= $cds_start and $residue_end >= $cds_start) {
#print DATA_7 "MORE_CCDS:$cds_strand: $residue_end < $cds_end and $residue_start >= $cds_start\n";

			push(@{$contained_cds},
				APPRIS::CDS->new
				(
					-start		=> $cds_start,
					-end		=> $cds_end,
					-strand		=> $cds_strand,
					-phase		=> $cds_phase
				)
			);
						
			for (my $j = $i+1; $j < scalar(@{$sort_cds_list}); $j++) {
				my ($next_cds_start) = $sort_cds_list->[$j]->start;
				my ($next_cds_end) = $sort_cds_list->[$j]->end;
				my ($next_cds_strand) = $sort_cds_list->[$j]->strand;
				my ($next_cds_phase) = $sort_cds_list->[$j]->phase;

#print DATA_7 "NEXT_CDS_COORD: cds_start: $next_cds_start cds_end: $next_cds_end\n";
#print DATA_7 "NEXT_CDS_COORD: $residue_end <= $next_cds_start\n";
				if ($residue_start <= $next_cds_end) {
#print DATA_7 "ENTRA\n";					
					push(@{$contained_cds},
						APPRIS::CDS->new
						(
							-start		=> $next_cds_start,
							-end		=> $next_cds_end,
							-strand		=> $next_cds_strand,
							-phase		=> $next_cds_phase
						)
					);					
				}
			}
			last;			
		}
	}
	
#print DATA_7 "CONTAINED_CDS:\n".Dumper($contained_cds)."\n";
#close(DATA_7);	
	
	return $contained_cds;
		
} # End get_contained_cds

=head2 get_cds_init_length

  Arg [1]    : Listref of APPRIS::CDS $cds_list
  Arg [2]    : Int - start postion of the transcript
  Arg [3]    : Int - end position of the transcript
  Arg [4]    : Char - '+','-' the strand the transcript is on  
  Example    : use APPRIS::Utils::ProCDS qw(get_cds_init_length);
               get_cds_init_length($cds, $trans_start, $trans_end, $trans_strand);
  Description: Get the init exon and its length depending orientation.
  Returntype : APPRIS::ProCDS or undef
  Exceptions : none

=cut

sub get_cds_init_length($$$$)
{
	my ($cds, $trans_start, $trans_end, $trans_strand) = @_;
	my ($init,$length);

	if ( $trans_strand eq '-' )
	{
		$init = abs($cds->end - $trans_end);
		$length = $cds->end - $cds->start +1;
	}
	else
	{
		$init = $cds->start - $trans_start;
		$length = $cds->end - $cds->start +1;
	}	
	return ($init,$length);
}

=head2 get_coords_from_residue

  Arg [1]    : APPRIS::Transcript
  Arg [2]    : Int - the aminoacid residue that for this protein
  Example    : use APPRIS::Utils::ProCDS qw(get_coords_from_residue);
               get_coords_from_residue($transcript,$residue);
  Description: Get genomic region from peptide position.
  Returntype : APPRIS::CDS or undef
  Exceptions : none

=cut

sub get_coords_from_residue($$)
{
	my ($transcript, $residue) = @_;
	my ($protein_cds);
	my ($residue_start);
	my ($residue_end);
	my ($trans_start) = $transcript->start;
	my ($trans_end) = $transcript->end;
	my ($trans_strand) = $transcript->strand;
	my ($trans_length) = length($transcript->sequence);
	my ($cds_list) = $transcript->translate->cds;

#print STDERR "TRANSCRIPT_ID:".$transcript->stable_id."\n";
#print STDERR "TRANS_LENGTH: $trans_length\n";	

	# Sort the cds depending orientation from transcript
	my ($sort_cds_list) = sort_cds($cds_list, $trans_strand);

	# Start transcript residue
	my ($transcript_nucleotide_relative_position) = 0;
	my ($j) = 1; # First residue
	while ( $j <= 3 )
	{
		$transcript_nucleotide_relative_position = ($residue-1)*3 + $j;

#print STDERR "J:$transcript_nucleotide_relative_position = ($residue-1)*3 + $j\n";
		my ($cds_length_accumulate) = 0;
		foreach my $cds (@{$sort_cds_list})
		{
			my ($cds_start) = $cds->start;
			my ($cds_end) = $cds->end;

			my ($cds_init, $cds_length) = get_cds_init_length($cds, $trans_start, $trans_end, $trans_strand);			
			$cds_length_accumulate += $cds_length;
#print STDERR "CSD: start: $cds_start end: $cds_end init: $cds_init length: $cds_length: accum: $cds_length_accumulate\n";	
#print STDERR "CSD: cond: ($cds_length_accumulate / $transcript_nucleotide_relative_position) >= 1\n";	
			if ( ($cds_length_accumulate / $transcript_nucleotide_relative_position) >= 1 )
			{
#print STDERR "CSD: cond: si\n";
				my ($relative_position_from_exon) = $cds_length_accumulate - $cds_length;
				if ($j == 1)
				{
					if ( $trans_strand eq '-' )
					{
						$residue_start = $cds_end - ($transcript_nucleotide_relative_position - $relative_position_from_exon) +1;
#print STDERR "RES_START:$residue:$residue_start = $cds_end - ($transcript_nucleotide_relative_position-$relative_position_from_exon) +1\n";
					} else {
						$residue_start = $cds_start + ($transcript_nucleotide_relative_position - $relative_position_from_exon) -1;
#print STDERR "RES_START:$residue:$residue_start = $cds_start + ($transcript_nucleotide_relative_position-$relative_position_from_exon) -1\n";
					}	
					last;				
				}
				elsif ($j == 3)
				{
					if ( $trans_strand eq '-' )
					{
						$residue_end = $cds_end - ($transcript_nucleotide_relative_position - $relative_position_from_exon) +1;
#print STDERR "RES_END:$residue:$residue_end = $cds_end - ($transcript_nucleotide_relative_position-$relative_position_from_exon) +1\n";						
					} else {
						$residue_end = $cds_start + ($transcript_nucleotide_relative_position - $relative_position_from_exon) -1;
#print STDERR "RES_END:$residue:$residue_end = $cds_start + ($transcript_nucleotide_relative_position-$relative_position_from_exon) -1\n";
					}			
					last;
				}
			}
		}
		$j = $j+2; # Third transcrip residue for aminoacid
	}
#print STDERR "FRAMESHIFT: $transcript_nucleotide_relative_position = ($trans_length +1) and !(defined $residue_end)\n";	
	# In the case that residue posotion is within "frameshift" (Start-End codon does not found)
	if ($transcript_nucleotide_relative_position = ($trans_length +1) and !(defined $residue_end)) {
		if ( $trans_strand eq '-' )
		{
			$residue_end = $sort_cds_list->[scalar(@{$sort_cds_list})-1]->start;
			$residue_start = $residue_end - 3;
#print STDERR "RES_END:$residue:$residue_end\n";
		} else {
			$residue_end = $sort_cds_list->[scalar(@{$sort_cds_list})-1]->end;
			$residue_start = $residue_end - 3;
#print STDERR "RES_END:$residue:$residue_end\n";
		}			
	}
#print STDERR "RES_START:$residue_start RES_END:$residue:$residue_end\n";	
	if ( defined $residue_start and defined $residue_end )
	{
		$protein_cds = APPRIS::CDS->new(
						-start		=> $residue_start,
						-end		=> $residue_end,
						-strand		=> $trans_strand,
		);
	}
	return $protein_cds;
	
} # End get_coords_from_residue

=head2 get_cds_from_residue

  Arg [1]    : APPRIS::Transcript
  Arg [2]    : Int - the aminoacid residue that for this protein
  Example    : use APPRIS::Utils::ProCDS qw(get_cds_from_residue);
               get_cds_from_residue($transcript,$residue);
  Description: Get cds coordinates from peptide position.
  Returntype : APPRIS::CDS or undef
  Exceptions : none

=cut

sub get_cds_from_residue($$)
{
	my ($transcript, $residue) = @_;
	my ($protein_cds);
	my ($residue_start);
	my ($residue_end);
	my ($trans_start) = $transcript->start;
	my ($trans_end) = $transcript->end;
	my ($trans_strand) = $transcript->strand;
	my ($trans_length) = length($transcript->sequence);
	my ($cds_list) = $transcript->translate->cds;
	
	# Sort the cds depending orientation from transcript
	my ($sort_cds_list) = sort_cds($cds_list, $trans_strand);

	# Start transcript residue
	my ($transcript_nucleotide_relative_position) = 0;
	my ($j) = 1; # First residue
	while ( $j <= 3 )
	{
		$transcript_nucleotide_relative_position = ($residue-1)*3 + $j;
		my ($cds_length_accumulate) = 0;
		foreach my $cds (@{$sort_cds_list})
		{
			my ($cds_start) = $cds->start;
			my ($cds_end) = $cds->end;

			my ($cds_init, $cds_length) = get_cds_init_length($cds, $trans_start, $trans_end, $trans_strand);			
			$cds_length_accumulate += $cds_length;			
			if ( ($cds_length_accumulate / $transcript_nucleotide_relative_position) >= 1 )
			{
				unless ( defined $protein_cds ) {
					push(@{$protein_cds}, APPRIS::CDS->new(
									-start		=> $cds_start,
									-end		=> $cds_end,
									-strand		=> $trans_strand
					));					
				}
			}
		}
		$j = $j+2; # Third transcrip residue for aminoacid
	}
	# In the case that residue posotion is within "frameshift" (Start-End codon does not found)
	if ($transcript_nucleotide_relative_position = ($trans_length +1) and !(defined $protein_cds)) {
		my ($cds) = $sort_cds_list->[scalar(@{$sort_cds_list})-1];
		my ($cds_start) = $cds->start;
		my ($cds_end) = $cds->end;		
		push(@{$protein_cds}, APPRIS::CDS->new(
						-start		=> $cds_start,
						-end		=> $cds_end,
						-strand		=> $trans_strand
		));
	}
	return $protein_cds;
	
} # End get_cds_from_residue

1;