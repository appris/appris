=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::CDS - Utility functions for info exporting

=head1 SYNOPSIS

  use APPRIS::Export::CDS
    qw(
       get_trans_annotations
     );

  or to get all methods just

  use APPRIS::Export::CDS;

  eval { get_trans_annotations($feature,$params) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Retrieves sequences of transcript as fasta format.

=head1 METHODS

=cut

package APPRIS::Export::CDS;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

use APPRIS::Utils::ProCDS qw(get_protein_cds_sequence);

use APPRIS::Utils::Exception qw(throw warning deprecate);

use vars qw(
	$NUM_RESIDUES
);

$NUM_RESIDUES = 70;

=head2 get_trans_annotations

  Arg [1]    : Listref of APPRIS::Gene or APPRIS::Transcript or undef
  Arg [2]    : String $type
               type of sequence ('na' or 'aa')
  Example    : $annot = $exporter->get_trans_annotations($feature,'aa');
  Description: Retrieves nucleotide o aminoacid sequence.
  Returntype : String or undef

=cut

sub get_trans_annotations {
    my ($feature,$type) = @_;
    my ($report);
	
    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {
		my ($id);
		if ($feature->stable_id) {
			$id = $feature->stable_id;
			if ( $feature->translate ) {
				my ($translate) = $feature->translate;
				my ($cds_rep);
				$report->{'id'} = $id;
				if ( $translate->cds and $feature->translate->sequence ) {
					my ($sequence) = $feature->translate->sequence;
					my ($exons) = $translate->cds;
					my ($cds_sequence) = get_protein_cds_sequence($exons, $sequence);
					
					for (my $icds = 0; $icds < scalar(@{$cds_sequence}); $icds++) {
						my ($cds) = $translate->cds->[$icds];
						my ($exon) = $exons->[$icds];
						my ($pro_cds) = $cds_sequence->[$icds];
	
						my ($cds_start) = $cds->start;
						my ($cds_end) = $cds->end;
						my ($cds_strand) = $cds->strand;
						my ($cds_phase) = $cds->phase;
						
						my ($pro_cds_start) = $pro_cds->start;
						my ($pro_cds_end) = $pro_cds->end;
						my ($pro_cds_end_phase) = $pro_cds->end_phase;
						my ($pro_cds_seq) = $pro_cds->sequence;
						
						push(@{$report->{'cds'}}, {
							'start'		=> $pro_cds->start,
							'end'		=> $pro_cds->end,
							'end_phase'	=> $pro_cds->end_phase,
							'seq'		=> $pro_cds->sequence,
							'length'	=> length($pro_cds->sequence),
						});

					}
				}
								
			}			
		}		
    }
    else {
		throw('Argument must be an APPRIS::Transcript');
   	}
	return $report;
}

1;
