=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::SEQ - Utility functions for info exporting

=head1 SYNOPSIS

  use APPRIS::Export::SEQ
    qw(
       get_trans_annotations
     );

  or to get all methods just

  use APPRIS::Export::SEQ;

  eval { get_trans_annotations($feature,$params) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

Retrieves sequences of transcript as fasta format.

=head1 METHODS

=cut

package APPRIS::Export::SEQ;

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

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
    my ($output) = '';
	
    if (ref($feature) and $feature->isa("APPRIS::Transcript")) {
		my ($id);
		my ($sequence);
		if ($feature->stable_id) {
			$id = $feature->stable_id;
			if ($type eq 'na') {
				if ($feature->sequence) {
					$sequence = $feature->sequence;
				}						
			} elsif ($type eq 'aa') {
				if ($feature->translate and $feature->translate->sequence) {
					$sequence = $feature->translate->sequence;
				}
			}

			my ($main_id) = ''; 
			if ($feature->xref_identify) {
				foreach my $xref_identify (@{$feature->xref_identify}) {
					if (($xref_identify->dbname eq 'Gene_Id') and $xref_identify->id) {
						$main_id = $xref_identify->id;
					}
				}				
			}
			my ($ext_name) = '';
			$ext_name = $feature->external_name if ($feature->external_name);
			
			if (defined $id and defined $sequence) {
				my($seqobj)=Bio::Seq->new(
					-display_id => $id,
					-seq => $sequence
				);
				my ($slength) = $seqobj->length;
				$output .= ">".$id."|".$main_id."|".$ext_name."|".$slength."\n";
				my ($index) = 1;
				while($index < $slength) {
					my ($send) = $NUM_RESIDUES+$index-1;
					$send = $slength if ($send >= $slength);
					$output .= $seqobj->subseq($index,$send)."\n";
					$index = $send+1;	
				}
			}
		}		
    }
    else {
		throw('Argument must be an APPRIS::Transcript');
   	}
	return $output;
}

1;
