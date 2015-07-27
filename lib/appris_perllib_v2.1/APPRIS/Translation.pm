=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Translation - Object representing a peptide

=head1 SYNOPSIS

  my $peptide = APPRIS::Translation->new(
    -stable_id	=> 'ENST000001',
    -protein_id	=> 'ENSP000001',
    -sequence => <AminoAcid sequence>,
  );

  # set some additional attributes
  $peptide->stable_id('ENSP000001');
  $peptide->description('This is the peptide description');

=head1 DESCRIPTION

A representation of a Translation within the APPRIS system.
A Translation object defines the CDS of a Transcript and the
aminoacid sequence of its protein.

=head1 METHODS

=cut

package APPRIS::Translation;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::ProCDS qw(get_protein_cds_sequence get_contained_cds);

=head2 new

  Arg [-stable_id] : (optional)
        string - the stable identifier of this protein
  Arg [-protein_id] : (optional)
        string - the stable identifier of this protein
  Arg [-description]: (optional)
        string - the peptide description
  Arg [-sequence] : (optional)
        string - the nucleotide sequence associated with this peptide
  Arg [-cds] : (optional)
        Listref of APPRIS::CDS - the set of CDS that for 
        this translation
  Arg [-codons] : (optional)
        Listref of APPRIS::Codon - the set of Codon that for 
        this translation
  Example    : $peptide = APPRIS::Translation->new(...);
  Description: Creates a new peptide object
  Returntype : APPRIS::Translation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
	my ($caller) = shift;

	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;
	
	my (
		$stable_id,		$protein_id,	
		$description,
		$sequence,		$cds,
		$cds_sequence,	$codons
    )
    = rearrange( [
		'stable_id',	'protein_id',
		'description',		
		'sequence',		'cds',
		'cds_sequence',	'codons'
	],
	@_
	);

 	$self->stable_id($stable_id) if(defined $stable_id);
 	$self->protein_id($protein_id) if(defined $protein_id);
	$self->description($description) if(defined $description);
	$self->sequence($sequence) if(defined $sequence);
	$self->cds($cds) if(defined $cds);
	$self->codons($codons) if(defined $codons);
	$self->cds_sequence($self) if(defined $cds and defined $sequence);
	$self->get_overlapping_cds();
	  
  return $self;
}

=head2 stable_id

  Arg [1]    : (optional) String - the stable ID to set
  Example    : $peptide->stable_id("ENSP0000000001");
  Description: Getter/setter for stable id for this peptide
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if(@_);
  return $self->{'stable_id'};
}

=head2 protein_id

  Arg [1]    : (optional) String - the stable ID to set
  Example    : $peptide->protein_id("ENSP0000000001");
  Description: Getter/setter for stable id for this peptide
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub protein_id {
  my $self = shift;
  $self->{'protein_id'} = shift if(@_);
  return $self->{'protein_id'};
}

=head2 description

  Arg [1]    : (optional) String - the description to set
  Example    : $peptide->description('This is the peptide\'s description');
  Description: Getter/setter for peptide description
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if( @_ );
    return $self->{'description'};
}

=head2 sequence

  Arg [1]    : (optional) String - the aminoacid sequence 
               that for this peptide
  Example    : $peptide->sequence();
  Description: Getter/setter for aminoacid sequence for this 
               peptide
  Returntype : String or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sequence {
  my  $self  = shift;

  $self->{'sequence'} = shift if(@_);
  return $self->{'sequence'};
}

=head2 cds

  Arg [1]    : (optional) Listref of APPRIS::CDS - the set of 
               CDS that for this translation
  Example    : $peptide->cds();
  Description: Getter/setter for the set of CDS that 
               for this translation
  Returntype : Listref of APPRIS::CDS or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub cds {
  my  $self  = shift;

  $self->{'cds'} = shift if(@_);
  return $self->{'cds'};
}

=head2 codons

  Arg [1]    : (optional) Listref of APPRIS::Codon - the set of 
               Codon that for this translation
  Example    : $peptide->codons();
  Description: Getter/setter for the set of codons that 
               for this translation
  Returntype : Listref of APPRIS::Codon or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub codons {
  my  $self  = shift;

  $self->{'codons'} = shift if(@_);
  return $self->{'codons'};
}

=head2 cds_sequence

  Arg [1]    : (optional) APPRIS::Translation - peptide object
  Example    : $peptide->cds_sequence();
  Description: Getter for the set of sequence of
               CDS that for this translation
  Returntype : Listref of APPRIS::ProCDS or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub cds_sequence {
  my  $self  = shift;

  if(@_) {
    $self->{'cds_sequence'} = get_protein_cds_sequence($self->cds, $self->sequence);
  }
  #else {
  #  $self->{'cds_sequence'} = get_protein_cds_sequence($self->cds, $self->sequence);
  #}
  return $self->{'cds_sequence'};
}

=head2 get_overlapping_cds

  Arg [1]    : Int - the start location of looking region
  Arg [2]    : Int - the end location of looking region
  Example    : $peptide->get_overlapping_cds($start,$end);
  Description: Getter for the set of CDS that are within
               given region
  Returntype : Listref of APPRIS::CDS or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_overlapping_cds {
  my  $self  = shift;

  if(@_) {
    my ($start) = shift;
    my ($end) = shift;
    $self->{'overlapping_cds'} = get_contained_cds($self->cds, $start, $end);
  }
  else {
  	$self->{'overlapping_cds'} = undef;
  }
  return $self->{'overlapping_cds'};
}

sub DESTROY {}

1;
