=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::INERTIA - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::INERTIA->new(
    -result      => <Analysis result>,
    -functional_residue  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of INERTIA within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::INERTIA;

use strict;
use warnings;

use APPRIS::Analysis::INERTIARegion;
use APPRIS::Analysis::Omega;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-unusual_evolution]:
        string - the appris annotation for this analysis
  Arg [-mafft_alignment]: (optional)
        APPRIS::Analysis::Omega - the omega object of MAFFT
        alignment
  Arg [-prank_alignment]: (optional)
        APPRIS::Analysis::Omega - the omega object of Prank
        alignment
  Arg [-kalign_alignment]: (optional)
        APPRIS::Analysis::Omega - the omega object of KAlign
        alignment
  Arg [-compara_alignment]: (optional)
        APPRIS::Analysis::Omega - the omega object of Compara
        alignment
  Arg [-regions]: (optional)
        Listref of APPRIS::Analysis::INERTIARegion - the 
        set of regions that was analysed        
  Example    : $analysis = APPRIS::Analysis::INERTIA->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::INERTIA
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
		$result,				$unusual_evolution,
		$mafft_alignment,		$prank_alignment,		$kalign_alignment,		$compara_alignment,
		$regions
	)
	= rearrange( [
		'result',				'unusual_evolution',
		'mafft_alignment',		'prank_alignment',		'kalign_alignment',		'compara_alignment',
		'regions'
	],
	@_
	);

	$self->result($result);
 	$self->unusual_evolution($unusual_evolution);
	$self->mafft_alignment($mafft_alignment) if(defined $mafft_alignment);
	$self->prank_alignment($prank_alignment) if(defined $prank_alignment);
	$self->kalign_alignment($kalign_alignment) if(defined $kalign_alignment);
	$self->compara_alignment($compara_alignment) if(defined $compara_alignment);
	$self->regions($regions) if(defined $regions);
	
	return $self;
}

=head2 result

  Arg [1]    : (optional) String - the result to set
  Example    : $analysis->result(123);
  Description: Getter/setter for the results for this analysis
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub result {
	my ($self) = shift;
	$self->{'result'} = shift if(@_);
	return $self->{'result'};
}

=head2 unusual_evolution

  Arg [1]    : (optional) String - the analysed unusual_evolution to set
  Example    : $analysis->unusual_evolution(123);
  Description: Getter/setter for the analysed unusual_evolution
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub unusual_evolution {
	my ($self) = shift;
	$self->{'unusual_evolution'} = shift if(@_);
	return $self->{'unusual_evolution'};
}

=head2 mafft_alignment

  Arg [1]    : (optional) APPRIS::Analysis::Omega - the omega object 
               of MAFFT alignment that for this analysis 
  Example    : $analysis->mafft_alignment($regions);
  Description: Getter/setter for the omega object of MAFFT alignment 
  Returntype : APPRIS::Analysis::Omega or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub mafft_alignment {
	my ($self) = shift;
	$self->{'mafft_alignment'} = shift if(@_);
	return $self->{'mafft_alignment'};
}

=head2 prank_alignment

  Arg [1]    : (optional) APPRIS::Analysis::Omega - the omega object 
               of Prank alignment that for this analysis 
  Example    : $analysis->prank_alignment($regions);
  Description: Getter/setter for the omega object of Prank alignment 
  Returntype : APPRIS::Analysis::Omega or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub prank_alignment {
	my ($self) = shift;
	$self->{'prank_alignment'} = shift if(@_);
	return $self->{'prank_alignment'};
}

=head2 kalign_alignment

  Arg [1]    : (optional) APPRIS::Analysis::Omega - the omega object 
               of KAlign alignment that for this analysis 
  Example    : $analysis->kalign_alignment($regions);
  Description: Getter/setter for the omega object of KAlign alignment 
  Returntype : APPRIS::Analysis::Omega or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub kalign_alignment {
	my ($self) = shift;
	$self->{'kalign_alignment'} = shift if(@_);
	return $self->{'kalign_alignment'};
}

=head2 compara_alignment

  Arg [1]    : (optional) APPRIS::Analysis::Omega - the omega object 
               of KAlign alignment that for this analysis 
  Example    : $analysis->compara_alignment($regions);
  Description: Getter/setter for the omega object of Compara alignment 
  Returntype : APPRIS::Analysis::Omega or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub compara_alignment {
	my ($self) = shift;
	$self->{'compara_alignment'} = shift if(@_);
	return $self->{'compara_alignment'};
}

=head2 regions

  Arg [1]    : (optional) Listref of APPRIS::Analysis::INERTIARegion - 
               the set of regions that for this analysis 
  Example    : $analysis->regions($regions);
  Description: Getter/setter for the set of analysed regions 
  Returntype : Listref of APPRIS::Analysis::INERTIARegion or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub regions {
	my ($self) = shift;
	$self->{'regions'} = shift if(@_);
	return $self->{'regions'};
}

sub DESTROY {}

1;
