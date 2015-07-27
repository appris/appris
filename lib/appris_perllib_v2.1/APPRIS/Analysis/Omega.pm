=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::Omega - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::Omega->new(
    -result      => <Analysis result>,
    -functional_residue  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of Omega within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::Omega;

use strict;
use warnings;

use APPRIS::Analysis::OmegaRegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-average]  : 
       string - the omega average for this analysis
  Arg [-st_desviation]:
        string - the standard desviation for this analysis
  Arg [-unusual_evolution] :
        string - the standard desviation for this analysis
  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-regions]: (optional)
        Listref of APPRIS::Analysis::OmegaRegion - the 
        set of residues that was analysed
  Example    : $analysis = APPRIS::Analysis::Omega->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::Omega
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
		$average,				$st_desviation,
		$unusual_evolution,		$result,
		$regions
	)
	= rearrange( [
		'average',				'st_desviation',
		'unusual_evolution',	'result',
		'regions'				
	],
	@_
	);

 	$self->average($average);
 	$self->st_desviation($st_desviation);
 	$self->unusual_evolution($unusual_evolution);
 	$self->result($result);
	$self->regions($regions) if(defined $regions);
		
	return $self;
}

=head2 average

  Arg [1]    : (optional) String - the omega average to set
  Example    : $analysis->average(123);
  Description: Getter/setter for the omega average for this analysis
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub average {
	my ($self) = shift;
	$self->{'average'} = shift if(@_);
	return $self->{'average'};
}

=head2 st_desviation

  Arg [1]    : (optional) String - the omega standard 
               desviation to set
  Example    : $analysis->st_desviation(123);
  Description: Getter/setter for the omega standard 
               desviation for this analysis
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub st_desviation {
	my ($self) = shift;
	$self->{'st_desviation'} = shift if(@_);
	return $self->{'st_desviation'};
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


=head2 regions

  Arg [1]    : (optional) Listref of APPRIS::Analysis::OmegaRegion - 
               the set of residues that for this analysis 
  Example    : $analysis->regions($residues);
  Description: Getter/setter for the set of analysed residues 
  Returntype : Listref of APPRIS::Analysis::OmegaRegion or undef
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
