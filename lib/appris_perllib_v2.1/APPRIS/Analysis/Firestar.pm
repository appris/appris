=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::Firestar - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::Firestar->new(
    -result      => <Analysis result>,
    -functional_residue  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of Firestar within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::Firestar;

use strict;
use warnings;

use APPRIS::Analysis::FirestarRegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-functional_residue]:
        string - the appris annotation for this analysis
  Arg [-num_residues] : (optional)
        int - the number of analysed residues
  Arg [-residues]: (optional)
        Listref of APPRIS::Analysis::FirestarRegion - the 
        set of residues that was analysed
  Example    : $analysis = APPRIS::Analysis::Firestar->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::Firestar
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
		$result,			$functional_residue,
		$num_residues,		$residues
	)
	= rearrange( [
		'result',			'functional_residue',
		'num_residues',		'residues'
	],
	@_
	);

 	$self->result($result);
 	$self->functional_residue($functional_residue);
	$self->num_residues($num_residues) if(defined $num_residues);
	$self->residues($residues) if(defined $residues);
		
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

=head2 functional_residue

  Arg [1]    : (optional) String - the analysed functional_residue to set
  Example    : $analysis->functional_residue(123);
  Description: Getter/setter for the analysed functional_residue
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub functional_residue {
	my ($self) = shift;
	$self->{'functional_residue'} = shift if(@_);
	return $self->{'functional_residue'};
}

=head2 num_residues

  Arg [1]    : (optional) Int - the number of residues for 
               the analysis
  Example    : $analysis->num_residues(5);
  Description: Getter/setter for the number of residues 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_residues {
	my ($self) = shift;
	$self->{'num_residues'} = shift if(@_);
	return $self->{'num_residues'};
}

=head2 residues

  Arg [1]    : (optional) Listref of APPRIS::Analysis::FirestarRegion - 
               the set of residues that for this analysis 
  Example    : $analysis->residues($residues);
  Description: Getter/setter for the set of analysed residues 
  Returntype : Listref of APPRIS::Analysis::FirestarRegion or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub residues {
	my ($self) = shift;
	$self->{'residues'} = shift if(@_);
	return $self->{'residues'};
}

sub DESTROY {}

1;
