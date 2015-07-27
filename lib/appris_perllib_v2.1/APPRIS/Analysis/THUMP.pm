=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::THUMP - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::THUMP->new(
    -result                => <Analysis result>,
    -transmembrane_signal  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of THUMP within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::THUMP;

use strict;
use warnings;

use APPRIS::Analysis::THUMPRegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-transmembrane_signal]:
        string - the appris annotation for this analysis
  Arg [-num_tmh] : (optional)
        int - the number of analysed transmembrane helixes
  Arg [-num_damaged_tmh] : (optional)
        int - the number of damaged transmembrane helixes
  Arg [-regions]: (optional)
        Listref of APPRIS::Analysis::THUMPRegion - the 
        set of regions that was analysed
  Example    : $analysis = APPRIS::Analysis::THUMP->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::THUMP
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
		$result,		$transmembrane_signal,
		$num_tmh,		$num_damaged_tmh,
		$regions
	)
	= rearrange( [
		'result',		'transmembrane_signal',
		'num_tmh',		'num_damaged_tmh',
		'regions'
	],
	@_
	);

 	$self->result($result);
 	$self->transmembrane_signal($transmembrane_signal);
	$self->num_tmh($num_tmh) if(defined $num_tmh);
	$self->num_damaged_tmh($num_damaged_tmh) if(defined $num_damaged_tmh);
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

=head2 transmembrane_signal

  Arg [1]    : (optional) String - the analysed transmembrane_signal to set
  Example    : $analysis->transmembrane_signal(123);
  Description: Getter/setter for the analysed transmembrane_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transmembrane_signal {
	my ($self) = shift;
	$self->{'transmembrane_signal'} = shift if(@_);
	return $self->{'transmembrane_signal'};
}

=head2 num_tmh

  Arg [1]    : (optional) Int - the number of transmembrane helixes for 
               the analysis
  Example    : $analysis->num_tmh(5);
  Description: Getter/setter for the number of transmembrane helixes 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_tmh {
	my ($self) = shift;
	$self->{'num_tmh'} = shift if(@_);
	return $self->{'num_tmh'};
}

=head2 num_damaged_tmh

  Arg [1]    : (optional) Int - the number of damaged transmembrane 
               helixes for the analysis
  Example    : $analysis->num_damaged_tmh(5);
  Description: Getter/setter for the number of damaged transmembrane 
               helixes that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_damaged_tmh {
	my ($self) = shift;
	$self->{'num_damaged_tmh'} = shift if(@_);
	return $self->{'num_damaged_tmh'};
}

=head2 regions

  Arg [1]    : (optional) Listref of APPRIS::Analysis::THUMPRegion - 
               the set of regions that for this analysis 
  Example    : $analysis->regions($regions);
  Description: Getter/setter for the set of analysed regions 
  Returntype : Listref of APPRIS::Analysis::THUMPRegion or undef
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
