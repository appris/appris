=head1 CONTACT

  For contact details see the L<APPRIS website|http://appris-tools.org>.

=cut

=head1 NAME

APPRIS::Analysis::TRIFID - Object representing a TRIFID analysis result

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::TRIFID->new(
    -trifid_score  => <Annotation analysed>
    ...
  );

=head1 DESCRIPTION

A representation of a TRIFID analysis result incorporated within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::TRIFID;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);

=head2 new

  Arg [-trifid_score]  :
       int - the TRIFID score of the transcript
  Example    : $analysis = APPRIS::Analysis::TRIFID->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::TRIFID
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
	my ($caller) = shift;

	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;

	my (
		$trifid_score
	)
	= rearrange( [
		'trifid_score'
	],
	@_
	);

	$self->trifid_score($trifid_score);

	return $self;
}

=head2 trifid_score

  Arg [1]    : Float - the trifid_score to set
  Example    : $analysis->trifid_score(0.9);
  Description: Getter/setter for the analysed trifid_score
  Returntype : Float
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub trifid_score {
	my ($self) = shift;
	$self->{'trifid_score'} = shift if(@_);
	return $self->{'trifid_score'};
}

sub DESTROY {}

1;
