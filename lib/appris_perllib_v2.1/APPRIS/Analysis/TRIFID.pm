=head1 CONTACT

  For contact details see the L<APPRIS website|https://appris.bioinfo.cnio.es/>.

=cut

=head1 NAME

APPRIS::Analysis::TRIFID - Object representing a TRIFID analysis result

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::TRIFID->new(
    -trifid_score  => <Raw TRIFID score>,
    -norm_trifid_score  => <Normalized TRIFID score>
    ...
  );

=head1 DESCRIPTION

A representation of a TRIFID analysis result incorporated within
the APPRIS system. Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::TRIFID;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);

=head2 new

  Arg [-trifid_score]  :
       float - the raw TRIFID score of the isoform
  Arg [-norm_trifid_score]  :
       float - the normalized TRIFID score of the isoform
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
		$trifid_score,
		$norm_trifid_score
	)
	= rearrange( [
		'trifid_score',
		'norm_trifid_score'
	],
	@_
	);

	$self->trifid_score($trifid_score);
	$self->norm_trifid_score($norm_trifid_score);

	return $self;
}

=head2 trifid_score

  Arg [1]    : Float - the raw TRIFID score to set
  Example    : $analysis->trifid_score(0.9);
  Description: Getter/setter for the raw TRIFID score
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

=head2 norm_trifid_score

  Arg [1]    : Float - the normalized TRIFID score to set
  Example    : $analysis->norm_trifid_score(0.9);
  Description: Getter/setter for the normalized TRIFID score
  Returntype : Float
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub norm_trifid_score {
	my ($self) = shift;
	$self->{'norm_trifid_score'} = shift if(@_);
	return $self->{'norm_trifid_score'};
}


sub DESTROY {}

1;
