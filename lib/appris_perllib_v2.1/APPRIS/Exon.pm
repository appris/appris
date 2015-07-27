=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Exon - Object representing a exon

=head1 SYNOPSIS

  my $exon = APPRIS::Exon->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print exon information
  print("exon start:end:strand is "
      . join( ":", map { $exon->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $exon->stable_id('ENSE000001');

=head1 DESCRIPTION

A representation of a Exon within the APPRIS system.
This is a class which represents an exon which is part of a transcript.

=head1 METHODS

=cut

package APPRIS::Exon;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-start]  : 
       int - start postion of the exon
  Arg [-end]    : 
       int - end position of the exon
  Arg [-strand] : 
       char - '+','-' the strand the exon is on
  Arg [-stable_id] :
        string - the stable identifier of this exon
  Example    : $exon = APPRIS::Exon->new(...);
  Description: Creates a new exon object
  Returntype : APPRIS::Exon
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
		$start,			$end,
		$strand,		$stable_id
    )
    = rearrange( [
		'start',		'end',
		'strand',		'stable_id'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->stable_id($stable_id);
  
	return $self;
}

=head2 start

  Arg [1]    : (optional) Int - the start location to set
  Example    : $exon->start(123);
  Description: Getter/setter for start location for this exon
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my $self = shift;
  $self->{'start'} = shift if(@_);
  return $self->{'start'};
}

=head2 end

  Arg [1]    : (optional) Int - the end location to set
  Example    : $exon->end(123);
  Description: Getter/setter for end location for this exon
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my $self = shift;
  $self->{'end'} = shift if(@_);
  return $self->{'end'};
}

=head2 strand

  Arg [1]    : (optional) Char - the strand the exon is on
  Example    : $exon->strand('-');
  Description: Getter/setter for strand for this exon
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strand {
  my $self = shift;
  $self->{'strand'} = shift if(@_);
  return $self->{'strand'};
}

=head2 stable_id

  Arg [1]    : (optional) String - the stable ID to set
  Example    : $exon->stable_id("ENSE0000000001");
  Description: Getter/setter for stable id for this exon
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

sub DESTROY {}

1;
