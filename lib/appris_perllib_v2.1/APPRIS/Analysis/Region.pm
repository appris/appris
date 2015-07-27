=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::Region - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::Region->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a Region within the APPRIS system.
A Region object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::Region;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-start]  : 
       int - start postion of the slice
  Arg [-end]    : 
       int - end position of the slice
  Arg [-strand] : 
       char - '+','-' the strand the slice is on
  Example    : $region = APPRIS::Analysis::Region->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::Region
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
		$start,		$end,
		$strand
	)
	= rearrange( [
		'start',	'end',
		'strand'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
  
	return $self;
}

=head2 start

  Arg [1]    : (optional) Int - the start location to set
  Example    : $slice->start(123);
  Description: Getter/setter for start location for this slice
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
	my ($self) = shift;
	$self->{'start'} = shift if(@_);
	return $self->{'start'};
}

=head2 end

  Arg [1]    : (optional) Int - the end location to set
  Example    : $slice->end(123);
  Description: Getter/setter for end location for this slice
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
	my ($self) = shift;
	$self->{'end'} = shift if(@_);
	return $self->{'end'};
}

=head2 strand

  Arg [1]    : (optional) Char - the strand the slice is on
  Example    : $slice->strand('-');
  Description: Getter/setter for strand for this slice
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strand {
	my ($self) = shift;
	$self->{'strand'} = shift if(@_);
	return $self->{'strand'};
}

sub DESTROY {}

1;
