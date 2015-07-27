=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::THUMPRegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::THUMPRegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a THUMPRegion within the APPRIS system.
A THUMPRegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::THUMPRegion;

use strict;
use warnings;
use vars qw(@ISA);

use APPRIS::Analysis::Region;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

@ISA = qw(APPRIS::Analysis::Region);

=head2 new

  Arg [-start]  : 
       int - start postion of the slice
  Arg [-end]    : 
       int - end position of the slice
  Arg [-strand] : 
       char - '+','-' the strand the slice is on
  Arg [-pstart]  : 
       int - start postion of the peptide region
  Arg [-pend]    : 
       int - end position of the peptide region
  Arg [-damaged] : 
       int - '0','1' if the regions is damaged transmembrane helix
  Example    : $region = APPRIS::Analysis::THUMPRegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::THUMPRegion
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
		$start,		$end,		$strand,
		$pstart,	$pend,
		$damaged
	)
	= rearrange( [
		'start',	'end',		'strand',
		'pstart',	'pend',
		'damaged'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->pstart($pstart);
 	$self->pend($pend);
 	$self->damaged($damaged);
  
	return $self;
}

=head2 pstart

  Arg [1]    : (optional) Int - the start of peptide location to set
  Example    : $region->pstart(123);
  Description: Getter/setter for start of hit location
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pstart {
	my ($self) = shift;
	$self->{'pstart'} = shift if(@_);
	return $self->{'pstart'};
}

=head2 pend

  Arg [1]    : (optional) Int - the end of peptide location to set
  Example    : $region->pend(123);
  Description: Getter/setter for end of hit location
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pend {
	my ($self) = shift;
	$self->{'pend'} = shift if(@_);
	return $self->{'pend'};
}

=head2 damaged

  Arg [1]    : (optional) Int - the flag that says if it is 
               damaged this region
  Example    : $region->damaged(1);
  Description: Getter/setter for the damaged flag
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub damaged {
	my ($self) = shift;
	$self->{'damaged'} = shift if(@_);
	return $self->{'damaged'};
}

sub DESTROY {}

1;
