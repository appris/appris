=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::CDS - Object representing a cds

=head1 SYNOPSIS

  my $cds = APPRIS::CDS->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print cds information
  print("cds start:end:strand is "
      . join( ":", map { $cds->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $cds->stable_id('ENSE000001');

=head1 DESCRIPTION

A representation of a CDS within the APPRIS system.
This is a class which represents a coding region of 
a transcript (translation). 

=head1 METHODS

=cut

package APPRIS::CDS;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-start]  : 
       int - start postion of the cds
  Arg [-end]    : 
       int - end position of the cds
  Arg [-strand] : 
       char - '+','-' the strand the exon is on
  Arg [-phase] :
        string - the phase of the exon
  Example    : $cds = APPRIS::CDS->new(...);
  Description: Creates a new cds object
  Returntype : APPRIS::CDS
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
		$strand,		$phase
    )
    = rearrange( [
		'start',		'end',
		'strand',		'phase'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->phase($phase);
  
	return $self;
}

=head2 start

  Arg [1]    : (optional) Int - the start location to set
  Example    : $cds->start(123);
  Description: Getter/setter for start location for this cds
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
  Example    : $cds->end(123);
  Description: Getter/setter for end location for this cds
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
  Example    : $cds->strand('-');
  Description: Getter/setter for strand for this exon
  Returntype : Char
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strand {
  my $self = shift;
  $self->{'strand'} = shift if(@_);
  return $self->{'strand'};
}

=head2 phase

  Arg [1]    : (optional) Int - the phase of the exon
  Example    : $cds->phase(1);
  Description: Getter/setter for the phase of the exon
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub phase {
  my $self = shift;
  $self->{'phase'} = shift if(@_);
  return $self->{'phase'};
}

sub DESTROY {}

1;
