=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Codon - Object representing a exon

=head1 SYNOPSIS

  my $exon = APPRIS::Codon->new(
    -type  => start,
    -start  => 123,
    -end    => 1045,
    -strand => '+',
    -phase => 0,
  );

  # print exon information
  print("exon start:end:strand is "
      . join( ":", map { $exon->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $exon->stable_id('ENSE000001');

=head1 DESCRIPTION

A representation of a Codon within the APPRIS system.
This is a class which represents an codon which is part of a transcript.

=head1 METHODS

=cut

package APPRIS::Codon;

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
  Arg [-phase] :
        string - the phase of the exon
  Arg [-type] :
        string - type of codon (start or end)
  Example    : $exon = APPRIS::Codon->new(...);
  Description: Creates a new codon object
  Returntype : APPRIS::Codon
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
		$strand,		$phase,
		$type
    )
    = rearrange( [
		'start',		'end',
		'strand',		'phase',
		'type'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->phase($phase); 	
 	$self->type($type);
  
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

=head2 type

  Arg [1]    : (optional) String - the type of codon to set
  Example    : $exon->type("start");
  Description: Getter/setter for the type for this codon
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  $self->{'type'} = shift if(@_);
  return $self->{'type'};
}

sub DESTROY {}

1;
