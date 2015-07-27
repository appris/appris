=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::ProCDS - Object representing a cds

=head1 SYNOPSIS

  my $cds = APPRIS::ProCDS->new(
    -start  => 123,
    -end    => 1045,  
    -start_phase  => 1,
    -end_phase    => 2,
    -sequence => <AminoAcid sequence>,    
  );

  # print cds information
  print("cds start:end:strand is "
      . join( ":", map { $cds->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a protein from the CDS within 
the APPRIS system.
This is a class which represents a coding region of 
a transcript (translation). 

=head1 METHODS

=cut

package APPRIS::ProCDS;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-start]  : 
       int - start postion of the protein
  Arg [-end]    : 
       int - end position of the protein
  Arg [-start_phase]  : (optional)
       int - start phase of the cds
  Arg [-end_phase]    : (optional)
       int - end phase of the cds
  Arg [-sequence] : (optional)
       string - sequence for this cds
  Example    : $cds = APPRIS::ProCDS->new(...);
  Description: Creates a new peptide cds object
  Returntype : APPRIS::ProCDS
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
		$start,				$end,
		$start_phase,		$end_phase,
		$sequence
    )
    = rearrange( [
		'start',			'end',
		'start_phase',		'end_phase',
		'sequence'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->start_phase($start_phase) if(defined $start_phase);
 	$self->end_phase($end_phase) if(defined $end_phase);
 	$self->sequence($sequence) if(defined $sequence);
  
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

=head2 start_phase

  Arg [1]    : (optional) Int - the start phase to set
  Example    : $cds->start_phase(1);
  Description: Getter/setter for start phase for this cds
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start_phase {
  my $self = shift;
  $self->{'start_phase'} = shift if(@_);
  return $self->{'start_phase'};
}

=head2 end_phase

  Arg [1]    : (optional) Int - the end phase to set
  Example    : $cds->end_phase(1);
  Description: Getter/setter for end phase for this cds
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end_phase {
  my $self = shift;
  $self->{'end_phase'} = shift if(@_);
  return $self->{'end_phase'};
}

=head2 sequence

  Arg [1]    : (optional) String - the sequence of the cds
  Example    : $cds->sequence();
  Description: Getter/setter for sequence for this cds
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sequence {
  my $self = shift;
  $self->{'sequence'} = shift if(@_);
  return $self->{'sequence'};
}

sub DESTROY {}

1;
