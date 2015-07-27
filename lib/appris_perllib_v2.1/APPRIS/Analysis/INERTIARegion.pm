=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::INERTIARegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::INERTIARegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a INERTIARegion within the APPRIS system.
A INERTIARegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::INERTIARegion;

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
  Arg [-unusual_evolution]:
        string - the appris annotation for this analysis
  Example    : $region = APPRIS::Analysis::INERTIARegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::INERTIARegion
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
		$unusual_evolution
	)
	= rearrange( [
		'start',	'end',		'strand',	
		'unusual_evolution'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->unusual_evolution($unusual_evolution);
  
	return $self;
}

=head2 unusual_evolution

  Arg [1]    : (optional) String - the analysed unusual_evolution to set
  Example    : $analysis->unusual_evolution(123);
  Description: Getter/setter for the analysed unusual_evolution
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub unusual_evolution {
	my ($self) = shift;
	$self->{'unusual_evolution'} = shift if(@_);
	return $self->{'unusual_evolution'};
}

sub DESTROY {}

1;
