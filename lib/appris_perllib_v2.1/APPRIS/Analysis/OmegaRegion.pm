=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::OmegaRegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::OmegaRegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a OmegaRegion within the APPRIS system.
A OmegaRegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::OmegaRegion;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-start]  : 
       int - start postion of the CDS
  Arg [-end]    : 
       int - end position of the CDS
  Arg [-omega_mean] : 
       string - the omega mean of the region
  Arg [-st_deviation] : 
       string - the standard desviation of the region
  Arg [-p_value] : 
       string - the p_value of the region
  Arg [-difference_value] : 
       string - the difference value of the region
  Arg [-unusual_evolution]:
        string - the appris annotation for this analysis
  Example    : $region = APPRIS::Analysis::OmegaRegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::OmegaRegion
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
		$omega_mean,	$st_deviation,
		$p_value,		$difference_value,
		$unusual_evolution
	)
	= rearrange( [
		'start',		'end',
		'omega_mean',	'st_deviation',
		'p_value',		'difference_value',
		'unusual_evolution'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->omega_mean($omega_mean);
 	$self->st_deviation($st_deviation);
 	$self->p_value($p_value);
 	$self->difference_value($difference_value);
 	$self->unusual_evolution($unusual_evolution);
  
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

=head2 omega_mean

  Arg [1]    : (optional) String - the omega_mean of the region is on
  Example    : $region->omega_mean('-');
  Description: Getter/setter for the omega_mean that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub omega_mean {
	my ($self) = shift;
	$self->{'omega_mean'} = shift if(@_);
	return $self->{'omega_mean'};
}

=head2 st_deviation

  Arg [1]    : (optional) String - the st_deviation of the region is on
  Example    : $region->st_deviation('-');
  Description: Getter/setter for the st_deviation that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub st_deviation {
	my ($self) = shift;
	$self->{'st_deviation'} = shift if(@_);
	return $self->{'st_deviation'};
}

=head2 p_value

  Arg [1]    : (optional) String - the p_value of the region is on
  Example    : $region->p_value('-');
  Description: Getter/setter for the p_value that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub p_value {
	my ($self) = shift;
	$self->{'p_value'} = shift if(@_);
	return $self->{'p_value'};
}

=head2 difference_value

  Arg [1]    : (optional) String - the difference_value of the region is on
  Example    : $region->difference_value('-');
  Description: Getter/setter for the difference_value that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub difference_value {
	my ($self) = shift;
	$self->{'difference_value'} = shift if(@_);
	return $self->{'difference_value'};
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
