=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::CExonic - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::CExonic->new(
    -result                => <Analysis result>,
    -conservation_exon  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of CExonic within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::CExonic;

use strict;
use warnings;

use APPRIS::Analysis::CExonicRegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-conservation_exon]:
        string - the appris annotation for this analysis
  Arg [-num_introns] : (optional)
        int - the number of introns for this analysis
  Arg [-first_specie_num_exons] : (optional)
        int - the number of exons of the first specie 
        for this analysis
  Arg [-second_specie_num_exons] : (optional)
        int - the number of exons of the second specie 
        for this analysis
  Arg [-regions]: (optional)
        Listref of APPRIS::Analysis::CExonicRegion - the 
        set of regions that was analysed
  Example    : $analysis = APPRIS::Analysis::CExonic->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::CExonic
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
		$result,						$conservation_exon,
		$num_introns,					$first_specie_num_exons,
		$second_specie_num_exons,		$regions
	)
	= rearrange( [
		'result',						'conservation_exon',
		'num_introns',					'first_specie_num_exons',
		'second_specie_num_exons',		'regions'
	],
	@_
	);

 	$self->result($result);
 	$self->conservation_exon($conservation_exon);
	$self->num_introns($num_introns) if(defined $num_introns);
	$self->first_specie_num_exons($first_specie_num_exons) if(defined $first_specie_num_exons);
	$self->second_specie_num_exons($second_specie_num_exons) if(defined $second_specie_num_exons);
	$self->regions($regions) if(defined $regions);
		
	return $self;
}

=head2 result

  Arg [1]    : (optional) String - the result to set
  Example    : $analysis->result(123);
  Description: Getter/setter for the results for this analysis
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub result {
	my ($self) = shift;
	$self->{'result'} = shift if(@_);
	return $self->{'result'};
}

=head2 conservation_exon

  Arg [1]    : (optional) String - the analysed conservation_exon to set
  Example    : $analysis->conservation_exon(123);
  Description: Getter/setter for the analysed conservation_exon
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub conservation_exon {
	my ($self) = shift;
	$self->{'conservation_exon'} = shift if(@_);
	return $self->{'conservation_exon'};
}

=head2 num_introns

  Arg [1]    : (optional) Int - the num introns of 
               the transcript is on
  Example    : $region->num_introns(1);
  Description: Getter/setter for the num introns for this region
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_introns {
	my ($self) = shift;
	$self->{'num_introns'} = shift if(@_);
	return $self->{'num_introns'};
}

=head2 first_specie_num_exons

  Arg [1]    : (optional) Int - the num exons of 
               the first alignment is on
  Example    : $region->first_specie_num_exons(1);
  Description: Getter/setter for the num exons of 
               the first alignment
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub first_specie_num_exons {
	my ($self) = shift;
	$self->{'first_specie_num_exons'} = shift if(@_);
	return $self->{'first_specie_num_exons'};
}

=head2 second_specie_num_exons

  Arg [1]    : (optional) Int - the num exons of 
               the second alignment is on
  Example    : $region->second_specie_num_exons(1);
  Description: Getter/setter for the num exons of 
               the second alignment
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub second_specie_num_exons {
	my ($self) = shift;
	$self->{'second_specie_num_exons'} = shift if(@_);
	return $self->{'second_specie_num_exons'};
}

=head2 regions

  Arg [1]    : (optional) Listref of APPRIS::Analysis::CExonicRegion - 
               the set of regions that for this analysis 
  Example    : $analysis->regions($regions);
  Description: Getter/setter for the set of analysed regions 
  Returntype : Listref of APPRIS::Analysis::CExonicRegion or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub regions {
	my ($self) = shift;
	$self->{'regions'} = shift if(@_);
	return $self->{'regions'};
}

sub DESTROY {}

1;
