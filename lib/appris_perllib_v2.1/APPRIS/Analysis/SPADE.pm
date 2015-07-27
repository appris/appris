=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::SPADE - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::SPADE->new(
    -result      => <Analysis result>,
    -domain_signal  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of SPADE within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::SPADE;

use strict;
use warnings;

use APPRIS::Analysis::SPADERegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-domain_signal]:
        string - the appris annotation for this analysis
  Arg [-num_domains] : (optional)
        int - the number of analysed domains
  Arg [-num_possibly_damaged_domains] : (optional)
        int - the number of possibly damaged domains
  Arg [-num_damaged_domains] : (optional)
        int - the number of damaged domains
  Arg [-num_wrong_domains] : (optional)
        int - the number of wrong domains
  Arg [-regions]: (optional)
        Listref of APPRIS::Analysis::SPADERegion - the 
        set of regions that was analysed
  Example    : $analysis = APPRIS::Analysis::SPADE->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::SPADE
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
		$result,						$domain_signal,
		$num_domains,					$num_damaged_domains,
		$num_possibly_damaged_domains,	$num_wrong_domains,
		$regions
	)
	= rearrange( [
		'result',						'domain_signal',
		'num_domains',					'num_damaged_domains',
		'num_possibly_damaged_domains',	'num_wrong_domains',
		'regions'
	],
	@_
	);

 	$self->result($result);
 	$self->domain_signal($domain_signal);
	$self->num_domains($num_domains) if(defined $num_domains);
	$self->num_possibly_damaged_domains($num_possibly_damaged_domains) if(defined $num_possibly_damaged_domains);
	$self->num_damaged_domains($num_damaged_domains) if(defined $num_damaged_domains);
	$self->num_wrong_domains($num_wrong_domains) if(defined $num_wrong_domains);
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

=head2 domain_signal

  Arg [1]    : (optional) String - the analysed domain_signal to set
  Example    : $analysis->domain_signal(123);
  Description: Getter/setter for the analysed domain_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub domain_signal {
	my ($self) = shift;
	$self->{'domain_signal'} = shift if(@_);
	return $self->{'domain_signal'};
}

=head2 num_domains

  Arg [1]    : (optional) Int - the number of domains for 
               the analysis
  Example    : $analysis->num_domains(5);
  Description: Getter/setter for the number of domains 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_domains {
	my ($self) = shift;
	$self->{'num_domains'} = shift if(@_);
	return $self->{'num_domains'};
}

=head2 num_possibly_damaged_domains

  Arg [1]    : (optional) Int - the number of possibly damaged domains for 
               the analysis
  Example    : $analysis->num_possibly_damaged_domains(5);
  Description: Getter/setter for the number of possibly damaged domains 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_possibly_damaged_domains {
	my ($self) = shift;
	$self->{'num_possibly_damaged_domains'} = shift if(@_);
	return $self->{'num_possibly_damaged_domains'};
}

=head2 num_damaged_domains

  Arg [1]    : (optional) Int - the number of damaged domains for 
               the analysis
  Example    : $analysis->num_damaged_domains(5);
  Description: Getter/setter for the number of damaged domains 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_damaged_domains {
	my ($self) = shift;
	$self->{'num_damaged_domains'} = shift if(@_);
	return $self->{'num_damaged_domains'};
}

=head2 num_wrong_domains

  Arg [1]    : (optional) Int - the number of wrong domains for 
               the analysis
  Example    : $analysis->num_wrong_domains(5);
  Description: Getter/setter for the number of wrong domains 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_wrong_domains {
	my ($self) = shift;
	$self->{'num_wrong_domains'} = shift if(@_);
	return $self->{'num_wrong_domains'};
}

=head2 regions

  Arg [1]    : (optional) Listref of APPRIS::Analysis::SPADERegion - 
               the set of regions that for this analysis 
  Example    : $analysis->regions($regions);
  Description: Getter/setter for the set of analysed regions 
  Returntype : Listref of APPRIS::Analysis::SPADERegion or undef
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
