=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis - Object representing the set analysis

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis->new();

  # set some additional attributes
  $analysis->firestar($firestar);

=head1 DESCRIPTION

A representation of a Analysis within the APPRIS system.
A analysis object consists of a set of Analysis object that executed for a
analysis.

=head1 METHODS

=cut

package APPRIS::Analysis;

use strict;
use warnings;
use Data::Dumper;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-number]  : (optional)
       int - number of analysis
  Arg [-firestar]  : (optional)
       APPRIS::Analysis::Firestar - analysis object of Firestar method
  Arg [-matador3d]  : (optional)
       APPRIS::Analysis::Matador3D - analysis object of Matador3D method
  Arg [-spade]  : (optional)
       APPRIS::Analysis::SPADE - analysis object of SPADE method
  Arg [-inertia]  : (optional)
       APPRIS::Analysis::INERTIA - analysis object of INERTIA method
  Arg [-crash]  : (optional)
       APPRIS::Analysis::CRASH - analysis object of CRASH method
  Arg [-thump]  : (optional)
       APPRIS::Analysis::THUMP - analysis object of THUMP method
  Arg [-cexonic]  : (optional)
       APPRIS::Analysis::CExonic - analysis object of CExonic method
  Arg [-corsair]  : (optional)
       APPRIS::Analysis::CORSAIR - analysis object of CORSAIR method
  Arg [-proteo]  : (optional)
       APPRIS::Analysis::PROTEO - analysis object of PROTEO method
  Arg [-trifid]  : (optional)
       APPRIS::Analysis::TRIFID - analysis object of TRIFID method
  Arg [-appris]  : (optional)
       APPRIS::Analysis::APPRIS - analysis object of APPRIS method
  Example    : $analysis = APPRIS::Analysis->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis
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
		$number,
		$firestar,	
		$matador3d, $matador3d2,
		$spade,		$inertia,
		$crash,		$thump,
		$cexonic,	$corsair,
		$proteo,	$trifid,
		$appris
		
	)
	= rearrange( [
		'number',
		'firestar',	
		'matador3d', 'matador3d2',
		'spade',	'inertia',
		'crash',	'thump',
		'cexonic',	'corsair',
		'proteo',	'trifid',
		'appris'
	],
	@_
	);

	if(defined $number) { $self->number($number); }
	else { $self->number(0); }
	$self->firestar($firestar) if(defined $firestar);
	$self->matador3d($matador3d) if(defined $matador3d);
	$self->matador3d2($matador3d2) if(defined $matador3d2);
	$self->spade($spade) if(defined $spade);
	$self->inertia($inertia) if(defined $inertia);
	$self->crash($crash) if(defined $crash);
	$self->thump($thump) if(defined $thump);
	$self->cexonic($cexonic) if(defined $cexonic);
	$self->corsair($corsair) if(defined $corsair);
	$self->proteo($proteo) if(defined $proteo);
	$self->trifid($trifid) if(defined $trifid);
	$self->appris($appris) if(defined $appris);
	
	return $self;
}

=head2 number

  Arg [1]    : (optional) Int - the number of analysis to set
  Example    : $analysis->number(21);
  Description: Getter/setter for the number of analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub number {
	my ($self) = shift;
	$self->{'number'} = shift if(@_);
	return $self->{'number'};
}

=head2 firestar

  Arg [1]    : (optional) APPRIS::Analysis::Firestar - the firestar 
               object to set
  Example    : $analysis->firestar($method);
  Description: Getter/setter for the firestar object
  Returntype : APPRIS::Analysis::Firestar or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub firestar {
	my ($self) = shift;
	$self->{'firestar'} = shift if(@_);
	return $self->{'firestar'};
}

=head2 matador3d

  Arg [1]    : (optional) APPRIS::Analysis::Matador3D - the matador3d 
               object to set
  Example    : $analysis->matador3d($method);
  Description: Getter/setter for the matador3d object
  Returntype : APPRIS::Analysis::Matador3D or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub matador3d {
	my ($self) = shift;
	$self->{'matador3d'} = shift if(@_);
	return $self->{'matador3d'};
}

=head2 matador3d2

  Arg [1]    : (optional) APPRIS::Analysis::Matador3D - the matador3d 
               object to set
  Example    : $analysis->matador3d2($method);
  Description: Getter/setter for the matador3d object
  Returntype : APPRIS::Analysis::Matador3D or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub matador3d2 {
	my ($self) = shift;
	$self->{'matador3d2'} = shift if(@_);
	return $self->{'matador3d2'};
}

=head2 spade

  Arg [1]    : (optional) APPRIS::Analysis::SPADE - the spade 
               object to set
  Example    : $analysis->spade($method);
  Description: Getter/setter for the spade object
  Returntype : APPRIS::Analysis::SPADE or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub spade {
	my ($self) = shift;
	$self->{'spade'} = shift if(@_);
	return $self->{'spade'};
}

=head2 inertia

  Arg [1]    : (optional) APPRIS::Analysis::INERTIA - the inertia 
               object to set
  Example    : $analysis->inertia($method);
  Description: Getter/setter for the inertia object
  Returntype : APPRIS::Analysis::INERTIA or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub inertia {
	my ($self) = shift;
	$self->{'inertia'} = shift if(@_);
	return $self->{'inertia'};
}

=head2 crash

  Arg [1]    : (optional) APPRIS::Analysis::CRASH - the crash 
               object to set
  Example    : $analysis->crash($method);
  Description: Getter/setter for the crash object
  Returntype : APPRIS::Analysis::CRASH or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub crash {
	my ($self) = shift;
	$self->{'crash'} = shift if(@_);
	return $self->{'crash'};
}

=head2 thump

  Arg [1]    : (optional) APPRIS::Analysis::THUMP - the thump 
               object to set
  Example    : $analysis->thump($method);
  Description: Getter/setter for the thump object
  Returntype : APPRIS::Analysis::THUMP or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub thump {
	my ($self) = shift;
	$self->{'thump'} = shift if(@_);
	return $self->{'thump'};
}

=head2 cexonic

  Arg [1]    : (optional) APPRIS::Analysis::CExonic - the cexonic 
               object to set
  Example    : $analysis->cexonic($method);
  Description: Getter/setter for the cexonic object
  Returntype : APPRIS::Analysis::CExonic or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub cexonic {
	my ($self) = shift;
	$self->{'cexonic'} = shift if(@_);
	return $self->{'cexonic'};
}

=head2 corsair

  Arg [1]    : (optional) APPRIS::Analysis::CORSAIR - the corsair 
               object to set
  Example    : $analysis->corsair($method);
  Description: Getter/setter for the corsair object
  Returntype : APPRIS::Analysis::CORSAIR or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub corsair {
	my ($self) = shift;
	$self->{'corsair'} = shift if(@_);
	return $self->{'corsair'};
}

=head2 proteo

  Arg [1]    : (optional) APPRIS::Analysis::PROTEO - the proteo 
               object to set
  Example    : $analysis->proteo($method);
  Description: Getter/setter for the proteo object
  Returntype : APPRIS::Analysis::PROTEO or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub proteo {
	my ($self) = shift;
	$self->{'proteo'} = shift if(@_);
	return $self->{'proteo'};
}

=head2 trifid

  Arg [1]    : (optional) APPRIS::Analysis::TRIFID - the trifid
               object to set
  Example    : $analysis->trifid($method);
  Description: Getter/setter for the trifid object
  Returntype : APPRIS::Analysis::TRIFID or undef
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub trifid {
	my ($self) = shift;
	$self->{'trifid'} = shift if(@_);
	return $self->{'trifid'};
}

=head2 appris

  Arg [1]    : (optional) APPRIS::Analysis::APPRIS - the appris 
               object to set
  Example    : $analysis->appris($method);
  Description: Getter/setter for the appris object
  Returntype : APPRIS::Analysis::APPRIS or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub appris {
	my ($self) = shift;
	$self->{'appris'} = shift if(@_);
	return $self->{'appris'};
}

sub DESTROY {}

1;
