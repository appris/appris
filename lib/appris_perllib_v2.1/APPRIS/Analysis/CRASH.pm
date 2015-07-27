=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::CRASH - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::CRASH->new(
    -result      => <Analysis result>,
    -peptide_signal  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of CRASH within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::CRASH;

use strict;
use warnings;

use APPRIS::Analysis::CRASHRegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-peptide_signal]:
        string - the appris annotation for this analysis
  Arg [-mitochondrial_signal]:
        string - the appris annotation for this analysis
  Arg [-num_regions] : (optional)
        int - the number of analysed residues
  Arg [-peptide]: (optional)
        Listref of APPRIS::Analysis::CRASHRegion - the 
        set of regions that was analysed
  Example    : $analysis = APPRIS::Analysis::CRASH->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::CRASH
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
		$result,
		$peptide_signal,	$mitochondrial_signal,
		$sp_score,			$tp_score,
		$regions,
	)
	= rearrange( [
		'result',
		'peptide_signal',	'mitochondrial_signal',
		'sp_score',			'tp_score',
		'regions',
	],
	@_
	);

 	$self->result($result);
 	$self->peptide_signal($peptide_signal);
 	$self->mitochondrial_signal($mitochondrial_signal);
 	$self->sp_score($sp_score);
 	$self->tp_score($tp_score);
	$self->regions($regions) if(defined $regions);
		
	return $self;
}

=head2 result

  Arg [1]    : (optional) String - the result to set
  Example    : $analysis->result(123);
  Description: Getter/setter for the result for this analysis
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

=head2 peptide_signal

  Arg [1]    : (optional) String - the analysed peptide_signal to set
  Example    : $analysis->peptide_signal(123);
  Description: Getter/setter for the analysed peptide_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_signal {
	my ($self) = shift;
	$self->{'peptide_signal'} = shift if(@_);
	return $self->{'peptide_signal'};
}

=head2 mitochondrial_signal

  Arg [1]    : (optional) String - the analysed mitochondrial_signal to set
  Example    : $analysis->mitochondrial_signal(123);
  Description: Getter/setter for the analysed mitochondrial_signal
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub mitochondrial_signal {
	my ($self) = shift;
	$self->{'mitochondrial_signal'} = shift if(@_);
	return $self->{'mitochondrial_signal'};
}

=head2 sp_score

  Arg [1]    : (optional) String - the analysed sp_score to set
  Example    : $analysis->sp_score(123);
  Description: Getter/setter for the analysed sp_score
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sp_score {
	my ($self) = shift;
	$self->{'sp_score'} = shift if(@_);
	return $self->{'sp_score'};
}

=head2 tp_score

  Arg [1]    : (optional) String - the analysed tp_score to set
  Example    : $analysis->tp_score(123);
  Description: Getter/setter for the analysed tp_score
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub tp_score {
	my ($self) = shift;
	$self->{'tp_score'} = shift if(@_);
	return $self->{'tp_score'};
}

=head2 regions

  Arg [1]    : (optional) Listref of APPRIS::Analysis::CRASHRegion - 
               the set of peptide-mitochondrial regions that for this analysis 
  Example    : $analysis->regions($regions);
  Description: Getter/setter for the set of analysed regions 
  Returntype : Listref of APPRIS::Analysis::CRASHRegion or undef
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
