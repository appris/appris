=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::PROTEO - Object representing an analysis run

=head1 SYNOPSIS

  my $analysis = APPRIS::Analysis::PROTEO->new(
    -result      => <Analysis result>,
    -peptide_evidence  => <Annotation analysed>,
  );

=head1 DESCRIPTION

A representation of analysis of PROTEO within the APPRIS system.
Object to store details of an analysis run.

=head1 METHODS

=cut

package APPRIS::Analysis::PROTEO;

use strict;
use warnings;

use APPRIS::Analysis::PROTEORegion;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-result]  : 
       string - the anlysis result of the transcript
  Arg [-peptide_evidence]:
        string - the appris annotation for this analysis
  Arg [-num_peptides] : (optional)
        int - the number of analysed peptides
  Arg [-num_experiments] : (optional)
        int - the number of analysed peptides
  Arg [-peptides]: (optional)
        Listref of APPRIS::Analysis::PROTEORegion - the 
        set of peptides that was analysed
  Example    : $analysis = APPRIS::Analysis::PROTEO->new(...);
  Description: Creates a new analysis object
  Returntype : APPRIS::Analysis::PROTEO
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
		$result,			$peptide_evidence,
		$num_peptides,		$num_experiments,
		$peptides
	)
	= rearrange( [
		'result',			'peptide_evidence',
		'num_peptides',		'num_experiments',
		'peptides'
	],
	@_
	);

 	$self->result($result);
 	$self->peptide_evidence($peptide_evidence);
	$self->num_peptides($num_peptides) if(defined $num_peptides);
	$self->num_experiments($num_experiments) if(defined $num_experiments);
	$self->peptides($peptides) if(defined $peptides);
		
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

=head2 peptide_evidence

  Arg [1]    : (optional) String - the analysed peptide_evidence to set
  Example    : $analysis->peptide_evidence(123);
  Description: Getter/setter for the analysed peptide_evidence
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_evidence {
	my ($self) = shift;
	$self->{'peptide_evidence'} = shift if(@_);
	return $self->{'peptide_evidence'};
}

=head2 num_peptides

  Arg [1]    : (optional) Int - the number of peptides for 
               the analysis
  Example    : $analysis->num_peptides(5);
  Description: Getter/setter for the number of peptides 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_peptides {
	my ($self) = shift;
	$self->{'num_peptides'} = shift if(@_);
	return $self->{'num_peptides'};
}

=head2 num_experiments

  Arg [1]    : (optional) Int - the number of peptides for 
               the analysis
  Example    : $analysis->num_experiments(5);
  Description: Getter/setter for the number of peptides 
               that for this analysis
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_experiments {
	my ($self) = shift;
	$self->{'num_experiments'} = shift if(@_);
	return $self->{'num_experiments'};
}

=head2 peptides

  Arg [1]    : (optional) Listref of APPRIS::Analysis::PROTEORegion - 
               the set of peptides that for this analysis 
  Example    : $analysis->peptides($peptides);
  Description: Getter/setter for the set of analysed peptides 
  Returntype : Listref of APPRIS::Analysis::PROTEORegion or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptides {
	my ($self) = shift;
	$self->{'peptides'} = shift if(@_);
	return $self->{'peptides'};
}

sub DESTROY {}

1;
