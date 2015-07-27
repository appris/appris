=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::CRASHRegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::CRASHRegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a CRASHRegion within the APPRIS system.
A CRASHRegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::CRASHRegion;

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
  Arg [-pstart]  : 
       int - start postion of the peptide region
  Arg [-pend]    : 
       int - end position of the peptide region
  Arg [-score] : 
       string - the score of the region
  Example    : $region = APPRIS::Analysis::CRASHRegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::CRASHRegion
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
		$start,			$end,		$strand,	
		$pstart,		$pend,
		$s_mean,		$s_prob,
		$d_score,		$c_max,
		$reliability,	$localization
	)
	= rearrange( [
		'start',		'end',		'strand',	
		'pstart',		'pend',
		's_mean',		's_prob',
		'd_score',		'c_max',
		'reliability',	'localization',
		
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->pstart($pstart) if (defined $pstart);
 	$self->pend($pend) if (defined $pend); 	
 	$self->s_mean($s_mean) if (defined $s_mean);
 	$self->s_prob($s_prob) if (defined $s_prob);
 	$self->d_score($d_score) if (defined $d_score);
 	$self->c_max($c_max) if (defined $c_max);
 	$self->reliability($reliability) if (defined $reliability);
 	$self->localization($localization) if (defined $localization);
  
	return $self;
}

=head2 pstart

  Arg [1]    : (optional) Int - the start of peptide location to set
  Example    : $region->pstart(123);
  Description: Getter/setter for start of hit location
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pstart {
	my ($self) = shift;
	$self->{'pstart'} = shift if(@_);
	return $self->{'pstart'};
}

=head2 pend

  Arg [1]    : (optional) Int - the end of peptide location to set
  Example    : $region->pend(123);
  Description: Getter/setter for end of hit location
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pend {
	my ($self) = shift;
	$self->{'pend'} = shift if(@_);
	return $self->{'pend'};
}

=head2 s_mean

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->s_mean('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub s_mean {
	my ($self) = shift;
	$self->{'s_mean'} = shift if(@_);
	return $self->{'s_mean'};
}

=head2 s_prob

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->s_prob('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub s_prob {
	my ($self) = shift;
	$self->{'s_prob'} = shift if(@_);
	return $self->{'s_prob'};
}

=head2 d_score

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->d_score('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub d_score {
	my ($self) = shift;
	$self->{'d_score'} = shift if(@_);
	return $self->{'d_score'};
}

=head2 c_max

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->c_max('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub c_max {
	my ($self) = shift;
	$self->{'c_max'} = shift if(@_);
	return $self->{'c_max'};
}

=head2 reliability

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->reliability('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub reliability {
	my ($self) = shift;
	$self->{'reliability'} = shift if(@_);
	return $self->{'reliability'};
}

=head2 localization

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->localization('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub localization {
	my ($self) = shift;
	$self->{'localization'} = shift if(@_);
	return $self->{'localization'};
}

sub DESTROY {}

1;
