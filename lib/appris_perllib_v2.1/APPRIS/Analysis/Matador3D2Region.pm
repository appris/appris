=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::Matador3D2Region - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::Matador3D2Region->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a Matador3D2Region within the APPRIS system.
A Matador3D2Region object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::Matador3D2Region;

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
  Arg [-bias] : 
       string - the bias of the region
  Arg [-pdb_id]    : (optional)
       string - the pdb id founded for the peptide region
  Example    : $region = APPRIS::Analysis::Matador3D2Region->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::Matador3D2Region
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
		$pstart,	$pend,
		$score,		$bias,
		$pdb_id
	)
	= rearrange( [
		'start',	'end',		'strand',
		'pstart',	'pend',
		'score',	'bias',
		'pdb_id'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->pstart($pstart);
 	$self->pend($pend);
 	$self->score($score);
 	$self->bias($bias);
 	$self->pdb_id($pdb_id) if (defined $pdb_id);
  
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

=head2 score

  Arg [1]    : (optional) String - the score of the region is on
  Example    : $region->score('-');
  Description: Getter/setter for the score that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score {
	my ($self) = shift;
	$self->{'score'} = shift if(@_);
	return $self->{'score'};
}

=head2 bias

  Arg [1]    : (optional) String - the bias of the region is on
  Example    : $region->bias('-');
  Description: Getter/setter for the bias that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub bias {
	my ($self) = shift;
	$self->{'bias'} = shift if(@_);
	return $self->{'bias'};
}

=head2 pdb_id

  Arg [1]    : (optional) String - the pdb id founded to set
  Example    : $region->pdb_id('3LUI_A');
  Description: Getter/setter for the pdb id founded
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pdb_id {
	my ($self) = shift;
	$self->{'pdb_id'} = shift if(@_);
	return $self->{'pdb_id'};
}

sub DESTROY {}

1;
