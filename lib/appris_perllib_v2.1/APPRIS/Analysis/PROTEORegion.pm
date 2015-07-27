=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::PROTEORegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::PROTEORegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a PROTEORegion within the APPRIS system.
A PROTEORegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::PROTEORegion;

use strict;
use warnings;
use vars qw(@ISA);

use APPRIS::Analysis::Region;
use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

@ISA = qw(APPRIS::Analysis::Region);

=head2 new

  Arg [-peptide_id]  : (optional)
       string - peptide id
  Arg [-sequence] :
       string - the sequence of the region
  Arg [-num_experiments] :
       integer - num. experiments for this region
  Arg [-experiments] :
       string - the list of experiments for this region
  Arg [-pstart]  : 
       int - start postion of the peptide region
  Arg [-pend]    : 
       int - end position of the peptide region
  Arg [-start]  : (optional)
       int - start postion of the slice 
  Arg [-end]    : (optional)
       int - end position of the slice
  Arg [-strand] : (optional)
       char - '+','-' the strand the slice is on
  Example    : $region = APPRIS::Analysis::PROTEORegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::PROTEORegion
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
		$peptide_id,	$sequence,
		$num_experiments,	$experiments,
		$pstart,	$pend,
		$start,		$end,		$strand,
	)
	= rearrange( [
		'peptide_id',		'sequence',
		'num_experiments',	'experiments',
		'pstart',	'pend',
		'start',	'end',		'strand',
	],
	@_
	);

 	$self->peptide_id($peptide_id) if(defined $peptide_id);
 	$self->sequence($sequence);
 	$self->num_experiments($num_experiments);
 	$self->experiments($experiments);
 	$self->pstart($pstart);
 	$self->pend($pend);
 	$self->start($start) if(defined $start);
 	$self->end($end) if(defined $end);
 	$self->strand($strand) if(defined $strand);
  
	return $self;
}

=head2 peptide_id

  Arg [1]    : (optional) String - peptide id
  Example    : $region->id(123);
  Description: Getter/setter for peptide id
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub peptide_id {
	my ($self) = shift;
	$self->{'peptide_id'} = shift if(@_);
	return $self->{'peptide_id'};
}

=head2 sequence

  Arg [1]    : (optional) String - the sequence of the region is on
  Example    : $region->sequence('-');
  Description: Getter/setter for the sequence that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sequence {
	my ($self) = shift;
	$self->{'sequence'} = shift if(@_);
	return $self->{'sequence'};
}

=head2 num_experiments

  Arg [1]    : (optional) Integer - num. experiments for this region
  Example    : $region->num_experiments('-');
  Description: Getter/setter for the num. experiments that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub num_experiments {
	my ($self) = shift;
	$self->{'num_experiments'} = shift if(@_);
	return $self->{'num_experiments'};
}

=head2 experiments

  Arg [1]    : (optional) String - the list of experiments for this region
  Example    : $region->experiments('-');
  Description: Getter/setter for the experiments that for this region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub experiments {
	my ($self) = shift;
	$self->{'experiments'} = shift if(@_);
	return $self->{'experiments'};
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

sub DESTROY {}

1;
