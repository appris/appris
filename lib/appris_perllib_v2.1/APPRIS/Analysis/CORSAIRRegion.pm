=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::CORSAIRRegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::CORSAIRRegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a CORSAIRRegion within the APPRIS system.
A CORSAIRRegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::CORSAIRRegion;

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
  Arg [-cds_id]    : (optional)
       int - the order of CDS of the peptide region
  Arg [-maxscore]    : (optional)
       string - the pdb id founded for the peptide region
  Arg [-sp_report]    : (optional)
       string - sp_report of pdb in the region
  Arg [-type]    : (optional)
       string - type of region: exon or mini-exon
  Example    : $region = APPRIS::Analysis::CORSAIRRegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::CORSAIRRegion
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
		$score,		$cds_id,
		$maxscore,	$sp_report,
		$type
	)
	= rearrange( [
		'start',	'end',		'strand',
		'pstart',	'pend',
		'score',	'cds_id',
		'maxscore',	'sp_report',
		'type'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->pstart($pstart);
 	$self->pend($pend);
 	$self->score($score);
 	$self->cds_id($cds_id) if (defined $cds_id);
 	$self->maxscore($maxscore) if (defined $maxscore);
 	$self->sp_report($sp_report) if (defined $sp_report);
	$self->type($type) if (defined $type);
  
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

=head2 cds_id

  Arg [1]    : (optional) Int - the order of CDS of 
               peptide region
  Example    : $region->cds_id(1);
  Description: Getter/setter for the order of CDS
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub cds_id {
	my ($self) = shift;
	$self->{'cds_id'} = shift if(@_);
	return $self->{'cds_id'};
}

=head2 maxscore

  Arg [1]    : (optional) String - the maximum score founded to set
  Example    : $region->maxscore('3LUI_A');
  Description: Getter/setter for the maximum score founded
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub maxscore {
	my ($self) = shift;
	$self->{'maxscore'} = shift if(@_);
	return $self->{'maxscore'};
}

=head2 sp_report

  Arg [1]    : (optional) String - report of species
  Example    : $region->sp_report('33.3');
  Description: Getter/setter for the report of species
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sp_report {
	my ($self) = shift;
	$self->{'sp_report'} = shift if(@_);
	return $self->{'sp_report'};
}

=head2 type

  Arg [1]    : (optional) String - type of region: exon or mini-exon
  Example    : $region->type('exon');
  Description: Getter/setter for the type of region
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
	my ($self) = shift;
	$self->{'type'} = shift if(@_);
	return $self->{'type'};
}

sub DESTROY {}

1;
