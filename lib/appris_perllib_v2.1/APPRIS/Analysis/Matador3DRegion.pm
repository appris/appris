=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::Matador3DRegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::Matador3DRegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a Matador3DRegion within the APPRIS system.
A Matador3DRegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::Matador3DRegion;

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
  Arg [-alignment_start]    : (optional)
       string - start position of pdb alignment
  Arg [-alignment_end]    : (optional)
       string - end position of pdb alignment
  Arg [-pdb_id]    : (optional)
       string - the pdb id founded for the peptide region
  Arg [-identity]    : (optional)
       string - identity of pdb in the region
  Arg [-external_id]    : (optional)
       string - the trans id which obtain the alignment
  Arg [-type]    : (optional)
       string - type of region: exon or mini-exon
  Example    : $region = APPRIS::Analysis::Matador3DRegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::Matador3DRegion
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
		$alignment_start, $alignment_end,
		$pdb_id,	$identity,  $external_id,
		$type
	)
	= rearrange( [
		'start',	'end',		'strand',
		'pstart',	'pend',
		'score',	'cds_id',
		'alignment_start', 'alignment_end',
		'pdb_id',	'identity', 'external_id',
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
 	$self->alignment_start($alignment_start) if (defined $alignment_start);
 	$self->alignment_end($alignment_end) if (defined $alignment_end);
 	$self->pdb_id($pdb_id) if (defined $pdb_id);
 	$self->identity($identity) if (defined $identity);
 	$self->external_id($external_id) if (defined $external_id);
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

=head2 alignment_start

  Arg [1]    : (optional) String - start position of pdb alignment
  Example    : $region->alignment_start('76');
  Description: Getter/setter for the alignment postion
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub alignment_start {
	my ($self) = shift;
	$self->{'alignment_start'} = shift if(@_);
	return $self->{'alignment_start'};
}

=head2 alignment_end

  Arg [1]    : (optional) String - start position of pdb alignment
  Example    : $region->alignment_end('76');
  Description: Getter/setter for the alignment postion
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub alignment_end {
	my ($self) = shift;
	$self->{'alignment_end'} = shift if(@_);
	return $self->{'alignment_end'};
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

=head2 identity

  Arg [1]    : (optional) String - identity of pdb in the region
  Example    : $region->identity('33.3');
  Description: Getter/setter for the identity of pdb
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub identity {
	my ($self) = shift;
	$self->{'identity'} = shift if(@_);
	return $self->{'identity'};
}

=head2 external_id

  Arg [1]    : (optional) String - the trans id which obtain the alignment
  Example    : $region->external_id('ENST...');
  Description: Getter/setter for the external id founded
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub external_id {
	my ($self) = shift;
	$self->{'external_id'} = shift if(@_);
	return $self->{'external_id'};
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
