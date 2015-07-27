=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Analysis::SPADERegion - Object representing a region

=head1 SYNOPSIS

  my $region = APPRIS::Analysis::SPADERegion->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print region information
  print("region start:end:strand is "
      . join( ":", map { $region->$_ } qw(start end strand) )
      . "\n" );

=head1 DESCRIPTION

A representation of a SPADERegion within the APPRIS system.
A SPADERegion object represents a region of a genome from analysis
method.

=head1 METHODS

=cut

package APPRIS::Analysis::SPADERegion;

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
  Arg [-cds_id]    : 
       int - the order of CDS of the peptide region
  Arg [-pdb_id]    : 
       string - the pdb id founded for the peptide region
  Example    : $region = APPRIS::Analysis::SPADERegion->new(...);
  Description: Creates a new region object
  Returntype : APPRIS::Analysis::SPADERegion
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
		$score,				$type_domain,
		$alignment_start,	$alignment_end,
		$envelope_start,	$envelope_end,
		$hmm_start,			$hmm_end,
		$hmm_length,		$hmm_acc,
		$hmm_name,			$hmm_type,
		$bit_score,			$evalue,
		$significance,		$clan,
		$predicted_active_site_residues,
		$external_id
	)
	= rearrange( [
		'start',	'end',		'strand',
		'score',			'type_domain',
		'alignment_start',	'alignment_end',
		'envelope_start',	'envelope_end',
		'hmm_start',		'hmm_end',
		'hmm_length',		'hmm_acc',
		'hmm_name',			'hmm_type',
		'bit_score',		'evalue',
		'significance',		'clan',
		'predicted_active_site_residues',
		'external_id'
	],
	@_
	);

 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->score($score);
 	$self->type_domain($type_domain);
 	$self->alignment_start($alignment_start) if (defined $alignment_start);
 	$self->alignment_end($alignment_end) if (defined $alignment_end);
 	$self->envelope_start($envelope_start) if (defined $envelope_start);
 	$self->envelope_end($envelope_end) if (defined $envelope_end);
 	$self->hmm_start($hmm_start) if (defined $hmm_start);
 	$self->hmm_end($hmm_end) if (defined $hmm_end);
 	$self->hmm_length($hmm_length) if (defined $hmm_length);
 	$self->hmm_acc($hmm_acc) if (defined $hmm_acc);
 	$self->hmm_name($hmm_name) if (defined $hmm_name);
 	$self->hmm_type($hmm_type) if (defined $hmm_type);
 	$self->bit_score($bit_score) if (defined $bit_score);
 	$self->evalue($evalue) if (defined $evalue);
 	$self->significance($significance) if (defined $significance);
 	$self->clan($clan) if (defined $clan);
	if (defined $predicted_active_site_residues and $predicted_active_site_residues ne 'NULL') {
		$self->predicted_active_site_residues($predicted_active_site_residues);
	}
	if (defined $external_id and $external_id ne 'NULL') {
		$self->external_id($external_id);
	}
  
	return $self;
}

=head2 score

  Arg [1]    : (optional) Int - the score of the region is on
  Example    : $region->score(1);
  Description: Getter/setter for the score that for this region
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score {
	my ($self) = shift;
	$self->{'score'} = shift if(@_);
	return $self->{'score'};
}

=head2 type_domain

  Arg [1]    : (optional) String - the type of domain for this region
  Example    : $region->type_domain('domain');
  Description: Getter/setter for the type of domain
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type_domain {
	my ($self) = shift;
	$self->{'type_domain'} = shift if(@_);
	return $self->{'type_domain'};
}

=head2 alignment_start

  Arg [1]    : (optional) Int - the start position of alignment
  Example    : $region->alignment_start(1);
  Description: Getter/setter for the start position of alignment
  Returntype : Int
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

  Arg [1]    : (optional) Int - the end position of alignment
  Example    : $region->alignment_end(1);
  Description: Getter/setter for the end position of alignment
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub alignment_end {
	my ($self) = shift;
	$self->{'alignment_end'} = shift if(@_);
	return $self->{'alignment_end'};
}

=head2 envelope_start

  Arg [1]    : (optional) Int - the start position of envelopment
  Example    : $region->envelope_start(1);
  Description: Getter/setter for the start position of envelopment
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub envelope_start {
	my ($self) = shift;
	$self->{'envelope_start'} = shift if(@_);
	return $self->{'envelope_start'};
}

=head2 envelope_end

  Arg [1]    : (optional) Int - the end position of envelopment
  Example    : $region->envelope_end(1);
  Description: Getter/setter for the end position of envelopment
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub envelope_end {
	my ($self) = shift;
	$self->{'envelope_end'} = shift if(@_);
	return $self->{'envelope_end'};
}

=head2 hmm_start

  Arg [1]    : (optional) Int - the start position of hmm
  Example    : $region->hmm_start(1);
  Description: Getter/setter for the start position of hmm
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hmm_start {
	my ($self) = shift;
	$self->{'hmm_start'} = shift if(@_);
	return $self->{'hmm_start'};
}

=head2 hmm_end

  Arg [1]    : (optional) Int - the end position of hmm
  Example    : $region->hmm_end(1);
  Description: Getter/setter for the end position of hmm
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hmm_end {
	my ($self) = shift;
	$self->{'hmm_end'} = shift if(@_);
	return $self->{'hmm_end'};
}

=head2 hmm_length

  Arg [1]    : (optional) Int - the length of hmm
  Example    : $region->hmm_length(1);
  Description: Getter/setter for the length of hmm
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hmm_length {
	my ($self) = shift;
	$self->{'hmm_length'} = shift if(@_);
	return $self->{'hmm_length'};
}

=head2 hmm_acc

  Arg [1]    : (optional) String - the accession of hmm
  Example    : $region->hmm_acc(1);
  Description: Getter/setter for the accesion of hmm
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hmm_acc {
	my ($self) = shift;
	$self->{'hmm_acc'} = shift if(@_);
	return $self->{'hmm_acc'};
}

=head2 hmm_name

  Arg [1]    : (optional) String - the name of hmm
  Example    : $region->hmm_name(1);
  Description: Getter/setter for the name of hmm
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hmm_name {
	my ($self) = shift;
	$self->{'hmm_name'} = shift if(@_);
	return $self->{'hmm_name'};
}

=head2 hmm_type

  Arg [1]    : (optional) String - the type of hmm
  Example    : $region->hmm_type(1);
  Description: Getter/setter for the type of hmm
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub hmm_type {
	my ($self) = shift;
	$self->{'hmm_type'} = shift if(@_);
	return $self->{'hmm_type'};
}

=head2 bit_score

  Arg [1]    : (optional) String - the score of alignment
  Example    : $region->bit_score(1);
  Description: Getter/setter for the score of alignment
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub bit_score {
	my ($self) = shift;
	$self->{'bit_score'} = shift if(@_);
	return $self->{'bit_score'};
}

=head2 evalue

  Arg [1]    : (optional) String - the evalue of alignment
  Example    : $region->evalue(1);
  Description: Getter/setter for the evalue of alignment
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub evalue {
	my ($self) = shift;
	$self->{'evalue'} = shift if(@_);
	return $self->{'evalue'};
}

=head2 significance

  Arg [1]    : (optional) String - the significance of alignment
  Example    : $region->significance(1);
  Description: Getter/setter for the significance of alignment
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub significance {
	my ($self) = shift;
	$self->{'significance'} = shift if(@_);
	return $self->{'significance'};
}

=head2 clan

  Arg [1]    : (optional) String - the clan of alignment
  Example    : $region->clan(1);
  Description: Getter/setter for the clan of alignment
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub clan {
	my ($self) = shift;
	$self->{'clan'} = shift if(@_);
	return $self->{'clan'};
}

=head2 predicted_active_site_residues

  Arg [1]    : (optional) String - the predicted active site residues of alignment
  Example    : $region->predicted_active_site_residues(1);
  Description: Getter/setter for the predicted active site residues of alignment
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub predicted_active_site_residues {
	my ($self) = shift;
	$self->{'predicted_active_site_residues'} = shift if(@_);
	return $self->{'predicted_active_site_residues'};
}

=head2 external_id

  Arg [1]    : (optional) String - transcript if from other Pfam result
  Example    : $region->external_id(1);
  Description: Getter/setter for the domain coming from Pfam result of other transcript
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

sub DESTROY {}

1;
