=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Transcript - Object representing a transcript

=head1 SYNOPSIS

  my $transcript = APPRIS::Transcript->new(
    -chr	=> 1,
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print transcript information
  print("transcript start:end:strand is "
      . join( ":", map { $transcript->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $transcript->stable_id('ENST000001');
  $transcript->description('This is the transcript description');

=head1 DESCRIPTION

A representation of a Transcript within the APPRIS system.
A transcript consists of a set of Exons and (possibly) a 
Translation which defines the coding and non-coding regions of the exons.
Also, a transcript is a set of analysis results of some methos 
for this transcript.

=head1 METHODS

=cut

package APPRIS::Transcript;

use strict;
use warnings;
use Data::Dumper;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-chr]  : 
       int - chromosome of the transcript
  Arg [-start]  : 
       int - start postion of the transcript
  Arg [-end]    : 
       int - end position of the transcript
  Arg [-strand] : 
       char - '+','-' the strand the transcript is on
  Arg [-stable_id] :
        string - the stable identifier of this transcript
  Arg [-description]: (optional)
        string - the transcripts description
  Arg [-biotype]:
        string - the biotype e.g. "protein_coding"
  Arg [-status]:
        string - the transcript status i.e. "KNOWN","NOVEL"
  Arg [-source]:
        string - the transcripts source, e.g. "ensembl"        
  Arg [-level] :
        int - the level of this transcript
  Arg [-version] :
        int - the version of this transcript
  Arg [-sequence] : (optional)
        string - the nucleotide sequence associated with this transcript
  Arg [-external_name] : (optional)
        string - the external database name associated with this transcript
  Arg [-xref_identify]: (optional)
        Listref of APPRIS::XrefEntry - the external database entry that for this transcript
  Arg [-translate]: (optional)
        APPRIS::Peptide - the translate object that for this transcript
  Arg [-exons]: (optional)
        Listref of APPRIS::Exon - exons which make up this transcript        
  Arg [-analysis]: (optional)
        APPRIS::Analysis - set of analysis for this transcript
  Example    : $transcript = APPRIS::Transcript->new(...);
  Description: Creates a new transcript object
  Returntype : APPRIS::Transcript
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
		$chr,				$start,
		$end,				$strand,
		$stable_id,			$description,
		$biotype,			$status,
		$source,			$level,
		$version,			$external_name,
		$xref_identify,		$sequence,
		$translate,			$exons,
		$analysis
	)
	= rearrange( [
		'chr',				'start',
		'end',				'strand',
		'stable_id',		'description',		
		'biotype',			'status',
		'source',			'level',
		'version',			'external_name',
		'xref_identify',	'sequence',
		'translate',		'exons',
		'analysis'
	],
	@_
	);

 	$self->chromosome($chr);
 	$self->start($start);
 	$self->end($end);
 	$self->strand($strand);
 	$self->stable_id($stable_id);	
	$self->description($description) if(defined $description);
	$self->biotype($biotype);
	$self->status($status);
	$self->source($source);
	$self->level($level);
	$self->version($version);
	$self->sequence($sequence) if(defined $sequence);
	$self->external_name($external_name) if(defined $external_name);
	$self->xref_identify($xref_identify) if(defined $xref_identify);
	$self->translate($translate) if(defined $translate);
	$self->exons($exons) if(defined $exons);
	$self->analysis($analysis) if(defined $analysis);
	
	return $self;
}

=head2 chromosome

  Arg [1]    : (optional) String - the chromosome to set
  Example    : $transcript->chromosome("21");
  Description: Getter/setter for chromosome for this transcript
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub chromosome {
	my ($self) = shift;
	$self->{'chr'} = shift if(@_);
	return $self->{'chr'};
}

=head2 start

  Arg [1]    : (optional) Int - the start location to set
  Example    : $transcript->start(123);
  Description: Getter/setter for start location for this transcript
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
	my ($self) = shift;
	$self->{'start'} = shift if(@_);
	return $self->{'start'};
}

=head2 end

  Arg [1]    : (optional) Int - the end location to set
  Example    : $transcript->end(123);
  Description: Getter/setter for end location for this transcript
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
	my ($self) = shift;
	$self->{'end'} = shift if(@_);
	return $self->{'end'};
}

=head2 strand

  Arg [1]    : (optional) Char - the strand the transcript is on
  Example    : $transcript->strand('-');
  Description: Getter/setter for strand for this transcript
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strand {
	my ($self) = shift;
	$self->{'strand'} = shift if(@_);
	return $self->{'strand'};
}

=head2 stable_id

  Arg [1]    : (optional) String - the stable ID to set
  Example    : $transcript->stable_id("ENST0000000001");
  Description: Getter/setter for stable id for this transcript
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
	my ($self) = shift;
	$self->{'stable_id'} = shift if(@_);
	return $self->{'stable_id'};
}

=head2 description

  Arg [1]    : (optional) String - the description to set
  Example    : $transcript->description('This is the transcript\'s description');
  Description: Getter/setter for transcript description
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
	my ($self) = shift;
	$self->{'description'} = shift if( @_ );
	return $self->{'description'};
}

=head2 biotype

  Arg [1]    : (optional) String - the biotype to set
  Example    : $transcript->biotype("protein_coding");
  Description: Getter/setter for the attribute biotype
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub biotype {
	my ($self) = shift;
	$self->{'biotype'} = shift if( @_ );
	return $self->{'biotype'};
}

=head2 status

  Arg [1]    : (optional) String - status to set
  Example    : $transcript->status('KNOWN');
  Description: Getter/setter for attribute status
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Medium Risk

=cut

sub status {
	my ($self) = shift;
	$self->{'status'} = shift if( @_ );
	return $self->{'status'};
}

=head2 source

  Arg [1]    : (optional) String - the source to set
  Example    : $transcript->source('ensembl');
  Description: Getter/setter for attribute source
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source {
	my ($self) = shift;
	$self->{'source'} = shift if( @_ );
	return $self->{'source'};
}

=head2 level

  Arg [1]    : (optional) Int - the level to set
  Example    : $transcript->level(2);
  Description: Getter/setter for attribute level
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub level {
	my ($self) = shift;
	$self->{'level'} = shift if( @_ );
	return $self->{'level'};
}

=head2 version

  Arg [1]    : (optional) Int - the version to set
  Example    : $transcript->version(2);
  Description: Getter/setter for attribute level
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
	my ($self) = shift;
	$self->{'version'} = shift if( @_ );
	return $self->{'version'};
}

=head2 sequence

  Arg [1]    : (optional) String - the nucleotide sequence 
               that for this transcript
  Example    : $transcript->sequence();
  Description: Getter/setter for nucleotide sequence for this 
               transcript
  Returntype : String or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sequence {
	my ($self) = shift;
	$self->{'sequence'} = shift if(@_);
	return $self->{'sequence'};
}

=head2 external_name

  Arg [1]    : (optional) String - the external name to set
  Example    : $transcript->external_name('BRCA2');
  Description: Getter/setter for attribute external_name
  Returntype : String or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub external_name {
	my ($self) = shift;
	$self->{'external_name'} = shift if(@_);
	return $self->{'external_name'};
}

=head2 xref_identify

  Arg [1]    : (optional) Listref - the external database entry 
               that for this transcript
  Example    : $transcript->xref_identify();
  Description: Getter/setter for the external database entry
  Returntype : Listref of APPRIS::XrefEntry or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub xref_identify {
	my ($self) = shift;
	$self->{'xref_identify'} = shift if(@_);
	return $self->{'xref_identify'};
}

=head2 translate

  Arg [1]    : (optional) APPRIS::Translation - the translate object 
               that for this transcript
  Example    : $transcript->translate();
  Description: Getter/setter for the external database entry
  Returntype : APPRIS::Translation or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translate {
	my ($self) = shift;
	$self->{'translate'} = shift if(@_);
	return $self->{'translate'};
}

=head2 exons

  Arg [1]    : (optional) Listref of APPRIS::Exon - reference to list 
               of APPRIS::Exon objects - exons which make up 
               this transcript
  Example    : $transcript->exons();
  Description: Getter/setter for the exons which make up 
               this transcript
  Returntype : Listref of APPRIS::Exon or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub exons {
	my ($self) = shift;
	$self->{'exons'} = shift if(@_);
	return $self->{'exons'};
}

=head2 analysis

  Arg [1]    : (optional) APPRIS::Analysis - reference to 
               APPRIS::Analysis object - set of analysis 
               which executed for this transcript
  Example    : $transcript->analysis();
  Description: Getter/setter for the analysis which executed for 
               this transcript
  Returntype : APPRIS::Analysis or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub analysis {
	my ($self) = shift;
	$self->{'analysis'} = shift if(@_);
	return $self->{'analysis'};
}

sub DESTROY {}

1;
