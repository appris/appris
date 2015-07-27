=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Gene - Object representing a gene

=head1 SYNOPSIS

  my $gene = APPRIS::Gene->new(
    -chr	=> 1,
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print gene information
  print("gene start:end:strand is "
      . join( ":", map { $gene->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $gene->stable_id('ENSG000001');
  $gene->description('This is the gene description');

=head1 DESCRIPTION

A representation of a Gene within the APPRIS system.
A gene is a set of one or more alternative transcripts.

=head1 METHODS

=cut

package APPRIS::Gene;

use strict;
use warnings;
use Data::Dumper;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-chr]  : 
       int - chromosome of the gene
  Arg [-start]  : 
       int - start postion of the gene
  Arg [-end]    : 
       int - end position of the gene
  Arg [-strand] : 
       char - '+','-' the strand the gene is on
  Arg [-stable_id] :
        string - the stable identifier of this gene
  Arg [-description]:
        string - the genes description
  Arg [-biotype]:
        string - the biotype e.g. "protein_coding"
  Arg [-status]:
        string - the gene status i.e. "KNOWN","NOVEL"
  Arg [-source]:
        string - the genes source, e.g. "ensembl"        
  Arg [-level] :
        int - the level of this gene
  Arg [-version] :
        int - the version of this gene
  Arg [-external_name] :
        string - the external database name associated with this gene
  Arg [-xref_identify]:
        Listref of APPRIS::XrefEntry - the external database entry that for this gene
  Arg [-transcripts]:
        Listref of APPRIS::Transcripts - this gene's transcripts
  Arg [-analysis]: (optional)
        APPRIS::Analysis - set of analysis for this gene        
  Example    : $gene = APPRIS::Gene->new(...);
  Description: Creates a new gene object
  Returntype : APPRIS::Gene
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
		$xref_identify,		$transcripts,
		$analysis
	)
	= rearrange( [
		'chr',				'start',
		'end',				'strand',
		'stable_id',		'description',		
		'biotype',			'status',
		'source',			'level',
		'version',			'external_name',
		'xref_identify',	'transcripts',
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
	$self->external_name($external_name) if(defined $external_name);
	$self->xref_identify($xref_identify) if(defined $xref_identify);
	$self->transcripts($transcripts) if(defined $transcripts);
	$self->analysis($analysis) if(defined $analysis);
	
	return $self;
}

=head2 chromosome

  Arg [1]    : (optional) String - the chromosome to set
  Example    : $gene->chromosome("21");
  Description: Getter/setter for chromosome for this gene
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
  Example    : $gene->start(123);
  Description: Getter/setter for start location for this gene
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
  Example    : $gene->end(123);
  Description: Getter/setter for end location for this gene
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

  Arg [1]    : (optional) Char - the strand the gene is on
  Example    : $gene->strand('-');
  Description: Getter/setter for strand for this gene
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
  Example    : $gene->stable_id("ENSG0000000001");
  Description: Getter/setter for stable id for this gene.
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
  Example    : $gene->description('This is the gene\'s description');
  Description: Getter/setter for gene description
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
  Example    : $gene->biotype("protein_coding");
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
  Example    : $gene->status('KNOWN');
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
  Example    : $gene->source('ensembl');
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
  Example    : $gene->level(2);
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
  Example    : $gene->version(2);
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

=head2 external_name

  Arg [1]    : (optional) String - the external name to set
  Example    : $gene->external_name('BRCA2');
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
               that for this gene.
  Example    : $gene->xref_identify();
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

=head2 transcripts

  Arg [1]    : (optional) Listref of APPRIS::Transcript
  Arg [2]    : (optional) Listref of Index of APPRIS::Transcript  
  Example    : $gene->transcripts();
  Description: Getter/setter for this gene's transcripts
  Returntype : Listref of APPRIS::Transcript or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcripts {
	my ($self) = shift;
	$self->{'transcripts'} = shift if(@_);	
	$self->{'_index_transcripts'} = shift if(@_); # index the list of transcripts	
	return $self->{'transcripts'};
}

=head2 transcript

  Arg [1]    : String - the stable identifier of this transcript.
  Example    : $gene->transcript();
  Description: Getter for this gene's transcript
  Returntype : APPRIS::Transcript or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
	my ($self) = shift;
	my ($id) = shift;
	my ($index) = $self->_index_transcripts->{$id};
	return $self->{'transcripts'}->[$index];
}

=head2 analysis

  Arg [1]    : (optional) APPRIS::Analysis - reference to 
               APPRIS::Analysis object - set of analysis 
               which executed for this gene
  Example    : $gene->analysis();
  Description: Getter/setter for the analysis which executed for 
               this gene
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
