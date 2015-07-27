=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::XrefEntry -
Object representing an external reference (xref)

=head1 SYNOPSIS

=head1 DESCRIPTION

This object holds information about external references (xrefs) to
Ensembl objects.

=head1 METHODS

=cut

package APPRIS::XrefEntry;

use strict;
use warnings;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

=head2 new

  Arg [-id]  : 
       int - chromosome of the gene
  Arg [-dbname]  : 
       string - start postion of the gene
  Arg [-description]:
        string - the genes description       
  Example    : $gene = APPRIS::XrefEntry->new(...);
  Description: Creates a new XrefEntry object
  Returntype : APPRIS::XrefEntry
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
		$id,	$dbname,
		$description
    )
    = rearrange( [
		'id',	'dbname',
		'description'
	],
	@_
	);

 	$self->id($id);
 	$self->dbname($dbname);
 	$self->description($description) if(defined $description);
  
	return $self;
}

=head2 id

  Arg [1]    : (optional) String - the ID to set
  Example    : $gene->id("Q8N6Q8");
  Description: Getter/setter for id for this entry
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub id {
  my $self = shift;
  $self->{'id'} = shift if(@_);
  return $self->{'id'};
}

=head2 dbname

  Arg [1]    : (optional) String - the database name
  Example    : $gene->dbname('UniProtKB_SwissProt');
  Description: Getter/setter for the name of database
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dbname {
    my $self = shift;
    $self->{'dbname'} = shift if( @_ );
    return $self->{'dbname'};
}

=head2 description

  Arg [1]    : (optional) String - the description to xref identifier
  Example    : $gene->description('This is the xref id description');
  Description: Getter/setter for the description to xref identifier
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if( @_ );
    return $self->{'description'};
}

sub DESTROY {}

1;
