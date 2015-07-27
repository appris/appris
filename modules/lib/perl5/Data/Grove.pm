#
# Copyright (C) 1999 Ken MacLeod
# Data::Grove is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.
#
# $Id: Grove.pm,v 1.6 1999/12/22 21:15:00 kmacleod Exp $
#

###
### For a similar package, see also:
###
###   Graph::Element -- elements for a directed graph
###     Neil Bowers <neilb@cre.canon.co.uk> (NIELB)
### 

package Data::Grove;

use vars qw{ $VERSION };

# will be substituted by make-rel script
$VERSION = "0.08";

sub new {
    my $type = shift;
    my $self = ($#_ == 0) ? { %{ (shift) } } : { @_ };

    if (defined $self->{Raw}) {
	# clone the raw object
	$self = { %{ $self->{Raw} } };
    }

    return bless $self, $type;
}

package Data::Grove::Characters;
use vars qw{ @ISA $type_name };
@ISA = qw{Data::Grove};
$type_name = 'characters';

1;

__END__

=head1 NAME

Data::Grove -- support for deeply nested structures

=head1 SYNOPSIS

 use Data::Grove;

 $object = MyPackage->new;

 package MyPackage;
 @ISA = qw{Data::Grove};

=head1 DESCRIPTION

C<Data::Grove> provides support for deeply nested tree or graph
structures.  C<Data::Grove> is intended primarily for Perl module
authors writing modules with many types or classes of objects that
need to be manipulated and extended in a consistent and flexible way.

C<Data::Grove> is best used by creating a core set of ``data'' classes
and then incrementally adding functionality to the core data classes
by using ``extension'' modules.  One reason for this design is so that
the data classes can be swapped out and the extension modules can work
with new data sources.  For example, these other data sources could be
disk-based, network-based or built on top of a relational database.

Two extension modules that come with C<Data::Grove> are
C<Data::Grove::Parent> and C<Data::Grove::Visitor>.
C<Data::Grove::Parent> adds a `C<Parent>' property to grove objects
and implements a `C<root>' method to grove objects to return the root
node of the tree from anywhere in the tree and a `C<rootpath>' method
to return a list of nodes between the root node and ``this'' node.
C<Data::Grove::Visitor> adds callback methods `C<accept>' and
`C<accept_name>' that call your handler or receiver module back by
object type name or the object's name.

C<Data::Grove> objects do not contain parent references, Perl garbage
collection will delete them when no longer referenced and
sub-structures can be shared among several structures.
C<Data::Grove::Parent> is used to create temporary objects with parent
pointers.  

Properties of data classes are accessed directly using Perl's hash
functions (i.e. `C<$object-E<gt>{Property}>').  Extension modules may
also define properties that they support or use, for example
Data::Grove::Parent adds `C<Parent>' and `C<Raw>' properties and
Visitor depends on `C<Name>' and `C<Content>' properties.

See the module C<XML::Grove> for an example implementation of
C<Data::Grove>.

=head1 METHODS

=over 4

=item new( PROPERTIES )

Return a new object blessed into the SubClass, with the given
properties.  PROPERTIES may either be a list of key/value pairs, a
single hash containing key/value pairs, or an existing C<Data::Grove>
object.  If an existing C<Data::Grove> is passed to `C<new()>', a
shallow copy of that object will be returned.  A shallow copy means
that you are returned a new object, but all of the objects underneath
still refer to the original objects.

=back

=head1 AUTHOR

Ken MacLeod, ken@bitsko.slc.ut.us

=head1 SEE ALSO

perl(1)

=cut
