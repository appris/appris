#
# Copyright (C) 1998,1999 Ken MacLeod
# Data::Grove::Visitor is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.
#
# $Id: Visitor.pm,v 1.6 2000/03/20 23:06:45 kmacleod Exp $
#

use strict;
use 5.005;

package Data::Grove::Visitor;

use vars qw{ $VERSION };

# will be substituted by make-rel script
$VERSION = "0.08";

# The following methods extend Data::Grove
package Data::Grove;

sub accept {
    my $self = shift;
    my $visitor = shift;

    my $type_name;
    my $package = ref($self);
    eval "\$type_name = \$${package}::type_name";
    if (!defined $type_name) {
	return (); # no action
    }

    my $method_name = 'visit_' . $type_name;
    if ($visitor->can($method_name)) {
	return $visitor->$method_name ($self, @_);
    } else {
	return (); # no action
    }
}

sub accept_name {
    my $self = shift;

    if (!defined $self->{Name}) {
	return $self->accept (@_);
    }

    my $visitor = shift;

    my $name = $self->{Name};
    $name =~ s/\W/_/g;
    my $name_method = "visit_name_$name";

    if (!$self->{'has'}{$name_method}) {
	return if (defined $self->{'has'}{$name_method});
	$self->{'has'}{$name_method} = $visitor->can($name_method);
	return $self->accept($visitor, @_) if (!$self->{'has'}{$name_method});
    }

    return $visitor->$name_method ($self, @_);
}

sub attr_accept {
    my $self = shift; my $attr = shift; my $visitor = shift;

    if (!defined $self->{Attributes}) {
	return (); # no action
    }

    my $attrs = $self->{Attributes}{$attr};
    if (ref($attrs) eq 'ARRAY') {
	return $self->_children_accept ($attrs, $visitor, @_);
    } else {

	if (!$self->{has_visit_characters}) {
	    return if (defined $self->{has_visit_characters});
	    $self->{has_visit_characters} = $visitor->can('visit_characters');
	    return if (!$self->{has_visit_characters});
	}
	# FIXME should be some other generic than XML::Grove::Characters
	return $visitor->visit_characters (XML::Grove::Characters->new(Data => $attrs), @_);
    }
}

sub children_accept {
    my $self = shift;

    if (defined $self->{Contents}) {
	return $self->_children_accept ($self->{Contents}, @_);
    } else {
	return (); # no action
    }
}

sub children_accept_name {
    my $self = shift;

    if (defined $self->{Contents}) {
	return $self->_children_accept_name ($self->{Contents}, @_);
    } else {
	return (); # no action
    }
}

sub _children_accept {
    my $self = shift; my $array = shift; my $visitor = shift;

    my @return;
    my $ii;
    for ($ii = 0; $ii <= $#$array; $ii ++) {
	push @return, $array->[$ii]->accept ($visitor, @_);
    }

    return @return;
}

sub _children_accept_name {
    my $self = shift; my $array = shift; my $visitor = shift;

    my @return;
    my $ii;
    for ($ii = 0; $ii <= $#$array; $ii ++) {
	push @return, $array->[$ii]->accept_name ($visitor, @_);
    }

    return @return;
}

1;

__END__

=head1 NAME

Data::Grove::Visitor - add visitor/callback methods to Data::Grove objects

=head1 SYNOPSIS

 use Data::Grove::Visitor;

 @results = $object->accept ($visitor, ...);
 @results = $object->accept_name ($visitor, ...);
 @results = $object->children_accept ($visitor, ...);
 @results = $object->children_accept_name ($visitor, ...);

=head1 DESCRIPTION

Data::Grove::Visitor adds visitor methods (callbacks) to Data::Grove
objects.  A ``visitor'' is a class (a package) you write that has
methods (subs) corresponding to the objects in the classes being
visited.  You use the visitor methods by creating an instance of your
visitor class, and then calling `C<accept($my_visitor)>' on the
top-most object you want to visit, that object will in turn call your
visitor back with `C<visit_I<OBJECT>>', where I<OBJECT> is the type of
object.

There are several forms of `C<accept>'.  Simply calling `C<accept>'
calls your package back using the object type of the object you are
visiting.  Calling `C<accept_name>' on an element object calls you
back with `C<visit_name_I<NAME>>' where I<NAME> is the tag name of the
element, on all other objects it's as if you called `C<accept>'.

All of the forms of `C<accept>' return a concatenated list of the
result of all `C<visit>' methods.

`C<children_accept>' calls `C<accept>' on each of the children of the
element.  This is generally used in element callbacks to recurse down
into the element's children, you don't need to get the element's
contents and call `C<accept>' on each item.  `C<children_accept_name>'
does the same but calling `C<accept_name>' on each of the children.
`C<attr_accept>' calls `C<accept>' on each of the objects in the named
attribute.

Refer to the documentation of the classes you are visiting
(XML::Grove, etc.) for the type names (`C<element>', `C<document>',
etc.) of the objects it implements.

=head1 RESERVED NAMES

The hash keys `C<Contents>' and `C<Name>' are used to indicate objects
with children (for `C<children_accept>') and named objects (for
`C<accept_name>').

=head1 NOTES

These are random ideas that haven't been implemented yet:

=over 4

=item *

Several objects fall into subclasses, or you may want to be able to
subclass a visited object and still be able to tell the difference.
In SGML::Grove I had used the package name in the callback
(`C<visit_SGML_Element>') instead of a generic name
(`C<visit_element>').  The idea here would be to try calling
`C<visit_I<PACKAGE>>' with the most specific class first, then try
superclasses, and lastly to try the generic.

=back

=head1 AUTHOR

Ken MacLeod, ken@bitsko.slc.ut.us

=head1 SEE ALSO

perl(1), Data::Grove

Extensible Markup Language (XML) <http://www.w3c.org/XML>

=cut
