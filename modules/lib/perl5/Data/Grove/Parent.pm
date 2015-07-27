#
# Copyright (C) 1998,1999 Ken MacLeod
# Data::Grove::Parent is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.
#
# $Id: Parent.pm,v 1.2 1999/12/22 21:15:00 kmacleod Exp $
#

###
### WARNING
###
###
### This code has a bug in it that renders it useless.  In the FETCH
### routines, the new object created should have a reference to the
### the tied object that has $self as the underlying value.  As of
### this version, I don't know of a way to get to the tied object.
###

# Search for places marked `VALIDATE' to see where validation hooks
# may be added in the future.

use strict;

#--------------------------------------------------------------------------
# Data::Grove::Parent
#--------------------------------------------------------------------------

package Data::Grove::Parent;

use UNIVERSAL;
use Carp;

use vars qw{ $VERSION };

# will be substituted by make-rel script
$VERSION = "0.08";

sub new {
    my $type = shift;
    my $raw = shift;
    my $parent = shift;

    if (UNIVERSAL::isa($raw, 'Data::Grove::Parent')) {
        return $raw;
    }

    my @properties = ( Raw => $raw );

    if (defined $parent) {
        push @properties, Parent => $parent;
    }

    my $dummy = bless {}, ref($raw);
    tie %$dummy, $type, @properties;
    return $dummy;
}

sub TIEHASH  {
    my $type = shift;

    return bless { @_ }, $type;
}

sub STORE    {
    my $self = shift;
    my $key = shift;
    my $value = shift;

    if (exists $self->{$key}) {
	$self->{$key} = $value;
    } else {
	# VALIDATE
	if (UNIVERSAL::isa($value, 'Data::Grove::Parent')) {
	    $value = $value->{Raw};
	} elsif (UNIVERSAL::isa($value, 'Data::Grove::ParentList')) {
	    $value = $value->[0];
	}
	$self->{Raw}{$key} = $value;
    }
}

sub FETCH {
    my $self = shift;
    my $key = shift;

    if (exists $self->{$key}) {
	return $self->{$key};
    } else {
	my $value = $self->{Raw}{$key};
	if (ref($value) eq 'ARRAY') {
	    $value = Data::Grove::ParentList->new($value, $self);
	}
	return $value;
    }
}

sub FIRSTKEY {
    my $self = shift;
    my $raw = $self->{Raw};

    $self->{'__each_in_raw'} = 1;
    my $a = scalar keys %$raw;
    each %$raw;
}

sub NEXTKEY  {
    my $self = shift;
    my $raw = $self->{Raw};

    my ($key, $value);
    if ($self->{'__each_in_raw'}) {
	if (($key, $value) = each %$raw) {
	    return $key;
	}
	delete $self->{'__each_in_raw'};
	my $a = scalar keys %$self;
    }

    return each %$self;
}

sub EXISTS {
    my $self = shift;
    my $key = shift;

    return (exists $self->{Raw}{$key})
	|| (exists $self->{$key});
}


sub DELETE {
    my $self = shift;
    my $key = shift;

    if (exists $self->{$key}) {
	croak "can't delete \`Parent' or \`Raw' properties\n"
	    if ($key eq 'Parent' || $key eq 'Raw');
	delete $self->{$key};
    } else {
	delete $self->{'Raw'}{$key};
    }
}

sub CLEAR {
    my $self = shift;

    %{ $self->{Raw} } = ();
}

#--------------------------------------------------------------------------
# Data::Grove::ParentList
#--------------------------------------------------------------------------

package Data::Grove::ParentList;

use UNIVERSAL;

sub new {
    my $type = shift;
    my $raw = shift;
    my $parent = shift;

    if (UNIVERSAL::isa($raw, 'Data::Grove::ParentList')) {
        return $raw;
    }

    my $dummy = [];
    tie @$dummy, $type, $raw, $parent;
    return $dummy;
}

sub TIEARRAY {
    my $type = shift;

    return bless [ @_ ], $type;
}

sub FETCHSIZE {
    scalar @{$_[0][0]};
}

sub STORESIZE {
    $#{$_[0][0]} = $_[1]-1;
}

sub STORE {
    my $self = shift;
    my $index = shift;
    my $value = shift;

    # VALIDATE
    if (UNIVERSAL::isa($value, 'Data::Grove::Parent')) {
	$value = $value->{Raw};
    } elsif (UNIVERSAL::isa($value, 'Data::Grove::ParentList')) {
	$value = $value->[0];
    }
    $self->[0][$index] = $value;
}

sub FETCH {
    my $self = shift;
    my $index = shift;

    my $value = $self->[0][$index];
    if (defined $value) {
	if (ref($value)) {
	    return Data::Grove::Parent->new($value, $self->[1]);
	} else {
	    return Data::Grove::Parent->new({ Data => $value }, $self->[1]);
	}
    }

    return $value;
}

sub CLEAR {
    @{$_[0][0]} = ();
}

sub POP {
    pop(@{$_[0][0]});
}

sub PUSH {
    my $o = shift;

    foreach my $value (@_) {
	# VALIDATE
	if (UNIVERSAL::isa($value, 'Data::Grove::Parent')) {
	    $value = $value->{Raw};
	} elsif (UNIVERSAL::isa($value, 'Data::Grove::ParentList')) {
	    $value = $value->[0];
	}
    }	
    push(@{$o->[0]},@_);
}

sub SHIFT {
    shift(@{$_[0][0]});
}

sub UNSHIFT {
    my $o = shift;

    foreach my $value (@_) {
	# VALIDATE
	if (UNIVERSAL::isa($value, 'Data::Grove::Parent')) {
	    $value = $value->{Raw};
	} elsif (UNIVERSAL::isa($value, 'Data::Grove::ParentList')) {
	    $value = $value->[0];
	}
    }	
    unshift(@{$o->[0]},@_);
} 

sub SPLICE
{
    my $ob  = shift;                    
    my $sz  = $ob->FETCHSIZE;
    my $off = @_ ? shift : 0;
    $off   += $sz if $off < 0;
    my $len = @_ ? shift : $sz-$off;

    foreach my $value (@_) {
	# VALIDATE
	if (UNIVERSAL::isa($value, 'Data::Grove::Parent')) {
	    $value = $value->{Raw};
	} elsif (UNIVERSAL::isa($value, 'Data::Grove::ParentList')) {
	    $value = $value->[0];
	}
    }	
    return splice(@{$ob->[0]},$off,$len,@_);
}

#--------------------------------------------------------------------------
# Data::Grove
#--------------------------------------------------------------------------

package Data::Grove;

sub root {
    my $self = shift;

    return $self
	if !defined $self->{Parent};

    return $self->{Parent}->root(@_);
}

sub rootpath {
    my $self = shift;

    if (defined $self->{Parent}) {
        return ($self->{Parent}->rootpath, $self);
    } else {
        return ($self);
    }
}

sub add_magic {
    my $self = shift;
    my $parent = shift;

    return Data::Grove::Parent->new($self, $parent);
}

1;

__END__

=head1 NAME

Data::Grove::Parent - provide parent properties to Data::Grove objects

=head1 SYNOPSIS

 use Data::Grove::Parent;

 $root = $object->root;
 $rootpath = $object->rootpath;
 $tied = $object->add_magic([ $parent ]);

 $node = Data::Grove::Parent->new($hash [, $parent]);
 $node_list = Data::Grove::ParentList->new($array [, $parent]);

=head1 DESCRIPTION

Data::Grove::Parent is an extension to Data::Grove that adds
`C<Parent>' and `C<Raw>' properties to Data::Grove objects and methods
for returning the root node of a grove, a list of nodes between and
including the root node and the current node, and a method that
creates parented nodes.

Data::Grove::Parent works by creating a Perl ``tied'' object that
contains a parent reference (`C<Parent>') and a reference to the
original Data::Grove object (`C<Raw>').  Tying-magic is used so that
every time you reference the Data::Grove::Parent object it actually
references the underlying raw object.

When you retrieve a list or a property of the Raw object,
Data::Grove::Parent automatically adds magic to the returned list or
node.  This means you only call `add_magic()' once to create the first
Data::Grove::Parent object and then use the grove objects like you
normally would.

The most obvious use of this is so you don't have to call a
`C<delete>' method when you want to release a grove or part of a
grove; since Data::Grove and Data::Grove::Parent objects have no
cyclic references, Perl can garbage collect them normally.

A secondary use is to allow you to reuse grove or property set
fragments in multiple trees.  WARNING: Data::Grove currently does not
protect you from creating your B<own> cyclic references!  This could
lead to infinite loops if you don't take care to avoid them.

=head1 METHODS

=over 4

=item $object->root()

=item $object->rootpath()

`C<root()>' returns the root node if `C<$object>' is a
`C<Data::Grove::Parent>' object.  `C<rootpath()>' returns an array of
all the nodes between and including the root node and `C<$object>'.

=item $tied = $object->add_magic([ $parent ])

`C<add_magic()>' returns a C<Data::Grove::Parent> object with
`C<$object>' as it's `C<Raw>' object.  If `C<$parent>' is given, that
becomes the tied object's parent object.

=back

=head1 AUTHOR

Ken MacLeod, ken@bitsko.slc.ut.us

=head1 SEE ALSO

perl(1), Data::Grove(3)

=cut
