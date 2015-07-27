package Module::Install::Admin;

use strict 'vars';
use File::Path           ();
use inc::Module::Install ();

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.10';
	@ISA     = 'Module::Install';
}

=pod

=head1 NAME

Module::Install::Admin - Author-side manager for Module::Install

=head1 SYNOPSIS

In a B<Module::Install> extension module:

    sub extension_method {
        my $self = shift;
        $self->admin->some_method(@args);
    }

As an one-liner:

    % perl "-MModule::Install::Admin" -e'&some_method(@args);'

The two snippets above are really shorthands for

    $some_obj->some_method(@args)

where C<$some_obj> is the singleton object of a class under the
C<Module::Install::Admin::*> namespace that provides the method
C<some_method>.  See L</METHODS> for a list of built-in methods.

=head1 DESCRIPTION

This module implements the internal mechanism for initializing,
including and managing extensions, and should only be of interest to
extension developers; it is I<never> included under a distribution's
F<inc/> directory, nor are any of the B<Module::Install::Admin::*>
extensions.

For normal usage of B<Module::Install>, please see L<Module::Install>
and L<Module::Install/"COOKBOOK / EXAMPLES"> instead.

=head2 Bootstrapping

When someone runs a F<Makefile.PL> that has C<use inc::Module::Install>,
and there is no F<inc/> in the current directory, B<Module::Install>
will load this module bootstrap itself, through the steps below:

=over 4

=item *

First, F<Module/Install.pm> is POD-stripped and copied from C<@INC> to
F<inc/>.  This should only happen on the author's side, never on the
end-user side.

=item *

Reload F<inc/Module/Install.pm> if the current file is somewhere else.
This ensures that the included version of F<inc/Module/Install.pm> is
always preferred over the installed version.

=item *

Look at F<inc/Module/Install/*.pm> and load all of them.

=item *

Set up a C<main::AUTOLOAD> function to delegate missing function calls
to C<Module::Install::Admin::load> -- again, this should only happen
at the author's side.

=item *

Provide a C<Module::Install::purge_self> function for removing included
files under F<inc/>.

=back

=head1 METHODS

=cut

sub import {
	my $class = shift;
	my $self  = $class->new( _top => Module::Install->new, @_ );
	local $^W;
	*{caller(0) . "::AUTOLOAD"} = sub {
		no strict 'vars';
		$AUTOLOAD =~ /([^:]+)$/ or die "Cannot load";
		return if uc($1) eq $1;
		my $obj = $self->load($1) or return;
		unshift @_, $obj;
		goto &{$obj->can($1)};
	};
}

sub new {
	my ($class, %args) = @_;
	return $class->SUPER::new(
		%{$args{_top}}, %args,
		extensions  => undef,
		pathnames   => undef,
	);
}

sub init {
	my $self = shift;
	$self->copy($INC{"$self->{path}.pm"} => $self->{file});

	unless ( grep { $_ eq $self->{prefix} } @INC ) {
		unshift @INC, $self->{prefix};
	}
 	delete $INC{"$self->{path}.pm"};

	local $^W;
	do "$self->{path}.pm";
}

sub copy {
	my ($self, $from, $to) = @_;

	my @parts = split('/', $to);
	File::Path::mkpath([ join('/', @parts[ 0 .. $#parts-1 ])]);

	chomp $to;

	local (*FROM, *TO, $_);
	open FROM, "< $from" or die "Can't open $from for input:\n$!";
	open TO,   "> $to"   or die "Can't open $to for output:\n$!";
	print TO "#line 1\n";

	my $content;
	my $in_pod;

	while ( <FROM> ) {
		if ( /^=(?:b(?:egin|ack)|head\d|(?:po|en)d|item|(?:ove|fo)r)/ ) {
			$in_pod = 1;
		} elsif ( /^=cut\s*\z/ and $in_pod ) {
			$in_pod = 0;
			print TO "#line $.\n";
		} elsif ( ! $in_pod ) {
			print TO $_;
		}
	}

	close FROM;
	close TO;

	print "include $to\n";
}

# scan through our target to find
sub load_all_extensions {
	my $self = shift;
	unless ($self->{extensions}) {
		$self->{extensions} = [];
		foreach my $inc (@INC) {
			next if ref($inc) or $inc eq $self->{prefix};
			$self->load_extensions("$inc/$self->{path}", $self->{_top});
		}
	}
	return @{$self->{extensions}};
}

sub load {
	my ($self, $method, $copy) = @_;

	my @extobj;
	foreach my $obj ($self->load_all_extensions) {
		next unless defined &{ref($obj)."::$method"};
		my $is_admin = (ref($obj) =~ /^\Q$self->{name}::$self->{dispatch}::/);
		# Don't ever include admin modules, and vice versa.
		# $copy = 0 if $XXX and $is_admin;
		push @extobj, $obj if $copy xor $is_admin;
	}
	unless ( @extobj ) {
		die "Cannot find an extension with method '$method'";
	}

	# XXX - do we need to reload $obj from the new location?
	my $obj = $self->pick($method, \@extobj);
	$self->copy_package(ref($obj)) if $copy;

	return $obj;
}

# Copy a package to inc/, with its @ISA tree. $pathname is optional.
sub copy_package {
	my ($self, $pkg, $pathname) = @_;
	return unless ($pathname ||= $self->{pathnames}{$pkg});

	my $file = $pkg; $file =~ s!::!/!g;
	$file = "$self->{prefix}/$file.pm";
	return if -f $file; # prevents infinite recursion

	$self->copy($pathname => $file);
	foreach my $pkg (@{"$pkg\::ISA"}) {
		$self->copy_package($pkg);
	}
}

sub pick {
	# determine which name to load
	my ($self, $method, $objects) = @_;

	# XXX this whole thing needs to be discussed
	return $objects->[0] unless $#{$objects} > 0 and -t STDIN;

	# sort by last modified time
	@$objects = map { $_->[0] }
	            sort { $a->[1] <=> $b->[1] }
	            map { [ $_ => -M $self->{pathnames}{ref($_)} ] } @$objects;

	print "Multiple extensions found for method '$method':\n";
	foreach my $i ( 1 .. @$objects ) {
		print "\t$i. ", ref($objects->[$i-1]), "\n";
	}

	while ( 1 ) {
		print "Please select one [1]: ";
		chomp(my $choice = <STDIN>);
		$choice ||= 1;
		return $objects->[$choice-1] if $choice > 0 and $choice <= @$objects;
		print "Invalid choice.  ";
	}
}

sub delete_package {
	my ($self, $pkg) = @_;

	# expand to full symbol table name if needed
	unless ( $pkg =~ /^main::.*::$/ ) {
		$pkg = "main$pkg"   if     $pkg =~ /^::/;
		$pkg = "main::$pkg" unless $pkg =~ /^main::/;
		$pkg .= '::'        unless $pkg =~ /::$/;
	}

	my($stem, $leaf) = $pkg =~ m/(.*::)(\w+::)$/;
	my $stem_symtab = *{$stem}{HASH};
	return unless defined $stem_symtab and exists $stem_symtab->{$leaf};

	# free all the symbols in the package
	my $leaf_symtab = *{$stem_symtab->{$leaf}}{HASH};
	foreach my $name (keys %$leaf_symtab) {
		next if $name eq "$self->{dispatch}::";
		undef *{$pkg . $name};
	}

	# delete the symbol table
	foreach my $name (keys %$leaf_symtab) {
		next if $name eq "$self->{dispatch}::";
		delete $leaf_symtab->{$name};
	}
}

sub AUTOLOAD {
	goto &{shift->autoload};
}

sub DESTROY { }

1;

__END__

=pod

=head1 SEE ALSO

L<Module::Install>

=head1 AUTHORS

Audrey Tang E<lt>autrijus@autrijus.orgE<gt>

=head1 COPYRIGHT

Copyright 2003, 2004 by Audrey Tang E<lt>autrijus@autrijus.orgE<gt>.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut
