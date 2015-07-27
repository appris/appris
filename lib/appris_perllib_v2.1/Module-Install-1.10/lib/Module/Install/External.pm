package Module::Install::External;

# Provides dependency declarations for external non-Perl things

use strict;
use Module::Install::Base ();

use vars qw{$VERSION $ISCORE @ISA};
BEGIN {
	$VERSION = '1.10';
	$ISCORE  = 1;
	@ISA     = qw{Module::Install::Base};
}

sub requires_xs {
	my $self = shift;

	# First check for the basic C compiler
	$self->requires_external_cc;

	# We need a C compiler that can build XS files
	unless ( $self->can_xs ) {
		print "Unresolvable missing external dependency.\n";
		print "This package requires perl's header files.\n";
		print STDERR "NA: Unable to build distribution on this platform.\n";
		exit(0);
	}

	1;
}

sub requires_external_cc {
	my $self = shift;

	# We need a C compiler, use the can_cc method for this
	unless ( $self->can_cc ) {
		print "Unresolvable missing external dependency.\n";
		print "This package requires a C compiler.\n";
		print STDERR "NA: Unable to build distribution on this platform.\n";
		exit(0);
	}

	# Unlike some of the other modules, while we need to specify a
	# C compiler as a dep, it needs to be a build-time dependency.

	1;
}

sub requires_external_bin {
	my ($self, $bin, $version) = @_;
	if ( $version ) {
		die "requires_external_bin does not support versions yet";
	}

	# Load the package containing can_run early,
	# to avoid breaking the message below.
	$self->load('can_run');

	# Locate the bin
	print "Locating bin:$bin...";
	my $found_bin = $self->can_run( $bin );
	if ( $found_bin ) {
		print " found at $found_bin.\n";
	} else {
		print " missing.\n";
		print "Unresolvable missing external dependency.\n";
		print "Please install '$bin' seperately and try again.\n";
		print STDERR "NA: Unable to build distribution on this platform.\n";
		exit(0);
	}

	# Once we have some way to specify external deps, do it here.
	# In the mean time, continue as normal.

	1;
}

1;

__END__

=pod

=head1 NAME

Module::Install::External - Specify dependencies on external non-Perl things

=head1 DESCRIPTION

C<Module::Install::External> provides command that allow you to
declaratively specify a dependency on a program or system that is not
Perl.

The commands it provides are similar to those in L<Module::Install::Can>,
except that they implement an explicit dependency, in addition to just
working out if the particular thing is available.

=head1 COMMANDS

=head2 requires_xs

  requires_xs;

The C<requires_xs> command explicitly specifies that a C compiler and the
perl header files are required in order to build (at F<make>-time) the
distribution (specifically XS files).

It does not take any params, and aborts the F<Makefile.PL> execution in a
way that an automated installation or testing system will interpret as a
C<NA> ("not applicable to this platform") result.

This may be changed to an alternative abort result at a later time.

Returns true as a convenience.

=head2 requires_external_cc

  requires_external_cc;

The C<requires_external_cc> command explicitly specifies that a C compiler
is required in order to build (at F<make>-time) the distribution.

It does not take any params, and aborts the F<Makefile.PL> execution
in a way that an automated installation or testing system will interpret
as a C<NA> ("not applicable to this platform") result.

This may be changed to an alternative abort result at a later time.

Returns true as a convenience.

=head2 requires_external_bin

  requires_external_bin 'cvs';

The C<requires_external_bin> command takes the name of a system command
or program, similar to the C<can_run> command, except that
C<requires_external_bin> checks in a way that is a declarative explicit
dependency.

The takes a single param of the command/program name, and aborts the
C<Makefile.PL> execution in a way that an automated installation or
testing system will interpret as a C<NA> ("not applicable to this
platform") result.

Returns true as a convenience.

=head1 TO DO

Work out how to save the external dependency metadata, in agreement with
the larger Perl community.

Implement the agreed external dependency metadata solution.

=head1 AUTHORS

Adam Kennedy E<lt>adamk@cpan.orgE<gt>

=head1 SEE ALSO

L<Module::Install>

=head1 COPYRIGHT

Copyright 2006 Adam Kennedy.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut
