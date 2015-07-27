package Module::Install::With;

# See POD at end for docs

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.10';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

=pod

=head1 NAME

Module::Install::With - find environnement for Module::Install

=cut

#####################################################################
# Installer Target

# Are we targeting ExtUtils::MakeMaker (running as Makefile.PL)
sub eumm {
	!! ($0 =~ /Makefile.PL$/i);
}

# You should not be using this, but we'll keep the hook anyways
sub mb {
	!! ($0 =~ /Build.PL$/i);
}





#####################################################################
# Testing and Configuration Contexts

=pod

=head2 interactive

The C<interactive> function tests for an install that has a user present
(or at least, one in which it is reasonable for us to present prompts
and other similar types of things).

Returns true if in an interactive environment, or false otherwise.

=cut

sub interactive {
	# Treat things interactively ONLY based on input
	!! (-t STDIN and ! automated_testing());
}

=pod

=head2 automated_testing

Are we currently running in an automated testing environment, such as
CPAN Testers.

This is primarily a cleaner and more human-readable equivalent of
checking $ENV{AUTOMATED_TESTING} yourself, but may be improved in line
with best practices at a later date.

=cut

sub automated_testing {
	!! $ENV{AUTOMATED_TESTING};
}

=pod

=head2 release_testing

Are we currently running in an release testing environment. That is,
are we in the process of running in a potential highly-intensive and
high dependency bloat testing process prior to packaging a module for
release.

This is primarily a cleaner and more human-readable equivalent of
checking $ENV{RELEASE_TESTING} yourself, but may be improved in line
with best practices at a later date.

=cut

sub release_testing {
	!! $ENV{RELEASE_TESTING};
}

sub author_context {
	!! $Module::Install::AUTHOR;
}





#####################################################################
# Operating System Convenience

=pod

=head2 win32

The C<win32> function tests if the Makefile.PL is currently running in a
native Microsoft Windows Perl, such as ActivePerl or Strawberry Perl.

This is primarily a cleaner and more human-readable equivalent of
checking C<$^O eq 'MSWin32'> yourself, but may be improved in line
with best practices at a later date.

=cut

sub win32 {
	!! ($^O eq 'MSWin32');
}

=pod

=head2 winlike

The C<winlike> function tests if the Makefile.PL is currently running
in a Microsoft Windows Perl, under either cygwin or a native Win32 Perl.

This is primarily a cleaner and more human-readable equivalent of
checking C<$^O eq 'MSWin32' or $^O eq 'cygwin'>yourself, but may be
improved in line with best practices at a later date.

=cut

sub winlike {
	!! ($^O eq 'MSWin32' or $^O eq 'cygwin');
}

1;

=pod

=head1 SEE ALSO

L<Module::Install>

=head1 AUTHORS

Adam Kennedy E<lt>adamk@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright 2007 - 2012 Adam Kennedy.

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=cut
