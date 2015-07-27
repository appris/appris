package Module::Install::Deprecated;

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.10';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}





#####################################################################
# Previous API for Module::Install::Compoler

sub c_files {
	warn "c_files has been changed to cc_files to reduce confusion and keep all compiler commands as cc_";
	shift()->cc_files(@_);
}

sub inc_paths {
	warn "inc_paths has been changed to cc_inc_paths due to confusion between Perl and C";
	shift()->cc_inc_paths(@_);
}

sub lib_paths {
	warn "lib_paths has been changed to cc_lib_paths due to confusion between Perl and C";
	shift()->cc_lib_paths(@_);
}

sub lib_links {
	warn "lib_links has been changed to cc_lib_links due to confusion between Perl and C";
	shift()->cc_lib_links(@_);
}

sub optimize_flags {
	warn "optimize_flags has been changed to cc_optimize_flags for consistency reasons";
	shift()->cc_optimize_flags(@_);
}

1;

__END__

=pod

=head1 NAME

Module::Install::Deprecated - Warnings and help for deprecated commands

=head1 DESCRIPTION

One of the nicest features of L<Module::Install> is that as it improves,
there is no need to take into account user compatibility, because users do
not need to install L<Module::Install> itself.

As a result, the L<Module::Install> API changes at a faster rate than usual,
and this results in deprecated commands.

C<Module::Install::Deprecated> provides implementations of the deprecated
commands, so that when an author is upgrading their L<Module::Install> and
they are using a deprecated command they will be told that the command has
been deprecated, and what the author should use instead.

This extension should NEVER end up bundled into the distribution tarball.

=head1 COMMANDS

=head2 c_files

The C<c_files> command has been changed to C<cc_files> to reduce confusion
and keep all compiler commands within a consistent C<ff_foo> naming scheme.

=head2 inc_paths

The C<inc_paths> command has been changed to C<cc_inc_paths> due to
confusion between Perl and C.

=head2 lib_paths

The C<lib_paths> command has been changed to C<cc_lib_paths> due to confusion
between Perl and C.

=head2 lib_links

The C<lib_links> command has been changed to C<cc_lib_links> due to confusion
between Perl and C.

=head2 optimize_flags

The C<optimize_flags> command has been changed to C<cc_optimize_flags> for
consistency reasons.

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
