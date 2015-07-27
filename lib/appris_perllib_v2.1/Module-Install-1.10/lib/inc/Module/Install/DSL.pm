package inc::Module::Install::DSL;

# This module ONLY loads if the user has manually installed their own
# installation of Module::Install, and are some form of MI author.
#
# It runs from the installed location, and is never bundled
# along with the other bundled modules.
#
# So because the version of this differs from the version that will
# be bundled almost every time, it doesn't have it's own version and
# isn't part of the synchronisation-checking.

use strict;
use vars qw{$VERSION};
BEGIN {
	# While this version will be overwritten when Module::Install
	# loads, it remains so Module::Install itself can detect which
	# version an author currently has installed.
	# This allows it to implement any back-compatibility features
	# it may want or need to.
	$VERSION = '1.10';
}

if ( -d './inc' ) {
	my $author = $^O eq 'VMS' ? './inc/_author' : './inc/.author';
	if ( -d $author ) {
		$Module::Install::AUTHOR = 1;
		require File::Path;
		File::Path::rmtree('inc');
	}
} else {
	$Module::Install::AUTHOR = 1;
}

unshift @INC, 'inc' unless $INC[0] eq 'inc';
require inc::Module::Install;
require Module::Install::DSL;

# Tie our import to the main one
sub import {
	goto &Module::Install::DSL::import;
}

1;

=pod

=head1 NAME

inc::Module::Install::DSL - Domain Specific Language for Module::Install

=head1 SYNOPSIS

  use inc::Module::Install::DSL 0.80;
  
  all_from       lib/ADAMK/Repository.pm
  requires       File::Spec            3.29
  requires       File::pushd           1.00
  requires       File::Find::Rule      0.30
  requires       File::Find::Rule::VCS 1.05
  requires       File::Flat            0
  requires       File::Remove          1.42
  requires       IPC::Run3             0.034
  requires       Object::Tiny          1.06
  requires       Params::Util          0.35
  requires       CPAN::Version         5.5
  test_requires  Test::More            0.86
  test_requires  Test::Script          1.03
  install_script adamk
  
  requires_external_bin svn

=head1 DESCRIPTION

One of the primary design goals of L<Module::Install> is to simplify
the creation of F<Makefile.PL> scripts.

Part of this involves the gradual reduction of any and all superfluous
characters, with the ultimate goal of requiring no non-critical
information in the file.

L<Module::Install::DSL> is a simple B<Domain Specific Language> based
on the already-lightweight L<Module::Install> command syntax.

The DSL takes one command on each line, and then wraps the command
(and its parameters) with the normal quotes and semi-colons etc to
turn it into Perl code.

=head1 SUPPORT

Bugs should be reported via the CPAN bug tracker at

L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Module-Install>

For other issues contact the author.

=head1 AUTHORS

Adam Kennedy E<lt>adamk@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright 2008 - 2014 Adam Kennedy.

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut
