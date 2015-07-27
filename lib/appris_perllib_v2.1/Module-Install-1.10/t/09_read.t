#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 7;

SCOPE: {
	ok( create_dist('Foo'), 'create_dist' );
	ok( build_dist(), 'built dist' );
	require_ok( file('inc/Module/Install.pm') );

	my $file = file('lib/Foo.pm');
	ok( -f $file, "Found test file '$file'" );
	my $pod = Module::Install::_readpod($file);
	is($pod, <<'END_POD', "_readpod($file)" );
=head1 NAME

Foo - A test module

=head1 AUTHORS

Foo Bar

=head1 COPYRIGHT

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

END_POD

	my $perl = Module::Install::_readperl($file);
	is($perl, <<'END_PERL', "_readperl($file)" );
package Foo;

use 5.005;
use strict;

$VERSION = '3.21';

use File::Spec 0.80;

1;
END_PERL

	ok( kill_dist(), 'kill dist' );
}
