#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 4;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use Module::Install 0.81;  # should use "use inc::Module::Install"!
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
WriteAll;
END_DSL

	if ( supports_capture() ) {
		my $error = capture_build_dist();
		ok $?, 'build failed';
		ok $error =~ /Please invoke Module::Install with/, 'correct error';
		diag $error if $ENV{TEST_VERBOSE};
	}
	else {
		ok !build_dist(), "build_dist failed";
		SKIP: {
			skip 'this platform does not support 2>&1', 1;
		}
	}

	ok( kill_dist(), 'kill_dist' );
}
