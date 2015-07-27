#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 18;

SCOPE: {  # runtime error
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use strict;
use warnings;
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';

unknown_func('args');

warn "after unknown func";

WriteAll;
END_DSL
	if ( supports_capture() ) {
		my $error = capture_build_dist();
		ok $?, 'build fails';
		ok( $error =~ /Unknown function is found/, 'correct error');
		ok( $error !~ /after unknown func/, 'no bogus warning');
		diag $error if $ENV{TEST_VERBOSE};
	}
	else {
		ok( !build_dist, 'build fails' );
		SKIP : {
			skip 'this platform does not support 2>&1', 2;
		}
	}
	my $file = makefile();
	ok(!-f $file, 'Makefile is not created');
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {  # Bareword not allowed while "strict subs" in use
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use strict;
use warnings;
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';

unknown_func;

warn "after unknown func";

WriteAll;
END_DSL
	if ( supports_capture() ) {
		my $error = capture_build_dist();
		ok $?, 'build fails';
		ok( $error =~ /Bareword .+ not allowed/, 'correct error');
		ok( $error !~ /after unknown func/, 'no bogus warning');
		diag $error if $ENV{TEST_VERBOSE};
	}
	else {
		ok( !build_dist, 'build fails' );
		SKIP : {
			skip "this platform does not support 2>&1", 2;
		}
	}
	my $file = makefile();
	ok(!-f $file, 'Makefile is not created');
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {  # String found where operator expected
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use strict;
use warnings;
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';

unknown_func 'args';

warn "after unknown func";

WriteAll;
END_DSL
	if ( supports_capture() ) {
		my $error = capture_build_dist();
		ok $?, 'build fails';
		ok( $error =~ /String found where operator expected/, 'correct error');
		ok( $error !~ /after unknown func/, 'no bogus warning');
		diag $error if $ENV{TEST_VERBOSE};
	}
	else {
		ok( !build_dist, 'build fails' );
		SKIP : {
			skip "this platform does not support 2>&1", 2;
		}
	}
	my $file = makefile();
	ok(!-f $file, 'Makefile is not created');
	ok( kill_dist(), 'kill_dist' );
}
