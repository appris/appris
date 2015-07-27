#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

eval "use Module::Install::ExtraTests 0.007";
plan skip_all => "requires Module::Install::ExtraTests 0.007" if $@;
plan tests => 38;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
extra_tests;
WriteAll;
END_DSL

	ok( add_test(qw(t test.t)), 'added t' );
	ok( add_test(qw(xt author test.t)), 'added xt/author' );
	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	diag my ($testline) = $content =~ /^#\s*(test => .+)$/m if $ENV{TEST_VERBOSE};
	ok($content !~ /#\s*test => \{ TESTS=>.+xt\/\*\.t/, 'has no xt/*.t');
    my $res = make('test');
    ok( $res =~ /Result:\s*PASS/ );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
extra_tests;
WriteAll;
END_DSL

	ok( add_test(qw(t test.t)), 'added t' );
	ok( add_test(qw(xt author test.t)), 'added xt/author' );
	ok( build_dist(), 'build_dist' );
	rmdir dir(qw(inc .author)); # non-author-mode
	unlink makefile();
	ok( run_makefile_pl(), 'build_dist again' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	diag my ($testline) = $content =~ /^#\s*(test => .+)$/m if $ENV{TEST_VERBOSE};
	ok($content !~ /#\s*test => \{ TESTS=>.+xt\/\*\.t/, 'has no xt/*.t');
    my $res = make('test');
    ok( $res =~ /Result:\s*PASS/ );
	ok( kill_dist(), 'kill_dist' );
}

# cases with (undocumented) tests_recursive()

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
tests_recursive;
extra_tests;
WriteAll;
END_DSL

	ok( add_test(qw(t test.t)), 'added t' );
	ok( add_test(qw(xt author test.t)), 'added xt/author' );
	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	diag my ($testline) = $content =~ /^#\s*(test => .+)$/m if $ENV{TEST_VERBOSE};
	ok($content !~ /#\s*test => \{ TESTS=>.+xt\/\*\.t/, 'has no xt/*.t');
    my $res = make('test');
    ok( $res =~ /Result:\s*PASS/ );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
tests_recursive;
extra_tests;
WriteAll;
END_DSL

	ok( add_test(qw(t test.t)), 'added t' );
	ok( add_test(qw(xt author test.t)), 'added xt/author' );
	ok( build_dist(), 'build_dist' );
	rmdir dir(qw(inc .author)); # non-author-mode
	unlink makefile();
	ok( run_makefile_pl(), 'build_dist again' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	diag my ($testline) = $content =~ /^#\s*(test => .+)$/m if $ENV{TEST_VERBOSE};
	ok($content !~ /#\s*test => \{ TESTS=>.+xt\/\*\.t/, 'has no xt/*.t');
    my $res = make('test');
    ok( $res =~ /Result:\s*PASS/ );
	ok( kill_dist(), 'kill_dist' );
}
