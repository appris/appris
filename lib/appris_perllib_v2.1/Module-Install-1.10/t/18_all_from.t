#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;
require ExtUtils::MakeMaker;

my $eumm = eval $ExtUtils::MakeMaker::VERSION;

plan tests => 18;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ /#\s*ABSTRACT => q\[A test module\]/, 'has abstract');
	ok($content =~ author_makefile_re("Foo Bar"), 'has author');
	ok($content =~ /#\s*VERSION => q\[3\.21\]/, 'has version');
	SKIP: {
		skip "requires ExtUtils::MakeMaker 6.31", 1 unless $eumm < 6.31;
		ok($content =~ /#\s*LICENSE => q\[perl\]/, 'has license');
	}
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
abstract      'overriden abstract';
author        'Bar Foo';
version       '0.01';
license       'MIT';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ /#\s*ABSTRACT => q\[overriden abstract\]/, 'has abstract');
	ok($content =~ author_makefile_re("Bar Foo"), 'has author');
	ok($content =~ /#\s*VERSION => q\[0\.01\]/, 'has version');
	SKIP: {
		skip "requires ExtUtils::MakeMaker 6.31", 1 if $eumm < 6.31;
		ok($content =~ /#\s*LICENSE => q\[mit\]/, 'has license');
	}
	ok( kill_dist(), 'kill_dist' );
}
