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

plan skip_all => 'test requires ExtUtils::MakeMaker >= 6.17' if eval($ExtUtils::MakeMaker::VERSION) < 6.17;

plan tests => 12;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
WriteAll(sign => 0);
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content !~ /#\s*SIGN => q\[[01]\]/, 'has no sign');
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
WriteAll(sign => 1);
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ /#\s*SIGN => q\[1\]/, 'has sign');

	# XXX: might be better test `$Config{make} distsign` here
	# but it's neither safe nor portable...

	ok( kill_dist(), 'kill_dist' );
}
