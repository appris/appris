#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 6;

SCOPE: { # PREOP/POSTOP
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';

makemaker_args(dist => {
	PREOP  => 'my_preop',
	POSTOP => 'my_postop',
});

WriteAll;
END_DSL
	ok( build_dist(), 'build_dist' );
	my $makefile = makefile();
	ok(-f $makefile, 'has Makefile');
	my $content = _read($makefile);
	ok( $content =~ /my_preop/, 'has PREOP' );
	ok( $content =~ /my_postop/, 'has POSTOP' );
	ok( kill_dist(), 'kill_dist' );
}
