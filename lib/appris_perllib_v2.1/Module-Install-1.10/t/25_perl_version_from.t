#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;
use Parse::CPAN::Meta;

plan tests => 36;

my @statements = (
	'use 5.005',
	'use v5.5.0',
	'require 5.005',
	'require v5.5.0',
);
for my $statement (@statements) {
SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name     'Foo';
version  '0.01';
license  'perl';
perl_version_from 'lib/Foo.pm';
WriteAll;
END_DSL

	ok( add_file('lib/Foo.pm', <<"END_PM") );
package Foo;
use strict;
$statement;
1;
END_PM

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok($content, 'file is not empty');
SKIP: {
    skip 'requires ExtUtils::MakeMaker > 6.48', 1 if eval($ExtUtils::MakeMaker::VERSION) < 6.48;
	ok($content =~ /#\s*MIN_PERL_VERSION => q\[5\.005\]/, 'has perl_version');
}
	my $meta = file('META.yml');
	ok( -f $meta, 'META.yml exists' );
	my $yaml = Parse::CPAN::Meta::LoadFile($meta);
	is( $yaml->{requires}->{perl}, '5.005', 'META has correct perl version requirement' );
	ok( kill_dist(), 'kill_dist' );
}}
