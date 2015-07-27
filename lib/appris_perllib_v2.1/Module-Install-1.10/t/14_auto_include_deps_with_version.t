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
use vars qw{ $PREREQ_PM $MIN_PERL_VERSION $BUILD_REQUIRES };

plan skip_all => 'your perl is new enough to have File::Spec 3.30 in core' if $] > 5.010000;
plan skip_all => 'your File::Spec is not new enough for this test' if $File::Spec::VERSION < 3.30;

plan tests => 8;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
version       '0.01';
author        'Someone';
license       'perl';
perl_version  '5.005';
requires_from 'lib/Foo.pm';
test_requires 'File::Spec' => 0.6;
auto_include_deps;
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file_spec = file(qw(inc File Spec.pm));
	ok( !-f $file_spec, 'File::Spec is not included');
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
test_requires 'File::Spec' => 3.30;
auto_include_deps;
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file_spec = file(qw(inc File Spec.pm));
	ok( -f $file_spec, 'File::Spec is included');
	ok( kill_dist(), 'kill_dist' );
}
