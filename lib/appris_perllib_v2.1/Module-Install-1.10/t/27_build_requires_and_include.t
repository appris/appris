#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use YAML::Tiny;
use t::lib::Test;

plan tests => 14;

SCOPE: {  # see if test_requires works correctly
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use strict;
use warnings;
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';

test_requires 'Test::More' => 9999;

WriteAll;
END_DSL
	ok( build_dist, 'build dist' );
	my $meta_yml = file('META.yml');
	ok(-f $meta_yml, 'has META.yml');

	my $meta = YAML::Tiny::LoadFile($meta_yml);
	ok $meta->{build_requires}{'Test::More'}, 'Test::More is listed in build_requires';

	my $makefile = makefile();
	ok(-f $makefile, 'has Makefile');

	my $content = _read($makefile);
	ok($content =~ /^#\s+(PREREQ_PM|BUILD_REQUIRES)\s*=>\s*{[^}]+Test::More=>q\[9999\]/m, 'Test::More is listed in PREREQ_PM|BUILD_REQUIRES in Makefile');

	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {  # include removes Test::More from the build_requires in META.yml
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use strict;
use warnings;
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';

test_requires 'Test::More' => 9999;
include 'Test::More';

WriteAll;
END_DSL
	ok( build_dist, 'build dist' );
	my $meta_yml = file('META.yml');
	ok(-f $meta_yml, 'has META.yml');

	my $meta = YAML::Tiny::LoadFile($meta_yml);
	ok !$meta->{build_requires}{'Test::More'}, 'Test::More is not listed in build_requires';

	my $makefile = makefile();
	ok(-f $makefile, 'has Makefile');

	my $content = _read($makefile);
	ok($content !~ /^#\s+(PREREQ_PM|BUILD_REQUIRES)\s*=>\s*{[^}]+Test::More=>q\[9999\]/m, 'Test::More is not listed in PREREQ_PM|BUILD_REQUIRES in Makefile');
	ok( kill_dist(), 'kill_dist' );
}
