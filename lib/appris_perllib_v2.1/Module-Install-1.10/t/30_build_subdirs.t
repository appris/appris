#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use t::lib::Test;
use YAML::Tiny;

plan tests => 20 * 2;

foreach my $command ('', "build_subdirs 'bar';") {
	SCOPE: {
		ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist Foo');
use inc::Module::Install;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
$command
WriteAll;
END_DSL

		ok( mkdir(dir('bar'), 0777),     'created bar directory');
		ok( mkdir(dir('bar/lib'), 0777), 'created bar/lib directory');
		ok( mkdir(dir('bar/t'), 0777),   'created bar/t directory');
		ok( add_file('bar/MANIFEST', <<'END_MANIFEST'), 'created MANIFEST');
MANIFEST
Makefile.PL
lib/Bar.pm
t/load.t
END_MANIFEST

		ok( add_file('bar/Makefile.PL', <<'END_DSL'), 'created Makefile.PL');
use inc::Module::Install;
name          'Bar';
abstract      'bar';
author        'foobar';
perl_version  '5.005';
version_from  'lib/Bar.pm';
WriteAll;
END_DSL

		ok( add_file('bar/lib/Bar.pm', <<'END_PERL'), 'created Bar.pm');
package Bar;

$VERSION = '0.01';

1;
END_PERL

		ok( add_file('bar/t/load.t', <<'END_T'), 'created load.t');
use Test::More tests => 1;
use_ok('Bar');
END_T

		ok( supports_capture() ? capture_build_dist() : build_dist(), 'build_dist' );

		my $makefile_foo = makefile();
		ok( -f $makefile_foo, 'has Makefile for Foo' );
		my $content_foo = _read($makefile_foo);
		ok( $content_foo =~ /DISTNAME\s*=>\s*q\[Foo\]/, 'content is correct');
		my $meta_yml_foo = file('META.yml');
		ok( -f $meta_yml_foo, 'has META.yml for Foo');
		my $meta_foo = YAML::Tiny::LoadFile($meta_yml_foo);
		ok( $meta_foo->{name} eq 'Foo', 'META.yml is correct' );
		ok( $meta_foo->{version} eq '3.21', 'META.yml is correct' );

		my $makefile_bar = file('bar/Makefile');
		ok( -f $makefile_bar, 'has Makefile for Bar' );
		my $content_bar = _read($makefile_bar);
		ok( $content_bar =~ /DISTNAME\s*=>\s*q\[Bar\]/, 'content is correct');
		my $meta_yml_bar = file('bar/META.yml');
		ok( -f $meta_yml_bar, 'has META.yml for Bar');
		my $meta_bar = YAML::Tiny::LoadFile($meta_yml_bar);
		ok( $meta_bar->{name} eq 'Bar', 'META.yml is correct' );
		ok( $meta_bar->{version} eq '0.01', 'META.yml is correct' );

		ok( kill_dist(), 'kill_dist' );
	}
}
