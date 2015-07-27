#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use Parse::CPAN::Meta;
use t::lib::Test;

plan tests => 24;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
makemaker_args 'INSTALLDIRS' => '~/local/';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	my ($installdirs) = $content =~ /^INSTALLDIRS\s*=\s*(.+)$/m;
	diag "INSTALLDIRS: $installdirs" if $ENV{TEST_VERBOSE};
	ok( $installdirs eq '~/local/', 'correct INSTALLDIRS');
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
installdirs   '~/local/';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	my ($installdirs) = $content =~ /^INSTALLDIRS\s*=\s*(.+)$/m;
	diag "INSTALLDIRS: $installdirs" if $ENV{TEST_VERBOSE};
	ok( $installdirs eq '~/local/', 'correct INSTALLDIRS');
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
installdirs   '~/local/';
makemaker_args 'INSTALLDIRS' => '~/old/';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	my ($installdirs) = $content =~ /^INSTALLDIRS\s*=\s*(.+)$/m;
	diag "INSTALLDIRS: $installdirs" if $ENV{TEST_VERBOSE};
	ok( $installdirs eq '~/local/', 'correct INSTALLDIRS');
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
install_as_core;
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	my ($installdirs) = $content =~ /^INSTALLDIRS\s*=\s*(.+)$/m;
	diag "INSTALLDIRS: $installdirs" if $ENV{TEST_VERBOSE};
	ok( $installdirs eq 'perl', 'correct INSTALLDIRS');
	ok( kill_dist(), 'kill_dist' );
}
