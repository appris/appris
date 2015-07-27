#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 35;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL', <<"END_DSL" }), 'create_dist' );
use strict;
use ExtUtils::MakeMaker;
WriteMakefile( NAME => 'Foo' );
END_DSL

	ok( add_file('MyFoo.pm.PL', '1;'), 'created .pm.PL');
	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file);
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($target) = $content =~ /^MyFoo.pm :: (.+)$/m;
	diag $target if $ENV{TEST_VERBOSE};
	ok( $target =~ /^MyFoo.pm.PL pm_to_blib/, 'has MyFoo target' );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL', <<"END_DSL" }), 'create_dist' );
use strict;
use ExtUtils::MakeMaker;
WriteMakefile(
	NAME => 'Foo',
	PL_FILES => { 'MyFoo.pm.PL' => 'MyFoo.pm' },
);
END_DSL

	ok( add_file('MyFoo.pm.PL', '1;'), 'created .pm.PL');
	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file);
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($target) = $content =~ /^MyFoo.pm :: (.+)$/m;
	diag $target if $ENV{TEST_VERBOSE};
	ok( $target =~ /^MyFoo.pm.PL pm_to_blib/, 'has MyFoo target' );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo'), 'create_dist' );
	ok( add_file('MyFoo.pm.PL', '1;'), 'created .pm.PL');
	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file);
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($target) = $content =~ /^MyFoo.pm :: (.+)$/m;
	diag $target if $ENV{TEST_VERBOSE};
	ok( $target =~ /^MyFoo.pm.PL pm_to_blib/, 'has MyFoo target' );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL', <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
makemaker_args(
	PL_FILES => { 'MyFoo.pm.PL' => 'MyFoo.pm' },
);
WriteAll;
END_DSL

	ok( add_file('MyFoo.pm.PL', '1;'), 'created .pm.PL');
	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file);
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($target) = $content =~ /^MyFoo.pm :: (.+)$/m;
	diag $target if $ENV{TEST_VERBOSE};
	ok( $target =~ /^MyFoo.pm.PL pm_to_blib/, 'has MyFoo target' );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL', <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteMakefile(
	PL_FILES => { 'MyFoo.pm.PL' => 'MyFoo.pm' },
);
END_DSL

	ok( add_file('MyFoo.pm.PL', '1;'), 'created .pm.PL');
	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file);
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($target) = $content =~ /^MyFoo.pm :: (.+)$/m;
	diag $target if $ENV{TEST_VERBOSE};
	ok( $target =~ /^MyFoo.pm.PL pm_to_blib/, 'has MyFoo target' );
	ok( kill_dist(), 'kill_dist' );
}
