#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 45;

# Let's see how MakeMaker behaves first

# normal run
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use ExtUtils::MakeMaker;
WriteMakefile(
	NAME => 'Foo',
	INC  => '-I/usr/include/',
);
END_DSL

	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc =~ m{/usr/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}

# added as ARGV
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use ExtUtils::MakeMaker;
WriteMakefile(
	NAME => 'Foo',
);
END_DSL

	ok( run_makefile_pl(run_params => ['INC=-I/usr/opt/include/']), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc =~ m{/usr/opt/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}

# combined
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use ExtUtils::MakeMaker;
WriteMakefile(
	NAME => 'Foo',
	INC  => '-I/usr/include/',
);
END_DSL

	ok( run_makefile_pl(run_params => ['INC=-I/usr/opt/include/']), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc !~ m{/usr/include/},     "INC is overriden";
	ok $inc && $inc =~ m{/usr/opt/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}

# Now test Module::Install

# normal run
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
cc_inc_paths  '/usr/include/';
WriteAll;
END_DSL

	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc =~ m{/usr/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}

# multiple inc paths
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
cc_inc_paths  '/usr/include/';
cc_inc_paths  '/usr/opt/include/';
WriteAll;
END_DSL

	ok( run_makefile_pl(), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc =~ m{/usr/include/},     "correct INC";
	ok $inc && $inc =~ m{/usr/opt/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}

# added as ARGV
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( run_makefile_pl(run_params => ['INC=-I/usr/opt/include/']), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc =~ m{/usr/opt/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}

# combined
SCOPE: {
	ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
cc_inc_paths  '/usr/include/';
WriteAll;
END_DSL

	ok( run_makefile_pl(run_params => ['INC=-I/usr/opt/include/']), 'build_dist' );
	my $file = makefile();
	ok( -f $file, 'Makefile exists' );
	my $content = _read($file);
	ok( $content,'file is not empty');
	my ($inc) = $content =~ /^INC\s*=\s*(.+)$/m;
	diag "INC: $inc" if $ENV{TEST_VERBOSE};
	ok $inc && $inc !~ m{/usr/include/},     "INC is overriden";
	ok $inc && $inc =~ m{/usr/opt/include/}, "correct INC";
	ok( kill_dist(), 'kill_dist' );
}
