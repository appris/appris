#!/usr/bin/perl

# Tests for Module::Install::DSL

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More tests => 8;
use t::lib::Test;

# Load the DSL module
require_ok( 'inc::Module::Install::DSL' );

# Generate code from a simple dsl block
my $code = Module::Install::DSL::dsl2code(<<'END_DSL');
all_from lib/My/Module.pm
requires perl 5.008
requires Carp 0
requires Win32 if win32
test_requires Test::More
install_share
END_DSL

is( $code, <<'END_PERL', 'dsl2code generates the expected code' );
all_from 'lib/My/Module.pm';
requires 'perl', '5.008';
requires 'Carp', '0';
requires 'Win32' if win32;
test_requires 'Test::More';
install_share;
END_PERL





######################################################################
# Automatic dynamic vs static detection

# Automatically set static_config if there are no conditionals
my $static = Module::Install::DSL::dsl2code(<<'END_DSL');
all_from lib/My/Module.pm
requires perl 5.008
requires Carp 0
requires Win32
test_requires Test::More
install_share
END_DSL

is( $static, <<'END_PERL', 'dsl2code generates the expected code' );
all_from 'lib/My/Module.pm';
requires 'perl', '5.008';
requires 'Carp', '0';
requires 'Win32';
test_requires 'Test::More';
install_share;
static_config;
END_PERL





#####################################################################
# Full scan dist run

ok( create_dist( 'Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install::DSL 0.81;
name          Foo
version       0.01
license       perl
requires_from lib/Foo.pm
requires      File::Spec   0.79
END_DSL
ok( build_dist(), 'build_dist' );
ok( -f makefile() );
ok( -f file('META.yml') );
ok( kill_dist(), 'kill_dist' );
