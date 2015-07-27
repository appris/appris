#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More tests => 24;
require_ok( 'inc::Module::Install' );

my @version = qw{
	0       0
	1       1
	1.1     1.1
	1234    1234
	1.2_01  1.20001
	1.2.3   1.002003
	1.2.3_1 1.0020031
	5.8.1   5.008001
	5.8.10  5.00801
	5.10.0  5.01
};

while ( @version ) {
	my $in  = shift @version;
	my $out = shift @version;
	my $ver = Module::Install::_version($in);
	my $two = Module::Install::_version($ver);
	is( $ver, $out, "$in => $out pass 1 ok" );
	is( $two, $out, "$in => $out pass 2 ok" );
}

my @cmp = qw{
	0    1.2.3    1.002003
	-1   1.2.3    1.002004
	1    1.2.3    1.002002
};

while ( @cmp ) {
	my $want  = shift @cmp;
	my $left  = shift @cmp;
	my $right  = shift @cmp;
	my $have = Module::Install::_cmp(undef, $left, $right);
	is( $have, $want, "_cmp($left, $right) ok" );
}
