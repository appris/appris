#!/usr/bin/perl

use strict;
BEGIN {
        $|  = 1;
        $^W = 1;
}

use Test::More tests => 14;

require_ok( 'Module::Install::Metadata' );
my $metadata = Module::Install::Metadata->new;





#####################################################################
# Simple Checks

my @tests = qw{
	5.1        5.001
	5.6        5.006
	5.8        5.008
	5.10       5.010
	5.11       5.011
	5.8.8      5.008008
	5.10.0     5.010
	5.008005   5.008005
};
while ( @tests ) {
	my $in  = shift @tests;
	my $out = shift @tests;
	is(
		$metadata->_perl_version($in),
		$out,
		"->_perl_version($in)",
	);
}





#####################################################################
# Practical Approach

$metadata->perl_version('5.008');
is($metadata->perl_version, 5.008);

$metadata->perl_version('5.8.1');
is($metadata->perl_version, 5.008001);

$metadata->perl_version('5.008001');
is($metadata->perl_version, 5.008001);

$metadata->perl_version('5.10.1');
is($metadata->perl_version, 5.010001);

$metadata->perl_version('5.8');
is($metadata->perl_version, 5.008);
