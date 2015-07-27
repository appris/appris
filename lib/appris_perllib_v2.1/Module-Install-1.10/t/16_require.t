#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use File::Spec;
use t::lib::Test;

plan tests => 1;

SCOPE: {
	eval "require Module::Install";
	ok !$@, "import succeeds: \$Module::Install::VERSION $Module::Install::VERSION";
	diag $@ if $@;
}
