#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use t::lib::Test;
require ExtUtils::MakeMaker;
use vars qw{ $PREREQ_PM $MIN_PERL_VERSION $BUILD_REQUIRES };

plan tests => 19;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
requires_from 'lib/Foo.pm';
WriteAll;
END_DSL

	ok( run_makefile_pl(run_params=>['PREREQ_PRINT >test']), 'build_dist' );
	my $file = file('test');
	ok( -f $file);
	my $content = _read($file);
	ok( $content, 'file is not empty');
	ok( $content =~ s/^.*\$PREREQ_PM = \{/\$PREREQ_PM = {/s,'PREREQ_PM found');
	eval ($content);
	ok( !$@,'correct content');
	ok( $PREREQ_PM->{'File::Spec'} eq '0.80', 'correct requirement' );
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
perl_version  '5.005';
all_from      'lib/Foo.pm';
test_requires_from 't/test.t';
WriteAll;
END_DSL

	ok( mkdir(dir('t')), 'created t/' );
	ok( add_file('t/test.t' => <<'END_TEST'), 'added test');
use strict;
use warnings;
use Test::More tests => 1;
use File::Spec 0.80;
ok("ok");
END_TEST

	ok( run_makefile_pl(run_params=>['PREREQ_PRINT >test']), 'build_dist' );
	my $file = file('test');
	ok( -f $file);
	my $content = _read($file);
	ok( $content, 'file is not empty');
	ok( $content =~ s/^.*\$PREREQ_PM = \{/\$PREREQ_PM = {/s,'PREREQ_PM found');
	eval ($content);
	ok( !$@,'correct content');
    if ( eval($ExtUtils::MakeMaker::VERSION) < 6.55_03 ) {
        ok( exists $PREREQ_PM->{'File::Spec'});
        ok( !exists $BUILD_REQUIRES->{'File::Spec'});
    } else { #best to check both because user can have any version
        ok( exists $BUILD_REQUIRES->{'File::Spec'});
        ok( !exists $PREREQ_PM->{'File::Spec'});
    }
	ok( kill_dist(), 'kill_dist' );
}
