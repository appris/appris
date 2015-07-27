#!/usr/bin/perl

# Check that install_share 

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use t::lib::Test;

plan tests => 22;

SCOPE: { # normal share dir/file
	ok( create_dist('Foo', { 'Makefile.PL' => <<'END_DSL' }), 'create_dist' );
use inc::Module::Install 0.82;
name 'Foo';
all_from 'lib/Foo.pm';
install_share;
WriteAll;
END_DSL

	ok( mkdir(dir('share')), 'mkdir share' );
	add_file('share/dist_file.txt', 'sample file');
	add_file('MANIFEST', <<'END_MANIFEST');
MANIFEST
Makefile.PL
lib/Foo.pm
share/dist_file.txt
END_MANIFEST

	ok( build_dist(), 'build_dist' );
	make();
	ok( !$?, 'make' );

	my $dir = dir(qw{ blib lib auto share dist Foo });
	ok( -d $dir, 'Found install_share in correct dist_dir location' );

	my $file = file(qw{ blib lib auto share dist Foo dist_file.txt });
	ok( -f $file, 'Found expected file in dist_dir location' );

	my $content = _read($file);
	ok( $content eq 'sample file', 'correct content' );

	ok( kill_dist(), 'kill_dist' );
}

SCOPE: { # normal share dir, but a file not listed in MANIFEST
	ok( create_dist('Foo', { 'Makefile.PL' => <<'END_DSL' }), 'create_dist' );
use inc::Module::Install 0.82;
name 'Foo';
all_from 'lib/Foo.pm';
install_share;
WriteAll;
END_DSL

	ok( mkdir(dir('share')), 'mkdir share' );
	add_file('share/dist_file.txt', 'sample file');
	add_file('MANIFEST', <<'END_MANIFEST');
MANIFEST
Makefile.PL
lib/Foo.pm
#share/dist_file.txt
END_MANIFEST

	ok( build_dist(), 'build_dist' );
	make();
	ok( !$?, 'make' );

	my $dir = dir(qw{ blib lib auto share dist Foo });
	ok( -d $dir, 'Found install_share in correct dist_dir location' );

	my $file = file(qw{ blib lib auto share dist Foo dist_file.txt });
	ok( !-f $file, 'Found expected file in dist_dir location' );

	ok( kill_dist(), 'kill_dist' );
}

SCOPE: { # file is listed in MANIFEST, but also in MANIFEST.SKIP
	ok( create_dist('Foo', { 'Makefile.PL' => <<'END_DSL' }), 'create_dist' );
use inc::Module::Install 0.82;
name 'Foo';
all_from 'lib/Foo.pm';
install_share;
WriteAll;
END_DSL

	ok( mkdir(dir('share')), 'mkdir share' );
	add_file('share/dist_file.txt', 'sample file');
	add_file('MANIFEST', <<'END_MANIFEST');
MANIFEST
Makefile.PL
lib/Foo.pm
share/dist_file.txt
END_MANIFEST

	add_file('MANIFEST.SKIP', <<'END_MANIFEST_SKIP');
share/dist_file.txt
END_MANIFEST_SKIP

	ok( build_dist(), 'build_dist' );
	make();
	ok( !$?, 'make' );

	my $dir = dir(qw{ blib lib auto share dist Foo });
	ok( -d $dir, 'Found install_share in correct dist_dir location' );

	my $file = file(qw{ blib lib auto share dist Foo dist_file.txt });
	ok( !-f $file, 'Found expected file in dist_dir location' );

	ok( kill_dist(), 'kill_dist' );
}
