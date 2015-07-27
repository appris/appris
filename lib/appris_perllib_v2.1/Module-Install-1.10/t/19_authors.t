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

plan tests => 32;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
author        'ishigaki';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ author_makefile_re("ishigaki"), 'has one author');
	my $metafile = file('META.yml');
	ok(-f $metafile);
	my $meta = Parse::CPAN::Meta::LoadFile($metafile);
	is_deeply($meta->{author}, [qw(ishigaki)]);
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
author        'ishigaki', 'charsbar';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ author_makefile_re("ishigaki, charsbar"), 'has two authors');
	my $metafile = file('META.yml');
	ok(-f $metafile);
	my $meta = Parse::CPAN::Meta::LoadFile($metafile);
	is_deeply($meta->{author}, [qw(ishigaki charsbar)]);
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
authors       'ishigaki', 'charsbar';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ author_makefile_re("ishigaki, charsbar"), 'has two authors');
	my $metafile = file('META.yml');
	ok(-f $metafile);
	my $meta = Parse::CPAN::Meta::LoadFile($metafile);
	is_deeply($meta->{author}, [qw(ishigaki charsbar)]);
	ok( kill_dist(), 'kill_dist' );
}

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<"END_DSL" }), 'create_dist' );
use inc::Module::Install 0.81;
name          'Foo';
author        'ishigaki';
author        'charsbar';
perl_version  '5.005';
all_from      'lib/Foo.pm';
WriteAll;
END_DSL

	ok( build_dist(), 'build_dist' );
	my $file = makefile();
	ok(-f $file);
	my $content = _read($file);
	ok($content, 'file is not empty');
	ok($content =~ author_makefile_re("ishigaki, charsbar"), 'has two authors');
	my $metafile = file('META.yml');
	ok(-f $metafile);
	my $meta = Parse::CPAN::Meta::LoadFile($metafile);
	is_deeply($meta->{author}, [qw(ishigaki charsbar)]);
	ok( kill_dist(), 'kill_dist' );
}
