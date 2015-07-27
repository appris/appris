#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use t::lib::Test;
use YAML::Tiny ();

plan tests => 14;

SCOPE: {
	ok( create_dist('Foo', { 'Makefile.PL' => <<'END_DSL' }), 'create_dist' );
use inc::Module::Install 0.82;

name          'Foo';
license       'perl';
author        'Foo Bar <foo@bar.com>';
all_from      'lib/Foo.pm';
requires      'perl'       => '5.008000';
test_requires 'Test::More' => '0.47';
no_index      'directory'  => qw{ t xt share inc };
install_share 'eg';
keywords      'kw1','kw 2';
keywords      'kw3';
license       'apache';

WriteAll;
END_DSL

	unlink file('META.yml');
	unlink file('MYMETA.yml');
	ok( mkdir(dir('eg')), 'created eg/' );
	ok( add_file('eg/sample', 'This is a sample'), 'added sample' );
	ok( mkdir(dir('t')), 'created t/' );
	ok( add_file('t/01_comile.t', <<'END_TEST'), 'added test' );
#!/usr/bin/perl

BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More tests => 2;

ok( $] >= 5.005, 'Perl version is new enough' );

use_ok( 'Foo', 'Loaded Foo.pm' );
END_TEST

	ok( build_dist(), 'build dist' );

	my $metafile = file('META.yml');
	ok( -f $metafile, 'META.yml created' );

	my $meta = YAML::Tiny::LoadFile($metafile);

	is_deeply(
		[ sort @{ $meta->{no_index}->{directory} } ],
		[ qw{ eg inc t } ],
		'no_index is ok',
	) or diag(
		"no_index: @{ $meta->{no_index}->{directory} }"
	);
	is_deeply(
		$meta->{keywords},
		[ 'kw1','kw 2','kw3'],
		'no_index is ok',
	) or diag(
		"no_index: @{ $meta->{no_index}->{directory} }"
	);

	is($meta->{license},'apache','license');
	is($meta->{resources}->{license},'http://www.apache.org/licenses/LICENSE-2.0.txt','license URL');

	my $makefile = makefile();
	ok( -f $makefile, 'Makefile created' );

	my $content = _read($makefile);
	ok( $content =~ /^#\s+PREREQ_PM/m, 'PREREQ_PM found' );

	ok( kill_dist(), 'kill dist' );
}
