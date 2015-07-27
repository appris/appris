#!perl

use strict;

BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More tests => 8;

my @package_name = (
	'package Foo::Bar;',
	'package Foo::Bar 1.23;',
	'package Foo::Bar { ... }',
	'package Foo::Bar 1.23 { ... }',
);

foreach (@package_name) {
	if (
		$_ =~ m/
		^ \s*
		package \s*
		([\w:]+)
		[\s|;]*
		/ixms
		)
	{
		my ($name, $module_name) = ($1, $1);
		$name =~ s{::}{-}g;
		is($module_name, 'Foo::Bar', "found module_name Foo::Bar in $_");
		is($name,        'Foo-Bar',  "found name Foo-Bar in $_");
	}
}


__END__

Perl 5.12 introduced:

	package Foo::Bar 1.23;
	...;

Perl 5.14 introduced:

	package Foo::Bar { ... }

and they can be combined as:

	package Foo::Bar 1.23 { ... }

The name_from regex doesn't support any of the above. It expects the package name to be followed by optional whitespace then a semicolon.

###

approx line 342 in package Module::Install::Metadata::name_from;












