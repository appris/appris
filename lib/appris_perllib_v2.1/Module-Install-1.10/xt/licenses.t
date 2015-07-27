#!perl

use strict;

BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More;
use Module::Install::Metadata;
use Module::Runtime qw( require_module );
use Test::Requires qw(
	Software::License
	Module::Find
);

my @licenses = Module::Find::findsubmod('Software::License');

plan tests => 1 * @licenses;

foreach my $license (sort @licenses) {

SKIP: {
		local $@;
		eval { require_module($license) };
		if ($@) {
			skip "Can't load $license: $@", 1;
			next;
		}

		# Custom is not defined hence skip here
		if ($license eq 'Software::License::Custom') {
			skip 'Software::License::Custom is not predefined', 1;
		}

		my $name = $license->name;
		my $meta = $license->meta_name;

		unless ($meta) {
			skip "$license has no meta_name", 1;
			next;
		}
		if ($meta =~ m/open_source|restrictive|unrestricted|unknown/) {
			skip "$license meta_name is $meta", 1;
			next;
		}

		# $meta =~ s/_\d+$//;

		my $got = Module::Install::Metadata::__extract_license($name);
		ok $got =~ /^$meta/, $name;

		# should also test license urls?
		my $url = $license->url;
	}
}

done_testing();

__END__

