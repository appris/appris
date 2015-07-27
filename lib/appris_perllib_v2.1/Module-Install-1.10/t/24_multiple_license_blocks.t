#!/usr/bin/perl

use strict;
BEGIN {
	$|  = 1;
	$^W = 1;
}

use Test::More tests => 3;
use Module::Install::Metadata;

my %p = _setup();

for (
	"$p{header}$p{author}$p{copyright}$p{license}$p{footer}",
	"$p{header}$p{author}$p{license}$p{copyright}$p{footer}",
	"$p{header}$p{license}$p{author}$p{copyright}$p{footer}",
) {
	my $license = Module::Install::Metadata::_extract_license($_);
	ok defined $license && $license eq 'perl', "my license is $license";
}

sub _setup {
	my %parts;

	$parts{header} =<<'POD';
=head1 NAME

Win32::UTCFileTime - Get/set UTC file times with stat/utime on Win32

=head1 SYNOPSIS

... (snip) ...

POD

	$parts{author} = <<'POD';
=head1 AUTHOR

Steve Hay E<lt>shay@cpan.orgE<gt>

POD

	$parts{copyright} = <<'POD';
=head1 COPYRIGHT

Copyright (C) 2003-2007 Steve Hay. All rights reserved.

Portions Copyright (C) 2001 Jonathan M Gilligan. Used with permission.

Portions Copyright (C) 2001 Tony M Hoyle. Used with permission.

POD

	$parts{copyright} = <<'POD';
=head1 COPYRIGHT

Copyright (C) 2003-2007 Steve Hay. All rights reserved.

Portions Copyright (C) 2001 Jonathan M Gilligan. Used with permission.

Portions Copyright (C) 2001 Tony M Hoyle. Used with permission.

POD

	$parts{license} = <<'POD';
=head1 LICENCE

This module is free software; you can redistribute it and/or modify it
under the
same terms as Perl itself, i.e. under the terms of either the GNU
General Public
License or the Artistic License, as specified in the F<LICENCE> file.

POD

	$parts{footer} = <<'POD';
=cut
POD

	return %parts;
}
