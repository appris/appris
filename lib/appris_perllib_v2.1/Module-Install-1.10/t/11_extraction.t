#!/usr/bin/perl

use strict;
BEGIN {
        $|  = 1;
        $^W = 1;
}

use Test::More tests => 16;

require_ok( 'Module::Install::Metadata' );

SCOPE: {
	my @links=Module::Install::Metadata::_extract_bugtracker('L<http://rt.cpan.org/test>');
	is_deeply(
		\@links,
		[ 'http://rt.cpan.org/test' ],
		'1 bugtracker extracted',
	) or diag(
		"bugtrackers: @links"
	);
}

SCOPE: {
	my @links=Module::Install::Metadata::_extract_bugtracker('L<http://rt.cpan.org/test1> L<http://rt.cpan.org/test1>');
	is_deeply(
		\@links,
		[ 'http://rt.cpan.org/test1' ],
		'1 bugtracker extracted (2 links)',
	) or diag(
		"bugtrackers: @links"
	);
}

SCOPE: {
	my @links=Module::Install::Metadata::_extract_bugtracker('L<http://rt.cpan.org/test1> L<http://rt.cpan.org/test2>');
	is_deeply(
		[ sort @links ],
		[ 'http://rt.cpan.org/test1', 'http://rt.cpan.org/test2' ],
		'2 bugtrackers extracted',
	) or diag(
		"bugtrackers: @links"
	);
}

SCOPE: {
	my @links=Module::Install::Metadata::_extract_bugtracker('L<http://search.cpan.org/test1>');
	is_deeply(
		\@links,
		[ ],
		'0 bugtrackers extracted',
	) or diag(
		"bugtrackers: @links"
	);
}

SCOPE: {
	my @links=Module::Install::Metadata::_extract_bugtracker('L<http://github.com/marcusramberg/mojomojo/issues>');
	is_deeply(
		\@links,
		[ 'http://github.com/marcusramberg/mojomojo/issues' ],
		'1 bugtracker (github.com) extracted',
	) or diag(
		"bugtrackers: @links"
	);
}

SCOPE: {
	my @links=Module::Install::Metadata::_extract_bugtracker('L<http://code.google.com/p/www-mechanize/issues/list>');
	is_deeply(
		\@links,
		[ 'http://code.google.com/p/www-mechanize/issues/list' ],
		'1 bugtracker (code.google.com) extracted',
	) or diag(
		"bugtrackers: @links"
	);
}




SCOPE: {
	my $l=Module::Install::Metadata::_extract_license("=head1 Copyright\nunder the same terms as the perl programming language\n=cut\n");
		is($l, 'perl', 'Perl license detected',
	);
}

SCOPE: {
        my $text="=head1 LICENSE

This is free software, you may use it and distribute it under
the same terms as Perl itself.

=head1 SEE ALSO

test

=cut
";
	my $l=Module::Install::Metadata::_extract_license($text);
		is($l, 'perl', 'Perl license detected',
	);
}

SCOPE: {
        my $text="=head1 COPYRIGHTS

This module is distributed under the same terms as Perl itself.

=cut
";
	my $l=Module::Install::Metadata::_extract_license($text);
		is($l, 'perl', 'Perl license detected',
	);
}

SCOPE: {
	my $l=Module::Install::Metadata::_extract_license("=head1 COPYRIGHT\nAs LGPL license\n=cut\n");
		is($l, 'lgpl', 'LGPL detected',
	);
}

SCOPE: {
        my $text=<<'EOT';
=head1 COPYRIGHT AND LICENCE

... is free software; you can redistribute it and/or modify it under
the terms of Perl itself, that is to say, under the terms of either:

=over 4

=item *

The GNU General Public License as published by the Free Software Foundation;
either version 2, or (at your option) any later version, or

=item *

The "Artistic License" which comes with Perl.

=back

=cut
EOT
	my $l=Module::Install::Metadata::_extract_license($text);
		is($l, 'perl', 'Perl license detected',
	);
}


SCOPE: {
        my $text=<<'EOT';
=head1 COPYRIGHT AND LICENCE

Copyright (C) 2010

This library is free software; you can redistribute it and/or modify it under the terms of the Artistic License 2.0. For details, see the full text of the license at http://opensource.org/licenses/artistic-license-2.0.php.

=cut
EOT
	my $l=Module::Install::Metadata::_extract_license($text);
		is($l, 'artistic_2', 'Artistic 2.0 license detected',
	);
}



SCOPE: {
	my $version=Module::Install::Metadata::_extract_perl_version("use 5.10.0;");
		is($version, '5.10.0', 'perl 5.10.0 detected',
	);
}

SCOPE: {
	my $version=Module::Install::Metadata::_extract_perl_version("   use 5.010;");
		is($version, '5.010', 'perl 5.10.0 detected',
	);
}

SCOPE: {
	my $version=Module::Install::Metadata::_extract_perl_version("use strict;");
		is($version, undef, 'no perl prereq',
	);
}
