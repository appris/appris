package Module::Install::Run;

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.10';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

# eventually move the ipc::run / open3 stuff here.

1;

__END__

=pod

=encoding UTF-8

=head1 COPYRIGHT

Copyright 2008 - 2014 Adam Kennedy.

=head1 LICENSE

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
