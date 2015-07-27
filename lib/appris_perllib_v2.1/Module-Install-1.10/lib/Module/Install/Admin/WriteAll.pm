package Module::Install::Admin::WriteAll;

use strict;
use Module::Install::Base;

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.10';
	@ISA     = qw{Module::Install::Base};
}

sub WriteAll {
	my ($self, %args) = @_;
	$self->load('Makefile');
	if ( $args{check_nmake} ) {
		$self->load($_) for qw(Makefile check_nmake can_run get_file);
	}
}

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
