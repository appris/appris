package Module::Install::Admin::WriteAll;

use strict;
use Module::Install::Base;

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.16';
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
