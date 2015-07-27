package Module::Install::Admin::Compiler;

use strict;
use Module::Install::Base;
use File::Remove  ();
use Devel::PPPort ();

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.10';
	@ISA     = qw{Module::Install::Base};
}

sub ppport {
	my $self   = shift;
	my $file   = shift || 'ppport.h';
	if ( -f $file ) {
		# Update the file to a newer version
		File::Remove::remove($file);
	}

	# Install the file (and remove on realclean)
	Devel::PPPort::WriteFile( $file ) or die "Failed to write $file";
	$self->realclean_files( $file );
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
