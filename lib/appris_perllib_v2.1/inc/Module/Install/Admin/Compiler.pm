package Module::Install::Admin::Compiler;

use strict;
use Module::Install::Base;
use File::Remove  ();
use Devel::PPPort ();

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.16';
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
