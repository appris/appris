package Module::Install::Run;

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.16';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

# eventually move the ipc::run / open3 stuff here.

1;
