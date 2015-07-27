#line 1
package Module::Install::With;

# See POD at end for docs

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.10';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

#line 21

#####################################################################
# Installer Target

# Are we targeting ExtUtils::MakeMaker (running as Makefile.PL)
sub eumm {
	!! ($0 =~ /Makefile.PL$/i);
}

# You should not be using this, but we'll keep the hook anyways
sub mb {
	!! ($0 =~ /Build.PL$/i);
}





#####################################################################
# Testing and Configuration Contexts

#line 53

sub interactive {
	# Treat things interactively ONLY based on input
	!! (-t STDIN and ! automated_testing());
}

#line 71

sub automated_testing {
	!! $ENV{AUTOMATED_TESTING};
}

#line 90

sub release_testing {
	!! $ENV{RELEASE_TESTING};
}

sub author_context {
	!! $Module::Install::AUTHOR;
}





#####################################################################
# Operating System Convenience

#line 118

sub win32 {
	!! ($^O eq 'MSWin32');
}

#line 135

sub winlike {
	!! ($^O eq 'MSWin32' or $^O eq 'cygwin');
}

1;

#line 163
