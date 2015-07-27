package t::lib::Test;

use strict;
use File::Spec   ();
use File::Remove ();
use File::Path ();
use Cwd;
use Config;

use vars qw{$VERSION @ISA @EXPORT $DIST};
BEGIN {
	$VERSION = '1.10';
	@ISA     = 'Exporter';
	@EXPORT  = qw{
		create_dist
		build_dist
		kill_dist
		run_makefile_pl
		add_file
		add_test
		_read
		file dir
		makefile
		make
		supports_capture
		capture_build_dist
		author_makefile_re
	};
	$DIST = '';
}

# Done in evals to avoid confusing Perl::MinimumVersion
eval( $] >= 5.006 ? <<'END_NEW' : <<'END_OLD' ); die $@ if $@;
sub _read {
	local *FH;
	open( FH, '<', $_[0] ) or die "open($_[0]): $!";
	my $string = do { local $/; <FH> };
	close FH or die "close($_[0]): $!";
	return $string;
}
END_NEW
sub _read {
	local *FH;
	open( FH, "< $_[0]"  ) or die "open($_[0]): $!";
	my $string = do { local $/; <FH> };
	close FH or die "close($_[0]): $!";
	return $string;
}
END_OLD

sub create_dist {
	$DIST = shift;
	my $opt  = shift || {};

	# Clear out any existing directory
	kill_dist( $DIST );

	my $home      = cwd;
	my $dist_path = dir();
	my $dist_lib  = dir('lib');
	mkdir($dist_path, 0777) or return 0;
	mkdir($dist_lib,  0777) or return 0;
	chdir($dist_path      ) or return 0;

	# Write the MANIFEST
	open( MANIFEST, '>MANIFEST' ) or return 0;
	print MANIFEST $opt->{MANIFEST} || <<"END_MANIFEST";
MANIFEST
Makefile.PL
lib/$DIST.pm
END_MANIFEST
	close MANIFEST;

	# Write the configure script
	open MAKEFILE_PL, '>Makefile.PL' or return 0;
	print MAKEFILE_PL $opt->{'Makefile.PL'} || <<"END_MAKEFILE_PL";
use inc::Module::Install 0.81;
name          '$DIST';
version       '0.01';
license       'perl';
requires_from 'lib/$DIST.pm';
requires      'File::Spec' => '0.79';
WriteAll;
END_MAKEFILE_PL
	close MAKEFILE_PL;

	# Write the module file
	open MODULE, ">lib/$DIST.pm" or return 0;
	print MODULE $opt->{"lib/$DIST.pm"} || <<"END_MODULE";
package $DIST;

=pod

=head1 NAME

$DIST - A test module

=cut

use 5.005;
use strict;

\$VERSION = '3.21';

use File::Spec 0.80;

=pod

=head1 AUTHORS

Foo Bar

=cut

1;

__END__

=head1 COPYRIGHT

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut
END_MODULE
	close MODULE;

	chdir $home or return 0;
	return 1;
}

sub file { File::Spec->catfile('t', $DIST . $$, @_) }
sub dir  { File::Spec->catdir('t', $DIST . $$, @_) }
sub makefile { file(@_, $^O eq 'VMS' ? 'Descrip.MMS' : 'Makefile' ) }

sub add_file {
	my $dist_path = dir();
	return 0 unless -d $dist_path;

	my $content = pop;
	my $file    = pop;
	my @subdir  = @_;
	my $dist_subdir = dir(@subdir);
	my $dist_file   = file(@subdir, $file);
	unless (-d $dist_subdir) {
		File::Path::mkpath($dist_subdir, 0, 0777) or return 0;
	}

	open FILE, "> $dist_file" or return 0;
	print FILE $content;
	close FILE;

	return 1;
}

sub add_test { add_file(@_, qq{print "1..1\nok 1\n";}) }

sub build_dist {
	my %params = @_;
	my $dist_path = dir();
	return 0 unless -d $dist_path;
	my $home = cwd;
	chdir $dist_path or return 0;
	my $X_MYMETA = $params{MYMETA} || '';
	local $ENV{X_MYMETA} = $X_MYMETA;

	my @run_params=@{ $params{run_params} || [] };
	my $ret = system($^X, "-I../../lib", "-I../../blib/lib", "Makefile.PL",@run_params);
	chdir $home or return 0;
	return $ret ? 0 : 1;
}

sub run_makefile_pl {
	my %params = @_;
	my $dist_path = dir();
	return 0 unless -d $dist_path;
	my $home = cwd;
	chdir $dist_path or return 1;
	my $X_MYMETA = $params{MYMETA} || '';
	local $ENV{X_MYMETA} = $X_MYMETA;

	my $run_params=join(' ',@{ $params{run_params} || [] });
	my $ret = system("$^X -I../../lib -I../../blib/lib Makefile.PL $run_params");
	#my $result=qx();
	chdir $home or return 0;
	return $ret ? 0 : 1;
}

sub kill_dist {
	my $dir = dir();
	return 1 unless -d $dir;
	windows_delay();
	File::Remove::remove( \1, $dir );
	windows_delay();
	return -d $dir ? 0 : 1;
}

sub windows_delay {
	return if $^O ne 'MSWin32';
	select undef, undef, undef, 0.1;
}

sub supports_capture {
	# stolen from ExtUtils::MakeMaker's test
	use ExtUtils::MM;

	# Unix, modern Windows and OS/2 from 5.005_54 up can handle 2>&1 
	# This makes our failure diagnostics nicer to read.
	return 1
		if (MM->os_flavor_is('Unix') or
			(MM->os_flavor_is('Win32') and !MM->os_flavor_is('Win9x')) or
			($] > 5.00554 and MM->os_flavor_is('OS/2')));
}

sub capture_build_dist {
	my %params = @_;
	my $dist_path = dir();
	return '' unless -d $dist_path;
	my $home = cwd;
	chdir $dist_path or return '';
	my $X_MYMETA = $params{MYMETA} || '';
	local $ENV{X_MYMETA} = $X_MYMETA;

	my @run_params=@{ $params{run_params} || [] };
	my $command = join ' ', $^X, "-I../../lib", "-I../../blib/lib", "Makefile.PL", @run_params;
	my $ret = `$command 2>&1`;
	chdir $home;
	return $ret;
}

sub extract_target {
	my $target  = shift;
	my $makefile = makefile();
	return '' unless -f $makefile;
	my $content = _read($makefile) or return '';
	my @lines;
	my $flag;
	foreach (split /\n/, $content) {
		if (/^$target\s*:/) { $flag++ }
		elsif (/^\S+\s*:/) { $flag = 0 }
		push @lines, $_ if $flag;
	}
	return wantarray ? @lines : join "\n", @lines;
}

sub make {
	my $target = shift || '';

	my $dist_path = dir();
	return '' unless -d $dist_path;
	my $home = cwd;
	chdir $dist_path or return '';

	my $make = $Config{make};
	my $ret = supports_capture()
		? `$make $target 2>&1`
		: `$make $target`;
	chdir $home;
	return $ret;
}

require ExtUtils::MakeMaker;
my $eumm = eval $ExtUtils::MakeMaker::VERSION;

sub author_makefile_re {
	my $author=shift;
	if ($eumm>=6.5702) {
		return qr/#\s*AUTHOR => \[q\[$author\]\]/;
	} else {
		return qr/#\s*AUTHOR => q\[$author\]/;
	}
}

1;
