package Module::Install::Admin::Manifest;

use strict;
use Module::Install::Base;

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.16';
	@ISA     = qw{Module::Install::Base};
}

use Cwd ();
use File::Spec;

# XXX I really want this method in Module::Install::Admin::Makefile
# But you can't call across Admin modules. Audrey??
sub dist_preop {
	my ($self, $distdir) = @_;
	return if $self->check_manifest;
	print <<"END_MESSAGE";

It appears that your MANIFEST does not contain the same components that
are currently in the 'inc' directory. 

Please try running 'make manifest' and then run 'make dist' again.

Remember to use the MANIFEST.SKIP file to control things that should not
end up in your MANIFEST. See 'perldoc ExtUtils::Manifest' for details.

END_MESSAGE
	return if $self->prompt(
		'Do you *really* want to continue making a distribution?', 'n'
	) =~ /^[Yy]/;

	if ( -d $distdir ) {
		require File::Path;
		File::Path::rmtree($distdir);
	}

	exit(1);
}

# XXX Needs a refactoring.
sub check_manifest {
	my $self = shift;
	my $prefix = $self->_top->{prefix};
	my ($manifest, $manifest_path, $relative_path) = $self->_read_manifest or return;
	my $manifest_skip = "$manifest_path.SKIP";
	my @skip;

	if ( -f "$manifest_path.SKIP" ) {
		open SKIP, $manifest_skip 
			or die "Can't open $manifest_skip for input:\n$!";
		@skip = map {chomp; $_} <SKIP>;
		close SKIP;
	}

	my %manifest;
	for ( my $i = 0; $i < @$manifest; $i++ ) {
		my $path = $manifest->[$i];
		$path =~ s/\s.*//;
		$path =~ s/^\.[\\\/]//;
		$path =~ s/[\\\/]/\//g;
		next unless $path =~ m/^\Q$prefix\E\b/i;
		$manifest{$path} = \$manifest->[$i];
	}

	ADDLOOP:
	for my $pathname ( sort $self->_find_files($prefix) ) {
		$pathname = "$relative_path/$pathname" if length($relative_path);
		$pathname =~ s!//+!/!g;
		next unless -f $pathname;
		if ( defined $manifest{$pathname} ) {
			delete $manifest{$pathname};
		} else {
			for ( @skip ) {
				next ADDLOOP if $pathname =~ /$_/;
			}
			return 0;
		}
	}
	if ( keys %manifest ) {
		foreach ( keys %manifest ) {
			print "Found extra file $_\n";
		}
		return 0;
	}
	return 1;
}

sub _read_manifest {
	my $manifest = [];
	my $manifest_path = '';
	my $relative_path = '';
	my @relative_dirs = ();
	my $cwd = Cwd::getcwd();
	my @cwd_dirs = File::Spec->splitdir($cwd);
	while ( @cwd_dirs ) {
		last unless -f File::Spec->catfile(@cwd_dirs, 'Makefile.PL');
		my $path = File::Spec->catfile(@cwd_dirs, 'MANIFEST');
		if ( -f $path ) {
			$manifest_path = $path;
			last;
		}
		unshift @relative_dirs, pop(@cwd_dirs);
	}

	unless ( length($manifest_path) ) {
		warn "Can't locate the MANIFEST file for '$cwd'\n";
		return;
	}

	$relative_path = join '/', @relative_dirs if @relative_dirs;

	local *MANIFEST;
	open MANIFEST, $manifest_path 
		or die "Can't open $manifest_path for input:\n$!";
	@$manifest = map { chomp; $_ } <MANIFEST>;
	close MANIFEST;

	return ($manifest, $manifest_path, $relative_path);
}

# XXX I copied this from M::I::A::Find because I can't call that one. Please
# refactor/fix.
sub _find_files {
	my ($self, $file, $path) = @_;
	$path = '' if not defined $path;
	$file = "$path/$file" if length($path);
	if ( -f $file ) {
		return ( $file );
	} elsif ( -d $file ) {
		my @files = ();
		local *DIR;
		opendir(DIR, $file) or die "Can't opendir $file";
		while (defined(my $new_file = readdir(DIR))) {
			next if $new_file =~ /^(\.|\.\.)$/;
			push @files, $self->_find_files($new_file, $file);
		}
		return @files;
	}
	return ();
}

1;

__END__

=pod

=head1 COPYRIGHT

Copyright 2003, 2004 by
Audrey Tang E<lt>autrijus@autrijus.orgE<gt>,
Brian Ingerson E<lt>ingy@cpan.orgE<gt>

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut
