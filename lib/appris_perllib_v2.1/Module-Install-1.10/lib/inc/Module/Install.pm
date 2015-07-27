package inc::Module::Install;

# This module ONLY loads if the user has manually installed their own
# installation of Module::Install, and are some form of MI author.
#
# It runs from the installed location, and is never bundled
# along with the other bundled modules.
#
# So because the version of this differs from the version that will
# be bundled almost every time, it doesn't have it's own version and
# isn't part of the synchronisation-checking.

use strict;
use vars qw{$VERSION};
BEGIN {
	# While this version will be overwritten when Module::Install
	# loads, it remains so Module::Install itself can detect which
	# version an author currently has installed.
	# This allows it to implement any back-compatibility features
	# it may want or need to.
	$VERSION = '1.10';
}

if ( -d './inc' ) {
	my $author = $^O eq 'VMS' ? './inc/_author' : './inc/.author';
	if ( -d $author ) {
		my $modified_at = (stat($author))[9];
		if ((time - $modified_at) > 24 * 60 * 60) {
			# inc is a bit stale; there may be a newer Module::Install
			_check_update($modified_at);
		}
		$Module::Install::AUTHOR = 1;
		require File::Path;
		File::Path::rmtree('inc');
	}
} else {
	$Module::Install::AUTHOR = 1;
}

unshift @INC, 'inc' unless $INC[0] eq 'inc';
local $^W;
require Module::Install;

sub _check_update {
	my $modified_at = shift;

	# XXX: We have several online services to get update information
	# including search.cpan.org. They are more reliable than the
	# 02packages.details.txt.gz on the local machine. We might be
	# better to depend on those services... but on which?

	my $cpan_version = 0;
	if (0) {  # XXX: should be configurable?
		my $url = "http://search.cpan.org/dist/Module-Install/META.yml";
		eval "require YAML::Tiny; 1" or return;

		if (eval "require LWP::UserAgent; 1") {
			my $ua = LWP::UserAgent->new(
				timeout   => 10,
				env_proxy => 1,
			);
			my $res = $ua->get($url);
			return unless $res->is_success;
			my $yaml = eval { YAML::Tiny::Load($res->content) } or return;
			$cpan_version = $yaml->{version};
		}
	}
	else {
		# If you don't want to rely on the net...
		require File::Spec;
		$cpan_version = _check_update_local($modified_at) or return;
	}

	# XXX: should die instead of warn?
	warn <<"WARN" if $cpan_version > $VERSION;
Newer version of Module::Install is available on CPAN.
CPAN:  $cpan_version
LOCAL: $VERSION
Please upgrade.
WARN
}

sub _check_update_local {
	my $modified_at = shift;

	return unless eval "require Compress::Zlib; 1";
	_require_myconfig_or_config() or return;
	my $file = File::Spec->catfile(
		$CPAN::Config->{keep_source_where},
		'modules',
		'02packages.details.txt.gz'
	);
	return unless -f $file;
#	return if (stat($file))[9] < $modified_at;

	my $gz = Compress::Zlib::gzopen($file, 'r') or return;
	my $line;
	while($gz->gzreadline($line)) {
		my ($cpan_version) = $line =~ /^Module::Install\s+(\S+)/ or next;
		return $cpan_version;
	}
	return;
}

# adapted from CPAN::HandleConfig
sub _require_myconfig_or_config {
	return 1 if $INC{"CPAN/MyConfig.pm"};
	local @INC = @INC;
	my $home = _home() or return;
	my $cpan_dir = File::Spec->catdir($home,'.cpan');
	return unless -d $cpan_dir;
	unshift @INC, $cpan_dir;
	eval { require CPAN::MyConfig };
	if ($@ and $@ !~ m#locate CPAN/MyConfig\.pm#) {
		warn "Error while requiring CPAN::MyConfig:\n$@\n";
		return;
	}
	return 1 if $INC{"CPAN/MyConfig.pm"};

	eval { require CPAN::Config; };
	if ($@ and $@ !~ m#locate CPAN/Config\.pm#) {
		warn "Error while requiring CPAN::Config:\n$@\n";
		return;
	}
	return 1 if $INC{"CPAN/Config.pm"};
	return;
}

# adapted from CPAN::HandleConfig
sub _home {
	my $home;
	if (eval {require File::HomeDir; 1}) {
		$home = File::HomeDir->can('my_dot_config')
			? File::HomeDir->my_dot_config
			: File::HomeDir->my_data;
		unless (defined $home) {
			$home = File::HomeDir->my_home
		}
	}
	unless (defined $home) {
		$home = $ENV{HOME};
	}
	$home;
}

1;

__END__

=pod

=head1 NAME

inc::Module::Install - Module::Install configuration system

=head1 SYNOPSIS

  use inc::Module::Install;

=head1 DESCRIPTION

This module first checks whether the F<inc/.author> directory exists,
and removes the whole F<inc/> directory if it does, so the module author
always get a fresh F<inc> every time they run F<Makefile.PL>.  Next, it
unshifts C<inc> into C<@INC>, then loads B<Module::Install> from there.

Below is an explanation of the reason for using a I<loader module>:

The original implementation of B<CPAN::MakeMaker> introduces subtle
problems for distributions ending with C<CPAN> (e.g. B<CPAN.pm>,
B<WAIT::Format::CPAN>), because its placement in F<./CPAN/> duplicates
the real libraries that will get installed; also, the directory name
F<./CPAN/> may confuse users.

On the other hand, putting included, for-build-time-only libraries in
F<./inc/> is a normal practice, and there is little chance that a
CPAN distribution will be called C<Something::inc>, so it's much safer
to use.

Also, it allows for other helper modules like B<Module::AutoInstall>
to reside also in F<inc/>, and to make use of them.

=head1 AUTHORS

Audrey Tang E<lt>autrijus@autrijus.orgE<gt>

=head1 COPYRIGHT

Copyright 2003, 2004 Audrey Tang E<lt>autrijus@autrijus.orgE<gt>.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut
