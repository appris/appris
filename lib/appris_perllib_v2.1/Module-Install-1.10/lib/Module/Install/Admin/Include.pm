package Module::Install::Admin::Include;

use strict;
use Module::Install::Base;

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.10';
	@ISA     = qw{Module::Install::Base};
}

sub include {
	my $self = shift;
	foreach my $rv ( $self->admin->glob_in_inc($_[0]) ) {
		$self->admin->copy_package(@$rv);
		my @build_requires;
		foreach (@{ $self->build_requires || [] }) {
			next if $_->[0] eq $rv->[0];
			push @build_requires, $_;
		}
		$self->Meta->{values}{build_requires} = \@build_requires;
	}
}

sub include_deps {
	my ($self, $module, $version) = @_;
	my $deps = $self->admin->scan_dependencies($module, $self->perl_version, $version) or return;
	foreach my $key ( sort keys %$deps ) {
		$self->include($key);
	}
}

sub auto_include {
	my $self = shift;
	foreach my $module (
		map  { $_->[0] }
		map  { @$_     }
		grep { $_      }
		$self->build_requires
	) {
		$self->include($module);
	}
}

sub auto_include_deps {
	my $self = shift;
	foreach my $module (
		map  { $_  }
		map  { @$_ }
		grep { $_  }
		$self->build_requires
	) {
		my ($name, $version) = @{$module};
		$self->include_deps($name, $version);
	}
}

=pod

=head1 NAME

Module::Install::Admin::Include - include methods for Module::Install

=head2 auto_include_dependent_dists

Grabs everything in this module's build_requires and attempts to
include everything (at the whole distribution level) recursively.

=cut

sub auto_include_dependent_dists {
	my $self = shift;
	foreach my $module (
		map  { $_->[0] }
		map  { @$_     }
		grep { $_      }
		$self->build_requires
	) {
		$self->include_dependent_dists($module);
	}
}

=pod

=head2 include_dependent_dists $package

Given a module package name, recursively include every package that
module needs.

=cut

sub include_dependent_dists {
	my $self = shift;
	my $pkg  = shift;
	return unless $pkg;
	return if $self->{including_dep_dist}->{ $self->_pkg_to_dist($pkg) }++;
	$self->include_one_dist($pkg);
	foreach my $mod ( @{ $self->_dist_to_mods( $self->_pkg_to_dist($pkg) ) } ) {
		my $deps = $self->admin->scan_dependencies($mod) or return;
		foreach my $key ( sort grep { $_ } keys %$deps ) {
			$self->include_dependent_dists($key);
		}
	}
}

=pod

=head2 include_one_dist $module

Given a module name, C<$module>, figures out which modules are in the
dist containing that module and copies all those files to ./inc. I bet
there's a way to harness smarter logic from L<PAR>.

=cut

sub include_one_dist {
	my $self = shift;
	my @mods = $self->_dist_to_mods( $self->_pkg_to_dist($_[0]) );
	foreach my $pattern ( grep { $_ } @mods ) {
		foreach my $rv ( $self->admin->glob_in_inc($pattern) ) {
			$self->admin->copy_package(@$rv);
			my @build_requires;
			foreach (@{ $self->build_requires || [] }) {
				next if $_->[0] eq $rv->[0];
				push @build_requires, $_;
			}
			$self->Meta->{values}{build_requires} = \@build_requires;
		}
	}
}

=pod

=for private _pkg_to_dist $modname

Given a module name, returns the file on CPAN containing
its latest version.

=cut

sub _pkg_to_dist {
	require CPAN;
	my $mod = CPAN::Shell->expand( Module => $_[1] ) or return;
	$mod->cpan_file;
}

=pod

=for private _dist_to_mods $distname

Takes the output of CPAN::Module->cpan_file and return all the modules
that CPAN.pm knows are in that dist. There's probably a better way using CPANPLUS

=cut

sub _dist_to_mods {
	CPAN::Shell->expand( Distribution => $_[1] )->containsmods;
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
