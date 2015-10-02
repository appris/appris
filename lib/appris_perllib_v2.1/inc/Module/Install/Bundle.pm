package Module::Install::Bundle;

use strict;
use File::Spec;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.16';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

sub auto_bundle {
    my $self = shift;

    return $self->_install_bundled_dists unless $self->is_admin;

    # Flatten array of arrays into a single array
    my @core = map @$_, map @$_, grep ref, $self->requires;

    $self->bundle(@core);
}

sub bundle {
    my $self = shift;

    return $self->_install_bundled_dists unless $self->is_admin;

    $self->admin->bundle(@_);
}

sub auto_bundle_deps {
    my $self = shift;

    return $self->_install_bundled_dists unless $self->is_admin;

    # Flatten array of arrays into a single array
    my @core = map @$_, map @$_, grep ref, $self->requires;
    while (my ($name, $version) = splice(@core, 0, 2)) {
        next unless $name;
         $self->bundle_deps($name, $version);
    }
}

sub bundle_deps {
    my ($self, $pkg, $version) = @_;

    return $self->_install_bundled_dists unless $self->is_admin;

    my $deps = $self->admin->scan_dependencies($pkg);
    if (scalar keys %$deps == 0) {
        # Probably a user trying to install the package, read the dependencies from META.yml
        %$deps = ( map { $$_[0] => undef } (@{$self->requires()}) );
    }
    foreach my $key (sort keys %$deps) {
        $self->bundle($key, ($key eq $pkg) ? $version : 0);
    }
}

sub _install_bundled_dists {
    my $self = shift;

    # process bundle only the first time this function is called
    return if $self->{bundle_processed};

    $self->makemaker_args->{DIR} ||= [];

    # process all dists bundled in inc/BUNDLES/
    my $bundle_dir = $self->_top->{bundle};
    foreach my $sub_dir (glob File::Spec->catfile($bundle_dir,"*")) {

        next if -f $sub_dir;

        # ignore dot dirs/files if any
        my $dot_file = File::Spec->catfile($bundle_dir,'\.');
        next if index($sub_dir, $dot_file) >= $[;

        # EU::MM can't handle Build.PL based distributions
        if (-f File::Spec->catfile($sub_dir, 'Build.PL')) {
            warn "Skipped: $sub_dir has Build.PL.";
            next;
        }

        # EU::MM can't handle distributions without Makefile.PL
        # (actually this is to cut blib in a wrong directory)
        if (!-f File::Spec->catfile($sub_dir, 'Makefile.PL')) {
            warn "Skipped: $sub_dir has no Makefile.PL.";
            next;
        }
        push @{ $self->makemaker_args->{DIR} }, $sub_dir;
    }

    $self->{bundle_processed} = 1;
}

1;

__END__

=pod

=head1 NAME

Module::Install::Bundle - Bundle distributions along with your distribution

=head1 SYNOPSIS

Have your Makefile.PL read as follows:

  use inc::Module::Install;
  
  name      'Foo-Bar';
  all_from  'lib/Foo/Bar.pm';
  requires  'Baz' => '1.60';
  
  # one of either:
  bundle    'Baz' => '1.60';
  # OR:
  auto_bundle;
  
  WriteAll;

=head1 DESCRIPTION

Module::Install::Bundle allows you to bundle a CPAN distribution within your
distribution. When your end-users install your distribution, the bundled
distribution will be installed along with yours, unless a newer version of
the bundled distribution already exists on their local filesystem.

While bundling will increase the size of your distribution, it has several
benefits:

  Allows installation of bundled distributions when CPAN is unavailable
  Allows installation of bundled distributions when networking is unavailable
  Allows everything your distribution needs to be packaged in one place

Bundling differs from auto-installation in that when it comes time to
install, a bundled distribution will be installed based on the distribution
bundled with your distribution, whereas with auto-installation the distribution
to be installed will be acquired from CPAN and then installed.

=head1 METHODS

=over 4

=item * auto_bundle()

Takes no arguments, will bundle every distribution specified by a C<requires()>.
When you, as a module author, do a C<perl Makefile.PL> the latest versions of
the distributions to be bundled will be acquired from CPAN and placed in
F<inc/BUNDLES/>.

=item * bundle($name, $version)

Takes a list of key/value pairs specifying a distribution name and version
number. When you, as a module author, do a perl Makefile.PL the distributions
that you specified with C<bundle()> will be acquired from CPAN and placed in
F<inc/BUNDLES/>.

=item * bundle_deps($name, $version)

Same as C<bundle>, except that all dependencies of the bundled modules are
also detected and bundled.  To use this function, you need to declare the
minimum supported perl version first, like this:

    perl_version( '5.005' );

=item * auto_bundle_deps

Same as C<auto_bundle>, except that all dependencies of the bundled
modules are also detected and bundled. This function has the same constraints as bundle_deps.

=back

=head1 BUGS

Please report any bugs to (patches welcome):

    http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Module-Install

=head1 AUTHORS

Audrey Tang E<lt>autrijus@autrijus.orgE<gt>

Documentation by Adam Foxson E<lt>afoxson@pobox.comE<gt>

=head1 COPYRIGHT

Copyright 2003, 2004, 2005 by Audrey Tang E<lt>autrijus@autrijus.orgE<gt>.

This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut
