package Module::Install::Admin::Bundle;

use strict;
use Module::Install::Base;

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.10';
	@ISA     = qw{Module::Install::Base};
}

sub bundle {
    my $self       = shift;
    my $bundle_dir = $self->_top->{bundle};

    require Cwd;
    require CPANPLUS::Backend;

    my $cwd = Cwd::getcwd();

    # This code is what we _should_ be doing, but CPANPLUS doesn't
    # let you have multiple Backends in one program.
    # my $cp   = CPANPLUS::Backend->new;
    #
    # Jos Boumans tells us that this is the best way to do what we want
    # It still scares me.
    my $cp      = CPANPLUS::Internals->_retrieve_id( CPANPLUS::Internals->_last_id )
               || CPANPLUS::Backend->new;
    my $conf    = $cp->configure_object;
    my $modtree = $cp->module_tree;

    $conf->set_conf( verbose   => 1 );
    $conf->set_conf( signature => 0 );
    $conf->set_conf( md5       => 0 );

    mkdir( $bundle_dir, 0777 );

    while ( my ( $name, $version ) = splice( @_, 0, 2 ) ) {
        my $mod = $cp->module_tree($name);
        if (not $mod) {
            warn "Warning: Could not find distribution for module $name on CPAN. Bundle it manually.\n";
            next;
        }

        if ( $mod->package_is_perl_core or $self->{already_bundled}{$mod->package} ) {
            next;
        }

        my $where = $mod->fetch( fetchdir => $bundle_dir, );
        unless ($where) {
            warn "Warning: Could not download distribution $bundle_dir. Download it manually.\n";
            next;
        }
        my $file = Cwd::abs_path($where);

        my $extract_result = $mod->extract(
            files      => [ $file ],
            extractdir => $bundle_dir,
        );

        unlink $file;
        unless ($extract_result) {
            warn "Warning: Could not extract distribution $bundle_dir. Extract it manually.\n";
            next;
        }

        $self->{already_bundled}{ $mod->package }++;
    }

    chdir $cwd;
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
