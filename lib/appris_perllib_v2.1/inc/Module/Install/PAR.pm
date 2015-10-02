package Module::Install::PAR;

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.16';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

=head1 NAME

Module::Install::PAR - Module::Install Support for PAR::Dist packages

=head1 SYNOPSIS

To offer your users the possibility to install binaries if we cannot
compile an XS version of the module, you could use this simplistic stub:

    use inc::Module::Install;
    
    name            'Foo';
    all_from        'lib/Foo.pm';
    
    # Which CPAN directory do we fetch binaries from?
    par_base        'SMUELLER';
    
    unless ( can_xs ) {
        my $okay = extract_par( fetch_par );
        if (not $okay) {
            die "No compiler and no binary package found. Aborting.\n";
        }
    }
    
    WriteAll;

=head1 DESCRIPTION

This module adds a couple of directives to Module::Install
related to installing and creating PAR::Dist distributions.

=head2 par_base

This directive sets the CPAN ID from whose CPAN directory to
fetch binaries from. For example, you can choose to download
binaries from http://www.cpan.org/authors/id/S/SM/SMUELLER/
or its ftp counterpart by writing:

  par_base 'SMUELLER';

By default, the name of the file to fetch is generated from
the distribution name, its version, your platform name and your
perl version concatenated with dashes.

The directive, however, takes an optional second
argument which specifies the name of the file to fetch.
(Though C<par_base> does not fetch files itself, see below.)

  par_base 'SMUELLER', 'foo';

Once C<fetch_par> is called, the file 'foo' will be downloaded
from SMUELLER's CPAN directory. (It doesn't exist.)

The second argument could be used to fetch platform-agnostic
binaries:

  par_base 'SMUELLER', "Some-Distribution-0.01.par";

(Documentation TODO: Use the previously defined distribution
name and version in example.)

=cut

sub par_base {
    my ($self, $base, $file) = @_;
    my $class     = ref($self);
    my $inc_class = join('::', @{$self->_top}{qw(prefix name)});
    my $ftp_base;

    if ( defined $base and length $base ) {
        if ( $base =~ m!^(([A-Z])[A-Z])[-_A-Z]+\Z! ) {
            $self->{mailto} = "$base\@cpan.org";
            $ftp_base = "ftp://ftp.cpan.org/pub/CPAN/authors/id/$2/$1/$base";
            $base     = "http://www.cpan.org/authors/id/$2/$1/$base";
        } elsif ( $base !~ m!^(\w+)://! ) {
            die "Cannot recognize path '$base'; please specify an URL or CPAN ID";
        }
        $base     .= '/' unless $base     =~ m!/\Z!;
        $ftp_base .= '/' unless $ftp_base =~ m!/\Z!;
    }

    require Config;
    my $suffix = "$Config::Config{archname}-$Config::Config{version}.par";

    unless ( $file ||= $self->{file} ) {
        my $name    = $self->name    or return;
        my $version = $self->version or return;
        $name =~ s!::!-!g;
        $self->{file} = $file = "$name-$version-$suffix";
    }

    my $perl = $^X;
    $perl = Win32::GetShortPathName($perl)
        if $perl =~ / / and defined &Win32::GetShortPathName;

    $self->preamble(<<"END_MAKEFILE") if $base;
# --- $class section:

all ::
\t\$(NOECHO) $perl "-M$inc_class" -e "extract_par(q($file))"

END_MAKEFILE

    $self->postamble(<<"END_MAKEFILE");
# --- $class section:

$file: all test
\t\$(NOECHO) \$(PERL) "-M$inc_class" -e "make_par(q($file))"

par :: $file
\t\$(NOECHO) \$(NOOP)

par-upload :: $file
\tcpan-upload -verbose $file

END_MAKEFILE

    $self->{url}     = $base;
    $self->{ftp_url} = $ftp_base;
    $self->{suffix}  = $suffix;

    return $self;
}

=head2 fetch_par

Fetches the .par file previously referenced in the documentation
of the C<par_base> directive.

C<fetch_par> can be used without arguments given the C<par_base>
directive was used before. It will return the name of the file it
fetched.

If the first argument is an URL or a CPAN user ID, the file is
fetched from that directory unless an URL has been previously set.
(Read that again.)

If the second argument is a file name
it is used as the name of the file to download.

If the file could not be fetched, a suitable error message
about no package being available, yada yada yada, is printed.
You can turn this off by specifying a true third argument.

  # Try to fetch the package (see par_base) but
  # don't be verbose about failures
  my $file = fetch_par('', '', undef);

=cut

sub fetch_par {
    my ($self, $url, $file, $quiet) = @_;
    $url = '' if not defined $url;
    $file = '' if not defined $file;
    
    $url = $self->{url} || $self->par_base($url)->{url};
    my $ftp_url = $self->{ftp_url};
    $file ||= $self->{file};

    return $file if -f $file or $self->get_file(
        url     => "$url$file",
        ftp_url => "$ftp_url$file"
    );

    require Config;
    print <<"END_MESSAGE" if $self->{mailto} and ! $quiet;
*** No installation package available for your architecture.
However, you may wish to generate one with '$Config::Config{make} par' and send
it to <$self->{mailto}>, so other people on the same platform
can benefit from it.
*** Proceeding with normal installation...
END_MESSAGE
    return;
}

=head2 extract_par

Takes the name of a PAR::Dist archive file as first argument. The 'blib/'
directory of this archive is extracted and the 'pm_to_blib' is created.

Typical shorthand usage:

  extract_par( fetch_par ) or die "Could not install PAR::Dist archive.";

=cut

sub extract_par {
    my ($self, $file) = @_;
    return unless -f $file;

    if ( eval { require Archive::Zip; 1 } ) {
        my $zip = Archive::Zip->new;
        return unless $zip->read($file) == Archive::Zip::AZ_OK()
                  and $zip->extractTree('', 'blib/') == Archive::Zip::AZ_OK();
    } elsif ( $self->can_run('unzip') ) {
        return if system( unzip => $file, qw(-d blib) );
    }
    else {
        die <<'HERE';
Could not extract .par archive because neither Archive::Zip nor a
working 'unzip' binary are available. Please consider installing
Archive::Zip.
HERE
    }

    local *PM_TO_BLIB;
    open PM_TO_BLIB, '> pm_to_blib' or die $!;
    close PM_TO_BLIB or die $!;

    return 1;
}

=head2 make_par

This directive requires PAR::Dist (version 0.03 or up) on your system.
(And checks that it is available before continuing.)

Creates a PAR::Dist archive from the 'blib/' subdirectory.

First argument must be the name of the PAR::Dist archive to create.

If your Makefile.PL has a C<par_base> directive, the C<make par>
make target will be available. It uses this C<make_par> directive
internally, so on your development system, you can do this to create
a .par binary archive for your platform:

  perl Makefile.PL
  make
  make par

=cut

sub make_par {
    my ($self, $file) = @_;
    unlink $file if -f $file;

    unless ( eval { require PAR::Dist; PAR::Dist->VERSION(0.03) } ) {
        warn "Please install PAR::Dist 0.03 or above first.";
        return;
    }

    return PAR::Dist::blib_to_par( dist => $file );
}

1;

=head1 AUTHOR

Audrey Tang <cpan@audreyt.org>

With documentation from Steffen Mueller <smueller@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2006. Audrey Tang.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut
