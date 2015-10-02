package Module::Install::Admin::ScanDeps;

use strict;
use Module::Install::Base ();
use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.16';
	@ISA = qw(Module::Install::Base);
}

sub scan_dependencies {
    my ($self, $pkg, $perl_version, $pkg_version) = @_;

    return if $pkg eq 'perl';

    $perl_version ||= $self->perl_version or die <<'END_MESSAGE';
Please first specify a required perl version, like this:
    perl_version('5.005');
END_MESSAGE
    $perl_version =~ s{^(\d+)\.(\d+)\.(\d+)}{$1 + $2/1_000 + $3/1_000_000}e;

    require Module::ScanDeps;
    require Module::CoreList;

    die "Module::CoreList has no information on perl $perl_version"
        unless exists $Module::CoreList::version{$perl_version};

    if (my $min_version = Module::CoreList->first_release($pkg, $pkg_version)) {
        return if $min_version <= $perl_version;
    }

    # We only need the first one in the @INC here
    my $file = $self->admin->find_in_inc($pkg)
        or die "Cannot find $pkg in \@INC";
    my %result = ($pkg => $file);

    my @files = ($file);
    while (@files) {
        my $deps = Module::ScanDeps::scan_deps(
            files   => \@files,
            recurse => 0,
        );

        @files = ();

        foreach my $key (keys %$deps) {
            if ($deps->{$key}{type} eq 'shared') {
                foreach my $used_by (@{$deps->{$key}{used_by}}) {
                    $used_by =~ s!/!::!g;
                    $used_by =~ s!\.pm\Z!!i or next;
                    next if exists $result{$used_by};
                    $result{$used_by} = undef;
                    my $min_version = Module::CoreList->first_release($used_by);
                    print "skipped $used_by (needs shared library)\n"
                      unless !$min_version || $min_version <= $perl_version;
                }
            }
        }

        foreach my $key (keys %$deps) {
            my $dep_pkg = $key;
            $dep_pkg =~ s!/!::!g;
            $dep_pkg =~ s!\.pm\Z!!i or next;

            if (my $min_version = Module::CoreList->first_release($dep_pkg)) {
                next if $min_version <= $perl_version;
            }
            next if $dep_pkg =~ /^(?:DB|(?:Auto|Dyna|XS)Loader|threads|warnings)\b/i;
            next if exists $result{$dep_pkg};

            $result{$dep_pkg} = $deps->{$key}{file};
            push @files, $deps->{$key}{file};
        }
    }

    while (my($k,$v) = each %result) {
        delete $result{$k} unless defined $v;
    }
    return \%result;
}

1;
