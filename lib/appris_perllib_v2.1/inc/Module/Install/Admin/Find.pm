package Module::Install::Admin::Find;

use strict;
use File::Find ();
use Module::Install::Base ();
use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.16';
	@ISA = qw(Module::Install::Base);
}

sub find_extensions {
    my $self = shift;
    $self->_top->find_extensions(@_);
}

sub find_in_inc {
    my ($self, $pkg) = @_;

    unless ($pkg =~ /\.pm$/) {
        $pkg =~ s!::!/!g;
        $pkg = "$pkg.pm";
    }

    my @found;
    foreach my $inc (@INC) {
        next if $inc eq $self->_top->{prefix} or ref($inc);
        push @found, "$inc/$pkg" if -f "$inc/$pkg";
    }

    wantarray ? @found : $found[0];
}

sub glob_in_inc {
    my ($self, $pkg) = @_;

    unless ($pkg =~ /\.pm$/) {
        $pkg =~ s!::!/!g;
        $pkg = "$pkg.pm";
    }

    my @found;
    foreach my $inc (@INC) {
        next if $inc eq $self->_top->{prefix} or ref($inc);
        push @found, [ do {
            my $p = $_;
            $p =~ s!^\Q$inc\E/!!;
            $p =~ s!/!::!g;
            $p =~ s!\.pm\Z!!gi;
            $p
        }, $_ ] for grep -e, glob("$inc/$pkg");
    }

    wantarray ? @found : $found[0];
}

sub find_files {
    my ($self, $file, $path) = @_;
    $path = '' if not defined $path;
    $file = "$path/$file" if length($path);
    if (-f $file) {
        return ($file);
    }
    elsif (-d $file) {
        my @files = ();
        local *DIR;
        opendir(DIR, $file) or die "Can't opendir $file";
        while (my $new_file = readdir(DIR)) {
            next if $new_file =~ /^(\.|\.\.)$/;
            push @files, $self->find_files($new_file, $file);
        }
        return @files;
    }
    return ();
}

1;
