package Module::Install::Inline;

use strict;
use Module::Install::Base ();

use vars qw{$VERSION @ISA $ISCORE};
BEGIN {
	$VERSION = '1.16';
	@ISA     = 'Module::Install::Base';
	$ISCORE  = 1;
}

sub Inline { $_[0] }

sub write {
    my $self = shift;
    my $name = $self->module_name || $self->name
        or die "Please set name() before calling &Inline->write\n";
    $name =~ s/-/::/g;
    my $object = (split(/::/, $name))[-1] or return;
    my $version = $self->version
        or die "Please set version() or version_from() before calling &Inline->write\n";

    $version =~ /^\d\.\d\d$/ or die <<"END_MESSAGE";
Invalid version '$version' for $name.
Must be of the form '#.##'. (For instance '1.23')
END_MESSAGE

    $self->clean_files('_Inline', "$object.inl");
    $self->build_requires('Inline' => 0.44); # XXX: check for existing? yagni?

    my $class = ref($self);
    my $prefix = $self->_top->{prefix};
    $self->postamble(<<"MAKEFILE");
# --- $class section:

.SUFFIXES: .pm .inl

.pm.inl:
\t\$(PERL) -I$prefix "-Mblib" "-MInline=NOISY,_INSTALL_" "-M$name" -e1 $version \$(INST_ARCHLIB)

pure_all :: $object.inl

MAKEFILE

    $self->Makefile->write;
}

1;
