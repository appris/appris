package Module::Install::Admin::Makefile;

use strict 'vars';
use Module::Install::Base;
use ExtUtils::MakeMaker ();

use vars qw{$VERSION @ISA};
BEGIN {
	$VERSION = '1.16';
	@ISA     = qw{Module::Install::Base};
}

sub postamble {
    my ($self, $text) = @_;
    my $class       = ref($self);
    my $top_class   = ref($self->_top);
    my $admin_class = join('::', @{$self->_top}{qw(name dispatch)});

    $self->{postamble} ||= << "END_MAKEFILE";
# --- $class section:

realclean purge ::
\t\$(RM_F) \$(DISTVNAME).tar\$(SUFFIX)
\t\$(RM_F) MANIFEST.bak _build
\t\$(PERL) "-Ilib" "-M$admin_class" -e "remove_meta()"
\t\$(RM_RF) inc

reset :: purge

upload :: test dist
\tcpan-upload -verbose \$(DISTVNAME).tar\$(SUFFIX)

grok ::
\tperldoc $top_class

distsign ::
\tcpansign -s

END_MAKEFILE

    $self->{postamble} .= $text if defined $text;
    return $self->{postamble};
}

sub preop {
    my $self = shift;
    my ($user_preop) = @_;
    my $admin_class = join('::', @{$self->_top}{qw(name dispatch)});
    $user_preop = qq{\$(PERL) -I. "-M$admin_class" -e "dist_preop(q(\$(DISTVNAME)))"} unless $user_preop;
    return { PREOP => $user_preop };
}

1;
