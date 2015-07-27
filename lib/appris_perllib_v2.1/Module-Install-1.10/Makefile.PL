# Normally we don't need strict, but for this one we'll do it anyway.
use strict;

# Warnings need to be enabled this way, because Makefile.PL is not
# an executable script (so we can't use a hashbang).
BEGIN {
	$|  = 1;
	$^W = 1;
}

# NOTE: Module::Install's Makefile.PL needs this "use lib 'lib'" to
# bootstrap itself. Other module distributions using Module::Install
# do not need it. It causes inc::Module::Install to load from the
# (unique to this distribution) file ./lib/inc/Module/Install.pm instead
# of the system inc::Module::Install or the bundled ./inc/Module/Install.pm
use lib 'lib';
use inc::Module::Install;

name          'Module-Install';
author        'Adam Kennedy <adamk@cpan.org>';
perl_version  '5.006';
all_from      'lib/Module/Install.pm';

requires      'Devel::PPPort'       => '3.16';
requires      'ExtUtils::Install'   => '1.52';
requires      'ExtUtils::MakeMaker' => '6.59';
requires      'ExtUtils::ParseXS'   => '2.19';
requires      'File::Path'          => 0; #needs version
requires      'File::Remove'        => '1.42';
requires      'File::Spec'          => '3.28';
requires      'Module::Build'       => '0.29';
requires      'Module::CoreList'    => '2.17';
requires      'Module::ScanDeps'    => '1.09'; #detects prereqs better
requires      'Parse::CPAN::Meta'   => '1.4413';
requires      'Win32::UTCFileTime'  => '1.56' if win32;
requires      'YAML::Tiny'          => '1.38';

test_requires 'Test::Harness'       => '3.13';
test_requires 'Test::More'          => '0.86';

recommends    'Archive::Zip'            => '1.37';
recommends    'File::HomeDir'           => '1';
recommends    'JSON'                    => '2.9';
recommends    'LWP::Simple'             => '6.00';
recommends    'LWP::UserAgent'          => '6.05';
recommends    'Module::Find'            => '0.12',
recommends    'Module::Runtime'         => '0.014',
recommends    'Software::License'       => '0.103009',
recommends    'Test::Requires'          => '0.07',
recommends    'Test::Software::License' => '0.002',
recommends    'PAR::Dist'               => '0.29';
 
# Remove some extra test files
clean_files( qw{ t/Foo } );

# Test configure_requires support (also needed by WriteMyMeta)
configure_requires 'YAML::Tiny' => '1.33';

# Makes the recursive test script run
# (to test recursive test support)
tests_recursive;

# This is a normal CPAN module
install_as_cpan;

# Normally this shouldn't be done in here, it should be done
# by release automation. This is mostly for testing reasons.
repository 'https://github.com/Perl-Toolchain-Gang/Module-Install.git';

# All of our dependencies are the same on all platforms
dynamic_config 0;

WriteAll;
