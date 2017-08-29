use inc::Module::Install;


# Define metadata
name			'APPRIS';
author			'Jose Manuel Rodriguez <jmrodriguez@cnio.es>';
version			'0.1';
abstract 		'This distribution is for communicating with a APPRIS annotations';
license			'perl';
perl_version	'5.18';

# Perl dependencies

# Module::Install
requires 'YAML::Tiny'				=> '1.69';
requires 'File::Remove'				=> '1.52';

# APPRIS
requires 'DB'							=> '1.04';
requires 'DBD::mysql'					=> '4.039';
requires 'POSIX'						=> '1.32';
requires 'FindBin'						=> '1.51';
requires 'Config::IniFiles'				=> '2.82';
requires 'Getopt::Long'					=> '2.39';
requires 'Data::Dumper'					=> '2.145';
requires 'Exporter'						=> '5.68';
requires 'Bio::Seq'						=> 0;
requires 'Bio::SeqIO'					=> 0;
requires 'Bio::SearchIO'				=> 0;
requires 'IPC::Run'						=> 0;
requires 'File::Temp'					=> '0.23';
requires 'File::Basename'				=> '2.84';
requires 'Time::localtime'				=> '1.02';
requires 'JSON'							=> '2.90';
requires 'SVG'							=> '2.59';
requires 'Statistics::Descriptive'		=> '3.0605';
requires 'Moose'						=> '2.1202';
requires 'Class::Load'					=> '0.20';
requires 'Data::OptList'				=> '0.109';
requires 'Module::Implementation'		=> '0.07';
requires 'Class::Load::XS'				=> '0.06';
requires 'MRO::Compat'					=> '0.12';
requires 'Data::Printer'				=> '0.35';
requires 'XML::DOM'						=> '1.44';
requires 'XML::Parser'					=> '2.44';
requires 'CGI'							=> '3.63';
requires 'HTTP::Status'					=> '6.03';
requires 'Email::Valid'					=> '1.202';
requires 'MIME::Lite'					=> '3.030';
requires 'Digest::MD5'					=> '2.52';
requires 'List::MoreUtils'				=> '0.33';

my (%arg_install);
foreach my $arg (@ARGV) {
    if ( $arg =~ /^INSTALL_BASE=(.*)$/ ) {    	
#    	$arg_install{makepl_arg} = "q[INSTALL_BASE=$1]";
#    	$arg_install{mbuild_arg} = "q[--install_base $1]";
#    	$arg_install{mbuild_install_arg} = "q[--install_path lib=$1]";
#    	$arg_install{mbuildpl_arg} = "q[--install_base $1]";
    }
}
    	
#auto_install( %arg_install );
#installdeps_target( %arg_install );
installdeps_target(
#	make_args   => "INSTALL_BASE=$install_base",   # option(s) for CPAN::Config
#	force       => 1,           # pseudo-option to force install
#	do_once     => 1
);

#install_script 'script/scriptname';

WriteAll;


