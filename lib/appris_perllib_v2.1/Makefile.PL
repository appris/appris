use inc::Module::Install;


# Define metadata
name			'APPRIS';
author			'Jose Manuel Rodriguez <jmrodriguez@cnio.es>';
version			'0.1';
abstract 		'This distribution is for communicating with a APPRIS annotations';
#license			'perl';
#perl_version	'5.16.2';

# Perl dependencies

# Module::Install
requires 'YAML::Tiny'				=> '1.69';
requires 'File::Remove'				=> '1.52';

# APPRIS
requires 'DB'						=> 0;
requires 'DBD::mysql'				=> 0;
requires 'POSIX'					=> 0;
requires 'FindBin'					=> 0;
requires 'Config::IniFiles'			=> 0;
requires 'Getopt::Long'				=> 0;
requires 'Data::Dumper'				=> 0;
requires 'Exporter'					=> 0;
requires 'Bio::Seq'					=> '1.06';
requires 'Bio::SeqIO'				=> '1.06';
requires 'Bio::SearchIO'			=> '1.06';
requires 'File::Temp'				=> 0;
requires 'File::Basename'			=> 0;
requires 'Time::localtime'			=> 0;
requires 'File::Basename'			=> 0;
requires 'JSON'						=> 0;
requires 'SVG'						=> 0;
requires 'Statistics::Descriptive'	=> 0;
requires 'Moose'					=> 0;
requires 'Class::Load'				=> 0;
requires 'Data::OptList'			=> 0;
requires 'Module::Implementation'	=> 0;
requires 'Class::Load::XS'			=> 0;
requires 'MRO::Compat'				=> 0;
requires 'Data::Printer'			=> 0;
requires 'XML::DOM'					=> 0;
requires 'CGI'						=> 0;
requires 'HTTP::Status'				=> 0;
requires 'Email::Valid'				=> 0;
requires 'MIME::Lite'				=> 0;
requires 'Digest::MD5'				=> 0;

requires 'List::MoreUtils'			=> 0;



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


