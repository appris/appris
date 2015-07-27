use inc::Module::Install;


# Define metadata
name			'APPRIS';
author			'Jose Manuel Rodriguez <jmrodriguez@cnio.es>';
version			'0.1';
abstract 		'This distribution is for communicating with a APPRIS annotations';
all_from		'APPRIS/Registry.pm';
license			'perl';
perl_version	'5.16.2';

# APPRIS::* depends
requires 'DB'						=> '1.04';
requires 'DBD::mysql'				=> 0;
requires 'POSIX'					=> 0;
requires 'FindBin'					=> 0;
requires 'Config::IniFiles'			=> 0;
requires 'Getopt::Long'				=> 0;
requires 'Data::Dumper'				=> 0;
requires 'Exporter'					=> 0;
requires 'Bio::Seq'					=> 0;
requires 'Bio::SeqIO'				=> 0;
requires 'Bio::SearchIO'			=> 0;
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
requires 'IPC::Run'					=> 0;
requires 'XML::DOM'					=> 0;
requires 'CGI'						=> 0;
requires 'HTTP::Status'				=> 0;
requires 'Email::Valid'				=> 0;
requires 'MIME::Lite'				=> 0;
requires 'Digest::MD5'				=> 0;

#auto_install();
#install_script 'script/scriptname';

WriteAll;


