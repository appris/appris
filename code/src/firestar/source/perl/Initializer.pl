#!/usr/bin/perl
# $Id: Initializer.pl
# Developed by: Paolo Maietta -pmaietta@cnio.es-
# Created: 24th of August 2011
# ________________________________________________________________
# This program automatically generates folders - subfolders and  
# files for a new check-out working directory of firestar-FireDB
#

use strict;
use FindBin;
use File::Copy;
use File::Path qw(make_path remove_tree);
use DBI;
use Cwd;
my $cwd=$FindBin::Bin;
use Config::IniFiles;

chdir("$cwd/..");
my $fire_local = getcwd;
my $path;

$path="$fire_local/Square";
unless (-d "$path/tmpout"){
	make_path "$path/tmpout", {verbose => 1};
}
unless (-d "$path/tmpchads"){
	make_path "$path/tmpchads", {verbose => 1};
}

$path="$fire_local/tmp/faatmp";
unless (-d "$fire_local/tmp/faatmp"){
	make_path "$path", {verbose => 1};
}

unless (-e "$path/FAA_LOG.txt"){
	open (LOG5,">$path/FAA_LOG.txt");
	close(LOG5);
	chmod 0666,"$path/FAA_LOG.txt";
}

print "\n\nPlease, introduce the date of your FireDB version (format ddmmyyyy ej. 01Jan2000): ";
my $release=<STDIN>;
chomp($release);

my $user=$ENV{"USER"};


open (VAR,'>',"$fire_local/CONFIG_fire_var.ini");

print VAR
"
[PATHS]
home=$fire_local
lib=$fire_local/lib
faatmp=$fire_local/tmp/faatmp
chads=$fire_local/database
[MYSQL]
dbHost=localhost
dbName=FireDB
user=$user
pass=
[DATABASES]
release=$release
nrdb=[put here the name of the blast database that will be used for the generation of the profile]
hhbdb=hhblits_
hhprof=nr20_12Aug11
blast_path=[put here the path to the blast databases]
hhdb_path=[put here the path to hhsuite databases]
[PROGRAMS]
square=$fire_local/Square/another_fire_mess_web
square_test=$fire_local/Square/another_fire_mess_test
blast=[put here the directory of blast binaries. Blast+ not supported]
";
close(VAR);

######################################################
######		Connection check-up		######
######################################################


# First of all Initializer gets information from the file FIRE.var generated, then
# checks for all the external resources FireDB-firestar need.
my $variables=Config::IniFiles->new(-file => "$fire_local/CONFIG_fire_var.ini");



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Prueba de conexion con base de datos
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my $dbHost=$variables->val('MYSQL','dbHost');
my $dbName=$variables->val('MYSQL','dbName');
my $user=$variables->val('MYSQL','user');
my $pass=$variables->val('MYSQL','pass');

# FireDB mysql database
print "\nChecking now the connection to FireDB Database ....";
my $ConnectError=0;
my $dataHandle = DBI->connect("DBI:mysql:database=$dbName;host=$dbHost","$user","$pass")
                        || {$ConnectError=1};
if($ConnectError==0)
{
        print "\t[OK]\n";
}
else
{
        print "\t[NO]\n";
}



# checking chads FireDB database
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "\nChecking the modules neeeded using the program module_checker.pl\n";
print "\n\tdirectory $fire_local/perl ....\n\n";
system ("perl $fire_local/perl/module_checker.pl $fire_local/perl");
print "\n\n\tdirectory $fire_local/Square ....\n\n";
system ("perl $fire_local/perl/module_checker.pl $fire_local/Square");
print "\n\n\tdirectory $fire_local/perl/lib ....\n\n";
system ("perl $fire_local/perl/module_checker.pl $fire_local/perl/lib");


__END__



