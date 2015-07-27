=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Utils::File - Utility functions for error handling

=head1 SYNOPSIS

  use APPRIS::Utils::File
    qw(
       getLocalTime
       printStringIntoTmpFile
       printStringIntoFile
       updateStringIntoFile
       updateStringIntoLockFile
       getStringFromFile
       getTotalStringFromFile
       add_last_character_directory
       prepare_workspace
       rm_dir
       rm_dir_inside
       parse_file
     );

  or to get all methods just

  use APPRIS::Utils::File;

  eval { updateStringIntoFile("text to file",file_path) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to handle files.

=head1 METHODS

=cut

package APPRIS::Utils::File;

use strict;
use warnings;
use Time::localtime;
use File::Basename;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	getLocalTime
	printStringIntoTmpFile
	printStringIntoFile
	updateStringIntoFile
	updateStringIntoLockFile
	getStringFromFile
	getTotalStringFromFile
	add_last_character_directory
	prepare_workspace
	open_dir
	rm_dir
	rm_dir_inside
	parse_file
);

=head2 getLocalTime

  Example    : use APPRIS::Utils::File qw(getLocalTime);
               getLocalTime();
  Description: Get local time as string format 
               "%04d%02d%02d%02d%02d%02d".
  Returntype : string $date_time
  Exceptions : none

=cut

sub getLocalTime()
{
	my($systime)=localtime();
	my($date_string)=sprintf("%04d%02d%02d%02d%02d%02d", $systime->year()+1900 ,$systime->mon()+1, $systime->mday(), $systime->hour(), $systime->min(), $systime->sec());
	return $date_string;
}

=head2 printStringIntoTmpFile

  Arg [1]    : string $msg
               text file
  Example    : use APPRIS::Utils::File qw(printStringIntoTmpFile);
               printStringIntoTmpFile("text to file");
  Description: Print string into temporal file within /tmp/ dir.
  Returntype : string $tmp_file or undef
  Exceptions : none

=cut

sub printStringIntoTmpFile($)
{
	my ($string) = @_;
	my($file)='/tmp/encode_'.getppid.'.faa';
	local(*HANDLE);
	open(HANDLE,">$file") or return undef;
	print HANDLE $string;
	close(HANDLE);
	return $file;
}

=head2 printStringIntoFile

  Arg [1]    : string $msg
               text file
  Arg [2]    : string $file
               file path
  Example    : use APPRIS::Utils::File qw(printStringIntoFile);
               printStringIntoFile("text to file",$file);
  Description: Print string into given file.
  Returntype : string $file or undef
  Exceptions : none

=cut

sub printStringIntoFile($$)
{
	my ($string,$file) = @_;
	local(*HANDLE);
	open(HANDLE,">$file") or return undef;
	print HANDLE $string;
	close(HANDLE);
	return $file;
}

=head2 updateStringIntoFile

  Arg [1]    : string $msg
               text file
  Arg [2]    : string $file
               file path
  Example    : use APPRIS::Utils::File qw(updateStringIntoFile);
               updateStringIntoFile("text to file",$file);
  Description: Modify file with given text.
  Returntype : string $file or undef
  Exceptions : none

=cut

sub updateStringIntoFile($$)
{
	my ($string,$file) = @_;
	local(*HANDLE);
	open(HANDLE,">>$file") or return undef;
	print HANDLE $string;
	close(HANDLE);
	return $file;
}

=head2 updateStringIntoLockFile

  Arg [1]    : string $msg
               string to insert into file
  Arg [2]    : string $file
               file path
  Example    : use APPRIS::Utils::File qw(updateStringIntoLockFile);
               updateStringIntoLockFile("text to file",$file);
  Description: Modify file with given text but locking the file.
  Returntype : string $file or undef
  Exceptions : none

=cut

sub updateStringIntoLockFile($$)
{
	my ($string,$file) = @_;
	local(*HANDLE);
	open(HANDLE,">>$file") or return undef;
	flock(HANDLE, 2);
	print HANDLE $string;		
	close(HANDLE);
	return $file;
}

=head2 getStringFromFile

  Arg [1]    : string $file
               file path
  Example    : use APPRIS::Utils::File qw(getStringFromFile);
               getStringFromFile($file);
  Description: Get string from file.
  Returntype : string $msg or undef
  Exceptions : none

=cut

#sub getStringFromFile($)
#{
#	my ($file) = @_;
#
#	$/=undef;
#	local(*FILE);
#	open(FILE,$file) or return undef;
#	my($string)=<FILE>;
#	close(FILE);
#	$/='\n';
#	return $string;
#}
sub getStringFromFile($)
{
	my ($file) = @_;

	local(*FILE);
	open(FILE,$file) or return undef;
	my(@array)=<FILE>;
	close(FILE);
	my($string)= join "", @array;

	return $string;
}

=head2 getTotalStringFromFile

  Arg [1]    : string $file
               file path
  Example    : use APPRIS::Utils::File qw(getTotalStringFromFile);
               getTotalStringFromFile($file);
  Description: Get string from file as arrayref of string.
  Returntype : Listref of string or undef
  Exceptions : none

=cut

sub getTotalStringFromFile($)
{
	my ($file) = @_;

	local(*FILE);
	open(FILE,$file) or return undef;
	my(@string)=<FILE>;
	close(FILE);

	return \@string;
}

=head2 add_last_character_directory (DEPRECATED)

  Arg [1]    : string $dir
               dir path
  Example    : use APPRIS::Utils::File qw(add_last_character_directory);
               add_last_character_directory($dir);
  Description: Add '/' at the end of dir if it has not it.
  Returntype : string $dir or undef
  Exceptions : none

=cut

sub add_last_character_directory($)
{
	my ($inputPath) = @_;
	my ($outputPath) = $inputPath;

	if (substr($outputPath, -1, 1) ne '/') { $outputPath = $inputPath.'/'; }
	return $outputPath;
}

=head2 prepare_workspace

  Arg [1]    : string $dir
               dir path
  Example    : use APPRIS::Utils::File qw(prepare_workspace);
               prepare_workspace($dir);
  Description: Create directories recursively.
  Returntype : string $dir or undef
  Exceptions : none

=cut

sub prepare_workspace($)
{
	my ($directory) = @_;

	my ($dir,$accum)=('','');
	
	foreach $dir (split(/\//, $directory))
	{
		$accum = "$accum$dir/";
		if($dir ne "")
		{
			if(! -d "$accum")
			{
				mkdir($accum) || return undef;
			}
		}
	}
	$directory = $directory.'/' if (substr($directory, -1, 1) ne '/');
	return $directory;
}

=head2 open_dir

  Arg [1]    : string $dir
               dir path
  Arg [2]    : (Optional) string $pat
               grep pattern NOT WORKS!!!
  Example    : use APPRIS::Utils::File qw(open_dir);
               open_dir($dir);
  Description: Open dir and retrieves the files.
  Returntype : string $files or undef
  Exceptions : none

=cut

sub open_dir($;$)
{
	my ($directory, $pattern) = @_;
	my ($files);

	if ( -e $directory and (opendir(INPUT_DIR, $directory)) )
	{
		if ( defined $pattern )
		{
			@{$files} = sort { $a cmp $b } grep { /$pattern/ } readdir(INPUT_DIR);
		}
		else
		{
			@{$files} = sort { $a cmp $b } readdir(INPUT_DIR);			
		}
	}
	
	return $files;
}
	
=head2 rm_dir

  Arg [1]    : string $dir
               dir path
  Example    : use APPRIS::Utils::File qw(rm_dir);
               rm_dir($dir);
  Description: Remove all directories.
  Returntype : string $dir or undef
  Exceptions : none

=cut

sub rm_dir($)
{
	my ($directory) = @_;

	if ( -e $directory )
	{
		system("rm -rf $directory") == 0 or return undef;
	}
	return $directory;
}

=head2 rm_dir_inside

  Arg [1]    : string $dir
               dir path
  Example    : use APPRIS::Utils::File qw(rm_dir_inside);
               rm_dir_inside($dir);
  Description: Remove all directories.
  Returntype : string $dir or undef
  Exceptions : none

=cut

sub rm_dir_inside($)
{
	my ($directory) = @_;

	if ( -e $directory )
	{
		system("rm -rf $directory/*") == 0 or return undef;
	}
	return $directory;
}

=head2 parse_file

  Arg [1]    : string $file
               file path
  Example    : use APPRIS::Utils::File qw(parse_file);
               parse_file($file);
  Description: Parse file paths into directory, filename and suffix.
  Returntype : array ($dirname,$basename)
  Exceptions : none

=cut

sub parse_file($;$)
{
	my($path,$suffix)=@_;
	my($dirname,$basename)=(undef,undef);
	
	if ( defined $suffix ) { $basename=basename($path,$suffix); }
	else { $basename=basename($path); }
	$dirname=dirname($path);
	
	return($dirname,$basename);	
}


1;