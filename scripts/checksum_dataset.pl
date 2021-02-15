#!/usr/bin/env perl

use 5.012;
use strict;
use warnings::register;
use Digest::MD5;
use File::Copy qw(move);
use File::Spec::Functions qw(catfile);
use File::Temp qw(tempdir);
use Getopt::Long;

use APPRIS::Utils::Logger;

####################
# Global variables #
####################

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($data_dir) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'data-dir|d=s' 			=> \$data_dir,
	'loglevel=s'			=> \$loglevel,
	'logfile=s'				=> \$logfile,
	'logpath=s'				=> \$logpath,
	'logappend'				=> \$logappend,
);

# Check required arguments
unless ( defined $data_dir ) {
	print `perldoc $0`;
	exit 1;
}

# Optional arguments

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);


# Main subroutine
sub main()
{
	my ($cksum_file_name) = 'md5checksums.txt';
	my ($cksum_file_path) = catfile($data_dir, $cksum_file_name);
	
	$logger->info("-- listing files in dataset directory: $data_dir\n");
	my (%file_to_path);
	opendir(my $dh, $data_dir) or $logger->error("failed to access dataset directory: '$data_dir'");
	while (readdir $dh) {
		my ($file_name) = $_;
		my ($file_path) = catfile($data_dir, $file_name);

		next unless( -f $file_path );
		next if ( $file_name eq $cksum_file_name );

		$file_to_path{$file_name} = $file_path;
	}
	closedir($dh) or $logger->error("failed to close dataset directory: '$data_dir'");

	$logger->info("-- calculating checksums for dataset files\n");
	my ($tmp_dir) = tempdir( CLEANUP => 1 );
	my ($tmp_file_path) = catfile($tmp_dir, $cksum_file_name);
	open(my $out_fh, '>', $tmp_file_path) or $logger->error("failed to open file: '$tmp_file_path'");
	foreach my $file_name ( sort(keys %file_to_path) ) {
		my ($file_path) = $file_to_path{$file_name};

		my ($ctx) = Digest::MD5->new;
		open(my $in_fh, $file_path) or $logger->error("failed to open file: '$file_path'");
		$ctx->addfile($in_fh);
		close($in_fh) or $logger->error("failed to close file: '$file_path'");

		print $out_fh $ctx->hexdigest."  $file_name\n";
	}
	close($out_fh) or $logger->error("failed to close file: '$tmp_file_path'");
	
	$logger->info("-- saving checksums to file: $cksum_file_path\n");
	move($tmp_file_path, $cksum_file_path);

	$logger->finish_log();
	
	exit 0;	

}

main();


1;

__END__

=head1 NAME

checksum_dataset

=head1 DESCRIPTION

Generate a checksum file for the specified APPRIS dataset

=head1 SYNOPSIS

checksum_dataset

=head2 Required arguments:

	--data-dir <Path of an APPRIS dataset directory>

=head2 Optional arguments (log arguments):

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>

	--logpath=PATH <Write logfile to PATH (default: .)>

	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

	perl checksum_dataset.pl --data-dir=data/2021_02.v39/homo_sapiens/e102v39

=head1 AUTHOR

Adapted by Thomas Walsh (CNIO) from APPRIS script code that
was originally written by Jose Manuel Rodriguez Carrasco.

For contact details see the L<APPRIS website|http://appris-tools.org>.

=cut

