=head1 CONTACT

For contact details see the L<APPRIS website|http://appris-tools.org>.

=cut

=head1 NAME

APPRIS::Utils::TRIFID

=head1 DESCRIPTION

A collection of utility functions related to TRIFID.

=cut

package APPRIS::Utils::TRIFID;

use strict;
use warnings;
use Exporter;

use APPRIS::Utils::Exception qw(throw);

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	get_trifid_pred_file_path
	get_trifid_pred_file_url
	parse_trifid_release
);


=head2 get_trifid_pred_file_path

  Arg[1]      : String $base_dir
                Base directory under which TRIFID
                prediction files are stored
  Arg[2]      : String $species_id
                Species ID corresponding to
                the specified TRIFID release
  Arg[3]      : String $trifid_release
                TRIFID release ID
  Example     : get_trifid_pred_file_path('appris/db/trifid',
                                         'homo_sapiens',
                                         'trifid_GRCh38_g33_20200710');
  Description : Get the path of a TRIFID prediction file,
                given its base URL, species ID and release ID.
  Return type : String containing the TRIFID prediction file path
  Exceptions  : if either base directory or species ID is an empty string,
                or if TRIFID release ID cannot be parsed
  Caller      : general
  Status      : At Risk

=cut

sub get_trifid_pred_file_path($$$) {
	my ($base_dir, $species_id, $trifid_release) = @_;

	if ( $base_dir eq '' ) {
		throw("invalid TRIFID base directory: '$base_dir'");
	}

	if ( $species_id eq '' ) {
		throw("invalid species ID: '$species_id'");
	}

	my ($trifid_meta) = parse_trifid_release($trifid_release);

	$base_dir =~ s/\/$//;
	my $pred_file_path = (
		$base_dir
		.'/'.$species_id
		.'/'.$trifid_meta->{'assembly'}
		.'/'.$trifid_meta->{'dataset'}
		.'/'.$trifid_release.'.tsv.gz'
	);

	return $pred_file_path
}

=head2 get_trifid_pred_file_url

  Arg[1]      : String $base_url
                Base URL under which TRIFID
                prediction files are stored
  Arg[2]      : String $species_id
                Species ID corresponding to
                the specified TRIFID release
  Arg[3]      : String $trifid_release
                TRIFID release ID
  Example     : get_trifid_pred_file_url('http://example.com/trifid',
                                         'homo_sapiens',
                                         'trifid_GRCh38_g33_20200710');
  Description : Get the URL of a TRIFID prediction file,
                given its base URL, species ID and release ID.
  Return type : String containing the TRIFID prediction file URL
  Exceptions  : if either base URL or species ID is an empty string,
                or if TRIFID release ID cannot be parsed
  Caller      : general
  Status      : At Risk

=cut

sub get_trifid_pred_file_url($$$) {
	my ($base_url, $species_id, $trifid_release) = @_;

	if ( $base_url eq '' ) {
		throw("invalid TRIFID base URL: '$base_url'");
	}

	if ( $species_id eq '' ) {
		throw("invalid species ID: '$species_id'");
	}

	my ($trifid_meta) = parse_trifid_release($trifid_release);

	$base_url =~ s/\/$//;
	my $pred_file_url = (
		$base_url
		.'/'.$species_id
		.'/'.$trifid_meta->{'assembly'}
		.'/'.$trifid_meta->{'dataset'}
		.'/'.$trifid_release.'.tsv.gz'
	);

	return $pred_file_url
}

=head2 parse_trifid_release

  Arg[1]      : String $trifid_release
                TRIFID release ID
  Example     : $parsed = parse_trifid_release('trifid_GRCh38_g33_20200710');
  Description : Parse the given TRIFID release ID.
  Return type : Hashref with info taken from the parsed TRIFID release ID
  Exceptions  : if TRIFID release ID cannot be parsed
  Caller      : general
  Status      : At Risk

=cut

sub parse_trifid_release($) {
	my ($trifid_release) = @_;

	my ($result);
	if ( $trifid_release =~ /^trifid_(?<assembly>.+)_(?<dataset>.+)_(?<version>.+)$/ ) {

		$result = {
			'assembly' => $+{'assembly'},
			'dataset' => $+{'dataset'},
			'version' => $+{'version'}
		};

	} else {
		throw("failed to parse TRIFID release: '$trifid_release'");
	}

	return $result;
}


1;
