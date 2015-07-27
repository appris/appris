#!/usr/bin/perl

use strict;
use warnings;
use CGI;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use HTTP qw( print_http_response );
use DBRetriever;
$|=1; # not use buffering

###################
# Global variable #
###################
use vars qw(
	$CONFIG_INI_FILE
	$FORMAT
	$TYPE_SEQ
	$ENCODING
);

$CONFIG_INI_FILE = $ENV{APPRIS_WSERVER_SCRIPTS_DB_INI};
$FORMAT = $ENV{APPRIS_WSERVER_OUTPUT_FORMAT};
$ENCODING = $ENV{APPRIS_WSERVER_OUTPUT_ENCODING};

#####################
# Method prototypes #
#####################
sub main();

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	# Get input parameters:
	# id/homo_sapiens/ENSG00000099904,ENST00000436518/firestar?format=raw
	# name/homo_sapiens/RNF215/firestar?format=raw
	# position/homo_sapiens/chr22:30773835-30821305/spade?format=tsv
	my ($cgi) = new CGI;	
	my ($url_path_info) = $cgi->path_info();
	my (@input_parameters) = split('/',$url_path_info);

	print_http_response(400)
		unless ( @input_parameters ); # defined path

	print_http_response(400)
		if ( scalar(@input_parameters) > 2 ); # only 1 parameters (+1 empty value)

	my ($inputs) = $input_parameters[1];

	print_http_response(400)
		unless ( defined $inputs );

	# Ensembl database version is checked
	my ($db) = $cgi->param('db') || undef ;
	
	my ($ens) = $cgi->param('ens') || undef ;
	if ( defined $ens and ($ens ne '') ) {
		my ($db) = new DBRetriever(
								-conf   => $CONFIG_INI_FILE,
								-specie => undef,
								-ens	=> $ens
		);
		unless ( $db->is_registry(undef,$ens) ) {
			print_http_response(405, "The parameter $ens is not allowed.")
		}
	} 
	
	# Get optional parameters
	my ($format) = $cgi->param('format') || undef ;
	if (defined $format and ($format ne '')) {
		$format = lc($format);
		if (($format ne '') and ( ($format eq 'json') or ($format eq 'xml') ) ) {			
			$FORMAT = $format;	
		}
		else {
			print_http_response(405, "The parameter $format is not allowed.")
		}
	}	
	
	# Get features from input
	my ($features);
	my ($retriever) = new WSRetriever(
									-jobid	=> $inputs
	);
	
	# if jobid does not exists means that is database input
	unless ( defined $retriever and $retriever->transl ) {
		$retriever = new DBRetriever(
									-conf   	=> $CONFIG_INI_FILE,
									-ens		=> $ens,
									-assembly	=> $db,
									-input		=> $inputs
		);		
	}
		
	# Export features (db) from format	
	my ($result) = $retriever->seek_features($FORMAT);
	
	# Everything was OK
	if ( defined $result and ($result ne '') ) {
		if ( $FORMAT eq 'json' ) {
			$ENCODING = 'application/json';
		}
		elsif ( $FORMAT eq 'xml' ) {
			$ENCODING = 'text/xml';
		}
		print_http_response(200, $result, $ENCODING);
	}
	else {
		print_http_response(404);
	}
}

main();


1;


__END__

=head1 NAME

seeker

=head1 DESCRIPTION

RESTful web services that find genes/transcripts from identifer, name or genome position.

=head1 SYNOPSIS

=head2 Required arguments:

=head1 EXAMPLES

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
