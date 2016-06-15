#!/usr/bin/perl

use strict;
use warnings;
use CGI;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use HTTP qw( print_http_response );
use DBRetriever;
use WSRetriever;
$|=1; # not use buffering

###################
# Global variable #
###################
use vars qw(
	$CONFIG_DB_FILE
	$CONFIG_SERVER_FILE
	$FORMAT
	$TYPE_SEQ
	$ENCODING
);

$CONFIG_DB_FILE 	= $ENV{APPRIS_WSERVER_CONF_DB_FILE};
$CONFIG_SERVER_FILE = $ENV{APPRIS_WSERVER_CONF_SR_FILE};
$FORMAT 			= $ENV{APPRIS_WSERVER_OUTPUT_FORMAT};
$ENCODING 			= $ENV{APPRIS_WSERVER_OUTPUT_ENCODING};

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
		
	# Optional paramteres
	my ($sp) = ( $cgi->param('sp') ) ? $cgi->param('sp') : undef; # species
	my ($as) = ( $cgi->param('as') ) ? $cgi->param('as') : undef; # assembly	
	my ($sc) = ( $cgi->param('sc') ) ? $cgi->param('sc') : undef; # source name
	my ($ds) = ( $cgi->param('ds') ) ? $cgi->param('ds') : undef; # dataset version
	
	# Chech if parameters exists
# TODO!! Fix!!!
#	if ( (defined $as and ($as ne '')) or (defined $ds and ($ds ne '')) or (defined $db and ($db ne '')) ) {
#		my ($db) = new DBRetriever(
#								-conf   => $CONFIG_INI_FILE,
#								-assembly	=> $as,
#								-dataset	=> $ds,
#								-database	=> $db
#		);
#		unless ( $db->is_registry(undef,$ds,$db) ) {
#			print_http_response(405, "The input parameters are not allowed.")
#		}
#	} 
	
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
									-dbconf   	=> $CONFIG_DB_FILE,
									-srvconf   	=> $CONFIG_SERVER_FILE,
									-species	=> $sp,
									-assembly	=> $as,
									-source		=> $sc,
									-dataset	=> $ds,
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
