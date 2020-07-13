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
	$CONFIG_DB_FILE
	$CONFIG_SERVER_FILE
	$FORMAT
	$METHODS
	$TYPE_BED
	$ENCODING
	$TRIFID_URL
	
);

$CONFIG_DB_FILE 	= $ENV{APPRIS_WSERVER_CONF_DB_FILE};
$CONFIG_SERVER_FILE = $ENV{APPRIS_WSERVER_CONF_SR_FILE};
$METHODS			= $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};
$FORMAT				= $ENV{APPRIS_WSERVER_OUTPUT_FORMAT};
$TYPE_BED			= $ENV{APPRIS_WSERVER_TYPEBED}; # NOT EXISTS
$ENCODING			= $ENV{APPRIS_WSERVER_OUTPUT_ENCODING};
$TRIFID_URL			= $ENV{APPRIS_WSERVER_TRIFID_URL};

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
	my ($cgi) = new CGI;
	my ($url_path_info) = $cgi->path_info();
	my (@input_parameters) = split('/',$url_path_info);

	print_http_response(400)
		unless ( @input_parameters ); # defined path

	print_http_response(400)
		if ( scalar(@input_parameters) > 4 );

	my ($type) = $input_parameters[1];
	my ($species) = $input_parameters[2];
	my ($inputs) = $input_parameters[3];

	print_http_response(400)
		unless ( defined $type and defined $species and defined $inputs );

	print_http_response(405, "The parameter $type is not allowed. The parameter must be: 'id', 'name', or 'position'.")
		if ( ($type ne 'id') and ($type ne 'name') and ($type ne 'position') );

	my ($methods) = $cgi->param('methods') || undef;
	if (defined $methods and ($methods ne '')) {
		$methods = lc($methods);		
		if (($methods ne '') and ( 
			($methods =~ /appris/) or
			($methods =~ /firestar/) or
			($methods =~ /matador3d/) or
			($methods =~ /spade/) or
			($methods =~ /corsair/) or
			($methods =~ /inertia/) or
			($methods =~ /crash/) or
			($methods =~ /thump/) or
			($methods =~ /proteo/) or
			($methods =~ /trifid/) )
		) {
			$METHODS = $methods;	
		}
		else {
			print_http_response(405, "The parameter $methods is not allowed. The parameter must be: $METHODS.")
		}
	}
	
	# Optional paramteres
	my ($as) = ( $cgi->param('as') ) ? $cgi->param('as') : undef; # assembly
	my ($sc) = ( $cgi->param('sc') ) ? $cgi->param('sc') : undef; # source name
	my ($ds) = ( $cgi->param('ds') ) ? $cgi->param('ds') : undef; # dataset version

	# Chech if parameters exists
# TODO!! Fix!!!	
#	my ($ens) = $cgi->param('ens') || undef ;
#	if ( defined $ens and ($ens ne '') ) {
#		my ($db) = new DBRetriever(
#								-conf   => $CONFIG_INI_FILE,
#								-specie => $species,
#								-ens	=> $ens
#		);
#		unless ( $db->is_registry($species,$ens) ) {
#			my ($species_ens) = '';
#			eval {
#				my (@info_list) = `grep '_DB' $CONFIG_INI_FILE`;
#				if ( scalar(@info_list) > 0 ) {
#					$species_ens = "The parameter is optional; otherwise it must be: ";
#					foreach my $info (@info_list) {
#						$info =~ s/\n*//g;
#						if ( lc($info) =~ /$species\_ens(\d*)\_db/) {
#							$species_ens .= $1.',';
#						}
#					}
#					$species_ens =~ s/\,$//;
#				}
#			};
#			print_http_response(405, "The parameters specie: $species + ensembl version: $ens is not allowed. $species_ens")
#		}
#	} 
	
	# Get optional parameters
	my ($format) = $cgi->param('format') || undef ;
	if (defined $format and ($format ne '')) {
		$format = lc($format);
		if (($format ne '') and ( ($format eq 'raw') or ($format eq 'tsv') or ($format eq 'gtf') or ($format eq 'json')  or ($format eq 'gff3') or ($format eq 'bed')  or ($format eq 'bed12') ) ) {
			$FORMAT = $format;	
		}
		else {
			print_http_response(405, "The parameter $format is not allowed.")
		}
	}	
	
	my ($ids) = $cgi->param('ids') || undef;
	
	my ($res) = $cgi->param('res') || undef;
	
	# Get features from input
	my ($features);
	my ($retriever) = new DBRetriever(
								-dbconf   	=> $CONFIG_DB_FILE,
								-srvconf   	=> $CONFIG_SERVER_FILE,
								-species 	=> $species,
								-assembly	=> $as,								
								-source		=> $sc,
								-dataset	=> $ds,
								-type   	=> $type,
								-input		=> $inputs,
								-trifid_url => $TRIFID_URL
	);
		
	# Export features (db) from format	
	my ($result) = $retriever->export_features($FORMAT, $METHODS, $ids, $res);
	
	# Everything was OK
	if ( defined $result and ($result ne '') ) {
		if ( $FORMAT eq 'json' ) {
			$ENCODING = 'application/json';
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

exporter

=head1 DESCRIPTION

RESTful web services that retrieves annotations from given gene/transcript or genome position.

=head1 SYNOPSIS

=head2 Required arguments:

=head1 EXAMPLES

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
