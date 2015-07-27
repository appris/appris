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
	$METHODS
	$TYPE_BED
	$ENCODING
);

$CONFIG_INI_FILE	= $ENV{APPRIS_WSERVER_SCRIPTS_DB_INI};
$METHODS			= $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};
$FORMAT				= $ENV{APPRIS_WSERVER_OUTPUT_FORMAT};
$TYPE_BED			= $ENV{APPRIS_WSERVER_TYPEBED}; # NOT EXISTS
$ENCODING			= $ENV{APPRIS_WSERVER_OUTPUT_ENCODING};

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
	my ($specie) = $input_parameters[2];
	my ($inputs) = $input_parameters[3];

	print_http_response(400)
		unless ( defined $type and defined $specie and defined $inputs );

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
			($methods =~ /proteo/) )
		) {
			$METHODS = $methods;	
		}
		else {
			print_http_response(405, "The parameter $methods is not allowed. The parameter must be: $METHODS.")
		}
	}
	
	# Ensembl database version is checked
	my ($db) = $cgi->param('db') || undef ;
	
	my ($ens) = $cgi->param('ens') || undef ;
	if ( defined $ens and ($ens ne '') ) {
		my ($db) = new DBRetriever(
								-conf   => $CONFIG_INI_FILE,
								-specie => $specie,
								-ens	=> $ens
		);
		unless ( $db->is_registry($specie,$ens) ) {
			my ($specie_ens) = '';
			eval {
				my (@info_list) = `grep '_DB' $CONFIG_INI_FILE`;
				if ( scalar(@info_list) > 0 ) {
					$specie_ens = "The parameter is optional; otherwise it must be: ";
					foreach my $info (@info_list) {
						$info =~ s/\n*//g;
						if ( lc($info) =~ /$specie\_ens(\d*)\_db/) {
							$specie_ens .= $1.',';
						}
					}
					$specie_ens =~ s/\,$//;
				}
			};
			print_http_response(405, "The parameters specie: $specie + ensembl version: $ens is not allowed. $specie_ens")
		}
	} 
	
	# Get optional parameters
	my ($format) = $cgi->param('format') || undef ;
	if (defined $format and ($format ne '')) {
		$format = lc($format);
		if (($format ne '') and ( ($format eq 'raw') or ($format eq 'tsv') or ($format eq 'gtf') or ($format eq 'json')  or ($format eq 'gff3') or ($format eq 'bed') ) ) {
			$FORMAT = $format;	
		}
		else {
			print_http_response(405, "The parameter $format is not allowed.")
		}
	}	
	
	my ($ids) = $cgi->param('ids') || undef;
	
	my ($res) = $cgi->param('res') || undef;
	
	my ($typebed) = $cgi->param('typebed') || undef;
	if (defined $typebed and ($typebed ne '')) {
		$typebed =~ s/\s*//g;
		$TYPE_BED = $typebed;
	}

	# Get features from input
	my ($features);
	my ($retriever) = new DBRetriever(
								-conf   	=> $CONFIG_INI_FILE,
								-specie 	=> $specie,
								-ens		=> $ens,
								-assembly	=> $db,
								-type   	=> $type,
								-input		=> $inputs
	);
		
	# Export features (db) from format	
	my ($result) = $retriever->export_features($FORMAT, $METHODS, $ids, $res, $TYPE_BED);
	
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
