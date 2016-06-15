#!/usr/bin/perl

use strict;
use warnings;
use CGI;
use File::Temp;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use HTTP qw( print_http_response );
use DBRetriever;
#use WSRetriever;
$|=1; # not use buffering

###################
# Global variable #
###################
use vars qw(
	$CONFIG_DB_FILE
	$CONFIG_SERVER_FILE
	$OPERATION
	$METHODS
	$FORMAT
	$TYPE_SEQ
	$ENCODING
);

$CONFIG_DB_FILE 	= $ENV{APPRIS_WSERVER_CONF_DB_FILE};
$CONFIG_SERVER_FILE = $ENV{APPRIS_WSERVER_CONF_SR_FILE};
$OPERATION			= $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE_VIEW};
$METHODS			= $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};
$FORMAT				= '';
$TYPE_SEQ			= $ENV{APPRIS_WSERVER_OUTPUT_TYPE_SEQ};
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

	my ($jobid) = undef;
	my ($type) = $input_parameters[1];
	my ($species) = $input_parameters[2];
	my ($inputs) = $input_parameters[3];

	if ( defined $type and !(defined $species) and !(defined $inputs) ) {
		$jobid = $type; $type = undef;
	}
	else {
		print_http_response(400)
			unless ( defined $type and defined $species and defined $inputs );

		print_http_response(405, "The parameter $type is not allowed. The parameter must be: 'id', 'name', or 'position'.")
			if ( ($type ne 'id') and ($type ne 'name') and ($type ne 'position') );	
	}

	my ($operation) = $cgi->param('operation') || undef;
	if (defined $operation and ($operation ne '')) {
		$operation = lc($operation);
		if ( $operation =~ /align/ or $operation =~ /sequences/ ) { 
			$OPERATION = $operation;
			$FORMAT = 'html';
		}
		elsif ( $operation =~ /genome/ ) { 
			$OPERATION = $operation;
			$FORMAT = 'html';
		}
		elsif ( $operation =~ /svg/ ) {
			$OPERATION = $operation;
			$FORMAT = 'svg';
		}
		else {
			print_http_response(405, "The parameter $operation is not allowed. The parameter must be: $OPERATION.")
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
#						if ( lc($info) =~ /$species\_ens(\d*)\_db/) {
#							$species_ens .= $1.',';
#						}
#					}
#					$species_ens =~ s/\,$//;
#				}
#			};
#			print_http_response(405, "The parameter $ens is not allowed. $species_ens")
#		}
#	}
	
	# Get optional parameters
	my ($methods) = $cgi->param('methods');# || undef;
	if ( defined $methods ) {
		$methods = lc($methods);		
		if ( ($methods ne '') and ( 
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
	my ($ids) = $cgi->param('ids');# || undef;
	if ( defined $ids and ($ids eq '') ) {
		print_http_response(405, "The parameter 'ids' must be defined by sequence identifier.")		
	}
	
	
	# Get result from input
	my ($result);
	my ($retriever);
	if ( defined $jobid and !(defined $type) and !(defined $inputs) and !(defined $species) ) {
		$retriever = new WSRetriever(
									-jobid	=> $jobid
		);
	}
	elsif ( defined $type and defined $inputs and defined $species and !(defined $jobid) ) {
		$retriever = new DBRetriever(
								-dbconf   	=> $CONFIG_DB_FILE,
								-srvconf   	=> $CONFIG_SERVER_FILE,
								-species 	=> $species,
								-assembly	=> $as,								
								-source		=> $sc,
								-dataset	=> $ds,
								-type   	=> $type,
								-input		=> $inputs
		);
	}
		
	if ( defined $retriever ) {
		# get sequences/align annotations
		if ( ($OPERATION eq 'align') or ($OPERATION eq 'sequences') ) {
			$result = $retriever->export_seq_aln_features($OPERATION, $FORMAT, $METHODS, $ids);
		}
		# get genome annotations
		elsif ( $OPERATION eq 'genome' ) {
			$result = $retriever->export_gen_features($METHODS, $ids);
		}
		# get svg annotations
		elsif ( $OPERATION eq 'svg' ) {
			$result = $retriever->export_seq_features($OPERATION, $FORMAT, undef, $METHODS, $ids);
		}
	}
		
	# Everything was OK
	if ( defined $result and ($result ne '') ) {
		if ( $FORMAT eq 'json' ) {
			$ENCODING = 'application/json';
		}
		elsif ( $FORMAT eq 'html' ) {
			$ENCODING = 'text/html';
		}
		elsif ( $FORMAT eq 'svg' ) {
			$ENCODING = 'image/svg+xml';
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

sequencer

=head1 DESCRIPTION

RESTful web services that retrieves sequences and more from given gene/transcript or genome position.

=head1 SYNOPSIS

=head2 Required arguments:

=head1 EXAMPLES

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
