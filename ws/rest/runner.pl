#!/usr/bin/perl

use strict;
use warnings;
use CGI;
use JSON;
use Email::Valid;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use HTTP qw( print_http_response );
use WSRunner;
#$|=1; # not use buffering

###################
# Global variable #
###################
use vars qw(
	$FORMAT
	$METHODS
	$HEAD_BED
	$ENCODING
);

$FORMAT = $ENV{APPRIS_WSERVER_OUTPUT_FORMAT};
$METHODS = $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};
$HEAD_BED = $ENV{APPRIS_WSERVER_HEADBEAD};
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
	my ($cgi) = new CGI;
	
	# Check input parameters:
	my ($url_path_info) = $cgi->path_info();
	my (@input_parameters) = split('/',$url_path_info);
		
	print_http_response(400)
		unless ( @input_parameters ); # defined path
		
	print_http_response(400)
		if ( scalar(@input_parameters) >= 4 );
	my ($type) = $input_parameters[1];
	my ($jobid) = $input_parameters[2];
	
	print_http_response(405, "The parameter $type is not allowed. The parameter must be: 'run', 'status', 'resulttypes', or 'result'.")
		if ( ($type ne 'run') and ($type ne 'status') and ($type ne 'resulttypes') and ($type ne 'result') );
	
	# get postdata if exits
	my ($pdata) = $cgi->param('POSTDATA');
	my ($postdata) = decode_json($pdata) if ( defined $pdata and ($pdata ne '') );
	
	# Collect parameters for execution:
	my ($params);	
	if ( $type eq 'run' )
	{					
		my ($species) = $cgi->param('species') || $postdata->{'species'} || undef;
		if ( defined $species ) {
			$params->{'species'} = $species;
		}
		else {
			print_http_response(400);
		}
			
		my ($id) = $cgi->param('id') || $postdata->{'id'} || undef;
		my ($e_version) = $cgi->param('e_version') || $postdata->{'e_version'} || undef;
		my ($sequences) = $cgi->param('sequences') || $postdata->{'sequences'} || undef;
	
		if ( defined $id and ($id ne '') and defined $e_version and ($e_version ne '') and !defined $sequences ) {
			$params->{'id'} = $id;
			$params->{'e_version'} = $e_version;
		}
		elsif ( defined $sequences and ($sequences ne '') and !defined $id and !defined $e_version ) {
			$params->{'sequences'} = $sequences;
		}
		else {
			print_http_response(400);
		}
		
		my ($methods) = $cgi->param('methods') || $postdata->{'methods'} || undef;
		if (defined $methods and ($methods ne '')) {
			$methods = lc($methods);		
			if ( 
				($methods =~ /appris/) or
				($methods =~ /firestar/) or
				($methods =~ /matador3d/) or
				($methods =~ /spade/) or
				($methods =~ /corsair/) or
				($methods =~ /crash/) or
				($methods =~ /thump/)
			) {
				$METHODS = $methods;	
			}
			else {
				print_http_response(405, "The parameter $methods is not allowed. The parameter must be: $METHODS.")
			}
		}
		$params->{'methods'} = $METHODS; # by default, we execute all methods
		
		my ($email) = $cgi->param('email') || $postdata->{'email'} || undef;
		if (defined $email and ($email ne '')) {
			if ( Email::Valid->address($email) ) {
				$params->{'email'} = $email;
			}
			else {
				print_http_response(405, "The parameter $email must be valid")
			}
		}	
	}
	elsif ( $type eq 'result' )
	{
		my ($methods) = $cgi->param('methods') || $postdata->{'methods'} || undef;
		if ( defined $methods ) {
			$params->{'methods'} = $methods;
		}
		
		my ($ids) = $cgi->param('ids') || $postdata->{'ids'} || undef;
		if ( defined $ids ) {
			$params->{'ids'} = $ids;
		}
		
		my ($res) = $cgi->param('res') || $postdata->{'res'} || undef;
		if ( defined $res ) {
			$params->{'res'} = $res;
		}

		my ($format) = $cgi->param('format') || $postdata->{'format'} || undef ;
		if (defined $format and ($format ne '')) {
			$format = lc($format);
			if (($format ne '') and ( ($format eq 'raw') or ($format eq 'tsv') or ($format eq 'gtf') or ($format eq 'bed') or ($format eq 'json') ) ) {
				$FORMAT = $format;	
			}
			else {
				print_http_response(405, "The parameter $format is not allowed.")
			}
		}
		$params->{'format'} = $FORMAT;
		
		my ($head) = $cgi->param('headbed') || $postdata->{'headbed'} || undef;
		if (defined $head and ($head ne '') and ($format eq 'bed') ) {
			$head =~ s/\s*//g;
			if ($head ne '' and ( ($head =~ /^yes/) or ($head =~ /^no/) or ($head =~ /^only/) ) ) {
				$HEAD_BED = $head;	
			}
			else {
				print_http_response(405, "The parameter $head is not allowed.")
			}
			$params->{'headbed'} = $HEAD_BED;
		}
		
		my ($email) = $cgi->param('email') || $postdata->{'email'} || undef;
		if (defined $email and ($email ne '')) {
			if ( Email::Valid->address($email) ) {
				$params->{'email'} = $email;
			}
			else {
				print_http_response(405, "The parameter $email must be valid")
			}
		}
	}
	elsif ( $type eq 'resulttypes' )
	{					
		my ($species) = $cgi->param('species') || $postdata->{'species'} || undef;
		if ( defined $species ) {
			$params->{'species'} = $species;
		}
	}	
	
	# Runner object
	my ($runner) = new WSRunner(
								-oper	=> $type,
								-params => $params,
								-jobid	=> $jobid,
	);

	# Print data
	my ($result);
	if ( $type eq 'run' ) {
		$result = $runner->apply_request($type);
	}
	elsif ( $type eq 'status' ) {
		$result = $runner->apply_request($type);
	}
	elsif ( $type eq 'resulttypes' ) {
		$result = $runner->apply_request($type);
	}
	elsif ( $type eq 'result' ) {
		$result = $runner->apply_request($type);
	}
	else {
		print_http_response(404);
	}
		
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

runner

=head1 DESCRIPTION

RESTful web services that execute APPRIS pipeline

=head1 SYNOPSIS

=head2 Required arguments:

=head1 EXAMPLES

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
