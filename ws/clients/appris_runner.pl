#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

###################
# Global variable #
###################

my ($SCRIPTNAME) = basename( $0, () ); # Get the script filename for use in usage messages
my ($BASE_URL) = 'http://apprisws.bioinfo.cnio.es'; # Base URL
my ($BASE_REST_URL) = $BASE_URL.'/rest'; # REST URL for service
my ($APPRIS_HEAD) = <<APPRIS_HEAD_DESC;
# ================================================================================= #
# ================== Thanks to use APPRIS RESTful web services ==================== #
#                                                                                   #
# APPRIS can be accessed through RESTful web services. Services can retrieve useful #
# information about genes/transcripts, or the results of individual APPRIS methods. #
#                                                                                   #
# APPRIS (http://appris-tools.org) is a web server for annotating alternative #
# splice isoforms in vertebrate genomes.                                            #
#                                                                                   #
#                                                                                   #
# If you have questions or comments, please write to:                               #
#   Jose Manuel Rodriguez, jmrodriguez\@cnio.es                                     #
#                                                                                   #
# ================================================================================= #
APPRIS_HEAD_DESC
my ($LOGLEVEL) = {
	'none'      => 1,
	'error'     => 2,
	'warning'   => 3,
	'info'      => 4,
	'debug'     => 5,
	'verbose'   => 6,
};
my ($INTERVALTIME) = 3; # Set interval for checking status
my ($DEFAULT_FORMAT) = 'json';


# Input parameters
my ($species) = undef;
my ($assembly) = undef;
my ($ensembl) = undef;
my ($gene) = undef;
my ($infile) = undef;
my ($outfile) = undef;
my ($inmethods) = undef;
my ($informat) = undef;
my ($inemail) = undef;
my ($jobid) = undef;
my ($async) = undef;
my ($status) = undef;
my ($resultTypes) = undef;
my ($polljob) = undef;
my ($loglevel) = 'info';
my ($logfile) = undef;
my ($logappend) = undef;

&GetOptions(
	'species|s=s'		=> \$species,
	'assembly|a=s'		=> \$assembly,	
	'gene|g=s'			=> \$gene,
	'infile|i=s'		=> \$infile,
	'outfile|o=s'		=> \$outfile,
	'methods|m=s'		=> \$inmethods,
	'format|f=s'		=> \$informat,
	'email|e=s'			=> \$inemail,
	'jobid=s'			=> \$jobid,
	'async'				=> \$async,
	'status'			=> \$status,
	'resultTypes'		=> \$resultTypes,
	'polljob'			=> \$polljob,
	'loglevel|l=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logappend'			=> \$logappend,
);

# Check input parameters
if ( !(defined $polljob or defined $resultTypes or defined $status) && !(defined $species) ) {
	print `perldoc $0`;
	print "\nBad option combination of input data\n\n";
	exit 1;
}
if ( !(defined $polljob or defined $resultTypes or defined $status) && !(defined $infile) && !(defined $gene and $assembly) ) {
	print `perldoc $0`;
	print "\nBad option combination of input data\n\n";
	exit 1;
}
# Lowercase for some input parameters
$species = lc($species) if ( defined $species );
$assembly = lc($assembly) if ( defined $assembly );
$loglevel = lc($loglevel) if ( defined $loglevel );

# Obtain the list of APPRIS methods we can execute
my ($METHODS) = undef;
my ($methods_json_str) = &appris_http::rest_request($BASE_URL.'/methods.json');
my ($methods_report) = JSON->new->decode($methods_json_str);
while ( my ($met_id,$met_report) = each(%{$methods_report}) ) {
	my ($met_name) = lc($met_report->{'name'});
	$METHODS->{$met_name} = $met_report->{'desc'};
}

# Obtain the list of species from server
# Get the common names of possible species with their assembly versions.
# The first value of array of assemblies is the default value.
my ($SPECIES_ASSEMBLIES)	= undef;
my ($SPECIES_ASSEMBLY)		= undef;
my ($species_json_str) = &appris_http::rest_request($BASE_URL.'/config.json');
my ($species_report) = JSON->new->decode($species_json_str);
while ( my ($spe_id,$spe_report) = each(%{$species_report->{'species'}}) ) {
	my ($spe_name) = lc($spe_report->{'common'});
	my ($spe_sci) = $spe_report->{'scientific'};
	my ($assemblies) = $spe_report->{'assemblies'};
	my ($assembly_default) = $spe_report->{'assemblies'}->[0];
	$SPECIES_ASSEMBLIES->{$spe_name} = $spe_report->{'assemblies'};
	$SPECIES_ASSEMBLY->{$spe_name} = lc($assembly_default->{'id'});
}

&check_parameters();

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	appris_http::print_log_message('appris_runner', 'Begin', 'info');
	
	# Job status
	if ( defined $status and defined $jobid ) {
		appris_http::print_job_status($jobid);
	}
	# Result types
	elsif ( defined $resultTypes and defined $jobid ) {
		appris_http::print_result_types($jobid);
	}
	# Poll job and get results
	elsif ( defined $polljob and defined $jobid ) {
		appris_http::get_results($jobid);
	}
	# Submit a job
	else {
		my ($params) = prepare_query_data($species, $infile, $gene, $assembly, $ensembl, $inmethods, $informat, $inemail);
		appris_http::submit_job($params);
	}
	
	appris_http::print_log_message('appris_runner', 'End', 'info');
}

# Check parameters
sub check_parameters()
{
	# Check if available specie is correct
	if ( defined $species ) {
		unless ( exists $SPECIES_ASSEMBLIES->{$species} ) {
			my ($l) = '';
			foreach my $assembly_rep (@{$SPECIES_ASSEMBLIES->{$species}}) {
				my ($ass) = $assembly_rep->{'name'}; $ass =~ s/\|/ or /g;
				$l .= $ass.' and ';
			}
			$l =~ s/ and $//;
			die "Species parameter is wrong. The correct values are: $l\n";
		}
		
		# Check if assembly is correct
		unless ( defined $assembly ) {
			$assembly = $SPECIES_ASSEMBLY->{$species};		
		}		
		# Obtain ensembl (gene dataset) version
		foreach my $assembly_rep (@{$SPECIES_ASSEMBLIES->{$species}}) {
			if ( $assembly_rep->{'id'} eq $assembly ) {
				$ensembl = $assembly_rep->{'datasets'}->[0]->{'source'}->{'version'};
			}
		}
		unless ( defined $ensembl ) {
			die "Gene dataset (Ensembl version) is not defined\n";
		}		
	}
	
	# Check if available methods are correct
	if ( defined $inmethods ) {
		foreach my $method ( split(',', $inmethods) ) {
			unless ( exists $METHODS->{$method} ) { die "Method $method do not exist\n" }
		}		
	}
		
} # end check_parameters

# Prepare query data
sub prepare_query_data($$$$$$$$)
{
	my ($species, $infile, $gene, $assembly, $ensembl, $inmethods, $informat, $inemail) = @_;
	my ($params);

	appris_http::print_log_message('prepare_query_data', 'Begin', 'info');
	appris_http::print_log_message('prepare_query_data', 'Species: ' . $species, 'debug') if ( defined $species );
	appris_http::print_log_message('prepare_query_data', 'Filename: ' . $infile, 'debug') if ( defined $infile );
	appris_http::print_log_message('prepare_query_data', 'Gene: ' . $gene, 'debug') if ( defined $gene );
	appris_http::print_log_message('prepare_query_data', 'Assembly: ' . $assembly . ' GeneDataset: ' . "Ensembl$ensembl", 'debug') if ( defined $assembly and defined $ensembl );
	appris_http::print_log_message('prepare_query_data', 'Methods: ' . $inmethods, 'debug') if ( defined $inmethods );
	appris_http::print_log_message('prepare_query_data', 'Format: ' . $informat, 'debug') if ( defined $informat );
	appris_http::print_log_message('prepare_query_data', 'Email: ' . $inemail, 'debug') if ( defined $inemail );
	
	# Extract query file
	my ($sequences);
	if ( defined $infile ) {
		if ( -f $infile || $infile eq '-' ) { # "-" for STDIN
			$sequences = appris_http::read_file( $infile );
		}
		else {
			die("Unable to open STDIN query file : $!");
		}		
	}
	
	# Add required parameters (exclusived)
	if ( defined $species and defined $sequences and (!(defined $gene) or !(defined $assembly)) ) {
		$params->{'species'} = $species;
		$params->{'sequences'} = $sequences;
	}
	elsif ( defined $species and defined $gene or defined $assembly and defined $ensembl and (!(defined $sequences) ) ) {
		$params->{'species'} = $species;
		$params->{'id'} = $gene;
		$params->{'e_version'} = $ensembl;		
	}
	
	# Add optional parameters
	if ( defined $inmethods ) { $params->{'methods'} = $inmethods }
	if ( defined $informat ) { $params->{'format'} = $informat }
	if ( defined $inemail ) { $params->{'email'} = $inemail }	
	
	appris_http::print_log_message('prepare_query_data', 'queryDat: ' . Dumper($params), 'verbose');
	appris_http::print_log_message('prepare_query_data', 'End', 'info');
	
	return $params;
	
} # end prepare_query_data


main();


##################
# COMMON PACKAGE #
##################
package appris_http;

use strict;
use warnings;
use English;
use Getopt::Long;
use LWP;
use JSON;
use Data::Dumper;

######################
# Common subroutines #
######################

# Print log message at specified debug level.
sub print_log_message($$$)
{
	my ($func, $msg, $trace) = @_;
	
	if ( exists $LOGLEVEL->{$loglevel} and exists $LOGLEVEL->{$trace} ) {
		if ( $LOGLEVEL->{$loglevel} >= $LOGLEVEL->{$trace} ) {
			my ($FH) = \*STDERR;
			if ( defined $logfile ) {
	    		my ($mode) = '>';
	    		$mode = '>>' if ( defined $logappend );			
				open ($FH, "$mode", $logfile) or die("Unable to open $logfile for writing: $!"); 
			}
			print $FH '[', $func, '] ', $msg, "\n";
		}
	}
	
} # end print_log_message

# Read a file into a scalar. The special filename '-' can be used to read from standard input (STDIN).
sub read_file($)
{
	my ($filename) = @_;
	my ($content, $buffer);
	
	print_log_message('read_file', 'Begin', 'info');
	print_log_message('read_file', 'Filename: ' . $filename, 'debug');
	
	if ( $filename eq '-' ) {
		while ( sysread( STDIN, $buffer, 1024 ) ) { $content .= $buffer }
	}
	else {
		open( my $FILE, '<', $filename ) or die "ERROR: unable to open input file $filename ($!)";
		while ( sysread( $FILE, $buffer, 1024 ) ) { $content .= $buffer }
		close($FILE);
	}
	
	print_log_message('read_file', 'End', 'info');
	
	return $content;
	
} # end read_file

# Write data to a file. The special filename '-' can be used to write to standard output (STDOUT).
sub write_file($$)
{
	my ($filename, $data) = @_;
	
	print_log_message('write_file', 'Begin', 'info');
	print_log_message('write_file', 'Filename: ' . $filename, 'debug');
	
	if ( $filename eq '-' ) {
		print STDOUT $data;
	}
	else {
		open( my $FILE, '>', $filename ) or die "ERROR: unable to open output file $filename ($!)";
		syswrite( $FILE, $data );
		close($FILE);
	}
	
	print_log_message('write_file', 'End', 'info');
	
} # end write_file


###########################
# Common HTTP subroutines #
###########################

# Get result data of a specified type for a finished job.
sub rest_get_exporter($;$;$)
{
	my ($queryid, $methods, $format) = @_;
	
	print_log_message('rest_get_exporter', 'Begin', 'debug');
	
	my ($url) = $BASE_REST_URL . '/exporter/' . $queryid;
	$url	 .=  '?' if ( defined $methods or defined $format or defined $ensembl );
	$url     .=  'methods=' . $methods . '&'  if ( defined $methods);
	$url     .=  'format=' . $format . '&' if ( defined $format);
	$url     .=  'ens=' . $ensembl . '&' if ( defined $ensembl); # required
	$url =~ s/&$//g;
	my ($result) = rest_request($url);
	
	print_log_message('rest_get_exporter', length($result) . ' characters', 'debug');
	print_log_message('rest_get_exporter', 'End', 'debug');
	
	return $result;
	
} # end rest_get_exporter

# Check the status of a job.
sub rest_get_status($)
{
	my ($jobid) = @_;
	my ($status_str,$log_str) = ('UNKNOWN','');
	
	print_log_message('rest_get_status', 'Begin', 'debug');
	
	my ($url) = $BASE_REST_URL . '/runner/status/' . $jobid;
	my ($status_json_str) = rest_request($url);	
	if ( defined $status_json_str ) {
		my ($status_json) = JSON->new->decode($status_json_str);
		$status_str = $status_json->{'status'};
		$log_str = $status_json->{'log'};
	}
		
	print_log_message('rest_get_status', 'Status_str:Log_str ' . $status_str . ':' . $log_str, 'debug');
	print_log_message('rest_get_status', 'End', 'debug');
	
	return ($status_str,$log_str);
	
} # end rest_get_status

# Get list of result types for finished job.
sub rest_get_result_types($)
{
	my ($jobid) = @_;
	my (@resultTypes);
		
	print_log_message('rest_get_result_types', 'Begin', 'debug');
	
	my ($url) = $BASE_REST_URL . '/runner/resulttypes/' . $jobid;
	my ($result_type_list_json_str) = rest_request($url);
	if ( defined $result_type_list_json_str ) {
		(@resultTypes) = @{JSON->new->decode($result_type_list_json_str)};		
	}
	
	print_log_message('rest_get_result_types', scalar(@resultTypes) . ' result types', 'debug');
	print_log_message('rest_get_result_types', 'End', 'debug');
	
	return (@resultTypes);
	
} # end rest_get_result_types

# Get result data of a specified type for a finished job.
sub rest_get_result($;$;$)
{
	my ($jobid, $type, $format) = @_;
	my (@resultTypes);
		
	print_log_message('rest_get_result', 'Begin', 'debug');
	print_log_message('rest_get_result', 'Type: ' . $type, 'debug') if ( defined $type);
	print_log_message('rest_get_result', 'Format: ' . $format, 'debug') if ( defined $format);

	my ($url) = $BASE_REST_URL . '/runner/result/' . $jobid;
	$url     .=  '?' if ( defined $type or defined $format );
	$url     .=  'methods=' . $type . '&'  if ( defined $type);
	$url     .=  'format=' . $format . '&' if ( defined $format);
	$url =~ s/&$//g;
	my ($result) = rest_request($url);
		
	print_log_message('rest_get_result', length($result) . ' characters', 'debug');
	print_log_message('rest_get_result', 'End', 'debug');	
	
	return $result;
	
} # end rest_get_result

# Perform a REST request (HTTP GET).
sub rest_request($)
{
	my ($requestUrl) = @_;
	
	print_log_message('rest_request', 'Begin', 'info');
	print_log_message('rest_request', 'URL: ' . $requestUrl, 'debug');

	# Get an LWP UserAgent.
	my ($ua) = rest_user_agent();
	
	# Available HTTP compression methods.
	my ($can_accept);
	eval {
	    $can_accept = HTTP::Message::decodable();
	};
	$can_accept = '' unless defined($can_accept);
	
	# Perform the request
	my ($response) = $ua->get( $requestUrl, 'Accept-Encoding' => $can_accept ); # HTTP compression
	print_log_message('rest_request', 'HTTP status: ' . $response->code, 'debug');
	print_log_message('rest_request', 'response length: ' . length($response->content()), 'debug');
	print_log_message('rest_request', 'request:' ."\n" . $response->request()->as_string(), 'verbose');
	print_log_message('rest_request', 'response: ' . "\n" . $response->as_string(), 'verbose');
	
	# Unpack possibly compressed response.
	my ($retVal);
	if ( defined($can_accept) && $can_accept ne '') {
	    $retVal = $response->decoded_content();
	}
	
	# If unable to decode use orginal content.
	$retVal = $response->content() unless defined($retVal);
	
	# Check for an error.
	rest_error($response, $retVal);
	
	print_log_message('rest_request', 'End', 'info');

	return $retVal;
	
} # end rest_request

# Submit a job.
sub post_run($)
{	
	my ($params) = @_;
		
	print_log_message('post_run', 'Begin', 'info');
	print_log_message('post_run', 'Params: ' . Dumper($params), 'verbose');
	
	# Get an LWP UserAgent.
	my ($ua) = rest_user_agent();
	
	# Submit the job as a POST
	my ($url) = $BASE_REST_URL . '/runner/run';	
	my ($response) = $ua->post( $url, $params );
	print_log_message('post_run', 'HTTP status: ' . $response->code, 'debug');
	print_log_message('post_run', 'response length: ' . length($response->content()), 'debug');
	print_log_message('post_run', 'request:' ."\n" . $response->request()->as_string(), 'verbose');
	print_log_message('post_run', 'response: ' . "\n" . $response->as_string(), 'verbose');
	
	# Check for an error.
	rest_error($response);

	# The job id is returned
	my ($jobid) = undef;
	my ($rest_run_out_json) = $response->content();	
	if ( defined $rest_run_out_json ) {
		$jobid = JSON->new->decode($rest_run_out_json)->{'jobid'};
	}
	
	print_log_message('post_run', 'End', 'info');
	
	return $jobid;
	
} # end post_run

# Get a LWP UserAgent to use to perform REST requests.
sub rest_user_agent()
{
	print_log_message('rest_user_agent', 'Begin', 'info');
	
	# Create an LWP UserAgent for making HTTP calls.
	my $ua = LWP::UserAgent->new();
	
	# Set 'User-Agent' HTTP header to identifiy the client.	
	$ua->agent("APPRIS-REST-Client/ ($SCRIPTNAME;$OSNAME) " . $ua->agent());
	
	# Configure HTTP proxy support from environment.
	$ua->env_proxy;
	
	print_log_message('rest_user_agent', 'End', 'info');
	
	return $ua;
	
} # end rest_user_agent

# Check a REST response for an error condition. An error is mapped to a die.
sub rest_error($;$)
{
	my ($response, $contentdata) = @_;
	
	print_log_message('rest_error', 'Begin', 'info');

	if ( !defined($contentdata) || $contentdata eq '') {
		$contentdata = $response->content();
	}
	
	# Check for HTTP error codes
	if ( $response->is_error ) {
		my ($message) = '';
		if ( $response->code == 500 ) {
			die 'http status: ' . $response->code . ' ' . $response->message . '  ' . $message;
		}
	}
	
	print_log_message('rest_error', 'End', 'info');
	
} # end rest_error

# Print status of a job.
sub print_job_status($)
{
	my ($jobid) = @_;
	
	print_log_message('print_job_status', 'Begin', 'info');
	print_log_message('print_job_status', 'Jobid: ' . $jobid, 'debug');
	
	my ($status,$statusLog) = rest_get_status($jobid);
	print STDOUT "$status\n$statusLog";
	if ( $status eq 'FINISHED' ) {
		print STDERR "\n", 'To get results:', "\n";
		print STDERR "  $SCRIPTNAME --polljob --jobid $jobid \n";
		print STDERR "  $SCRIPTNAME --polljob --jobid $jobid --format <type> \n";
	}
	
	print_log_message('print_job_status', 'End', 'info');
	
}

# Print available result types for a job.
sub print_result_types($)
{
	my ($jobid) = @_;
	
	print_log_message('print_result_types', 'Begin', 'info');
	print_log_message('print_result_types', 'Jobid: ' . $jobid, 'debug');
	
	my ($status,$statusLog) = rest_get_status($jobid);
	if ( $status eq 'PENDING' || $status eq 'RUNNING' ) {
		print STDERR 'ERROR: Job status is ', $status,'. To get result types the job must be finished.', "\n";
	}
	else {
		my (@resultTypes) = rest_get_result_types($jobid);
		foreach my $resultType (@resultTypes) {
			print STDOUT $resultType->{'name'}, "\n";
			if ( defined( $resultType->{'description'} ) ) {
				print STDOUT "\t", $resultType->{'description'}, "\n";
			}
			foreach my $rstType (@{$resultType->{'types'}}) {
				if ( defined( $rstType->{'fileSuffix'} ) and defined( $rstType->{'mediaType'} ) ) {
					print STDOUT "\t", $rstType->{'fileSuffix'}, ",",$rstType->{'mediaType'}, "\n";
				}
			}
		}
		if ( $status eq 'FINISHED' ) {
			print STDERR "\n", 'To get results:', "\n";
			print STDERR "  $SCRIPTNAME --polljob --jobid $jobid \n";
			print STDERR "  $SCRIPTNAME --polljob --jobid $jobid --format <type> \n";
		}
	}
	
	print_log_message('print_result_types', 'End', 'info');
	
}

# Client-side job polling.
sub client_poll
{
	my ($jobid) = @_;
	my ($status,$statusLog) = ('PENDING','The job is waiting to be processed');
	
	print_log_message('client_poll', 'Begin', 'debug');

	# Check status and wait if not finished. Terminate if three attempts get "ERROR".
	my ($errorCount) = 0;
	while ( $status eq 'RUNNING' || $status eq 'PENDING' || ( $status eq 'ERROR' && $errorCount < 2 ) ) {
		($status,$statusLog) = rest_get_status($jobid);
		if ( $status eq 'ERROR' ) {
			$errorCount++;
		}
		elsif ( $errorCount > 0 ) {
			$errorCount--;
		}
		# Wait before polling again.
		if ( $status eq 'RUNNING' || $status eq 'PENDING' || $status eq 'ERROR' ) {
			sleep $INTERVALTIME;
		}
	}
	
	print_log_message('client_poll', 'Begin', 'debug');
	
	return $status;
	
} # end client_poll

# Get the results for a job identifier.
sub get_results($)
{
	my ($jobid) = @_;
	
	print_log_message('get_results', 'Begin', 'info');
	print_log_message('get_results', 'Jobid: ' . $jobid, 'debug');

	# Check status, and wait if not finished
	client_poll($jobid);

	# Use JobId if output file name is not defined
	unless ( defined $outfile ) {
		$outfile = $jobid;
	}

	# Get results depending on methods
	my (@resultTypes) = rest_get_result_types($jobid);
	if ( scalar(@resultTypes) > 0 ) {
		foreach my $resultType (@resultTypes) {
			if ( $resultType->{'name'} ) {
				my ($name) = $resultType->{'name'};
				my ($format);
				if ( $informat ) {
					$format = $informat;
				}
				else {
					#my ($rstTypes) = $resultType->{'types'};
					#$format = $rstTypes->[0]->{'fileSuffix'};
					$format = $DEFAULT_FORMAT;
				}
				my ($result) = '';
				if ( defined $inmethods ) {
					if ( $inmethods =~ /$name/ ) {			
						$result = rest_get_result($jobid, $name, $format);
					}
				}
				else {
					$result = rest_get_result($jobid, $name, $format);
				}
				if ( $result ne '' ) {
					if ( $outfile eq '-' ) {
						write_file($outfile, $result);
					}
					else {
						write_file($outfile . '.' . $name . '.' . $format, $result);
					}						
				}
			}
		}
	}
	print_log_message('get_results', 'End', 'info');
	
} # end get_results

# Submit a job to the service.
sub submit_job($)
{
	my ($params) = @_;
	
	print_log_message('submit_job', 'Begin', 'info');
	print_log_message('submit_job', 'Params: ' . Dumper($params), 'verbose');	
	
	# Submit the job
	my ($jobid) = post_run($params);

	# Simulate sync/async mode
	if ( defined $async ) {
		print STDOUT $jobid, "\n";
		print STDERR "To check status: $SCRIPTNAME --status --jobid $jobid\n";
	}
	else {
		print STDERR "JobId: $jobid\n";
		sleep 1;
		get_results($jobid);
	}
	
	print_log_message('submit_job', 'End', 'info');
	
} # end submit_job


1;

__END__

=pod

=head1 NAME

appris_runner

=head1 DESCRIPTION

script that executes APPRIS RESTful services 

=head1 SYNOPSIS

appris_runner

=head2 Required arguments (data input):

  -s, --species {string} <Name of species: human, mouse, rat, zebrafish>

=head2 Required arguments (exclusived):

  -g, --gene {string} <Ensembl gene id or gene name>
  
  -a, --assembly {string} <Assemby of given specie>
    * human:     GRCh38 (Ensembl77 gene dataset) - default -     
                 GRCh37 (Ensembl74 gene dataset)
    * mouse:     GRCm38 (Ensembl77 gene dataset)	  
    * rat:       Rnor_5.0 (Ensembl77 gene dataset)	  
    * zebrafish: Zv9 (Ensembl77 gene dataset)
      
Or
  
  -i, --infile {file} <Input file that contains the queries. Tabular or CSV file>
	
=head2 Required arguments (outputs):
	
  -o, --outfile {file} <Output file>
  
=head2 Optional arguments:
		
  -m, --methods {vector} <List of retrieved methods separated by commas: 'firestar,matador3d,spade,corsair,thump,crash,appris' (default: all)>
	* Functionally important residues, firestar
	* Protein structural information, Matador3D
	* Presence of whole protein domains, SPADE
	* Conservation against vertebrates, CORSAIR
	* Presence of whole trans-membrane helices, THUMP
	* Prediction of signal peptide and sub-cellular location, CRASH
	* Proteomic evidence, PROTEO
  
  -f, --format {vector} <Output format: 'gtf', 'bed', or 'json' (default: JSON)>
  
  -e, --email {email} <E-mail address>
  
=head2 Runner arguments:
	
  --jobid {string} <Jobid that was returned when an asynchronous job>
  
  --async <Forces to make an asynchronous query>
  
  --status <Get job status>
  
  --resultTypes <Get available result types for job>
  
  --polljob <Poll for the status of a job was submitted>
    
=head2 Optional arguments (log arguments):
	
  -l, --loglevel {vector} <Define log level: INFO,DEBUG,NONE (default: NONE)>
  
  --logfile {file} <Log to FILE (default: *STDERR)>
  
  --logappend <Append to logfile (default: truncate)>
  
=head1 EXAMPLE:

=head2 Synchronous job:

The results/errors are returned as soon as the job is finished.
  
  > perl appris_runner.pl
    --species=human
    --assembly=hg19
    --gene=ENSG00000099904
    --outfile=ENSG00000099904.out.gtf
    --methods=firestar,spade
    --format=gtf
	
  > perl appris_runner.pl -s human -g ENSG00000099904 -a hg19 -o ENSG00000099904.out -m firestar,spade -f gtf
  
  > perl appris_runner.pl -s human -i seq1.fa -o seq1.out -m firestar,spade -f tsv
  
=head2 Asynchronous job:

  Use this if you want to retrieve the results at a later time. The results are stored for up to 1 week. 	

  > perl appris_runner.pl
    --async
    --species=human
    --assembly=hg19
    --gene=ENSG00000099904
    --outfile=ENSG00000099904.out.gtf
    --methods=firestar,spade
    --format=gtf
	
  > perl appris_runner.pl --async -s human -i seq1.fa -o seq1.out -m firestar,spade -f tsv
  
  > perl appris_runner.pl --status --jobid <jobId>
  Returns: string indicating the status of the job.
  
  > perl appris_runner.pl --resultTypes --jobid <jobId>
  Returns: available result types for job
  
  > perl appris_runner.pl --polljob --jobid <jobId> [--outfile string]
  Returns: string indicating the status of the job and if applicable, results as an attachment.


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB,CNIO)

=cut
