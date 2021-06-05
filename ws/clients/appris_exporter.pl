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
my ($BASE_URL) = 'https://apprisws.bioinfo.cnio.es'; # Base URL
my ($BASE_REST_URL) = $BASE_URL.'/rest'; # REST URL for service
my ($APPRIS_HEAD) = <<APPRIS_HEAD_DESC;
# ================================================================================= #
# ================== Thanks to use APPRIS RESTful web services ==================== #
#                                                                                   #
# APPRIS can be accessed through RESTful web services. Services can retrieve useful #
# information about genes/transcripts, or the results of individual APPRIS methods. #
#                                                                                   #
# APPRIS (https://appris.bioinfo.cnio.es/) is a web server for annotating           #
# alternative splice isoforms in vertebrate genomes.                                #
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

# Input parameters
my ($species) = undef;
my ($assembly) = undef;
my ($ensembl) = undef;
my ($infile) = undef;
my ($id_col) = undef;
my ($name_col) = undef;
my ($pos_col) = undef;
my ($inmethods) = undef;
my ($informat) = undef;
my ($outfile) = undef;
my ($loglevel) = 'info';
my ($logfile) = undef;
my ($logappend) = undef;
  
&GetOptions(
	'species|s=s'		=> \$species,
	'assembly|a=s'		=> \$assembly,
	'infile|i=s'		=> \$infile,
	'gene-id|id=s'		=> \$id_col,
	'gene-name|n=s'		=> \$name_col,
	'position|p=s'		=> \$pos_col,
	'methods|m=s'		=> \$inmethods,
	'format|f=s'		=> \$informat,
	'outfile|o=s'		=> \$outfile,
	'loglevel|l=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logappend'			=> \$logappend,
);

# Check input parameters
unless ( defined $species and defined $infile and defined $outfile ) {
	print `perldoc $0`;
	print "\nPlease, check parameters\n\n";
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
my ($SPECIES_IDS) 			= undef;
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
	appris_http::print_log_message('appris_exporter', 'Begin', 'info');
	
	my ($params) = prepare_query_data($infile, $id_col, $name_col, $pos_col, $inmethods, $informat);
	
	my ($result) = run_query($params, $outfile);
	
	appris_http::write_file($outfile, $result);	
	
	appris_http::print_log_message('appris_exporter', 'End', 'info');
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
				my ($datasets_rep) = [
					grep { ! exists($_->{'queryable'}) || $_->{'queryable'} }
					@{$assembly_rep->{'datasets'}}
				];
				$ensembl = $datasets_rep->[0]->{'source'}->{'version'};
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
	
	# Column data has to be integer
	if ( defined $id_col or defined $name_col or defined $pos_col ) {
		if ( defined $id_col and !($id_col =~ /^\d*$/) ) { die "Num. Column of gene ids has to be integer\n" }
		if ( defined $name_col and !($name_col =~ /^\d*$/) ) { die "Num. Column of gene names has to be integer\n" }
		if ( defined $pos_col and !($pos_col =~ /^\d*$/) ) { die "Num. Column of gene ids has to be integer\n" }		
	}
	else {
		die "Data column parameter has not been defined\n"
	}
		
} # end check_parameters

# Prepare query data
sub prepare_query_data($$$$$$)
{
	my ($infile, $id_col, $name_col, $pos_col, $inmethods, $informat) = @_;
	my ($params);

	appris_http::print_log_message('prepare_query_data', 'Begin', 'info');
	appris_http::print_log_message('prepare_query_data', 'Filename: ' . $infile, 'debug');
	appris_http::print_log_message('prepare_query_data', 'Col. ids: ' . $id_col, 'debug') if ( defined $id_col );
	appris_http::print_log_message('prepare_query_data', 'Col. names: ' . $name_col, 'debug') if ( defined $name_col );
	appris_http::print_log_message('prepare_query_data', 'Col. pos: ' . $pos_col, 'debug') if ( defined $pos_col );
	
	# Extract query file
	my ($query_data);
	if ( -f $infile || $infile eq '-' ) { # "-" for STDIN
		$query_data = appris_http::read_file( $infile );
	}
	else {
		die("Unable to open STDIN query file : $!");
	}
		
	# Extract the query data from given column
	foreach my $line ( split('\n',$query_data) ) {
		my ($query, $data) = ('','');
		my (@cols) = split(/\t|\,/, $line);
		if ( scalar(@cols) > 0 ) {
			if ( defined $id_col and defined $cols[$id_col] ) {
				$data .= $cols[$id_col];
				$query .= 'id/'.$SPECIES_IDS->{$species}.'/'.$data;
			} #else { die("Not defined column data\n") }
			if ( defined $name_col and defined $cols[$name_col] ) {
				$data .= $cols[$name_col];
				$query .= 'name/'.$SPECIES_IDS->{$species}.'/'.$data;
			} #else { die("Not defined column data\n") }
			if ( defined $pos_col and defined $cols[$pos_col] and ($cols[$pos_col] =~ /^([^\:]*)\:(\d*)-(\d*)$/) ) {
				$data .= $cols[$pos_col];
				$query .= 'position/'.$SPECIES_IDS->{$species}.'/'.$data;
			} #else { die("Not defined column data\n") }
		}
		if ( $query ne '' ) {
			push(@{$params->{'queries'}}, {				
				'raw'	=> $line,
				'data'	=> $data,
				'query'	=> $query,
			}); 
		}
	}
	
	# Add optional parameters
	if ( defined $inmethods ) { $params->{'methods'} = $inmethods }
	if ( defined $informat ) { $params->{'format'} = $informat }	
	
	appris_http::print_log_message('prepare_query_data', 'queryDat: ' . Dumper($params), 'verbose');
	appris_http::print_log_message('prepare_query_data', 'End', 'info');
	
	return $params;
	
} # end prepare_query_data

sub run_query($$)
{
	my ($params, $ouffile) = @_;
	my ($result) = $APPRIS_HEAD;
	
	appris_http::print_log_message('run_query', 'Begin', 'info');
	appris_http::print_log_message('prepare_query_data', 'Num. queries: ' . scalar(@{$params->{'queries'}}), 'debug');
	
	if ( exists $params->{'queries'} ) {
		my ($methods) = $params->{'methods'} if ( exists $params->{'methods'} );
		my ($format) = $params->{'format'} if ( exists $params->{'format'} );
		foreach my $parm ( @{$params->{'queries'}} ) {
			my ($raw) = $parm->{'raw'};
			my ($data) = $parm->{'data'};
			my ($query) = $parm->{'query'};		
			my ($q_result) = appris_http::rest_get_exporter($query,$methods,$format);
			if ( defined $q_result ) {
				$result .= "# QUERY_RAW:$raw\n";
				$result .= "# QUERY_DATA:$data\n";
				$result .= "# QUERY_REST:$query\n";
				$result .= "$q_result\n";
			}
		}	
	}
	
	return $result;
	
	appris_http::print_log_message('run_query', 'End', 'info');
}

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
	my ($response) = $ua->get($requestUrl,
		'Accept-Encoding' => $can_accept, # HTTP compression.
	);
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


1;

__END__

=pod

=head1 NAME

appris_exporter

=head1 DESCRIPTION

script that make a query to APPRIS database using RESTful services 

=head1 SYNOPSIS

appris_exporter

=head2 Required arguments (data input):

  -s, --species {string} <Name of species: human, mouse, rat, zebrafish>

  -a, --assembly {string} <Assemby of given specie>
    * human:     GRCh38 (Ensembl77 gene dataset) - default -     
                 GRCh37 (Ensembl74 gene dataset)
    * mouse:     GRCm38 (Ensembl77 gene dataset)	  
    * rat:       Rnor_5.0 (Ensembl77 gene dataset)	  
    * zebrafish: Zv9 (Ensembl77 gene dataset)
  
  -i, --infile {file} <Input file that contains the queries. Tabular or CSV file>

=head2 Required arguments (exclusived):

  -id, --gene-id {integer} <Num. column that contains the Ensembl gene ids to make the queries>
	
  -n, --gene-name {integer} <Num. column that contains the Ensembl gene ids to make the queries>
	
  -p, --position {integer} <Num. column that contains the query data>

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
  
=head2 Optional arguments (log arguments):
	
  -l, --loglevel {vector} <Define log level: INFO,DEBUG,NONE (default: NONE)>
  
  --logfile {file} <Log to FILE (default: *STDERR)>
  
  --logappend <Append to logfile (default: truncate)>
  
=head1 EXAMPLEs

  > perl appris_exporter.pl
    --species=human
    --assembly=hg19
    --infile=query1.csv
    --position=1
    --outfile=query1.out.gtf
    --methods=firestar,spade,proteo
    --format=gtf
    --logfile=appris_exporter.log
    --loglevel=debug
    --logappend
    
  > perl appris_exporter.pl
    -s human
    -a hg19
    -p 1
    -i query1.csv
    -o query1.out.gtf
    -m firestar,spade,proteo
    -f gtf
    -l debug
        	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB,CNIO)

=cut
