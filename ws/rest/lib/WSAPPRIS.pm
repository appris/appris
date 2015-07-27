=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

WSAPPRIS - Utility functions for error handling

=head1 SYNOPSIS

  use WSAPPRIS qw(get_annotations);
  
  or to get all methods just

  use WSAPPRIS;

  eval { get_annotations("text to file",file_path) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to export database values as JSON format.

=head1 METHODS

=cut

package WSAPPRIS;

use strict;
use warnings;
use Data::Dumper;
use JSON;

use WSRetriever;

use APPRIS::Utils::Exception qw(throw warning deprecate info);
use APPRIS::Utils::File qw( rm_dir prepare_workspace getStringFromFile printStringIntoFile updateStringIntoFile );
use APPRIS::Utils::Clusters;
use APPRIS::Parser;
use APPRIS::Exporter;

###################
# Global variable #
###################
use vars qw(
	$STATUS_LOGS
	$SPECIES
	$METHODS
	$METHOD_IDS
);

$STATUS_LOGS = {
	'PENDING'	=> 'The job is waiting to be processed',
	'RUNNING'	=> 'The job is currently being processed',
	'FINISHED'	=> 'Job has finished, and the results can then be retrieved',
	'ERROR'		=> 'An error occurred attempting to get the job status',
	'FAILURE'	=> 'The job failed',
	'NOT_FOUND'	=> 'The job cannot be found'
};

$SPECIES = JSON->new()->decode( getStringFromFile($ENV{APPRIS_WSERVER_HOME} . '/species.json') );
$METHODS = JSON->new()->decode( getStringFromFile($ENV{APPRIS_WSERVER_HOME} . '/methods.json') );
while ( my ($met_id,$met_report) = each(%{$METHODS}) ) {
	my ($met_name) = lc($met_report->{'name'});
	$METHOD_IDS->{$met_name} = $met_id;
}

sub run_appris {
	my ($jobid) = shift;	
	my ($params) = shift;
	my ($run);
	
	if ( defined $jobid and defined $params ) {
		
		# prepare workspace
		my ($workspace) = $ENV{APPRIS_WORKSPACE} . '/' . $jobid;
		$workspace = rm_dir($workspace); # clear old executions
		$workspace = prepare_workspace($workspace);

		# check parameters...
		my ($p_status, $p_log) = 1;
		my ($parameters) = '';
		
		# ...check species
		my ($species);
		if ( exists $params->{'species'} ) {
			$species = _create_param_species($params->{'species'});
			($p_status, $p_log) = (0, "'Species' parameter is invalid. Please, rephrase your query.") unless ( defined $species );
		}
		if ( $p_status == 1 ) {
			while ( my ($k,$v) = each(%{$params}) ) {
				if ( $k eq 'species' ) {
					$parameters .= " --species='$species' " if ( defined $species );
				}
				# For Web Server, we discard ID and Ensembl version parameters.
				# We only execute with Gencode mode or sequence
				elsif ( $k eq 'e_version' ) {
					$parameters .= " --e-version=$v ";
				}
				#elsif ( $k eq 'id' ) {
				#	$parameters .= " --id=$v ";
				#}
				elsif ( $k eq 'sequences' ) {
					my ($transl_file) = _create_param_seq($workspace, $jobid, $v);
					if ( defined $transl_file and ($transl_file ne '') ) {				
						$parameters .= " --transl=$transl_file ";
					}
					else { ($p_status, $p_log) = (0, "'Sequences' parameter is invalid. Please, rephrase your query.") }
				}
				elsif ( $k eq 'methods' ) {
					my ($mets) = _create_param_methods($v);
					if ( defined $mets and ($mets ne '') ) {
						$parameters .= " --methods=$mets ";	
					}
					else { ($p_status, $p_log) = (0, "'Methods' parameter is invalid. Please, rephrase your query.") }
				}
			}
			#unless ( exists $params->{'id'} ) {
			#	$parameters .= " --id=$jobid ";
			#}
			$parameters .= " --outpath=$workspace ";
			# create data infiles for gencode/ensembl case
			if ( exists $params->{'id'} and exists $params->{'e_version'} and exists $params->{'species'} ) {
				my ($data_file, $transc_file, $transl_file) = _create_param_ensembl($workspace, $params->{'id'}, $species, $params->{'e_version'});
				if ( defined $data_file and defined $transc_file and defined $transl_file ) {
					if ( -e ($data_file) and (-s $data_file > 0) and
						 -e ($transc_file) and (-s $transc_file > 0) and
						 -e ($transl_file) and (-s $transl_file > 0) 
					) {
						$parameters .= " --data=$data_file ";
						$parameters .= " --transc=$transc_file ";
						$parameters .= " --transl=$transl_file ";
					}
					else { ($p_status, $p_log) = (0, "'Gene' parameter is invalid. Please, rephrase your query.") }
				}
			}
		}
		# email
		if ( $params->{'email'} ) {
			#$retriever->add_control("RUNNER_EMAIL",$params->{'email'});
			$parameters .= " --email=".$params->{'email'}." ";
		}
		
		# Create Retriever object from job id
		my ($retriever) = new WSRetriever( -jobid => $jobid );
		# checking param
		if ( $p_status == 0 ) {
			$retriever->add_control("RUNNER_STATUS",'ERROR');
			$retriever->add_control("RUNNER_LOG",$p_log);
			$run = {
				'jobid'	=> $jobid
			};		
			return $run;
		}
								
		# run script!!!
		my ($logfile) = $ENV{APPRIS_WORKSPACE_LOG_DIR} . '/' . $jobid . '.log'; # log params
		my ($loglevel) = 'INFO';		
		if ( defined $ENV{APPRIS_SCRIPTS_CLUSTER_INI_WSERVER} ) { # in cluster
			my ($pid) = fork();
			if ( !defined $pid ) {
			    die "Cannot fork: $!";
			}
			elsif ( $pid == 0 ) {
				close STDOUT; close STDERR; # so parent can go on
			    # client process
				my ($cmd) =	"perl $ENV{APPRIS_SCRIPTS_DIR}/apprisall.pl ".
							" --wserver=$jobid ".
							" $parameters ".
							" --cluster-conf=$ENV{APPRIS_SCRIPTS_CLUSTER_INI_WSERVER} ".
							" --loglevel=$loglevel --logfile=$logfile --logappend ";
				my ($cmd_status) = system($cmd);
			    exit 0;
			}
		}
		else { # in local
			my ($pid) = fork();
			if ( !defined $pid ) {
			    die "Cannot fork: $!";
			}
			elsif ( $pid == 0 ) {
				close STDOUT; close STDERR; # so parent can go on
			    # client process
				my ($cmd) =	"perl $ENV{APPRIS_SCRIPTS_DIR}/apprisall.pl ".
							" --wserver=$jobid ".
							" $parameters ".
							" --loglevel=$loglevel --logfile=$logfile --logappend ";
				my ($cmd_status) = system($cmd);
			    exit 0;
			}			
		}
		$run = {
			'jobid'	=> $jobid
		};		
	}
	else {
		$run = {
			'status'	=> "ERROR",
			'log'		=> "Your query is malformed. Please rephrase your query"
		};
	}
				
	return $run;
	
} # end run_appris

sub status_appris {
	my ($jobid) = shift;	
	my ($jobstatus,$joblog,$jobprogress);

	if ( defined $jobid ) {
			
		# Create Retriever object from job id
		my ($retriever) = new WSRetriever( -jobid => $jobid );
		
		# Parse trace log from local or cluster server
		if ( defined $ENV{APPRIS_SCRIPTS_CLUSTER_INI_WSERVER} ) {
			if ( $retriever->host and $retriever->jobid and $retriever->wsbase ) {
				my ($clusters) = new APPRIS::Utils::Clusters( -conf => $ENV{APPRIS_SCRIPTS_CLUSTER_INI_WSERVER} );
				my ($clusterstatus,$clusterlog) = $clusters->jstatuslogs($retriever->host, $retriever->jobid, $retriever->wsbase);
				($jobstatus,$joblog,$jobprogress) = $retriever->status_log($clusterstatus, $clusterlog);				
			}
		}
		else {
			($jobstatus,$joblog,$jobprogress) = $retriever->status_log();
		}
	}
	unless ( defined $jobstatus and defined $joblog ) {
		($jobstatus, $joblog, $jobprogress) = ('NOT_FOUND', $STATUS_LOGS->{'NOT_FOUND'}, 0);
	}
	my ($status) = {
		'status'	=> $jobstatus,
		'log'		=> $joblog,
		'progress'	=> $jobprogress
	};
	
	return $status;
	
} # end status_appris

sub resulttypes_appris {
	my ($jobid) = shift;
	my ($params) = shift;
	my ($resulttypes);
	my (@methods);
	
	my ($run_mode, $run_methods);
	if ( defined $jobid ) { # runner mode
				
		# Create Retriever object from job id
		my ($retriever) = new WSRetriever( -jobid => $jobid );
		if ( $retriever->data and $retriever->transc and $retriever->transl ) { $run_mode = 'gencode' }
		elsif ( $retriever->transl ) { $run_mode = 'sequence' }
		
		if ( $retriever->methods ) { $run_methods = $retriever->methods }
	}
	else { # exporter mode
		$run_mode = 'gencode';
		$run_methods = $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};
		# Add proteo method in the case of human
		if ( defined $params ) {
			if ( exists $params->{'species'} ) {
				my ($species) = _create_param_species($params->{'species'});
				$run_methods .= ',proteo' if ( $species eq 'Homo sapiens' );
			}			
		}
	}
	# Retrieve type of results
	foreach my $method ( split(',', $run_methods) ) {
		my ($m_results) = [
			{
				'fileSuffix'	=> 'raw',
				'mediaType'		=> 'text/plain',		
			},
			{
				'fileSuffix'	=> 'tsv',
				'mediaType'		=> 'text/plain',		
			},
			{
				'fileSuffix'	=> 'json',
				'mediaType'		=> 'application/json',		
			}		
		];
		if ( defined $run_mode ) {
			if ( $run_mode eq 'gencode' ) {
				push(@{$m_results}, {
							'fileSuffix'	=> 'gtf',
							'mediaType'		=> 'text/plain' });
				push(@{$m_results}, {
							'fileSuffix'	=> 'bed',
							'mediaType'		=> 'text/plain' });
			}
		}
		my ($method_id) = $METHOD_IDS->{$method};
		my ($type) = {			
			'id'			=> $method_id,
			'name'			=> $METHODS->{$method_id}->{'name'},
			'label'			=> $METHODS->{$method_id}->{'label'},
			'desc'			=> $METHODS->{$method_id}->{'desc'},
			'mode'			=> $run_mode,
			'types'			=> $m_results
		};
		push(@{$resulttypes}, $type);
	}					
	
	return $resulttypes;
	
} # end resulttypes_appris

sub result_appris {
	my ($jobid) = shift;
	my ($params) = shift;
	my ($cutoffs);	
	my ($result);

	if ( defined $jobid and defined $params ) {
		
		# Create Retriever object from job id
		my ($retriever) = new WSRetriever( -jobid => $jobid );
		
		# get input methods, if we don't already know
		my ($format) = $params->{'format'};			
		my ($methods);
		if ( exists $params->{'methods'} ) {
			$methods = $params->{'methods'};
		}
		else {
			if ( $retriever->methods ) {
				$methods = $retriever->methods;
			}
			else {
				$methods = $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE};	
			}			
		}
		
		my ($ids); $ids = $params->{'ids'} if ( exists $params->{'ids'} );
		my ($res); $res = $params->{'res'} if ( exists $params->{'res'} );
		my ($headbed); $headbed = $params->{'headbed'} if ( exists $params->{'headbed'} );
		$result = $retriever->export_features($format, $methods, $ids, $res, $headbed);
		
		# Send email from exported result
		#my ($email);
		#if ( $params->{'email'} ) {
		#	$email = $params->{'email'};
		#}
		#elsif ( $retriever->email ) {
		#	$email = $retriever->email;
		#}
		#if ( defined $email ) {
		#	$retriever->send_features($email, $jobid, $result, $format);
		#} 
	}

	return $result;
		
} # end result_appris

sub _create_param_seq {
	my ($workspace, $id, $raw_seqs) = @_;
	my ($gtransl_file);	
	
	my ($gtransl_cont) = '';
	my ($transl_cont) = APPRIS::Parser::_parse_inseq_transl($raw_seqs) if (defined $raw_seqs);
	if ( defined $transl_cont ) {
		$gtransl_file = $workspace.'/'.$id.'.transl.fa';
		while ( my ($transcript_id, $transcript_features) = each(%{$transl_cont}) ) {
			if ( defined $transl_cont and exists $transl_cont->{$transcript_id} and defined $transl_cont->{$transcript_id} ) {
				my ($seq) = $transl_cont->{$transcript_id};
				my ($len) = length($seq);
				$gtransl_cont .= ">$transcript_id|$len\n";
				$gtransl_cont .= $seq."\n"
			}			
		}
		if ( $gtransl_cont ne '' ) {
			my ($l) = printStringIntoFile($gtransl_cont, $gtransl_file);
			throw("creating transl infile: $gtransl_file: $!\n") unless (defined $l);
		}
	}
	
	return $gtransl_file;
	
} # end _create_param_seq

sub _create_param_ensembl {
	my ($workspace, $id, $species, $e_version) = @_;
	my ($gdata_file, $gtransc_file, $gtransl_file);	
	
	# check if main data files exist
	my ($spe) = lc($species); $spe =~ s/\s/\_/g;
	my ($data_file) = $ENV{APPRIS_FEATURES_DIR}.'/'.$spe.'/'."ensembl$e_version".'/'.$spe.'.annot.gtf';
	my ($transc_file) = $ENV{APPRIS_FEATURES_DIR}.'/'.$spe.'/'."ensembl$e_version".'/'.$spe.'.transc.fa';
	my ($transl_file) = $ENV{APPRIS_FEATURES_DIR}.'/'.$spe.'/'."ensembl$e_version".'/'.$spe.'.transl.fa';
	if ( -e ($data_file) and (-s $data_file > 0) and
		 -e ($transc_file) and (-s $transc_file > 0) and
		 -e ($transl_file) and (-s $transl_file > 0) 
	) {
		$gdata_file = $workspace.'/'.$id.'.annot.gtf';
		$gtransc_file = $workspace.'/'.$id.'.transc.fa';
		$gtransl_file = $workspace.'/'.$id.'.transl.fa';
		if ( defined $id ) {
			# customize data
			my ($g_cond) = '';
			if ( $id =~ /^ENS/ ) {
				$g_cond .= ' $9=="gene_id" && $10 ~ /'.$id.'/ ';
			}
			else {
				$g_cond .= ' ($11=="gene_name" && $12 ~ /\"'.$id.'\"/) || ' .
							'($13=="gene_name" && $14 ~ /\"'.$id.'\"/) || '.
							'($15=="gene_name" && $16 ~ /\"'.$id.'\"/) || '.
							'($17=="gene_name" && $18 ~ /\"'.$id.'\"/) || '.
							'($19=="gene_name" && $20 ~ /\"'.$id.'\"/) ';
			}
			if ( $g_cond ne '' ) {					
				eval {
					my ($cmd) = "awk 'BEGIN{IGNORECASE=1} { if ( $g_cond ) {print \$0} }' $data_file >> $gdata_file";
					info("** script: $cmd\n");
					my (@cmd_out) = `$cmd`;
				};
				return undef if($@);			
			}
			# customize sequences
			my ($gtransc_cont) = '';
			my ($gtransl_cont) = '';			
			my ($data_cont) = APPRIS::Parser::_parse_indata($gdata_file);
			my ($transc_cont) = APPRIS::Parser::_parse_inseq_transc($transc_file) if (defined $transc_file);
			my ($transl_cont) = APPRIS::Parser::_parse_inseq_transl($transl_file) if (defined $transl_file);
			while ( my ($gene_id, $gene_features) = each(%{$data_cont}) ) {
				if ( exists $gene_features->{'transcripts'} ) {
					while (my ($transcript_id, $transcript_features) = each(%{$gene_features->{'transcripts'}}) ) {
						my ($translate_id) = $transcript_id;
						my ($transcript_name) = $transcript_id;
						if ( exists $transcript_features->{'protein_id'}) { $translate_id = $transcript_features->{'protein_id'} }
						if ( exists $transcript_features->{'external_id'}) { $transcript_name = $transcript_features->{'external_id'} }
						if ( defined $transc_cont and exists $transc_cont->{$transcript_id} and defined $transc_cont->{$transcript_id} ) {
							my ($seq) = $transc_cont->{$transcript_id};
							my ($len) = length($seq);
							$gtransc_cont  .= ">$transcript_id|$gene_id|$transcript_name|$len\n";
							$gtransc_cont .= $seq."\n"
						}			
						if ( defined $transl_cont and exists $transl_cont->{$transcript_id} and defined $transl_cont->{$transcript_id} ) {
							my ($seq) = $transl_cont->{$transcript_id};
							my ($len) = length($seq);
							$gtransl_cont .= ">$transcript_id|$translate_id|$gene_id|$transcript_name|$len\n";
							$gtransl_cont .= $seq."\n"
						}			
					}			
				}
			}
			if ( $gtransc_cont ne '' ) {
				my ($l) = printStringIntoFile($gtransc_cont, $gtransc_file);
				throw("creating transc infile: $gtransc_file: $!\n") unless (defined $l);
			}
			if ( $gtransl_cont ne '' ) {
				my ($l) = printStringIntoFile($gtransl_cont, $gtransl_file);
				throw("creating transl infile: $gtransl_file: $!\n") unless (defined $l);
			}
		}
	}
	
	return ($gdata_file, $gtransc_file, $gtransl_file);
		
} # end _create_param_ensembl

sub _create_param_species
{
	my ($param) = @_;
	my ($result);
	
	while ( my ($specie_id, $specie) = each(%{$SPECIES}) ) {
		if ( 	( lc($param) eq lc($specie_id) ) or
				( lc($param) eq lc($specie->{'scientific'}) ) or
				( lc($param) eq lc($specie->{'common'}) ) 
		) {
			$result = $specie->{'scientific'};
		}
	}	
	
	return $result;
		
} # end _create_param_species

# Check only runner methods
sub _create_param_methods
{
	my ($param) = @_;
	my ($result) = '';
	
	while ( my ($method_id, $method) = each(%{$METHODS}) ) {
		my ($name) = $method->{'name'}; 
		if ( $method->{'control'} =~ /runner/ ) { # only runner methods
			if ( $param =~ /$name/ ) {
				$result .= $name.',';
			}			
		}
	}
	$result =~ s/\,$//;
	
	return $result;
		
} # end _create_param_methods

1;