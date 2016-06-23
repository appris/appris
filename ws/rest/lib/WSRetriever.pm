=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

WSRetriever

=head1 DESCRIPTION

Package the retrieves the results of methods

=head1 SYNOPSIS

	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			conf			=>  undef,
			species			=>  undef,
		);

	#_____________________________________________________________
    
  my $features = WSRetriever->new(
					-species => 'mouse',
				    -conf   => 'conf/config.ini',
  );
  
  
=head1 METHODS

=cut

package WSRetriever;

use strict;
use warnings;
use FindBin;
use HTTP::Request;
use HTTP::Headers; 
use LWP::UserAgent;
use HTML::TreeBuilder;
use JSON;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;

use APPRIS::Utils::Argument qw( rearrange );
use APPRIS::Utils::Exception qw( throw warning deprecate );
use APPRIS::Utils::File qw( printStringIntoFile getStringFromFile updateStringIntoFile );
use APPRIS::Exporter;
use APPRIS::Parser qw(
	parse_gencode
	parse_firestar_rst
	parse_matador3d_rst
	parse_spade_rst
	parse_corsair_rst
	parse_thump_rst
	parse_crash_rst
	parse_appris_rst
	parse_appris_methods
	create_appris_entity
);
use APPRIS::Utils::Constant qw(
	$METHOD_DESC
);
my ($METHOD_DESC) = $APPRIS::Utils::Constant::METHOD_DESC;

use vars qw(
	$STATUS_LOGS
	$SERVER
	$METHODS
	$METHOD_IDS
	$SEQ_PANEL
);

$STATUS_LOGS = {
	'PENDING'	=> 'The job is waiting to be processed',
	'RUNNING'	=> 'The job is currently being processed',
	'FINISHED'	=> 'Job has finished, and the results can then be retrieved',
	'ERROR'		=> 'An error occurred attempting to get the job status',
	'FAILURE'	=> 'The job failed',
	'NOT_FOUND'	=> 'The job cannot be found'
};

my ($server_json) = JSON->new(); $SERVER = $server_json->decode( getStringFromFile($ENV{APPRIS_WSERVER_CONF_SR_FILE}) );
my ($methods_json) = JSON->new(); $METHODS = $methods_json->decode( getStringFromFile($ENV{APPRIS_WSERVER_CONF_ME_FILE}) );
while ( my ($met_id,$met_report) = each(%{$METHODS}) ) {
	my ($met_name) = lc($met_report->{'name'});
	$METHOD_IDS->{$met_name} = $met_id;
}

{
	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			status			=>  undef,
			tracelog		=>  undef,
			
			jobid			=>  undef,
			species			=>  undef,
			dataset			=>  undef,
			assembly		=>  undef,
			source			=>  undef,
			
			data			=>  undef,
			pdata			=>  undef,
			transc			=>  undef,
			transl			=>  undef,

			methods			=>  undef,
			firestar		=>  undef,
			matador3d		=>  undef,
			spade			=>  undef,
			corsair			=>  undef,
			thump			=>  undef,
			crash			=>  undef,
			inertia			=>  undef,
			proteo			=>  undef,
			appris			=>  undef,
			appris_lb		=>  undef,
			logfile			=>  undef,
			
			host			=>  undef,
			pid				=>  undef,
			wsbase			=>  undef,
			
			email			=>  undef,
		);

	#_____________________________________________________________

	# Classwide default value for a specified object attribute
	sub _default_for {
		my ($self, $attr) = @_;
		$_attr_data{$attr};
	}

	# List of names of all specified object attributes
	sub _standard_keys {
		keys %_attr_data;
	}
	
	# Substitute values within template
	sub status {
		my ($self, $arg) = @_;
		$self->{'status'} = $arg if defined $arg;
		return $self->{'status'};
	}
	sub tracelog {
		my ($self, $arg) = @_;
		$self->{'tracelog'} = $arg if defined $arg;
		return $self->{'tracelog'};
	}

	sub jobid {
		my ($self, $arg) = @_;
		$self->{'jobid'} = $arg if defined $arg;
		return $self->{'jobid'};
	}
	sub species {
		my ($self, $arg) = @_;
		$self->{'species'} = $arg if defined $arg;
		return $self->{'species'};
	}
	sub dataset {
		my ($self, $arg) = @_;
		$self->{'dataset'} = $arg if defined $arg;
		return $self->{'dataset'};
	}
	sub assembly {
		my ($self, $arg) = @_;
		$self->{'assembly'} = $arg if defined $arg;
		return $self->{'assembly'};
	}	
	sub source {
		my ($self, $arg) = @_;
		$self->{'source'} = $arg if defined $arg;
		return $self->{'source'};
	}
	
	sub transl {
		my ($self, $arg) = @_;
		$self->{'transl'} = $arg if defined $arg;
		return $self->{'transl'};
	}
	sub data {
		my ($self, $arg) = @_;
		$self->{'data'} = $arg if defined $arg;
		return $self->{'data'};
	}
	sub pdata {
		my ($self, $arg) = @_;
		$self->{'pdata'} = $arg if defined $arg;
		return $self->{'pdata'};
	}
	sub transc {
		my ($self, $arg) = @_;
		$self->{'transc'} = $arg if defined $arg;
		return $self->{'transc'};
	}

	sub methods {
		my ($self, $arg) = @_;
		$self->{'methods'} = $arg if defined $arg;
		return $self->{'methods'};
	}
	sub firestar {
		my ($self, $arg) = @_;
		$self->{'firestar'} = $arg if defined $arg;
		return $self->{'firestar'};
	}
	sub matador3d {
		my ($self, $arg) = @_;
		$self->{'matador3d'} = $arg if defined $arg;
		return $self->{'matador3d'};
	}
	sub spade {
		my ($self, $arg) = @_;
		$self->{'spade'} = $arg if defined $arg;
		return $self->{'spade'};
	}
	sub corsair {
		my ($self, $arg) = @_;
		$self->{'corsair'} = $arg if defined $arg;
		return $self->{'corsair'};
	}	
	sub thump {
		my ($self, $arg) = @_;
		$self->{'thump'} = $arg if defined $arg;
		return $self->{'thump'};
	}
	sub crash {
		my ($self, $arg) = @_;
		$self->{'crash'} = $arg if defined $arg;
		return $self->{'crash'};
	}
	sub inertia {
		my ($self, $arg) = @_;
		$self->{'inertia'} = $arg if defined $arg;
		return $self->{'inertia'};
	}
	sub proteo {
		my ($self, $arg) = @_;
		$self->{'proteo'} = $arg if defined $arg;
		return $self->{'proteo'};
	}
	sub appris {
		my ($self, $arg) = @_;
		$self->{'appris'} = $arg if defined $arg;
		return $self->{'appris'};
	}
	sub appris_lb {
		my ($self, $arg) = @_;
		$self->{'appris_lb'} = $arg if defined $arg;
		return $self->{'appris_lb'};
	}
	sub logfile {
		my ($self, $arg) = @_;
		$self->{'logfile'} = $arg if defined $arg;
		return $self->{'logfile'};
	}

	sub host {
		my ($self, $arg) = @_;
		$self->{'host'} = $arg if defined $arg;
		return $self->{'host'};
	}
	sub pid {
		my ($self, $arg) = @_;
		$self->{'pid'} = $arg if defined $arg;
		return $self->{'pid'};
	}
	sub wsbase {
		my ($self, $arg) = @_;
		$self->{'wsbase'} = $arg if defined $arg;
		return $self->{'wsbase'};
	}

	sub email {
		my ($self, $arg) = @_;
		$self->{'email'} = $arg if defined $arg;
		return $self->{'email'};
	}

}

=head2 new

  Arg [-transl]: (required)
        string - protein sequences
  Arg [-data]: (optional)
        string - gene annotations
  Arg [-pdata]: (optional)
        string - protein annotations
  Arg [-transc]: (optional)
        string - transcript sequences
  Example    : $features = WSRetriever->new(...);
  Description: Creates a new retrieve object
  Returntype : WSRetriever
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
	my ($caller) = shift;
	
	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;

	foreach my $attrname ($self->_standard_keys) {
		$self->{$attrname} = $self->_default_for($attrname);
	}
	
	my (
		$jobid,
	) = rearrange( [
		'jobid',
	], @_ );	
	
	# optional paramater	
	if ( defined $jobid ) {
		$self->jobid($jobid);
		$self->load_control($jobid);
		#return undef unless ( defined $self->transl ); # WARNING!!! You need this for 'sequencer'.pl REST
	}
		
	return $self;
}

=head2 load_control

  Example    : use WSRetriever qw(load_entity);
               my $rst = load_control();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut
sub load_control {
	my ($self) = shift;
	my ($jobid) = shift;
	
	# retrieve runner variables: outfiles, methods
	my ($workspace) = $ENV{APPRIS_WORKSPACE} . '/' . $jobid;
	my ($wsrunnerctrfile) = $workspace . '/wsrunner.ctr';
	my (%outfiles);
	if ( -e $wsrunnerctrfile and (-s $wsrunnerctrfile > 0) ) {
		my ($wsrunnerctrcont) = getStringFromFile($wsrunnerctrfile);
		if ( $wsrunnerctrcont =~ /RUNNER_OUTFILES\:([^\n]*)/ ) {
			%outfiles =  map {$_ => 1 } split(',', $1);
		}
		# save the runner method by order!! By default all
		if ( $wsrunnerctrcont =~ /--methods=([^\s]*)/ ) {
			my ($aux_run_methods) = $1;
			my ($run_methods) = '';
			foreach my $method ( split(',', $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE}) ) {
				$run_methods .= $method.',' if ( $aux_run_methods =~ /$method/ );
			}
			if ( $run_methods ne '' ) { $run_methods =~ s/\,$// }
			else { $run_methods = $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE} }
			$self->methods($run_methods);
		}
		if ( $wsrunnerctrcont =~ /--species='([^\']*)'/ ) {
			my ($species) = $1;
			$species =~ s/\'//;
			$self->species($species);					
		}
		if ( $wsrunnerctrcont =~ /--e-version=([^\s]*)/ ) {
			$self->dataset($1);					
		}
		if ( $wsrunnerctrcont =~ /RUNNER_HOST\:([^\n]*)\n*RUNNER_PID\:([^\n]*)\n*RUNNER_JOBID\:([^\n]*)\n*RUNNER_WSPACE\:([^\n]*)\n/ ) {
			$self->host($1);
			$self->pid($2);
			$self->wsbase($4);
		}		
		if ( $wsrunnerctrcont =~ /RUNNER_STATUS\:([^\n]*)/ ) {
			$self->status($1);
		}
		if ( $wsrunnerctrcont =~ /RUNNER_LOG\:([^\n]*)/ ) {
			$self->tracelog($1);
		}
		if ( $wsrunnerctrcont =~ /RUNNER_EMAIL\:([^\n]*)/ ) {
			$self->email($1);
		}
	}
	
	# get assembly and source
	if ( $self->species and $self->dataset ) {
		my ($species) = $self->species;
		my ($dataset) = $self->dataset;		
		my ($species_id) = lc($species); $species_id =~ s/\s/\_/g;
		my ($cfg_species) = $SERVER->{'species'}->{$species_id};
		if ( defined $cfg_species and $cfg_species->{'assemblies'} ) {
			foreach my $cfg_assembly (@{$cfg_species->{'assemblies'}}) {
				my ($assembly_id) = $cfg_assembly->{'id'};
				foreach my $cfg_dataset (@{$cfg_assembly->{'datasets'}}) {
					if ( $dataset  eq $cfg_dataset->{'source'}->{'version'} ) {
						$self->assembly($cfg_assembly->{'id'});
						$self->source($cfg_dataset->{'source'}->{'name'});
					}	
				}
			}
        }		
	}
	
	# get input files
	my ($data_file) = grep(/\/annot\.gtf$/, keys %outfiles);	
	my ($pdata_file) = grep(/\/pannot\.gtf$/, keys %outfiles);
	my ($transc_file) = grep(/\/transc\.fa$/, keys %outfiles);
	my ($transl_file) = grep(/\/transl\.fa$/, keys %outfiles);
	my ($log_file) = grep(/\/log$/, keys %outfiles);
	$self->data($data_file) if ( defined $data_file and -e $data_file and (-s $data_file > 0) );	
	$self->pdata($pdata_file) if ( defined $pdata_file and -e $pdata_file and (-s $pdata_file > 0) );	
	$self->transc($transc_file) if ( defined $transc_file and -e $transc_file and (-s $transc_file > 0) );	
	$self->transl($transl_file) if ( defined $transl_file and -e $transl_file and (-s $transl_file > 0) );
	$self->logfile($log_file) if ( defined $log_file and -e $log_file and (-s $log_file > 0) );	

	# assign file results OF ALL METHODS!!
	foreach my $method ( split(',', $ENV{APPRIS_WSERVER_PIPELINE_STRUCTURE}) ) {
		if ( $method eq 'firestar' ) {
			my ($metfile) = grep(/\/$method$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) ) {
				$self->firestar(getStringFromFile($metfile));
			}
		}
		if ( $method eq 'matador3d' ) {
			my ($metfile) = grep(/$method$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) ) {
				$self->matador3d(getStringFromFile($metfile));
			}
		}
		if ( $method eq 'corsair' ) {
			my ($metfile) = grep(/$method$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) ) {
				$self->corsair(getStringFromFile($metfile));
			}			
		}
		if ( $method eq 'spade' ) {
			my ($metfile) = grep(/$method$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) ) {
				$self->spade(getStringFromFile($metfile));
			}			
		}
		if ( $method eq 'thump' ) {
			my ($metfile) = grep(/$method$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) ) {
				$self->thump(getStringFromFile($metfile));
			}
		}
		if ( $method eq 'crash' ) {
			my ($metfile) = grep(/$method$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) ) {
				$self->crash(getStringFromFile($metfile));
			}
		}
		if ( $method eq 'appris' ) {
			my ($metfile) = grep(/\/$method$/, keys %outfiles);
			my ($metfile2) = grep(/\/$method\.label$/, keys %outfiles);
			if ( defined $metfile and -e $metfile and (-s $metfile > 0) and defined $metfile2 and -e $metfile2 and (-s $metfile2 > 0) ) {
				$self->appris(getStringFromFile($metfile));
				$self->appris_lb(getStringFromFile($metfile2));
			}
		}
	}		
		
} # end load_control

=head2 add_control

  Example    : use WSRetriever qw(load_entity);
               my $rst = add_control();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut
sub add_control {
	my ($self) = shift;
	my ($c_type) = shift if (@_);
	my ($c_value) = shift if (@_);
		
	my ($workspace) = $ENV{APPRIS_WORKSPACE} . '/' . $self->jobid;
	my ($wsrunnerctrfile) = $workspace.'/wsrunner.ctr';
	if ( defined $c_type and defined $c_value ) {
		my ($ctrcont) = $c_type . ':' . $c_value . "\n";
		my ($file) = updateStringIntoFile($ctrcont, $wsrunnerctrfile);
	}
			
} # end add_control

=head2 status_log

  Example    : use WSRetriever qw(load_entity);
               my $rst = status_log();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut
sub status_log {
	my ($self) = shift;
	my ($c_status) = shift if (@_);
	my ($c_tracelog) = shift if (@_);
	my ($status, $log, $progress) = ('NOT_FOUND', $STATUS_LOGS->{'NOT_FOUND'}, 0);
	
	my ($workspace) = $ENV{APPRIS_WORKSPACE} . '/' . $self->jobid;
	my ($wsrunnerctrfile) = $workspace.'/wsrunner.ctr';
	
	# logs before started to run
	if ( $self->status and $self->tracelog ) {
		($status, $log, $progress) = ($self->status, $self->tracelog, 0);
	}
	else { # logs after started to run
		if ( $self->methods ) {
			my ($jobstatus,$joblog);
			if ( defined $c_status and defined $c_tracelog ) {
				($status, $log, $progress) = $self->parse_logfile($c_tracelog, $self->methods, $c_status);						
			}
			else {
				if ( $self->logfile and -e $self->logfile and (-s $self->logfile > 0) ) {
					my ($tracelog) = getStringFromFile($self->logfile);
					($status, $log, $progress) = $self->parse_logfile($tracelog, $self->methods);
				}				
			}
		}
		#else {
		#	($status, $log, $progress) = ('PENDING', $STATUS_LOGS->{'PENDING'}, 0); # at least, is pending
		#}
	}
			
	return ($status, $log, $progress);
		
} # end status_log

=head2 get_reports

  Example    : use WSRetriever qw(get_reports);
               my $rst = get_reports($reports);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_reports
{
	my ($self) = shift;
	my ($ids) = shift if (@_);
	my ($analyses);
	
	if (
	(
		defined $self->data and -e $self->data and (-s $self->data > 0) and
		defined $self->transc and -e $self->transc and (-s $self->transc > 0) and
		defined $self->transl and -e $self->transl and (-s $self->transl > 0)
	)
	or
	(
		defined $self->transl and -e $self->transl and (-s $self->transl > 0)
	) 
	) {	
		$analyses = APPRIS::Parser::create_appris_entity(
										$self->data,
										$self->transc,
										$self->transl,
										$self->firestar,
										$self->matador3d,
										$self->spade,
										$self->corsair,
										$self->crash,
										$self->thump,
										$self->inertia,
										$self->proteo,
										$self->appris,
										$self->appris_lb,
										$ids
		);
	}
	else {
		throw("entity files are not correct");		
	}	
	return $analyses;
		
} # end get_reports

=head2 get_results

  Example    : use WSRetriever qw(get_results);
               my $rst = get_results();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_results
{
	my ($self) = shift;
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);	
	my ($analyses);
	
	my $grep_ids = sub {
		my ($cutoff, $ids) = @_;
		my ($analyses);	
		if ( defined $cutoff ) {
			if ( defined $ids ) {
				my ($cutoff_ids);
				while ( my ($id,$rep) = each(%{$cutoff}) ) {
					if ( $id eq 'result' ) { $cutoff_ids->{'result'} = $rep; }
					elsif ( $ids =~ /$id/ ) { $cutoff_ids->{$id} = $rep; }					
				}
				$analyses = $cutoff_ids if ( defined $cutoff_ids );
			}
			else {
				$analyses = $cutoff;
			}
		}
		return $analyses;
	};
	
	if ( $self->firestar ) {
		my ($cutoff) = APPRIS::Parser::parse_firestar_rst($self->firestar);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'firestar'} = $cutoff_ids;			
		}
	}	
	if ( $self->matador3d ) {
		my ($cutoff) = APPRIS::Parser::parse_matador3d_rst($self->matador3d);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'matador3d'} = $cutoff_ids;			
		}
	}
	if ( $self->corsair ) {
		my ($cutoff) = APPRIS::Parser::parse_corsair_rst($self->corsair);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'corsair'} = $cutoff_ids;			
		}
	}
	if ( $self->spade ) {
		my ($cutoff) = APPRIS::Parser::parse_spade_rst($self->spade);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'spade'} = $cutoff_ids;			
		}
	}
	if ( $self->thump ) {
		my ($cutoff) = APPRIS::Parser::parse_thump_rst($self->thump);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'thump'} = $cutoff_ids;			
		}
	}
	if ( $self->crash ) {
		my ($cutoff) = APPRIS::Parser::parse_crash_rst($self->crash);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'crash'} = $cutoff_ids;			
		}
	}
	if ( $self->appris and $self->appris_lb ) {
		my ($cutoff) = APPRIS::Parser::parse_appris_rst($self->appris, $self->appris_lb);
		my ($cutoff_ids) = $grep_ids->($cutoff, $ids);
		if ( defined $cutoff_ids ) {
			$analyses->{'appris'} = $cutoff_ids;			
		}
	}

	return $analyses;	
	
} # end get_results

=head2 get_features

  Arg [1]    : Hash $report
  Example    : use WSRetriever qw(get_features);
               my $rst = get_features($reports);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_features
{
	my ($self) = shift;
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($features);
	
	# load files from jobid
	$self->load_control($self->jobid);
	
	if ( defined $self->data and -e $self->data and (-s $self->data > 0) and
 		 defined $self->transc and -e $self->transc and (-s $self->transc > 0) and
		 defined $self->transl and -e $self->transl and (-s $self->transl > 0)
	) {		
		$features = $self->get_reports($ids);
	}
	elsif ( defined $self->transl and -e $self->transl and (-s $self->transl > 0) ) {
		$features = $self->get_reports($ids);
	}	
	return $features;
	
} # end get_features

=head2 export_features

  Arg [1]    : String $format
  Arg [3]    : String $methods
  Example    : use WSRetriever qw(export_features);
               my $rst = export_features($format);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut 

sub export_features
{
	my ($self) = shift;
	my ($format) = shift;
	my ($methods) = shift;
	my ($ids) = shift if (@_);
	my ($res) = shift if (@_);
	my ($exporter) = APPRIS::Exporter->new();
	my ($result) = '';	
	my ($features) = $self->get_features($methods, $ids); # methods is not defined => retrieve all features
	
	if ( defined $features ) {
		if (
			defined $self->data and -e $self->data and (-s $self->data > 0) and
			defined $self->transc and -e $self->transc and (-s $self->transc > 0) and
			defined $self->transl and -e $self->transl and (-s $self->transl > 0) 
		) {
			if ($format eq 'tsv') {
				$result = $exporter->get_tsv_annotations($features, $methods, $res);
			}
			elsif ( ($format eq 'bed') or ($format eq 'bed12') ) {
				$result = $exporter->get_bed_annotations($features, $methods, undef, $format);
		    }
			elsif ($format eq 'json') {
				$result = $exporter->get_json_annotations($features, $methods);
		    }
			elsif ($format eq 'gtf') {
				$result = $exporter->get_gtf_annotations($features, $methods);
		    }
			elsif ($format eq 'raw') {
				$result = $exporter->get_raw_annotations($features, $methods);
		    }
		}
		elsif ( defined $self->transl and -e $self->transl and (-s $self->transl > 0) ) {
			if ($format eq 'tsv') {
				$result = $exporter->get_tsv_annotations($features, $methods, $res);
			}
			elsif ($format eq 'json') {
				$features = undef; $features = $self->get_results($methods, $ids); # IMPORTANT, OVERWRTITE FEATURES TO PRINT RESULTS
				$result = $exporter->get_json_results($features, $methods);
		    }
			elsif ($format eq 'raw') {
				$result = $exporter->get_raw_annotations($features, $methods);
		    }
			else {
				$result = "The format is not available for this job";
				warning($result);
			}			    
		}
	}	
	return $result;
	
} # end export_features

##
## SEQUENCER
##

=head2 export_seq_features

  Arg [1]    : String $type
               type of input: 'id' or 'name'
  Arg [2]    : String $input
               gene or transcript identifier (or name)
  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(export_seq_features);
               my $rst = export_seq_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub export_seq_features
{
	my ($self) = shift;
	my ($operation) = shift;
	my ($format) = shift;
	my ($type) = shift if (@_);
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($res) = shift if (@_);
	my ($result) = '';	
	my ($features) = $self->get_features($methods, $ids); # if methods is not defined => retrieve all features
		
	if ( defined $features ) {
		my ($exporter) = APPRIS::Exporter->new();
		if ( $operation eq 'sequences' ) {
			$result = $exporter->get_seq_annotations($self->source, $features, $type, $format);
		}
		elsif ( $operation eq 'align' ) { # "A minimum of 2 sequences is required\n"
			my ($seq_result) = $exporter->get_seq_annotations($self->source, $features, $type, 'fasta');
			my ($num_seq) = $seq_result =~ tr/\>//;
			if ( $num_seq >= 2 ) {
				my ($wsretriever) = new WSRetriever();
				$result = $wsretriever->get_aln_annotations($seq_result, $format, $ids);
			}
			else {
				$result = $exporter->get_seq_annotations($self->source, $features, $type, 'json');				
			}			
		}		
		elsif ($operation eq 'residues') {
			$result = $exporter->get_res_annotations($self->source, $features, $methods, $res);
	    }
		elsif ($operation eq 'cds') {
			$result = $exporter->get_cds_annotations($self->source, $features, $methods, $res);
	    }
		elsif ($operation eq 'genome') {
			$result = $self->get_gen_annotations($self->source, $methods, $ids);
	    }
	}
	
	return $result;
	
} # end export_seq_features

=head2 get_aln_annotations

  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(export_seq_features);
               my $rst = get_aln_annotations($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_aln_annotations
{
	my ($self) = shift;
	my ($seq_result) = shift;
	my ($format) = shift;
	my ($ids) = shift if (@_);
	my ($string) = '';
	my ($result) = '';
	
	# retrieve align from REST service of Muscle(EBI)
	if ( defined $seq_result and ($seq_result ne '') ) {
		$string = $self->get_tool_result($seq_result, 'muscle', 'clwstrict'); # overwrite
	}
		
	if ( defined $format and ($format eq 'json') and ($string ne '') ) {
		my ($alns);
		my ($talns);
		my ($stringfh) = IO::String->new($string);
		my ($alnio) = Bio::AlignIO->new(-fh => $stringfh, -format => 'clustalw');
		my ($aln) = $alnio->next_aln();
		my ($i) = 0;
		for my $seq ( $aln->each_seq ) {
			my ($seq_id) = $seq->id; $seq_id =~ s/^([^\|]*)\|[^\$]*$/$1/g;				
			my ($seq_seq) = $seq->seq;
			push(@{$alns}, {
				'id'	=> $seq_id,
				'seq'	=> $seq->seq
			});
			$talns->{$seq_id} = $i;
			$i++;
		}
		my ($report);
		if ( defined $ids ) {
			foreach my $t ( split(',', $ids) ) {
				if ( defined $t and ($t ne '') ) {
					if ( exists $talns->{$t} ) {
						my ($i) = $talns->{$t};
						push(@{$report->{'aln'}}, $alns->[$i]);							
					}						
				}
			}
		}
		else { $report->{'aln'} = $alns }
		
		$report->{'match'} = $aln->match_line();
		$report->{'length'} = $aln->length;
		$report->{'num'} = scalar(@{$report->{'aln'}});	
		require JSON;
		my ($json) = new JSON;
		$result = $json->encode($report);
	}
	elsif ( defined $format and ($format eq 'clw') and ($string ne '') ) {
		$result = $string;
	}
	return $result;
	
} # end get_aln_annotations

=head2 get_gen_annotations

  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(export_seq_features);
               my $rst = export_seq_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub get_gen_annotations
{
	my ($self) = shift;
	my ($source) = shift if (@_);
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($result) = '';
	
	# retrieve images from UCSC
	my ($query_id) = 'runner/result/' . $self->jobid;		
	$result = $self->get_gen_features($query_id, $self->species, $self->assembly, $source, $self->dataset, $methods, $ids);		
	unless ( defined $result ) {
		$result = "Job query needs genome information\n";
	}
	
	return $result;
	
} # end get_gen_annotations




#=head2 export_seq_aln_features
#
#  Arg [3]    : String $method
#               method name
#  Example    : use DBRetriever qw(export_seq_features);
#               my $rst = export_seq_features($type, $input, $method);
#  Description: Throws an exception which if not caught by an eval will
#               provide a stack trace to STDERR and die.  If the verbosity level
#               is lower than the level of the throw, then no error message is
#               displayed but the program will still die (unless the exception
#               is caught).
#  Returntype : none
#  Exceptions : thrown every time
#  Caller     : generally on error
#
#=cut

#sub export_seq_aln_features
#{
#	my ($self) = shift;
#	my ($operation) = shift;	
#	my ($format) = shift;
#	my ($methods) = shift if (@_);
#	my ($ids) = shift if (@_);
#	my ($result) = '';
#
#	# get protein fasta
#	my ($seq_result) = $self->export_seq_features('sequences', 'fasta', 'aa', $methods, $ids);
#	my ($num_seq) = $seq_result =~ tr/\>//;
#		
#	# get annotated residues
#	my ($res_result) = $self->export_seq_features('residues', 'fasta', 'aa', $methods, $ids);
#	
#	# retrieve align from REST service of Muscle(EBI)
#	my ($wsretriever) = new WSRetriever();
#	if ( $operation eq 'align' ) {
#		if ( $num_seq <= 1 ) { # "A minimum of 2 sequences is required\n"
#			$operation = 'sequences'; # sequence operator
#		}
#		else { # overwrite seq report by align report
#			if ( defined $seq_result ) {
#				$seq_result = $wsretriever->get_tool_result($seq_result, 'muscle', 'clwstrict');		
#			}			
#		}	
#	}
#		
#	# get sequence/align annotations
#	if ( defined $res_result and defined $seq_result ) {
#		$result = $wsretriever->get_seq_aln_features($operation, $seq_result, $res_result, $ids, $format);
#	}		
#		
#	return $result;
#	
#} # end export_seq_aln_features

#=head2 export_gen_features
#
#  Arg [3]    : String $method
#               method name
#  Example    : use DBRetriever qw(export_seq_features);
#               my $rst = export_seq_features($type, $input, $method);
#  Description: Throws an exception which if not caught by an eval will
#               provide a stack trace to STDERR and die.  If the verbosity level
#               is lower than the level of the throw, then no error message is
#               displayed but the program will still die (unless the exception
#               is caught).
#  Returntype : none
#  Exceptions : thrown every time
#  Caller     : generally on error
#
#=cut
#
#sub export_gen_features
#{
#	my ($self) = shift;
#	my ($methods) = shift if (@_);
#	my ($ids) = shift if (@_);
#	my ($result) = '';
#	
#	# retrieve images from UCSC
#	if ( $self->species and $self->ens ) {
#		my ($query_id) = 'runner/result/' . $self->jobid;		
#		$result = $self->get_gen_features($query_id, $self->species, $self->ens, $methods, $ids);		
#	}
#	else {
#		$result = "Job query needs genome information\n";		
#	}
#	
#	return $result;
#	
#} # end export_gen_features





sub parse_logfile
{
	my ($self) = shift;
	my ($tracelog) = shift;
	my ($run_methods) = shift;
	my ($clusterstatus) = shift if (@_);
	my ($status, $log, $progress) = (undef, '', 0);
	
	if ( defined $clusterstatus ) {
		if ( ($clusterstatus eq 'FINISHED') and ($tracelog eq '') ) {
			$status = 'ERROR';
			$log = 'An error occurred attempting to get the job status';
		}
		elsif ( $clusterstatus eq 'PENDING' ) {
			$status = 'PENDING';
			$log = 'The job is waiting to be processed';
		}
	}
	unless ( defined $status ) {
		# create log report
		my ($met_log);
		my (@methods) = split(',', $run_methods);
		my ($num_run_methods) = scalar(@methods); # deleting appris
		for ( my $i = 0; $i < scalar(@methods); $i++ ) {
			my ($met) = $methods[$i];
			my ($met_pl) = "$met\.pl"; if ( $met eq 'appris' ) { $met_pl = "$met/$met\.pl"; }
			my ($finished);
			if ( index($tracelog, "All done for $met") != -1 ) {
				push(@{$met_log}, {
					'method'	=> $met,
					'status'	=> 'finished'
				});
				$finished = 1;
			}
			elsif ( index($tracelog, "Prematurely exit for $met") != -1 ) {
				push(@{$met_log}, {
					'method'	=> $met,
					'status'	=> 'error'
				});						
				$finished = 1;
			}
			elsif ( index($tracelog, $met_pl ) != -1 ) {
				push(@{$met_log}, {
					'method'	=> $met,
					'status'	=> 'running'
				});
			}
			else {
				push(@{$met_log}, {
					'method'	=> $met,
					'status'	=> 'pending'
				});
			}
			
			# exception when there is an error that has not been captched.
			if ( defined $finished and ($i >= 1) and ($met ne 'appris') ) {
				my ($lastmet) = $methods[$i-1];
				my ($lastmet_log) = $met_log->[$i-1];
				if ( ($lastmet_log->{'status'} eq 'running') or ($lastmet_log->{'status'} eq 'pending') ) {
					$lastmet_log->{'status'} = 'error';
				}
			}
			
		}
				
		# create the final values of status-log
		my ($failed) = 0;
		my ($finish) = 1;
		my ($run_mets) = 0;
		foreach my $l_rep (@{$met_log}) {
			if ( defined $l_rep->{'method'} and defined $l_rep->{'status'} ) {
				my ($m) = $l_rep->{'method'};
				my ($m_id) = $METHOD_IDS->{$m};
				my ($m_label) = $METHODS->{$m_id}->{'label'};					
				if ( $l_rep->{'status'} eq "running" ) {
					$log .= "'". $m_label . "'" . " is calculating...\n";
					#$log .= $l_rep->{'method'}." is running...\n";
					$finish = 0;
					$run_mets += 0.3;
				}
				elsif ( $l_rep->{'status'} eq "finished" ) {
					$log .= "'". $m_label . "'" . " was obtained successfully\n";			
					#$log .= $l_rep->{'method'}." was successfully finished\n";
					$run_mets += 0.7;		
				}
				elsif ( $l_rep->{'status'} eq "error" ) {
					# Hard-Core with MATADOR3D!!!!!
					if ( $m eq 'matador3d' ){
						$log .= "'". $m_label . "'" . " was acquired\n"
					}
					else {
						$log .= "'". $m_label . "'" . " was acquired with errors\n"
					}					
					# Hard-Core with MATADOR3D!!!!!
					#if ( $l_rep->{'method'} eq 'matador3d' ){
					#	$log .= $l_rep->{'method'}." finished\n";
					#}
					#else {
					#	$log .= $l_rep->{'method'}." finished with errors\n";	
					#}					
					$failed = 1;
					$run_mets += 0.7;
				}
			}
		}
		if ( ($finish == 1) and ($log ne '') ) {
			if ( $failed == 1 ) {
				$status = 'FAILURE';
				$progress = 1;
			}
			else {
				$status = 'FINISHED';
				$progress = 1;
			}
		}
		else {
			$status = 'RUNNING';
			$progress = ($run_mets/$num_run_methods);
		}
	}
	
	return ($status, $log, $progress);
	
} # end parse_logfile

sub get_tool_result {
	my ($self) = shift;
	my ($sequences, $type, $format) = @_;
	my ($result) = undef;
	
	if ( defined $sequences and ($sequences ne '') ) {
		# create tmp links
		my ($intmpfile) = File::Temp->new( UNLINK => 0, SUFFIX => ".appris_tool.$type.in" );
		my ($outtmpfile) = File::Temp->new( UNLINK => 0, SUFFIX => ".appris_tool.$type.out" );
		print $intmpfile $sequences;	
		
		# create command
		my ($cmd) = "perl $FindBin::Bin/";
		my ($suffix) = '';		
		if ( $type eq 'muscle' ) {
			$cmd .= "ebi/muscle_lwp.pl --email appris\@cnio.es --format $format $intmpfile --outfile $outtmpfile 1>&2 2> /dev/null ";
			if ( $format eq 'clw' ) { $suffix = "aln-clustalw.clw" }
			if ( $format eq 'clwstrict' ) { $suffix = "aln-clustalw_strict.clwstrict" }
			elsif ( $format eq 'fasta' ) { $suffix = "aln-fasta.fasta" }
			elsif ( $format eq 'html' ) { $suffix = "aln-html.html" }			
		}
		elsif ( $type eq 'blastp' ) {
			$cmd .= "ncbi/web_blast.pl $type nr $intmpfile 1>&2 2> /dev/null ";
		}
		else {
			$result = undef;
		}

		# execute command
		eval {
			system($cmd);
		};
		return undef if($@);
		
		# retrieve output
		my ($outfile) = $outtmpfile.'.'.$suffix;
		my ($outcontent) = getStringFromFile($outfile);
		if ( defined $outcontent ) {
			$result = $outcontent;
			$result =~ s/STYLE\=\"[^\"]*\"//mg;
			$result =~ s/BGCOLOR\=\"[^\"]*\"//mg;
			$result =~ s/CLUSTAL multiple sequence alignment by MUSCLE \(3\.8\)\n*//mg;			
		};
		
		# destroy tmp links
		$intmpfile->unlink_on_destroy(1);
		eval {
			my ($cmd) = "rm -rf $outtmpfile*";
			system($cmd);
		};
		return undef if($@);
	}
	
	return $result;
	
} # end get_tool_result

sub get_seq_features {
	my ($self) = shift;
	my ($seq_string) = shift;
	my ($residues) = shift;	
	my ($max_peptides) = shift; $max_peptides = 90 unless ( defined $max_peptides );
	my ($panel, $result) = (undef,undef);
	my ($seq_panel);
	my ($seq_report);
	my ($seq_gaps);		
	
	my ($stringfh) = IO::String->new($seq_string);
	my ($in) = Bio::SeqIO-> new(
							-fh     => $stringfh,
							-format => 'Fasta'
	);		
	while ( my $seq = $in->next_seq() ) {
		my ($seq_id) = $seq->id; $seq_id =~ s/^([^\|]*)\|[^\$]*$/$1/g;				
		my ($seq_seq) = $seq->seq;
		local $seq_panel->{$seq_id} = $SEQ_PANEL;
		push(@{$seq_report}, {
			'id'	=> $seq_id,
			'seq'	=> $seq->seq
		});
		$seq_gaps->{$seq_id} = 0;
	}
		
	my ($num_title) = '';
	my ($num_pep) = 1;
	for ( my $i = $num_pep; $i < $max_peptides; $i++) {
		if ( ($i % 10) == 0 ) { $num_pep .= '|' }
		else { $num_pep .= '.' }
	}
	$num_pep .= "$max_peptides";
	my ($seq_result) = 	$num_title . $num_pep . "\n\n";
		
	foreach my $seq_rep ( @{$seq_report} ) {
		my ($seq_id) =  $seq_rep->{'id'};
		my ($seq_seq) =  $seq_rep->{'seq'};
		my ($s_idx) = 0;
		my ($seq_length) = length($seq_seq);
		my ($seq_pep) = '';
		
		# print features with residues annotations
		do {
			my ($s_seq) = substr( $seq_seq, $s_idx, $max_peptides);				
	    	my ($pep_idx) = $s_idx - $seq_gaps->{$seq_id} + 1;
	    	my ($pep_len) = 0;			    	
			my ($pep_res) = '';
			for ( my $i = 0; $i <= length($s_seq); $i++ ) { 
			    my ($pep) = substr( $s_seq, $i, 1);
			    if ( $pep ne '-' ) {			    	
				    if ( exists $residues->{$seq_id} and exists $residues->{$seq_id}->{$pep_idx} ) {
				    	my ($r_clas) = '';
				    	my ($pep_met) = '';
				    	foreach my $residue (@{$residues->{$seq_id}->{$pep_idx}}) {
				    		my ($met) = $residue->{'method'};				    	
				    		if ( $met eq $METHOD_DESC->{'firestar'} ) {
				    			my ($idx) = 0;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
				    			$r_clas .= " $m ";
				    			$pep_met .= $met.',';
				    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
				    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
				    			}				    			
				    		}
				    		elsif ( $met eq $METHOD_DESC->{'matador3d'} ) {
				    			my ($idx) = 1;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
				    			$r_clas .= " $m ";
				    			$pep_met .= $met.',';
				    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
				    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
				    			}				    			
				    		}
				    		elsif ( $met eq $METHOD_DESC->{'spade'} ) {
				    			my ($idx) = 2;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    		#if ( exists $residue->{'damaged'} and defined $residue->{'damaged'} ) {
					    		#	$r_clas .= " sd ";
					    		#}
					    		#else {
					    		#	$r_clas .= " s ";
					    		#}
					    		$r_clas .= " $m ";
				    			$pep_met .= $met.',';
				    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
				    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
				    			}				    			
				    		}
				    		elsif ( $met eq $METHOD_DESC->{'corsair'} ) {
				    			my ($idx) = 3;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
				    			#$r_clas .= " $m ";
				    			#$pep_met .= $met.',';
				    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
				    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
				    			}				    			
				    		}				    		
				    		elsif ( $met eq $METHOD_DESC->{'thump'} ) {
				    			my ($idx) = 4;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
				    			#$r_clas .= " $m ";
				    			#$pep_met .= $met.','
				    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
				    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
				    			}				    			
				    		}
				    		elsif ( $met eq $METHOD_DESC->{'crash'} ) {
				    			my ($idx) = 5;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
				    			$r_clas .= " $m ";
				    			$pep_met .= $met.',';
				    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
				    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
				    			}				    			
				    		}
				    		elsif ( $met eq $METHOD_DESC->{'appris'} ) {
				    			my ($idx) = 6;
				    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
				    			if ( exists $residue->{'annot'} and defined $residue->{'annot'} and ($residue->{'annot'} =~ /PRINCIPAL/) ) {
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}
				    			}				    			
				    		}
				    		elsif ( $met eq $METHOD_DESC->{'proteo'} ) {
				    			$r_clas .= " p ";
				    			$pep_met .= $met.',';
				    		}				    		
				    	}
				    	if ( ($r_clas ne '') and ($r_clas) ) {
					    	$pep_met =~ s/\,$//;				    	
					    	$pep_res .= "<seq-viewer id='$seq_id' methods='$pep_met' residue='$pep_idx' class='seq-browser-label $r_clas' title='Click to show the annotations'>" . $pep . "</seq-viewer>";				    		
				    	}
				    	else {
				    		$pep_res .= $pep;
				    	}				    	
				    }
				    else {
				    	$pep_res .= $pep;
				    }
			    	$pep_idx++;
			    	$pep_len++;
			    }
			    else {
			    	$pep_res .= $pep;
			    	$seq_gaps->{$seq_id}++;
			    }					    	
			}
			$seq_pep .= $pep_res . "\n";
			
			$s_idx += $max_peptides;
			
		} until ($s_idx >= $seq_length);
		
		$seq_result .= ">" . $seq_id . "|" . $seq_length . "\n";
		$seq_result .= $seq_pep . "\n";			
		
	}
	
	if ( $seq_result ne '' ) {
		foreach my $seq_rep ( @{$seq_report} ) {
			my ($seq_id) =  $seq_rep->{'id'};
			my ($seq_pan) = '';
			foreach my $seq_met (@{$seq_panel->{$seq_id}}) {
				if ( defined $seq_met->{'span'} ) {
					$seq_pan .= $seq_met->{'span'};	
				}					
			}
			if ( $seq_pan ne '' ) {
				$panel .= $seq_id . "   " . $seq_pan. "\n";	
			}
		}
	}	
	$result = $seq_result;

	$result = $seq_result;
				  	
	return ($panel, $result);
	
} # end get_seq_features

sub get_aln_features {
	my ($self) = shift;
	my ($seq_string) = shift;
	my ($residues) = shift;	
	my ($max_peptides) = shift; $max_peptides = 90 unless ( defined $max_peptides );
	my ($max_char_title) = shift; $max_char_title = 30 unless ( defined $max_char_title );
	my ($panel, $result) = (undef,undef);
	my ($seq_panel);
	my ($seq_report);
	my ($seq_gaps);		
	
	require IO::String;
	require Bio::AlignIO;
	my ($stringfh) = new IO::String($seq_string);
	my ($alnio) = Bio::AlignIO->new(-fh => $stringfh, -format => 'clustalw'); # have to be clwstrict
	while( my $aln = $alnio->next_aln ) {	
		for my $seq ( $aln->each_seq ) {
			my ($seq_id) = $seq->id; $seq_id =~ s/^([^\|]*)\|[^\$]*$/$1/g;				
			my ($seq_seq) = $seq->seq;
			local $seq_panel->{$seq_id} = $SEQ_PANEL;
			push(@{$seq_report}, {
				'id'	=> $seq_id,
				'seq'	=> $seq->seq
			});			
			$seq_gaps->{$seq_id} = 0;
		}
		my ($aln_length) = $aln->length;
		my ($aln_match_line) = $aln->match_line()."\n";
		
		my ($num_title) = '';
		my ($n) = $max_char_title - 1;
		$num_title =~ s/^(.*)/$1 . ' ' x $n/mge;
		my ($num_pep) = 1;
		for ( my $i = $num_pep; $i < $max_peptides; $i++) {
			if ( ($i % 10) == 0 ) { $num_pep .= '|' }
			else { $num_pep .= '.' }
		}
		$num_pep .= "$max_peptides";
		my ($seq_result) = 	$num_title . $num_pep . "\n\n";

		my ($s_idx) = 0;
		do {
			foreach my $seq_rep ( @{$seq_report} ) {
				my ($seq_id) =  $seq_rep->{'id'};
				my ($seq_seq) =  $seq_rep->{'seq'};				
				my ($s_seq) = substr( $seq_seq, $s_idx, $max_peptides);				
				my ($pep_res) = '';				
		    	my ($pep_idx) = $s_idx - $seq_gaps->{$seq_id} + 1;
		    	my ($pep_len) = 0;			    	
				for ( my $i = 0; $i < length($s_seq); $i++ ) { 
				    my ($pep) = substr( $s_seq, $i, 1);
				    if ( $pep ne '-' ) {
					    if ( exists $residues->{$seq_id} and exists $residues->{$seq_id}->{$pep_idx} ) {
					    	my ($r_clas) = '';
					    	my ($pep_met) = '';
					    	foreach my $residue (@{$residues->{$seq_id}->{$pep_idx}}) {
					    		my ($met) = $residue->{'method'};				    		
					    		if ( $met eq $METHOD_DESC->{'firestar'} ) {
					    			my ($idx) = 0;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    			$r_clas .= " $m ";
					    			$pep_met .= $met.',';
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}
					    		}
					    		elsif ( $met eq $METHOD_DESC->{'matador3d'} ) {
					    			my ($idx) = 1;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    			$r_clas .= " $m ";
					    			$pep_met .= $met.',';
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}
					    		}
					    		elsif ( $met eq $METHOD_DESC->{'spade'} ) {
					    			my ($idx) = 2;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
						    		#if ( exists $residue->{'damaged'} and defined $residue->{'damaged'} ) {
						    		#	$r_clas .= " sd ";
						    		#}
						    		#else {
						    		#	$r_clas .= " s ";
						    		#}
						    		$r_clas .= " $m ";
					    			$pep_met .= $met.',';
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}
					    		}
					    		elsif ( $met eq $METHOD_DESC->{'corsair'} ) {
					    			my ($idx) = 3;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    			#$r_clas .= " $m ";
					    			#$pep_met .= $met.',';
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}				    			
					    		}
					    		elsif ( $met eq $METHOD_DESC->{'thump'} ) {
					    			my ($idx) = 4;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    			#$r_clas .= " $m ";
					    			#$pep_met .= $met.','
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}				    			
					    		}
					    		elsif ( $met eq $METHOD_DESC->{'crash'} ) {
					    			my ($idx) = 5;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    			$r_clas .= " $m ";
					    			$pep_met .= $met.',';
					    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
					    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
					    			}				    			
					    		}
					    		elsif ( $met eq $METHOD_DESC->{'appris'} ) {
					    			my ($idx) = 6;
					    			my ($m) = $SEQ_PANEL->[$idx]->{'class'};
					    			if ( exists $residue->{'annot'} and defined $residue->{'annot'} and ($residue->{'annot'} =~ /PRINCIPAL/) ) {
						    			unless ( defined $seq_panel->{$seq_id}->[$idx]->{'span'} ) {
						    				$seq_panel->{$seq_id}->[$idx]->{'span'} = "<span class='seq-browser-label $m'><span class='glyphicon glyphicon-check'></span></span>";
						    			}
					    			}				    			
					    		}					    		
					    		elsif ( $met eq $METHOD_DESC->{'proteo'} ) {
					    			$r_clas .= " p ";
					    			$pep_met .= $met.',';
					    		}					    		
					    	}
					    	if ( ($r_clas ne '') and ($r_clas) ) {
						    	$pep_met =~ s/\,$//;				    	
						    	$pep_res .= "<seq-viewer id='$seq_id' methods='$pep_met' residue='$pep_idx' class='seq-browser-label $r_clas' title='Click to show the annotations'>" . $pep . "</seq-viewer>";				    		
					    	}
					    	else {
					    		$pep_res .= $pep;
					    	}
					    }
					    else {
					    	$pep_res .= $pep;
					    }
				    	$pep_idx++;
				    	$pep_len++;
				    }
				    else {
				    	$pep_res .= $pep;
				    	$seq_gaps->{$seq_id}++;
				    }					    	
				}
				my ($pep_title) = $seq_id;
				my ($pep_coord_i) = 0;
				my ($pep_coord_l) = 0;				
				if ( $pep_len > 0 ) {
					$pep_coord_i = $s_idx - $seq_gaps->{$seq_id} + 1;
					$pep_coord_l = $pep_coord_i + $pep_len - 1;
				}
				#$pep_title .= " " . "[" . $pep_coord_i . '-' . $pep_coord_l . "]";
				my ($n) = $max_char_title - length($pep_title);
				$pep_title =~ s/^(.*)/$1 . ' ' x $n/mge;								
				$seq_result .= $pep_title . $pep_res . "\n";
			}
			my ($consensus_title) = '';
			$consensus_title =~ s/^(.*)/$1 . ' ' x $max_char_title/mge; 				
			my ($consensus_pep) = substr( $aln_match_line, $s_idx, $max_peptides);
			
			$seq_result .= 	$consensus_title . $consensus_pep . "\n";				
			
			$s_idx += $max_peptides;
			
		
		} until ($s_idx >= $aln_length);
		
		if ( $seq_result ne '' ) {
			foreach my $seq_rep ( @{$seq_report} ) {
				my ($seq_id) =  $seq_rep->{'id'};
				my ($seq_pan) = '';
				foreach my $seq_met (@{$seq_panel->{$seq_id}}) {
					if ( defined $seq_met->{'span'} ) {
						$seq_pan .= $seq_met->{'span'};	
					}					
				}
				if ( $seq_pan ne '' ) {
					$panel .= $seq_id . "   " . $seq_pan. "\n";	
				}
			}
		}
			
		$result = $seq_result;
		
	}
	return ($panel, $result);
	
} # end get_aln_features

sub get_gen_features {
	my ($self) = shift;
	my ($query_id) = shift;
	my ($species) = shift;
	my ($assembly) = shift;
	my ($source) = shift;
	my ($dataset) = shift;
	my ($methods) = shift if (@_);
	my ($ids) = shift if (@_);
	my ($result) = undef;
	
	# species has to be defined
	if ( defined $species ) {
		my ($species_id) = lc($species); $species_id =~ s/\s/\_/g;
		my ($cfg_species) = $SERVER->{'species'}->{$species_id};

		# If not defined: get the official assembly (the last one), and the first dataset
		unless ( defined $assembly ) { if ( exists $cfg_species->{'official'} ) { $assembly = $cfg_species->{'official'} } }
				
		# create exporter URL to get the BED tracks
		my ($required_params) = 'format=bed' . '%26' . "as=$assembly";
		my ($optional_params) = '';
		if ( defined $methods ) {
			$optional_params .= '%26'.'methods='.$methods;
			if ( $species_id eq 'homo_sapiens' ) { $optional_params .= ',proteo' }
		}
		if ( defined $ids ) { $optional_params .= '%26'.'ids='.$ids }
		if ( defined $source ) { $optional_params .= '%26'.'sc='.$source }
		if ( defined $dataset ) { $optional_params .= '%26'.'ds='.$dataset }
		my ($query) = 'http://' . CGI->new()->server_name() . '/rest/' . $query_id . '?' . $required_params . $optional_params ;
		
		# make a request to render tracks of UCSC
		my ($params) = 'db=' . $assembly;
		$params .= '&' . 'textSize=' . '12';
		$params .= '&' . 'hgt.labelWidth=' . '20';
		$params .= '&' . 'pix=' . '1100';
		#$params .= '&' . 'hideTracks=1';
		$params .= '&' . 'hgt.customText=' . $query;
		#$params .= '&' . 'ccdsGene=full';
		#$params .= '&' . 'ensGene=full';
		#$params .= '&' . 'knownGene=full';
		
		if ( defined $source and $source eq 'refseq' ) {
			$params .= '&' . 'refGene.label.acc=1';
			$params .= '&' . 'refGene.label.gene=0';
			$params .= '&' . 'refGene.label.omimhg38=0';
			$params .= '&' . 'refGene.hideNoncoding=1';
		}
		else {
			$params .= '&' . 'knownGene.label.gencodeId=1';
			$params .= '&' . 'knownGene.label.gene=0';
			$params .= '&' . 'knownGene.label.kgId=0';
			$params .= '&' . 'knownGene.label.prot=0';
			$params .= '&' . 'knownGene.label.omimhg38=0';
			$params .= '&' . 'knownGene.show.comprehensive=1';
			$params .= '&' . 'knownGene.show.spliceVariants=1';
			$params .= '&' . 'knownGene.show.noncoding=0';
			# for mouse and old species
			my ($gencode_db_for_mouse) = 'M4'; # HARD-CORE!!! due UCSC does not update correctly
			$params .= '&' . 'wgEncodeGencodeCompV'.$gencode_db_for_mouse.'_sel=' . '1';
			$params .= '&' . 'wgEncodeGencodeBasicV'.$gencode_db_for_mouse.'_sel=' . '0';
			$params .= '&' . 'wgEncodeGencodePseudoGeneV'.$gencode_db_for_mouse.'_sel=' . '0';
			$params .= '&' . 'wgEncodeGencodeCompV'.$gencode_db_for_mouse.'.label=' . 'accession';				
		}
		
		my ($ucsc_render_track_url) = $ENV{APPRIS_WSERVER_UCSC_RENDER_URL} . '?' . $params;
		my ($ucsc_query_link) = $ENV{APPRIS_WSERVER_UCSC_URL} . '?' . $params;

#return $query."<br><br>"."<br><br>".$ucsc_query_link."<br><br>"."<br><br>".$ucsc_render_track_url;

		$result = "<a class='imgUCSC' target='_blank' title='Click to alter the display of original UCSC Genome Browser' href='".$ucsc_query_link."'>".
					"<img class='imgTrackUCSC img-responsive' src='".$ucsc_render_track_url."'>".
				  "</a>";
	}
			  
	return $result;

} # end get_gen_features


##
## SEEKER
##

=head2 seek_features

  Arg [1]    : String $type
               type of input: 'id' or 'name'
  Arg [2]    : String $input
               gene or transcript identifier (or name)
  Arg [3]    : String $method
               method name
  Example    : use DBRetriever qw(seek_features);
               my $rst = seek_features($type, $input, $method);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub seek_features
{
	my ($self) = shift;
	my ($format) = shift;
	my ($result) = '';	
	my ($feat) = $self->get_features('none'); # methods is not defined => retrieve all features
		
	if ( defined $feat ) {		
		my ($specie_id) = $self->species; $specie_id = lc($specie_id); $specie_id =~ s/\s/_/; # get specie_id		
		my ($features);
		push(@{$features}, {
					'species'	=> $specie_id,
					'assembly'	=> 'none',
					'dataset'	=> 'none',
					'entity'	=> $feat			
		});
		my ($report) = $self->create_seeker_report($self->jobid, $features );

		if ( defined $report ) {
			if ($format eq 'json') {
				require JSON;
				my ($json) = new JSON;
				$result = $json->encode($report);
		    }
			elsif ($format eq 'xml') {
				require XML::LibXML;
				my ($xml_doc) = $self->get_xml_report($report);
				if ( defined $xml_doc ) {
					$result = $xml_doc->toString(1);
				}
		    }
		}		
	}
	
	return $result;
	
} # end seek_features


sub create_seeker_report($)
{
	my ($self) = shift;
	my ($query) = shift;
	my ($features) = shift;
	my ($report);
	$report->{'query'} = $query;
	
	foreach my $feature (@{$features})
	{
		my ($species) = $feature->{'species'};
		my ($assembly) = $feature->{'assembly'};
		my ($source) = $feature->{'source'};
		my ($dataset) = $feature->{'dataset'};
		my ($entity) = $feature->{'entity'};
		$assembly = undef if ( $assembly eq 'none' ); # undef assembly from WSRetriever jobs
		if ( $entity and ref($entity) ) {
			my ($match) = {
				'species'	=> $species,
				'assembly'	=> $assembly,
				'source'	=> $source,
				'dataset'	=> $dataset,
				'id'		=> $entity->stable_id,
				'version'	=> $entity->version,
				'label'		=> $entity->stable_id,
				'chr'		=> $entity->chromosome,
				'start'		=> $entity->start,
				'end'		=> $entity->end,
				'biotype'	=> $entity->biotype,
				'status'	=> $entity->status,
			};
			if ( $entity->chromosome and $entity->start and $entity->end ) {
				if ( $entity->isa("APPRIS::Gene") ) {
					$match->{'namespace'} = 'Gene_Id';	
				}
				elsif ( $entity->isa("APPRIS::Transcript") ) {
					$match->{'namespace'} = 'Transcript_Id';
				}					
			}
			if ( $entity->external_name ) {
				my ($dblink) = {
					'id'		=> $entity->external_name,
					'namespace'	=> 'External_Id'
				};
				push(@{$match->{'dblink'}}, $dblink);
			}
			if ( $entity->xref_identify ) {
				foreach my $xref_identify (@{$entity->xref_identify}) {
					my ($dblink) = {
						'id'		=> $xref_identify->id,
						'namespace'	=> $xref_identify->dbname
					};
					push(@{$match->{'dblink'}}, $dblink);
				}			
			}
			if ( ($entity->isa("APPRIS::Gene")) and $entity->transcripts ) {
				foreach my $transcript (@{$entity->transcripts}) {
					my ($dblink) = {
						'id'		=> $transcript->stable_id,
						'namespace'	=> 'Transcript_Id'
					};
					push(@{$match->{'dblink'}}, $dblink);
				}
			}
			push(@{$report->{'match'}}, $match);
		}
		else {
			return $report; # Empty return
		}
	}
	
	return $report;
	
} # end create_seeker_report

sub get_xml_report($)
{
	my ($self) = shift;
	my ($report) = shift;	
	my ($xml_doc) = XML::LibXML::Document->new('1.0','UTF-8');
	
	my ($e_query) = $xml_doc->createElement('query');
	$e_query->setNamespace("http://appris.bioinfo.cnio.es", "appris", 0);
	$xml_doc->setDocumentElement($e_query);
	my ($a_query) = $xml_doc->createAttribute('query', $report->{'query'});
	$e_query->setAttributeNode($a_query);

	foreach my $match (@{$report->{'match'}})
	{
		my ($e_match) = $xml_doc->createElement('match');
		my ($a_match) = $xml_doc->createAttribute('label', $match->{'label'});
		my ($a_match2) = $xml_doc->createAttribute('namespace', $match->{'namespace'});
		my ($a_match3) = $xml_doc->createAttribute('chr', $match->{'chr'});
		my ($a_match4) = $xml_doc->createAttribute('start', $match->{'start'});
		my ($a_match5) = $xml_doc->createAttribute('end', $match->{'end'});
		my ($a_match6) = $xml_doc->createAttribute('species', $match->{'species'});
		$e_match->setAttributeNode($a_match);
		$e_match->setAttributeNode($a_match2);
		$e_match->setAttributeNode($a_match3);
		$e_match->setAttributeNode($a_match4);
		$e_match->setAttributeNode($a_match5);
		$e_match->setAttributeNode($a_match6);

		my ($e_class) = $xml_doc->createElement('biotype');
		my ($e_class_txt) = $xml_doc->createCDATASection($match->{'biotype'});
		$e_class->appendChild($e_class_txt);
		$e_match->appendChild($e_class);

		my ($e_status) = $xml_doc->createElement('status');
		my ($e_status_txt) = $xml_doc->createCDATASection($match->{'status'});
		$e_status->appendChild($e_status_txt);		
		$e_match->appendChild($e_status);
		
		foreach my $dblink (@{$match->{'dblink'}}) {
			my ($e_dblink) = $xml_doc->createElement('dblink');
			my ($a_dblink) = $xml_doc->createAttribute('namespace', $dblink->{'namespace'});
			my ($a_dblink2) = $xml_doc->createAttribute('id', $dblink->{'id'});
			$e_dblink->setAttributeNode($a_dblink);
			$e_dblink->setAttributeNode($a_dblink2);
			$e_match->appendChild($e_dblink);
		}
		$e_query->appendChild($e_match);
	}

	return $xml_doc;
}

sub DESTROY {}

$SEQ_PANEL = [{			
	'firestar'	=> undef,
	'class'		=> "f",
	'span'		=> undef,
},{
	'matador3d'	=> undef,
	'class'		=> "m",
	'span'		=> undef,
},{
	'spade'		=> undef,
	'class'		=> "s",
	'span'		=> undef,
},{
	'corsair'	=> undef,
	'class'		=> "c",
	'span'		=> undef,
},{
	'thump'		=> undef,
	'class'		=> "t",
	'span'		=> undef,
},{
	'crash'		=> undef,
	'class'		=> "cr",
	'span'		=> undef,
},{
	'appris'	=> undef,	
	'class'		=> "a",
	'span'		=> undef,
}];



1;
