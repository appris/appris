=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

WSRunner

=head1 DESCRIPTION

Package the retrieves the results of methods

=head1 SYNOPSIS

	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			conf			=>  undef,
			specie			=>  undef,
		);

	#_____________________________________________________________
    
  my $features = WSRunner->new(
					-specie => 'mouse',
				    -conf   => 'conf/config.ini',
  );
  
  
=head1 METHODS

=cut

package WSRunner;

use strict;
use warnings;
use FindBin;
use Time::localtime;
use Digest::MD5;
use Bio::SeqIO;
use JSON;
use Data::Dumper;

use WSAPPRIS;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);

###################
# Global variable #
###################
use vars qw(
	
);

{
	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = # DEFAULT
		(
			oper			=>  undef,
			params			=>  undef,
			jobid			=>  undef,
			run				=>  undef,
			status			=>  undef,
			result			=>  undef, # tsv, json, raw
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
	sub oper {
		my ($self, $arg) = @_;
		$self->{'oper'} = $arg if defined $arg;
		return $self->{'oper'};
	}
	
	sub params {
		my ($self, $arg) = @_;
		$self->{'params'} = $arg if defined $arg;
		return $self->{'params'};
	}
	
	sub jobid {
		my ($self, $arg) = @_;
		$self->{'jobid'} = $arg if defined $arg;
		return $self->{'jobid'};
	}
}

=head2 new

  Arg [-oper]:
        string - type of operation (run,status,result)
  Arg [-params]: (optional)
        string - list of parameters that will uses for execution (by default, all)
  Arg [-jobid]: (optional)
        string - job identifier
  Example    : $features = WSRunner->new(...);
  Description: Creates a new retrieve object
  Returntype : WSRunner
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

	my ($oper, $params, $jobid) = rearrange( ['oper','params', 'jobid'], @_ );

	# require paramaters
	if ( defined $oper ) { $self->oper($oper); }
	else { return undef; }
	
	# optional parameters
	$self->jobid($jobid) if ( defined $jobid );
	$self->params($params) if ( defined $params );
	
	return $self;
}

=head2 create_jobid

  Arg [1]    : (optional) string $file
               override the default level
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->create_id();
  Description: Get the ticket id from fasta input.
  Returntype : string
  Exceptions : thrown every time

=cut

sub create_jobid
{
	my ($self) = shift;
	my ($id) = undef;
	
	if ( $self->params ) {
		# create md5
		my ($params) = $self->params; 
		my ($data);
		
		# ... from sequence
		if ( exists $params->{'sequences'} ) {			
			my ($stringfh);			
			eval {
				open($stringfh, "<", \$params->{'sequences'}) or die "Could not open string for reading: $!";
				my ($in) = Bio::SeqIO->new(
				 			-fh     => $stringfh,
							-format	=> 'Fasta'
				);
				while ( my $seqObj = $in->next_seq() ) {
					my ($seq_id) = $seqObj->id;
					my ($seq) = $seqObj->seq;			
					$data .= $seq_id.':'.$seq.'|';
				}
				$data =~ s/\|$//;
			};
			warning('Argument must be a correct fasta file') if ($@);			
		}
		# ... from id_ever
		elsif ( exists $params->{'species'} and exists $params->{'id'} and exists $params->{'e_version'} ) {
			$data = $params->{'species'}.'_'.$params->{'id'}.'_'.$params->{'e_version'};
			# create (id) from datefrom sequence
			#$id = sprintf("%04d%02d%02d%02d%02d%02d",
			#						localtime()->year()+1900,
			#						localtime()->mon()+1,
			#						localtime()->mday(),
			#						localtime()->hour(),
			#						localtime()->min(),
			#						localtime()->sec()
			#);
			#$id = $params->{'id'}.'_'.$id if ( exists $params->{'id'} );
		}
		if ( defined $data ) {
			eval {
				my ($ctx) = Digest::MD5->new;
				$ctx->add($data);
				$id = $ctx->hexdigest;		
			};
			warning('Creating md5') if ($@);		
			$self->jobid($id);			
		}
	}
	
	return $id;
}

=head2 apply_request

  Arg [1]    : (optional) String $type
               type of request
  Example    : use WSRunner;
               my $request = apply_request($type);
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub apply_request {
	my ($self) = shift;
	$self->oper(shift) if(@_);
	my ($request);
		
	if ( $self->oper eq 'run' ) {
		my ($run) = $self->run;
		if ( defined $run ) {
			$request = JSON->new->encode($run);				
		}
	}
	elsif ( $self->oper eq 'status' ) {
		if ( $self->jobid ) {
			my ($status) = $self->status($self->jobid);
			if ( defined $status ) {
				$request = JSON->new->encode($status);				
			}
		}
	}
	elsif ( $self->oper eq 'resulttypes' ) {
		my ($jobid) = $self->jobid;
		my ($resulttypes) = $self->resulttypes($jobid);
		if ( defined $resulttypes ) {
			$request = JSON->new->encode($resulttypes);
		}
	}
	elsif ( $self->oper eq 'result' ) {
		if ( $self->jobid ) {
			my ($result) = $self->result($self->jobid);
			if ( defined $result ) {
				$request = $result;
			}
		}
	}
	return $request;
}


=head2 run

  Example    : use WSRunner;
               my $request = run();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub run {
	my ($self) = shift;

	my ($params) = $self->params;
	my ($jobid) = $self->create_jobid;
	my ($run) = WSAPPRIS::run_appris($jobid, $params);
	
	return $run;
}

=head2 status

  Arg [1]    : String $jobid
               jobid
  Example    : use WSRunner;
               my $request = status();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub status {
	my ($self) = shift;
	my ($jobid) = shift if(@_);	

	my ($status) = WSAPPRIS::status_appris($jobid);

	return ($status);
}

=head2 resulttypes

  Arg [1]    : (optional) String $type
               type of request
  Example    : use WSRunner;
               my $request = resulttypes();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub resulttypes {
	my ($self) = shift;
	my ($jobid) = shift;
	my ($resulttypes);
		
	if ( defined $jobid ) {
		my ($status) = $self->status($jobid);
		if ( ($status->{'status'} eq 'FINISHED') or ($status->{'status'} eq 'FAILURE') ) {
			my ($params) = $self->params;
			$resulttypes = WSAPPRIS::resulttypes_appris($jobid, $params);
		}
	}
	else {
		my ($params) = $self->params;
		$resulttypes = WSAPPRIS::resulttypes_appris($jobid, $params);
	}
	
	return $resulttypes;
}

=head2 result

  Arg [1]    : (optional) String $type
               type of request
  Example    : use WSRunner;
               my $request = result();
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub result {
	my ($self) = shift;
	my ($jobid) = shift;
	my ($result);
		
	if ( defined $jobid ) {
		my ($status) = $self->status($jobid);
		if ( ($status->{'status'} eq 'FINISHED') or ($status->{'status'} eq 'FAILURE') ) {
			my ($params) = $self->params;
			$result = WSAPPRIS::result_appris($jobid, $params);
		}
	}
	
	return $result;
}

sub DESTROY {}

1;
