=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Utils::WSpace - Object representing a gene

=head1 SYNOPSIS

  my $gene = APPRIS::Utils::WSpace->new(
    -chr	=> 1,
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print gene information
  print("gene start:end:strand is "
      . join( ":", map { $gene->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $gene->stable_id('ENSG000001');
  $gene->description('This is the gene description');

=head1 DESCRIPTION

A representation of a Gene within the APPRIS system.
A gene is a set of one or more alternative transcripts.

=head1 METHODS

=cut

package APPRIS::Utils::WSpace;

use strict;
use warnings;
use Digest::MD5;
use Bio::SeqIO;

use APPRIS::Utils::Argument qw(rearrange);
use APPRIS::Utils::Exception qw(throw warning deprecate);
use APPRIS::Utils::File qw( prepare_workspace getStringFromFile getTotalStringFromFile printStringIntoFile updateStringIntoLockFile );


use Data::Dumper;

{
    # Encapsulated class data
    #___________________________________________________________
    my %_attr_data =
		(
			id				=>  undef,
			file			=>  undef,
			path			=>  undef,
			findex			=>  undef,
			ws				=>  undef,
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
	
	sub id {
		my ($self, $arg) = @_;
		$self->{id} = $arg if defined $arg;
		return $self->{id};
	}

	sub file {
		my ($self, $arg) = @_;
		$self->{file} = $arg if defined $arg;
		return $self->{file};
	}
	
	sub path {
		my ($self, $arg) = @_;
		$self->{path} = $arg if defined $arg;
		return $self->{path};
	}
	
	sub findex {
		my ($self, $arg) = @_;
		$self->{findex} = $arg if defined $arg;
		return $self->{findex};
	}
	
}


=head2 new

  Arg [-file]: (optional) 
       string - fasta file
  Arg [-ticket_id]: (optional) 
       string - ticket id in md5
  Example    : $gene = APPRIS::Utils::WSpace->new(...);
  Description: Creates a new ticket object
  Returntype : APPRIS::Utils::WSpace
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub new {
	my ($caller, %args) = @_;
	
	my ($caller_is_obj) = ref($caller);
	return $caller if $caller_is_obj;
	my ($class) = $caller_is_obj || $caller;
	my ($self) = bless {}, $class;

	foreach my $attrname ($self->_standard_keys) {
		my ($attr) = "-".$attrname;
		if (exists $args{$attr} && defined $args{$attr}) {
			$self->{$attrname} = $args{$attr};
		} else {
			$self->{$attrname} = $self->_default_for($attrname);
		}
	}
		
	return $self;
}

=head2 ws

  Arg [1]    : (optional) string $file
               override the default level
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->ws();
  Description: Get the ticket id from fasta input.
  Returntype : string
  Exceptions : thrown every time

=cut

sub ws
{
	my ($self) = shift;
	my ($ws) = undef;	

	if ( $self->{id} and $self->{path} ) {
		$ws = $self->{path}.'/'.$self->{id}.'/';
	}
	elsif ( $self->{file} and $self->{path} ) {
		
		# get md5 from seq file
		my ($idx) = $self->seq_idx($self->{file});				
		$ws = $self->{path}.'/'.$idx.'/' if ( defined $idx );
	}
		
	return $ws;
}

=head2 idx

  Arg [1]    : string $i
               gene/transcript id 
  Arg [2]    : string $data
               data file 
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->idx();
  Description: Create idx
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub idx
{
	my ($self, $id, $data) = @_;
	my ($idx) = undef;
		
	if ( defined $id and defined $data ) {
		# extract transc/cds_coords
		my (@out);
		eval {
			#my ($cmd) = 'grep '.$id.' '.$data.' | awk -F "\t" \'{if( $3=="CDS" && match($9,/transcript_id\s*([^ ]*)/, TRANS) ){ if( match(TRANS[1],/^\"ENS/) ) { gsub(/\.[0-9]*/,"",TRANS[1]) } print TRANS[1]">"$1":"$4"-"$5}}\' | sed -r \'s/chr|\"|\;//g\' | sort -t- -nk1 | tr \'\n\' \';\' ';
			my ($cmd) = 'grep '.$id.' '.$data.' | awk -F "\t" \'{if( $3=="CDS" && ( match($9,/transcript_id\s*([^ ]*)/, TRANS) || match($9,/Genbank:([^,]*)/, TRANS) ) ){ if( match(TRANS[1],/^\"ENS/) ) { gsub(/\.[0-9]*/,"",TRANS[1]) } print TRANS[1]">"$1":"$4"-"$5}}\' | sed -r \'s/chr|\"|\;//g\' | sort -t- -nk1 | tr \'\n\' \';\' ';
			@out = `$cmd`;
		};
		throw("getting coords: ".$!) if($@);
		
		# create id (md5)
		if ( scalar(@out) > 0 ) {
			my ($data) = $out[0];
			eval {
				my ($ctx) = Digest::MD5->new;
				$ctx->add($data);
				$idx = $ctx->hexdigest;		
			};
			throw('Creating md5') if ($@)
		}
		else { throw("creating idx") }
	}
	return $idx;
}

=head2 seq_idx

  Arg [1]    : string $file
               sequence file
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->seq_idx();
  Description: Get the ticket id from fasta input.
  Returntype : string
  Exceptions : thrown every time

=cut

sub seq_idx
{
	my ($self, $file) = @_;
	my ($idx) = undef;
	
	my ($data) = '';
	eval {		
		my ($in) = Bio::SeqIO->new(
					-file => $file,
					-format => 'Fasta'
		);
		while ( my $seqObj = $in->next_seq() ) {
			my ($seq_id) = $seqObj->id;
			my ($seq) = $seqObj->seq;			
			$data .= $seq_id.':'.$seq.'|';
		}
		$data =~ s/\|$//;
	};
	throw('Argument must be a correct fasta file') if ($@);
	
	# create id (md5)
	eval {
		my ($ctx) = Digest::MD5->new;
		$ctx->add($data);
		$idx = $ctx->hexdigest;		
	};
	throw('Creating md5') if ($@);
	
	return $idx;	
}

=head2 t_idxs

  Arg [1]    : string $i
               gene/transcript id 
  Arg [2]    : string $data
               data file 
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->t_idxs();
  Description: Create idx
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub t_idxs
{
	my ($self, $id, $data) = @_;
	my ($t_idxs) = undef;
	
	# get the list of trancs with their idx
	if ( defined $id and defined $data ) {	
		my (@tid_list);
		eval {
			my ($cmd) = 'grep '.$id.' '.$data.' | awk -F "\t" \'{if( $3=="CDS" && match($9,/transcript_id\s*([^ ]*)/, TRANS) ){ if( match(TRANS[1],/^\"ENS/) ) { gsub(/\.[0-9]*/,"",TRANS[1]) } print TRANS[1]}}\' | sed -r \'s/chr|\"|\;//g\' | sort -u ';
			@tid_list = `$cmd`;
		};
		throw("getting list of tid: ".$!) if($@);
		
		foreach my $tid (@tid_list) {
			$tid =~ s/\s*//g;
			my ($tidx) = $self->idx($tid, $data);
			if ( defined $tid and defined $tidx ) {
				$t_idxs->{$tidx} = $tid;
			}
		}		
	}
	
	return $t_idxs;
}

=head2 add_idx

  Arg [1]    : string $idx
               Index data
  Arg [2]    : string $id
               Id associated to Index data
  Arg [3]    : string $path
               Path where cached file is saved
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->add_idx();
  Description: Insert the filename and ticket_id into index file.
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub add_idx
{
	my ($self, $idx, $id, $path) = @_;
	my ($added) = undef;
	
	if ( defined $idx and defined $id and defined $path ) {
		my ($exists) = $self->exist_idx($idx, $id, $path);
		unless ( defined $exists ) {
			my ($findex) = $path.'/cached.txt';
			my ($idx_cont) = $idx."\t".$id."\n";
			my ($p) = APPRIS::Utils::File::updateStringIntoLockFile($idx_cont, $findex);
			throw('Updating index file') unless ( defined $p );
		}
		$added = 1;
	}
	
	return $added;
}

=head2 exist_idx

  Arg [1]    : string $idx
               Index data
  Arg [2]    : string $id
               Id associated to Index data
  Arg [3]    : string $path
               Path where cached file is saved
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->exist_idx();
  Description: Says if ticket exist_idx
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub exist_idx
{
	my ($self, $idx, $id, $path) = @_;
	my ($exist) = undef;
	
	if ( defined $idx and defined $id and defined $path ) {
		my ($findex) = $path.'/cached.txt';
		my ($report) = $self->extract_findex($findex);
		if ( defined $report and exists $report->{$idx} and defined $report->{$idx} ) {
			my ($cached_id) = $report->{$idx};
			if ( defined $cached_id and ($id eq $cached_id) ) {
				$exist = $cached_id;	
			}
		}			
	}
	
	return $exist;
}

=head2 extract_findex

  Arg [1]    : string $file
               Cached file
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->extract_findex();
  Description: Says if ticket extract_findex
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub extract_findex
{
	my ($self, $findex) = @_;
	my ($report) = undef;
	
	my ($cont) = APPRIS::Utils::File::getTotalStringFromFile($findex);
	foreach my $line (@{$cont}) {
		my (@cols) = split("\t", $line);
		if ( scalar(@cols) > 1 ) {
			my ($idx) = $cols[0]; $idx =~ s/\s*//mg;
			my ($id) = $cols[1]; $id =~ s/\s*//mg;
			$report->{$idx} = $id;
		}
	}
	
	return $report;
}

=head2 cdir

  Arg [1]    : (optional) List of string $dirs
               override the default value
  Arg [2]    : (optional) string $path
               base of dire where we cdir the workspace
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->cdir();
  Description: Create workspace from ticketid and paths.
  Returntype : base path or undef 
  Exceptions : thrown every time

=cut

sub cdir
{
	my ($self, $dirs, $path) = @_;
	
	eval {
		if ( UNIVERSAL::isa($dirs,'ARRAY') ) {
			if ( defined $path ) {
				foreach my $d (@{$dirs}) {
					my ($p) = APPRIS::Utils::File::prepare_workspace($path.'/'.$d);
					throw('Creating directory') unless ( defined $p );
				}				
			}
		}	
		else {
			my ($p) = APPRIS::Utils::File::prepare_workspace($dirs);
			throw('Creating directory') unless ( defined $p );
		}	
	};
	throw('Argument must be a correct') if ($@);
	
	return $dirs;
}

sub DESTROY {}

1;
