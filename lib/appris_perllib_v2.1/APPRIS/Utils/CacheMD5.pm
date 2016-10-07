=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Utils::CacheMD5 - Object representing a gene

=head1 SYNOPSIS

  my $gene = APPRIS::Utils::CacheMD5->new(
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

package APPRIS::Utils::CacheMD5;

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
			data			=>  undef,
			idx				=>  undef,			
			ws				=>  undef,
			md5_len			=>  32,
			subdir_len		=>  8,
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
	
	sub data {
		my ($self, $arg) = @_;
		$self->{data} = $arg if defined $arg;
		return $self->{data};
	}

	sub idx {
		my ($self, $arg) = @_;
		$self->{idx} = $arg if defined $arg;
		return $self->{idx};
	}

	sub ws {
		my ($self, $arg) = @_;
		$self->{ws} = $arg if defined $arg;
		return $self->{ws};
	}

}


=head2 new

  Arg [-file]: (optional) 
       string - fasta file
  Arg [-ticket_id]: (optional) 
       string - ticket id in md5
  Example    : $gene = APPRIS::Utils::CacheMD5->new(...);
  Description: Creates a new ticket object
  Returntype : APPRIS::Utils::CacheMD5
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

	my ($dat, $ws) = rearrange( ['dat','ws'], @_ );

	# require paramaters
	if ( defined $dat ) {
		$self->data($dat);
		$self->md5($dat);
	}
	else { return undef; }
	
	# optional parameters
	my ($d) = ( defined $ws ) ? $ws.'/'.$self->idx_s : $self->idx_s;
	$self->ws($d);
	
	return $self;
}

=head2 md5

  Arg [1]    : string $seq
               string to get MD5 value 
  Example    : use APPRIS::Utils::CacheMD5;
               $id = $ticket->idx();
  Description: Create idx
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub md5
{
	my ($self, $seq) = @_;
	my ($idx) = undef;
	
	if ( defined $seq ) {
		eval {
			my ($ctx) = Digest::MD5->new;
			$ctx->add($seq);
			my ($md5) = $ctx->hexdigest;
			if ( defined $md5 ) {
				$idx = $md5.'_'.length($seq);
				$self->idx($idx);
			}
		};
		throw('Creating md5') if ($@);
	}
	
	return $idx;
}

=head2 idx_s

  Arg [1]    : string $seq
               string to get MD5 value 
  Example    : use APPRIS::Utils::CacheMD5;
               $id = $ticket->idx();
  Description: Create idx
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub idx_s
{
	my ($self) = @_;
	my ($idx) = $self->idx;
	my ($idx_s) = '';	
	my ($md5_len) = $self->{md5_len};
	my ($subdir_len) = $self->{subdir_len};			
	for ( my $i=0; $i <= $md5_len-1; $i+=$subdir_len ) {
		my ($prefix_idx) = substr($idx, $i, $subdir_len);
		$idx_s .= $prefix_idx.'/';
	}
	my ($len) = $idx =~ /\_([^\$]*)/;
	$idx_s =~ s/\/$/\/$len/;	
	return $idx_s;
}

=head2 cdir

  Arg [1]    : (optional) tring $dir
               starting dir where to create cache
  Arg [2]    : (optional) string $ix
               MD5 value
  Example    : use APPRIS::Utils::CacheMD5;
               $id = $ticket->cdir();
  Description: Create workspace from ticketid and paths.
  Returntype : base path or undef 
  Exceptions : thrown every time

=cut

sub cdir
{
	my ($self) = @_;
	my ($ws) = $self->ws;
	my ($dat) = $self->data;
	
		# init control
	if ( ! -d $ws ) { # first time we create idx
		$self->mkidir($ws);
		$self->add_data($dat,'seq');
	}
	else {
		if ( ! $self->same_data($dat,'seq') ) { # idx exist with diffent data
			# create new idx
			eval {
				my ($cmd) = "ls -1d $ws\_*";
				my (@dirs) = `$cmd`;
				my ($inc) = scalar(@dirs) + 1;
				my ($idx) = $self->idx.'_'.$inc;
				$self->idx($idx);
				my ($idx_s) 
			};
			throw('Increasing idx') if ($@);
			
		}		
	}
	
	return $ws;
}

=head2 mkidir

  Arg [1]    : (optional) tring $dir
               starting dir where to create cache
  Arg [2]    : (optional) string $ix
               MD5 value
  Example    : use APPRIS::Utils::CacheMD5;
               $id = $ticket->cdir();
  Description: Create workspace from ticketid and paths.
  Returntype : base path or undef 
  Exceptions : thrown every time

=cut

sub mkidir
{
	my ($self, $dir) = @_;
	
	# create dirs
	eval {
		my ($cmd) = "mkdir -p $dir";
		system($cmd);
	};
	throw('Creating subdir idx') if ($@);
	
	return $dir;
}

=head2 add_data

  Arg [1]    : string $dat
               data
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->add_idx();
  Description: Insert the filename and ticket_id into index file.
  Returntype : boolean 
  Exceptions : thrown every time

=cut

#sub add_data
#{
#	my ($self, $dat, $fname) = @_;
#	my ($file) = $self->ws.'/'.$fname;
#	my ($p) = APPRIS::Utils::File::updateStringIntoLockFile($dat, $file);
#	throw('Add data') unless ( defined $p );
#}

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

sub add_data
{
	my ($self, $dat, $fname) = @_;
	my ($file) = $self->ws.'/'.$fname;
	if ( -e $file and -s $file > 0 ) {
		my ($rep) = $self->ext_meta($dat);
		my ($loc_dat) = APPRIS::Utils::File::getStringFromFile($file);
		my ($loc_rep) = $self->ext_meta($loc_dat);
		while ( my ($k,$v) = each(%{$rep}) ) {
			if ( !exists $loc_rep->{$k} ) {
				my ($p) = APPRIS::Utils::File::updateStringIntoLockFile($dat, $file);
				throw('Add data') unless ( defined $p );				
			}
		}
	}
	else {
		my ($p) = APPRIS::Utils::File::updateStringIntoLockFile($dat, $file);
		throw('Add data') unless ( defined $p );		
	}
}

=head2 ext_meta

  Arg [1]    : string $file
               Cached file
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->ext_meta();
  Description: Says if ticket ext_meta
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub ext_meta
{
	my ($self, $dat) = @_;
	my ($report) = undef;
	foreach my $line (split("\n", $dat)) {
		my (@cols) = split("\t", $line);
		if ( scalar(@cols) > 1 ) {
			my ($ds) = $cols[0]; $ds =~ s/\s*//mg;
			my ($ids) = $cols[1]; $ids =~ s/\s*//mg;
			my ($id) = ( $ids =~ /([^\|]*)/ ) ? $1 : $ids; 
			$report->{$id} = {
				'ds'  => $ds,
				'ids' => $ids
			};
		}
	}	
	return $report;
}

=head2 same_data

  Arg [1]    : string $dat
               data
  Example    : use APPRIS::Utils::WSpace;
               $id = $ticket->add_idx();
  Description: Insert the filename and ticket_id into index file.
  Returntype : boolean 
  Exceptions : thrown every time

=cut

sub same_data
{
	my ($self, $dat, $fname) = @_;	
	my ($file) = $self->ws.'/'.$fname;
	my ($loc_dat) = APPRIS::Utils::File::getStringFromFile($file);
	if ( defined $loc_dat and $loc_dat eq $dat ) { return 1 }
	else { return undef }
}

sub DESTROY {}

1;
