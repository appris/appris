=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::CDS - Object representing a cds

=head1 SYNOPSIS

  my $cds = APPRIS::CDS->new(
    -start  => 123,
    -end    => 1045,
    -strand => '+',
  );

  # print cds information
  print("cds start:end:strand is "
      . join( ":", map { $cds->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $cds->stable_id('ENSE000001');

=head1 DESCRIPTION

A representation of a CDS within the APPRIS system.
This is a class which represents a coding region of 
a transcript (translation). 

=head1 METHODS

=cut

package HTTP;

use strict;
use warnings;

use CGI;
use HTTP::Status qw(is_success status_message);

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	print_http_response
);

=head2 print_http_response

  Arg [1]    : int $num
               HTTP status code
  Arg [2]    : string $content
               output result
  Arg [3]    : string $type
               content type
  Example    : use HTTP qw(print_http_response);
               print_http_response('500', 'content result', 'application/json');
  Description: Print the output.
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut


sub print_http_response($;$;$)
{
	my ($http_error_num,$content,$type) = @_;
	
	my($http_response)=new CGI();
	my($status_line)=$http_error_num." ".status_message($http_error_num);
	if (is_success($http_error_num)) 
	{
		$type = 'text/plain' unless ( defined $type );
		if ( $type eq 'redirect' ) {
			print $http_response->redirect($content);
		}
		else {
			print $http_response->header(
					-status	=> $status_line,
					-type	=> $type
			);
			print $content if (defined $content and ($content ne ''));			
		}
	}
	else
	{
		print $http_response->header(
				-status => $status_line,
				-type => 'text/html'
		);
		print $status_line;
		print ": ".$content if (defined $content and ($content ne ''));
		#print $http_response->start_html($status_line);
		#print $http_response->h1($status_line);
		#print $http_response->p($content);
		#print $http_response->end_html;
	}
	$http_error_num==200? exit 0: exit $http_error_num;
}


1;