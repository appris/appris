=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::Export::JSON - Utility functions for error handling

=head1 SYNOPSIS

  use APPRIS::Export::JSON qw(get_annotations);
  
  or to get all methods just

  use APPRIS::Export::JSON;

  eval { get_annotations("text to file",file_path) };
  if ($@) {
    print "Caught exception:\n$@";
  }

=head1 DESCRIPTION

The functions exported by this package provide a set of useful methods 
to export database values as JSON format.

=head1 METHODS

=cut

package APPRIS::Export::JSON;

use strict;
use warnings;
use JSON;
use Data::Dumper;

use APPRIS::Utils::Exception qw(throw warning deprecate);


=head2 create

  Arg [1]    : String - $res
  Example    : create($res);  
  Description: Retrieves text as JSON format from given params.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub create {
    my ($report) = @_;
	my ($json) = new JSON;
    my ($output) = '';
    
	if (defined $report) {
		$output .= $json->encode($report);		
	}
	else {
		$output .= $json->encode([]);		
	}
	
	return $output;
}

=head2 create_array

  Arg [1]    : Object - $res
               GTF annotation
  Example    : create_array($res);  
  Description: Retrieves text as JSON format from given params.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub create_array {
	my ($report) = @_;
	my ($list);
	my ($json) = new JSON;
	my ($output) = '';
	
	push(@{$list}, $report);
	if (defined $list) {
		$output .= $json->encode($list);		
	}
	else {
		$output .= $json->encode([]);		
	}	
	
	return $output;
}


=head2 get_annotations

  Arg [1]    : String - $gtf
               GTF annotation
  Example    : get_annotations($gtf);  
  Description: Retrieves text as JSON format with the annotations.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_annotations {
    my ($gtf) = @_;
    my ($list);
	my ($json) = new JSON;    
    my ($output) = '';
    
    my (@lines) = split("\n", $gtf);
    foreach my $line (@lines) {
		next if($line =~ /^##/); #ignore header
		
		my ($seqname,$source,$type,$start,$end,$score,$strand,$frame,$attributes) = split("\t", $line);
		next unless(defined $seqname and 
					defined $source and
					defined $type and
					defined $start and
					defined $end and
					defined $score and
					defined $strand and
					defined $frame and 
					defined $attributes);
		#store nine columns in hash
		my ($fields) = {
				seqname        => $seqname,
				source     => $source,
				type       => $type,
				start      => $start,
				end        => $end,
				score      => $score,
				strand     => $strand,
				frame      => $frame
		};	
		#store ids and additional information in second hash
		my (@add_attributes) = split(";", $attributes);				
		for ( my $i=0; $i<scalar @add_attributes; $i++ )
		{
			$add_attributes[$i] =~ s/^\s//; $add_attributes[$i] =~ s/\s$//;			
			$add_attributes[$i] =~ /^([^\s]*)\s*([^\$]*)$/;
			my ($c_type) = $1;
			my ($c_value) = $2;
			if(	defined $c_type and !($c_type=~/^\s*$/) and
				defined $c_value and !($c_value=~/^\s*$/))
			{
				$c_type =~ s/^\s//;
				$c_value =~ s/"//g;
				if(!exists($fields->{$c_type})) {
					$fields->{$c_type} = $c_value;
				}
			}
		}
		push(@{$list}, $fields);
	}
	
	if (defined $list) {
		$output .= $json->encode($list);		
	}
	else {
		$output .= $json->encode([]);		
	}
	
	return $output;
}

=head2 get_residues

  Arg [1]    : Object - $res
               GTF annotation
  Example    : get_residues($res);  
  Description: Retrieves text as JSON format with the residues.
  Returntype : String or undef
  Exceptions : if we cant get the gene or transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub get_residues {
    my ($residues) = @_;
	my ($json) = new JSON;
    my ($output) = '';
    
	if (defined $residues) {
		$output .= $json->encode($residues);		
	}
	else {
		$output .= $json->encode([]);		
	}
	
	return $output;
}


=head2 get_results

  Arg [1]    : Hash report or undef
  Example    : $annot = get_results($report);
  Description: Retrieves tabular information of transcript.
  Returntype : String or undef

=cut

sub get_results {
	my ($report) = @_;
	my ($list);
	my ($json) = new JSON;
	my ($output) = '';
	
	push(@{$list}, $report);
	if (defined $list) {
		$output .= $json->encode($list);		
	}
	else {
		$output .= $json->encode([]);		
	}	
	
	return $output;
}


1;