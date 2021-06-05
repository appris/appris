#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;


###################
# Global variable #
###################
my ($LOGLEVEL) = {
	'none'      => 1,
	'error'     => 2,
	'warning'   => 3,
	'info'      => 4,
	'debug'     => 5,
	'verbose'   => 6,
};
my ($HEAD_COMMENT_LINE) = <<HEAD_COMMENT;
gene_stable_id	exon_coord_hg19	annot_label	annot_id	annot_transc_list
HEAD_COMMENT

# Input parameters
my ($infile) = undef;
my ($informat) = undef;
my ($outfile) = undef;
my ($loglevel) = 'info';
my ($logfile) = undef;
my ($logappend) = undef;
  
&GetOptions(
	'infile|i=s'		=> \$infile,
	'format|f=s'		=> \$informat,
	'outfile|o=s'		=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logappend'			=> \$logappend,
);

# Check input parameters
unless ( defined $infile and defined $informat and defined $outfile ) {
	print `perldoc $0`;
	print "\nPlease, check parameters\n\n";
	exit 1;
}
$loglevel = lc($loglevel) if ( defined $loglevel );

#################
# Method bodies #
#################

# Main subroutine
sub main()
{
	appris_http::print_log_message('appris_parser', 'Begin', 'info');
	
	my ($result) = parse_infile_data($infile, $informat);
	
	appris_http::write_file($outfile, $result);	
	
	appris_http::print_log_message('appris_parser', 'End', 'info');
}

# Parser input file 
sub parse_infile_data($)
{
	my ($infile, $informat) = @_;
	my ($result) = $HEAD_COMMENT_LINE;

	appris_http::print_log_message('parse_infile_data', 'Begin', 'info');
	appris_http::print_log_message('parse_infile_data', 'Filename: ' . $infile, 'debug');
	
	# Extract query file
	my ($query_data);
	if ( -f $infile || $infile eq '-' ) { # "-" for STDIN
		$query_data = appris_http::read_file( $infile );
	}
	else {
		die("Unable to open STDIN query file : $!");
	}
		
	# Split the input data by query
	foreach my $line ( split('# QUERY_RAW:',$query_data) ) {
		my ($raw,$data,$query,$q_result) = (undef,undef,undef,undef);	
		if ( $line =~ /^([^\n]*)\n+/m ) {
			$raw = $1;
		}
		if ( $line =~ /^# QUERY_DATA:([^\n]*)\n+/m ) {
			$data = $1;
		}
		if ( $line =~ /^# QUERY_REST:([^\n]*)\n+([^\$]*)$/m ) {
			$query = $1;
			$q_result = $2;
		}
		appris_http::print_log_message('parse_infile_data', 'queRaw: ' . $raw, 'verbose');
		appris_http::print_log_message('parse_infile_data', 'queData: ' . $data, 'verbose');
		appris_http::print_log_message('parse_infile_data', 'queQuery: ' . $query, 'verbose');
		if ( defined $raw and defined $data and defined $query and defined $q_result ) {
			my ($q_rst_report) = parse_result($q_result, $informat);
			if ( defined $q_rst_report ) {
				$result .= print_result($raw, $q_rst_report);
			}			
		}
	}
	
	appris_http::print_log_message('parse_infile_data', 'rstDat: ' . Dumper($result), 'verbose');
	appris_http::print_log_message('parse_infile_data', 'End', 'info');
	
	return $result;
	
} # end parse_infile_data

# Parser result depending on given format
sub parse_result($$)
{
	my ($indata, $informat) = @_;
	my ($outdata);

	appris_http::print_log_message('parse_result', 'Begin', 'info');
	appris_http::print_log_message('parse_result', "inData: \n" . $indata, 'verbose');
	
	
	if ( $informat eq 'gtf' ) {
		$outdata = parse_result_gtf($indata);		
	}
	elsif ( $informat eq 'json' ) {
		# TODO
		#$outdata = parse_result_json($indata);		
	}
	elsif ( $informat eq 'bed' ) {
		# TODO
		#$outdata = parse_result_bed($indata);		
	}
		
	appris_http::print_log_message('parse_result', 'outDat: ' . Dumper($outdata), 'verbose');
	appris_http::print_log_message('parse_result', 'End', 'info');
	
	return $outdata;
	
} # end parse_result

# Parser result depending on given format
sub parse_result_gtf($)
{
	my ($indata) = @_;
	my ($outdata);

	appris_http::print_log_message('parse_result_gtf', 'Begin', 'info');
	
	# Split the input data by line
	foreach my $line ( split('\n',$indata) ) {
		next if ( $line =~ /^##/ );
		my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$attributes) = split("\t", $line);
		next unless(defined $chr and 
					defined $source and
					defined $type and
					defined $start and
					defined $end and
					defined $score and
					defined $strand and
					defined $phase and 
					defined $attributes);
		my ($fields) = {
				chr        => $chr,
				source     => $source,
				type       => $type,
				start      => $start,
				end        => $end,
				score      => $score,
				strand     => $strand,
				phase      => $phase,
				attributes => $attributes,
		};
		if (defined $chr and $chr=~/chr(\w*)/) { $chr = $1 if(defined $1) }
				
		my ($attribs);
		my (@add_attributes) = split(";", $attributes);				
		for ( my $i=0; $i<scalar @add_attributes; $i++ )
		{
			$add_attributes[$i] =~ /^(.+)\s(.+)$/;
			my ($c_type) = $1;
			my ($c_value) = $2;
			if(	defined $c_type and !($c_type=~/^\s*$/) and
				defined $c_value and !($c_value=~/^\s*$/))
			{
				$c_type =~ s/^\s//;
				$c_value =~ s/"//g;
				if(!exists($attribs->{$c_type}))
				{
					$attribs->{$c_type} = $c_value;
				}
			}
		}
		if(	exists $attribs->{'gene_id'} and defined $attribs->{'transcript_id'} )
		{
			my ($gene_id) = $attribs->{'gene_id'};
			my ($transcript_id) = $attribs->{'transcript_id'};
			my ($region) = $chr.':'.$start.'-'.$end.'['.$strand.']';
			if ( defined $type and defined $region ) {
				my ($annot_id);
				my ($rpt) = {
					#'gene_id'	=> $gene_id,
					'transc_id'	=> $transcript_id,
					'region'	=> $region
				};
				if ( exists $attribs->{'note'} ) {
					my ($note) = $attribs->{'note'};
					$rpt->{'note'} = $note;
					if ( $note =~ /hmm_name:([^\,]*)/ ) { # Pfam domain
						$annot_id = $1;
					}
					if ( $note =~ /ligands:([^\,]*)/ ) { # Ligand
						$annot_id = $1;
					}
					if ( $note =~ /pep_seq:([^\,]*)/ ) { # Peptide seq
						$annot_id = $1;
					}					
				}
				# group the SPADE labels
				if ( $type eq 'domain_possibly_damaged' ) { $type = 'domain' }
				if ( $type eq 'domain_wrong' ) { $type = 'domain_damaged' }					
				push(@{$outdata->{$type}->{$annot_id}->{$gene_id}}, $rpt);
			}			
		}
	}	
	
	appris_http::print_log_message('parse_result_gtf', 'End', 'info');
	
	return $outdata;
	
} # end parse_result_gtf

# Print result
sub print_result($$)
{
	my ($init, $report) = @_;
	my ($result) = '';
	
	# domain damaged
	if ( exists $report->{'domain'} ) {
		my ($type) = 'domain';
		my ($annot) = $report->{$type};		
		while ( my ($annot_id,$annot_rep) = each(%{$annot}) ) {
			my ($list) = '';			
			while ( my ($gene_id,$g_rep) = each(%{$annot_rep}) ) {
				my ($t_list) = '';
				foreach my $rep ( @{$g_rep} ) {
					$t_list .= $rep->{'transc_id'}.",";				
				}
				$t_list =~ s/\,$//;
				$list .= $gene_id.":".$t_list.";";		
			}
			$list =~ s/\;$//;
			$result .= $init."\t".$type."\t"."https://pfam.xfam.org/family/$annot_id"."\t".$list."\n";
		}
	}
	# domain damaged
	if ( exists $report->{'domain_damaged'} ) {
		my ($type) = 'domain_damaged';
		my ($annot) = $report->{$type};		
		while ( my ($annot_id,$annot_rep) = each(%{$annot}) ) {
			my ($list) = '';			
			while ( my ($gene_id,$g_rep) = each(%{$annot_rep}) ) {
				my ($t_list) = '';
				foreach my $rep ( @{$g_rep} ) {
					$t_list .= $rep->{'transc_id'}.",";				
				}
				$t_list =~ s/\,$//;
				$list .= $gene_id.":".$t_list.";";		
			}
			$list =~ s/\;$//;
			$result .= $init."\t".$type."\t"."https://pfam.xfam.org/family/$annot_id"."\t".$list."\n";
		}
	}	
	# func residue section
	if ( exists $report->{'functional_residue'} ) {
		my ($type) = 'functional_residue';
		my ($annot) = $report->{$type};		
		while ( my ($annot_id,$annot_rep) = each(%{$annot}) ) {
			my ($list) = '';			
			while ( my ($gene_id,$g_rep) = each(%{$annot_rep}) ) {
				my ($t_list) = '';
				foreach my $rep ( @{$g_rep} ) {
					$t_list .= $rep->{'transc_id'}.",";				
				}
				$t_list =~ s/\,$//;
				$list .= $gene_id.":".$t_list.";";		
			}
			$list =~ s/\;$//;
			$result .= $init."\t".$type."\t".''."\t".$list."\n";			
		}
	}
	# proteo section
	if ( exists $report->{'proteomic_evidence'} ) {
		my ($type) = 'proteomic_evidence';
		my ($annot) = $report->{$type};		
		while ( my ($annot_id,$annot_rep) = each(%{$annot}) ) {
			my ($list) = '';			
			while ( my ($gene_id,$g_rep) = each(%{$annot_rep}) ) {
				my ($t_list) = '';
				foreach my $rep ( @{$g_rep} ) {
					$t_list .= $rep->{'transc_id'}.",";				
				}
				$t_list =~ s/\,$//;
				$list .= $gene_id.":".$t_list.";";		
			}
			$list =~ s/\;$//;
			$result .= $init."\t".$type."\t".''."\t".$list."\n";
		}
	}
	return $result;
		
} # end print_result


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
	    		my ($mode) = '>>';
	    		#$mode = '>>' if ( defined $logappend );			
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



1;

__END__

=pod

=head1 NAME

appris_parser

=head1 DESCRIPTION

script that make a query to APPRIS database using RESTful services 

=head1 SYNOPSIS

appris_access

=head2 Required arguments (inputs):

  -i, --infile= <Input file that contains APPRIS annotations>

  -f, --format <Format of APPRIS annotations: 'gtf', 'bed', or 'json'>
  
=head2 Required arguments (outputs):

  -o, --outfile= <Output file>

=head2 Optional arguments:
		
  
=head2 Optional arguments (log arguments):
	
  --loglevel= <Define log level: INFO,DEBUG,NONE (default: NONE)>
  
  --logfile=FILE <Log to FILE (default: *STDERR)>
  
  --logappend= <Append to logfile (default: truncate)>
  
=head1 EXAMPLEs

appris_access

	--infile=query1.out
	
	--format=gtf

	--outfile=query1.out.csv	

	--logfile=appris_parser.log
	
	--loglevel=debug
	
	--logappend
		
	
=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB,CNIO)

=cut


