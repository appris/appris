#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;

###################
# Global variable #
###################

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($fasta_file) = undef;
my ($data_file) = undef;
my ($out_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'faa|f=s'			=> \$fasta_file,
	'dat|d=s'			=> \$data_file,
	'out|o=s'			=> \$out_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $fasta_file and defined $data_file and defined $out_file )
{
	print `perldoc $0`;
	exit 1;
}


# Optional arguments

# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log($str_params);

#####################
# Method prototypes #
#####################
sub create_xreference($$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->info("-- create cross-reference report -------\n");
	my ($xref) = create_xreference($fasta_file, $data_file);
	$logger->debug("XREF\n".Dumper($xref)."\n");

	$logger->info("-- scan fasta sequence -------\n");
	my ($output) = "";
	if (-e $fasta_file and (-s $fasta_file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $fasta_file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			my ($s_id) = $seq->id; 
			my ($s_desc) = $seq->desc;
			my ($s_seq) = $seq->seq;
			my ($s_len) = length($s_seq);
			my ($db);
			my ($isof_id);
			my ($transc_id) = '-';
			my ($gene_id) = '-';
			my ($gene_name) = '-';
			my ($ccds_id) = '-';
			
			if ( $s_id =~ /^(gi)\|[^|]*\|[^|]*\|([^|]*)/ ) { # RefSeq sequences
				$db = $1;
				$isof_id = $2;
				$transc_id = $xref->{'var'}->{$isof_id} if ( exists $xref->{'var'}->{$isof_id} );
				$gene_id = $xref->{'gene_id'}->{$isof_id} if ( exists $xref->{'gene_id'}->{$isof_id} );
				$gene_name = $xref->{'gene'}->{$isof_id} if ( exists $xref->{'gene'}->{$isof_id} );
				$ccds_id = $xref->{'ccds'}->{$isof_id} if ( exists $xref->{'ccds'}->{$isof_id} );
				$s_id =~ s/^gi\|//g;
				$s_id =~ s/\|$//g;
				$output .= ">".$db.'_a'.'|'.$s_id.'|'.$transc_id.'|'.$gene_id.'|'.$gene_name.'|'.$ccds_id.'|'.$s_len.' '.$s_desc."\n".
							$s_seq."\n";
			}
		}
	}
	
	# Print output
	if ($output ne '') {
		$logger->info("-- print output -------\n");
		my ($printing_file_log) = printStringIntoFile($output, $out_file);
		$logger->error("printing") unless(defined $printing_file_log);		
	}
	
	$logger->finish_log();
	
	exit 0;	
}

sub create_xreference($$) {
	my ($fasta_file, $data_file) = @_;
	my ($report);

	my ($in2) = Bio::SeqIO->new(
						-file	=> $data_file,
						-format	=> 'GenBank'
	);
	while ( my $seq = $in2->next_seq() ) {

		my ($isof_id) = $seq->accession_number().'.'.$seq->version();
		
		my $ann = $seq->annotation();
		foreach my $ref ( $ann->get_Annotations('dblink') ) {
			my ($transc_id) = $ref->primary_id();
			if ( $transc_id =~ /^[NM_|XM_]/ ) {
				if ( exists $report->{'var'}->{$isof_id} ) { $report->{'var'}->{$isof_id} .= ';'.$transc_id }
				else { $report->{'var'}->{$isof_id} = $transc_id }
			}
		}
		    
		for my $feat_object ($seq->get_SeqFeatures('CDS')) {
			for my $tag ($feat_object->get_all_tags) {
				if ( $tag eq 'coded_by' ) {
#			        for my $value ($feat_object->get_tag_values($tag)) {
#			        	my (@vals) = split(':', $value);
#			        	my ($transc_id) = $vals[0];
#						if ( exists $report->{'var'}->{$isof_id} ) { $report->{'var'}->{$isof_id} .= ';'.$transc_id }
#						else { $report->{'var'}->{$isof_id} = $transc_id }	        	
#			        }			
				}
				elsif ( $tag eq 'db_xref' ) {
			        for my $value ($feat_object->get_tag_values($tag)) {
			        	if ( $value =~ /CCDS:([^\$]*)$/ ) {
			        		my ($ccds_id) = $1;
							if ( exists $report->{'ccds'}->{$isof_id} ) { $report->{'ccds'}->{$isof_id} .= ';'.$ccds_id }
							else { $report->{'ccds'}->{$isof_id} = $ccds_id }	        	
			        	}
			        	elsif ( $value =~ /GeneID:([^\$]*)$/ ) {
			        		my ($gene_id) = $1;
							if ( exists $report->{'gene_id'}->{$isof_id} ) { $report->{'gene_id'}->{$isof_id} .= ';'.$gene_id }
							else { $report->{'gene_id'}->{$isof_id} = $gene_id }	        	
			        	}
			        }			
				}
				elsif ( $tag eq 'gene' ) {
			        for my $value ($feat_object->get_tag_values($tag)) {
						if ( exists $report->{'gene'}->{$isof_id} ) { $report->{'gene'}->{$isof_id} .= ';'.$value }
						else { $report->{'gene'}->{$isof_id} = $value }	        	
			        }			
				}
				elsif ( $tag eq 'gene_synonym' ) {
			        for my $value ($feat_object->get_tag_values($tag)) {
			        	$value =~ s/\s*//g;
						if ( exists $report->{'gene_syn'}->{$isof_id} ) { $report->{'gene_syn'}->{$isof_id} .= ';'.$value }
						else { $report->{'gene_syn'}->{$isof_id} = $value }	        	
			        }			
				}
		    }
		}    
	}
	
	return $report;
}
main();


1;

__END__

=head1 NAME

add_extraVals_into_RSseq

=head1 DESCRIPTION

Script that add the CCDS ids, Gene names, etc. into the comment of RefSeq fasta sequence.
 
=head1 SYNOPSIS

add_extraVals_into_UPseq

=head2 Required arguments:

	-f,--faa=  <RefSeq FASTA sequence>
	
	-d,--dat=  <Protein GenBak report>
	
	-o,--out=    <Output FASTA file with extra values>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl add_extraVals_into_UPseq.pl
	
	-f,--faa=protein.fa
	
	-d,--dat=protein.gbk

	-o,--out=protein.extra.fa


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
