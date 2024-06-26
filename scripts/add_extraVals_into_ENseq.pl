#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Data::Dumper;

use APPRIS::Parser;
use APPRIS::Utils::File qw( getTotalStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;

use lib "$FindBin::Bin/lib";
use appris;

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
	$logger->info("-- get gene report from GTF file -------\n");
	my ($genedata) = APPRIS::Parser::_parse_indata($data_file);
	$logger->debug(Dumper($genedata)."\n");
	

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
			my ($isof_id);
			my ($transl_id) = '-';
			my ($transc_id) = '-';
			my ($gene_id) = '-';
			my ($gene_name) = '-';
			my ($ccds_id) = '-';
			
			#ENSP00000452494.1 pep:known chromosome:GRCh38:14:22449113:22449125:1 gene:ENSG00000228985.1 transcript:ENST00000448914.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRDD3 description:T cell receptor delta dive
			if ( $s_desc =~ /gene:([^\s]+)\s*transcript:([^\s]+)\s*/ ) { # Ensembl sequences
				$gene_id = $1;
				$transc_id = $2;
				$gene_name = $4;
				( my $tid = $transc_id) =~ s/\.[0-9]*$//g;
				( my $gid = $gene_id) =~ s/\.[0-9]*$//g;
				
				if ( exists $genedata->{$gid} and exists $genedata->{$gid}->{'transcripts'}->{$tid} and exists $genedata->{$gid}->{'transcripts'}->{$tid}->{'ccdsid'} ) {
					$ccds_id = $genedata->{$gid}->{'transcripts'}->{$tid}->{'ccdsid'};
				}
				
				$output .= '>en_a'.'|'.$s_id.' '.$s_desc.' ccds:'.$ccds_id."\n".
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

main();


1;

__END__

=head1 NAME

add_extraVals_into_ENseq

=head1 DESCRIPTION

Script that add the CCDS ids, Gene names, etc. into the comment of ENSEMBL fasta sequence.
 
=head1 SYNOPSIS

add_extraVals_into_ENseq

=head2 Required arguments:

	-f,--faa=  <FASTA sequence>
	
	-d,--dat=  <ENSEMBL report>
	
	-o,--out=    <Output FASTA file with extra values>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl add_extraVals_into_ENseq.pl
	
	-f,--faa=Homo_sapiens.GRCh38.pep.all.fa
	
	-d,--dat=Homo_sapiens.GRCh38.84.gtf

	-o,--out=Homo_sapiens.GRCh38.pep.all.extra.fa


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
