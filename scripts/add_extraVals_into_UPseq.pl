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
use vars qw(
	$LOCAL_PWD
	$GENCODE_TSL_FILE
);

$LOCAL_PWD					= $FindBin::Bin; $LOCAL_PWD =~ s/bin//;

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($fasta_file) = undef;
my ($tab_file) = undef;
my ($data_file) = undef;
my ($out_file) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'faa|f=s'			=> \$fasta_file,
	'tab|t=s'			=> \$tab_file,
	'dat|d=s'			=> \$data_file,
	'out|o=s'			=> \$out_file,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $fasta_file and defined $tab_file and defined $data_file and defined $out_file )
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
sub create_xreference($$$);

#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	$logger->info("-- create cross-reference report -------\n");
	my ($xref) = create_xreference($fasta_file, $tab_file, $data_file);
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
			my ($id);
			my ($name);
			my ($isof_id);
			my ($gene_name) = '-';
			my ($ccds_id) = '-';
			
			if ( $s_id =~ /^(sp|tr)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences
				$db = $1;
				$isof_id = $2;
				$name = $3;
				$id = $isof_id; $id =~ s/\-[0-9]*$//g;
				$gene_name = $xref->{'gene'}->{$id} if ( exists $xref->{'gene'}->{$id} );
				$ccds_id = $xref->{'ccds'}->{$isof_id} if ( exists $xref->{'ccds'}->{$isof_id} );
				$output .= ">".$db.'_a'.'|'.$isof_id.'|'.$name.'|'.$id.'|'.$gene_name.'|'.$ccds_id.'|'.$s_len.' '.$s_desc."\n".
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

sub create_xreference($$$) {
	my ($fasta_file, $tab_file, $data_file) = @_;
	my ($report);
	
	# get the index of proteins -> isoforms id
	my ($index);
	my ($in) = Bio::SeqIO->new(
						-file => $fasta_file,
						-format => 'Fasta'
	);
	while ( my $seq = $in->next_seq() ) {
			my ($s_id) = $seq->id;
			my ($s_desc) = $seq->desc;
			my ($s_seq) = $seq->seq;
			my ($db);
			my ($id);
			my ($name);
			my ($isof_id);
			my ($gene_name) = '-';
			my ($ccds_id) = '-';
			
			if ( $s_id =~ /^(sp|tr)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences
				$db = $1;
				$isof_id = $2;
				$name = $3;
				$id = $isof_id; $id =~ s/\-[0-9]*$//g;
				$index->{$id}->{$isof_id} = 1;
			}
	}
#print STDERR "INDEX\n".Dumper($index)."\n";
	
	# get the CCDS -> protein
	# Add genename into output report
	my ($prot_ccds);
	my ($xref_list) = getTotalStringFromFile($tab_file);
	for ( my $i=1; $i < scalar(@{$xref_list}); $i++ ) { # first line is comment
		my (@cols) = split('\t', $xref_list->[$i]);
		my ($id) = $cols[0];
		my ($name) = $cols[1];
		my ($gene_name) = $cols[2];
		my ($ccds_str) = $cols[3];
		if ( defined $ccds_str and $ccds_str ne '' ) {
			my (@ccds_list) = split(';', $ccds_str);
			foreach my $ccds_id ( @ccds_list ) {
				$ccds_id =~ s/\s*//g;
				$ccds_id =~ s/\.[0-9]*//g;	
				$prot_ccds->{$ccds_id} = $id;
			}
		}
		if ( defined $gene_name and $gene_name ne '' ) {
			$gene_name =~ s/\s*//g;
			$report->{'gene'}->{$id} = $gene_name;
		}
	}
#print STDERR "CCDS\n".Dumper($prot_ccds)."\n";
#print STDERR "GENE\n".Dumper($prot_gn)."\n";
	
	
	# get the CCDS -> isoform id
	my (@ccds_isof_list) = `grep -e "CCDS;" $data_file`;
	foreach my $ccds_isof_line (@ccds_isof_list) {
		my (@cols) = split(';', $ccds_isof_line);
		if ( scalar(@cols) >= 3 ) {
			my ($ccds_id) = $cols[1];
			$ccds_id =~ s/\s*//g;
			$ccds_id =~ s/\.[0-9]*//g;
			my ($isof_id);			
			if ( $cols[2] =~ /\[([^\]]*)\]/ ) {
				$isof_id = $1;
			}
			if ( defined $ccds_id and defined $isof_id ) {
				my ($id) = $isof_id; $id =~ s/\-[0-9]*$//g;
				if ( !exists $index->{$id}->{$isof_id} ) {
					$isof_id =~ s/\-[0-9]*$//g;
				}
				$report->{'ccds'}->{$isof_id} = $ccds_id;
			}
			elsif ( defined $ccds_id ) {
				if ( exists $prot_ccds->{$ccds_id} ) {
					my ($id) = $prot_ccds->{$ccds_id};
					$report->{'ccds'}->{$id} = $ccds_id;
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

add_extraVals_into_UPseq

=head1 DESCRIPTION

Script that add the CCDS ids, Gene names, etc. into the comment of UniProt fasta sequence.
 
=head1 SYNOPSIS

add_extraVals_into_UPseq

=head2 Required arguments:

	-f,--faa=  <UniProt FASTA sequence>
	
	-t,--tab=    <Extra MetaData from UniProt>
	
	-d,--dat=   <UniProt data report>
	
	-o,--out=    <Output FASTA file with extra values>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>

=head1 EXAMPLE

perl add_extraVals_into_UPseq.pl
	
	-f,--faa=uniprot-proteome.fasta
	
	-t,--tab=uniprot-proteome.tab
	
	-d,--dat=uniprot-proteome.txt

	-o,--out=uniprot-proteome.extra.fasta


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
