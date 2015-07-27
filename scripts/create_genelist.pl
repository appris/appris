#!/usr/bin/perl -W

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use FindBin;
use Config::IniFiles;
use Data::Dumper;

use lib "$FindBin::Bin/lib";
use appris qw( retrieve_gene_list );
use APPRIS::Utils::File qw( prepare_workspace printStringIntoFile getTotalStringFromFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw );

###################
# Global variable #
###################
use vars qw(
	$LOCAL_PWD
);

$LOCAL_PWD					= $FindBin::Bin;

# Input parameters
my ($numlist) = undef;
my ($chrlist) = undef;
my ($reverse) = undef;
my ($data_file) = undef;
my ($translfile) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;

&GetOptions(
	'num=s'				=> \$numlist,
	'chr=s'				=> \$chrlist,
	'rever'				=> \$reverse,
	'data=s'			=> \$data_file,
	'transl=s'			=> \$translfile,
	'outfile=s'			=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless ( defined $numlist and defined $translfile and defined $outfile ){
	print `perldoc $0`;
	exit 1;
}
	
# Get log filehandle and print heading and parameters to logfile
my ($logger) = new APPRIS::Utils::Logger(
	-LOGFILE      => $logfile,
	-LOGPATH      => $logpath,
	-LOGAPPEND    => $logappend,
	-LOGLEVEL     => $loglevel,
);
$logger->init_log();

#################
# Method bodies #
#################
sub _parse_len_seq($;$);
sub _avg_len_seq($);

# Main subroutine
sub main()
{	
	# if we 
	my ($filtered_genes);
	if ( defined $chrlist and defined $data_file ) {
		my ($gdata);
		$filtered_genes = appris::retrieve_gene_list($data_file, $chrlist);
	}
		
	# parse sequence file
	$logger->info("-- parsing sequence file\n");
	my ($lenseqreport) = _parse_len_seq($translfile,$filtered_genes);
#$logger->debug("SEQ_REPORT:\n".Dumper($lenseqreport)."\n");
	
	# get the avg of length of trans sequences
	$logger->info("-- getting the average of length\n");
	my ($avglenreport) = _avg_len_seq($lenseqreport);
#$logger->debug("AVG_REPORT:\n".Dumper($avglenreport)."\n");

	# sort genes by avg length
	my (@sortedgenes);
	$logger->info("-- sorting genes by the average length...");
	if ( defined $reverse ) {
		$logger->info("reversed sort\n");
		@sortedgenes = sort ({$avglenreport->{$a} <=> $avglenreport->{$b}} keys %{$avglenreport});		
	}
	else {
		$logger->info("normal sort\n");
		@sortedgenes = sort ({$avglenreport->{$b} <=> $avglenreport->{$a}} keys %{$avglenreport});
	}
#$logger->debug("SORTED_GENES:\n".Dumper(@sortedgenes)."\n");

	# create lists sorting by avg leng
	$logger->info("-- creating the list of genes\n");
	my ($outreport);
	my ($num) = 1;
	foreach my $gene_id (@sortedgenes) {
		if ( $num > $numlist ) { $num = 1 }
		push(@{$outreport->{$num}}, $gene_id);
		$num++;
		
	}
	
	# print gene list
	$logger->info("-- printing the list of genes\n");
	while (my ($i, $genes) = each(%{$outreport}) ) {
		my ($outfile_i) =  $outfile.".$i";
		my ($outcontext) = join("\n", @{$genes});
		my ($printing_file_log) = printStringIntoFile($outcontext, $outfile_i);
		$logger->error("printing output\n") unless(defined $printing_file_log);		
	}

	$logger->finish_log();
	
	exit 0;	
}

sub _parse_len_seq($;$)
{
	my ($file, $filter_list) = @_;
	my ($data);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() )
		{
			if ( $seq->id=~/([^|]*)\|([^|]*)/ )
			{
				my ($transc_id) = $1;
				my ($gene_id) = $2;
				$transc_id =~ s/\.\d*$//; # delete ensembl version
				$gene_id =~ s/\.\d*$//; # delete ensembl version
				if ( defined $filter_list and exists $filter_list->{$gene_id} ) { # apply filter
					if(exists $data->{$gene_id}->{$transc_id}) {
						throw("Duplicated sequence: $transc_id");
					}
					else {
						my ($sequence) = $seq->seq; # control short sequences
						my ($seq_len) = length($sequence);
						$data->{$gene_id}->{$transc_id} = $seq_len;
						unless ( exists $data->{$gene_id}->{'longest'} ) {
							$data->{$gene_id}->{'longest'}->{'id'} = $transc_id;
							$data->{$gene_id}->{'longest'}->{'len'} = $seq_len;
							$data->{$gene_id}->{'longest'}->{'num'} = 1;
						}
						else {
							if ( $seq_len > $data->{$gene_id}->{'longest'}->{'len'} ) {
								$data->{$gene_id}->{'longest'}->{'id'} = $transc_id;
								$data->{$gene_id}->{'longest'}->{'len'} = $seq_len;
							}
							$data->{$gene_id}->{'longest'}->{'num'}++;
						}
						if ( $seq_len <= 2 ) {
							$logger->warning("Short sequence: $transc_id\n");
						}
					}					
				}
			}
		}		
	}
	return $data;	
}

sub _avg_len_seq($)
{
	my ($report) = @_;
	my ($data);
	while ( my ($gene_id,$g_report) = each (%{$report}) ) {
		my ($num_transc) = 0;
		my ($sumleng_tranc) = 0;
		while ( my ($transc_id,$t_report) = each (%{$g_report}) ) {
			if ( $transc_id eq 'longest' ) {
				$num_transc = $t_report->{'num'};
			}
			else {
				$sumleng_tranc += $t_report;
			}
		}
		if ( ($num_transc != 0) and ($sumleng_tranc != 0) ) {
			$data->{$gene_id} = ($sumleng_tranc/$num_transc);	
		}
	}
	return $data;
}

main();


1;

__END__

=head1 NAME

create_genelist

=head1 DESCRIPTION

create some list of genes depending on length of protein sequences 

=head1 SYNOPSIS

apprisall

=head2 Input arguments:

	--num=  <Number of list>

	--transl= <Translation sequences file>

	--outfile= <Output file>
		
=head2 Optional arguments (log arguments):
	
	--chr <list of chromosomes you want to use>
	
	--data=  <Gene annotation file>
	
	--rever <Flag if we want reversed list of genes>

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
		
	--logpath=PATH <Write logfile to PATH (default: .)>
		
	--logappend= <Append to logfile (default: truncate)>

=head1 EXAMPLE

create_genelist

	--num=5	

	--chr=chr19,chr21
	
	--data={ABSOLUTE_PATH}/gencode.v20.annotation.GRCh38.gtf
	
	--transl={ABSOLUTE_PATH}/gencode.v20.pc_translations.GRCh38.fa
	
	--outfile={ABSOLUTE_PATH}/genelist/gencode.v20.chr19-21.genelist.txt

	--loglevel=INFO
	
	--logappend
	
	--logpath=logs/
			
	--logfile=create_genelist.log


=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
