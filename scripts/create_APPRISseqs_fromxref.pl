#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Bio::SeqIO;
use Digest::MD5;
use Data::Dumper;

use APPRIS::Utils::CacheMD5;
use APPRIS::Utils::File qw( getTotalStringFromFile getStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw warning deprecate );

use lib "$FindBin::Bin/lib";
use appris;

###################
# Global variable #
###################
use vars qw(
	$SOURCE_DB_LIST
);

$SOURCE_DB_LIST = ['ensembl','refseq','uniprot'];

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($xref_file) = undef;
my ($ensembl_seqfile) = undef;
my ($refseq_seqfile) = undef;
my ($uniprot_seqfile) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'xref|x=s'			=> \$xref_file,
	'ens|e=s'			=> \$ensembl_seqfile,
	'ref|r=s'			=> \$refseq_seqfile,
	'uni|u=s'			=> \$uniprot_seqfile,
	'outfile|o=s'		=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined $xref_file and 
	defined $ensembl_seqfile and 
	defined $refseq_seqfile and 
	defined $uniprot_seqfile and 
	defined $outfile
) {
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


#################
# Method bodies #
#################
# Main subroutine
sub main()
{
	my ($outcontent) = 	"";
	
	$logger->info("-- create seqdata from ENSEMBL -------\n");
	my ($ens_seqdata) = appris::create_seqdata($ensembl_seqfile);
#	$logger->debug("ENS_SEQDATA:\n".Dumper($ens_seqdata)."\n");
	
	$logger->info("-- create seqdata from REFSEQ -------\n");
	my ($ref_seqdata) = appris::create_seqdata($refseq_seqfile);
#	$logger->debug("REF_SEQDATA:\n".Dumper($ref_seqdata)."\n");

	$logger->info("-- create seqdata from UNIPROT -------\n");
	my ($uni_seqdata) = appris::create_seqdata($uniprot_seqfile);
#	$logger->debug("UNI_SEQDATA:\n".Dumper($uni_seqdata)."\n");
	
	#Ensemb_ID	RefSeq_ID	UniProt_IDs
	my ($flines) = getTotalStringFromFile($xref_file);
	for ( my $i=1; $i < scalar(@{$flines}); $i++ ) { # first line is the header
		my ($fline) = $flines->[$i];
		chomp($fline);
		my (@cols) = split('\t', $fline);
		my ($hgnc_en) = $cols[0];
		my ($ensembl_id) = $cols[1];
		my ($entrez_id) = $cols[2];
		my ($uniprot_list_ids) = $cols[3];
		my ($uniprot_id) = '-';
		$ensembl_id =~ s/\s*//g; $ensembl_id =~ s/\.[0-9]$//g;
		$entrez_id =~ s/\s*//g;
		$uniprot_list_ids =~ s/\s*//g;
		
		#ÊCreate report with all the protein sequences from XREF gene
		my ($seqreport);
		if ( exists $ens_seqdata->{$ensembl_id} ) {
			my ($found) = create_seqrep('ensembl', $ens_seqdata->{$ensembl_id}, \$seqreport);
		}
		if ( exists $ref_seqdata->{$entrez_id} ) {
			my ($found) = create_seqrep('refseq', $ref_seqdata->{$entrez_id}, \$seqreport);
		}
		foreach my $uniprot_ids (split('\|', $uniprot_list_ids) ) {
			if ( $uniprot_ids =~ /([^\>]*)\>([^\$]*)/ ) {
				foreach my $up_id (split(',', $2) ) {
					if ( exists $uni_seqdata->{$up_id} ) {
						my ($found) = create_seqrep('uniprot', $uni_seqdata->{$up_id}, \$seqreport);
					}
				}
			}
		}
		
		#ÊCreate MD5 idx from XREF genes
		my ($seqreport_idx) = '-';		
		my ($seqreport_id) = join('/', @cols);			
		eval {
			my ($ctx) = Digest::MD5->new;
			$ctx->add($seqreport_id);
			($seqreport_idx) = $ctx->hexdigest;
		};
		throw('Creating md5') if ($@);		
		
		# Extract the Fasta sequences
		$outcontent .= extract_seqfasta($seqreport_idx, $seqreport_id, $seqreport);
	}
	
	my ($p) = printStringIntoFile($outcontent, $outfile);
	
	$logger->finish_log();
	
	exit 0;
}

sub create_seqrep($$\$)
{
	my ($db, $g_report, $ref_seqrep) = @_;
	my ($gene_id) = $g_report->{'id'};
	my ($gene_name) = ( exists $g_report->{'name'} ) ? $g_report->{'name'} : '-';
	
	while (my ($transc_id, $t_report) = each(%{$g_report->{'varsplic'}}) ) {
		if ( exists $t_report->{'seq'} ) {
			# create cache sequence idx
			my ($seq_s) = $t_report->{'seq'};
			my ($cache) = APPRIS::Utils::CacheMD5->new( -dat => $seq_s );		
			my ($seq_idx) = $cache->idx;
			#my ($seq_sidx) = $cache->sidx;
			my ($ccds_id) = ( exists $t_report->{'ccds'} ) ? $t_report->{'ccds'} : '-';
			unless ( exists $$ref_seqrep->{$seq_idx} ) {
				$$ref_seqrep->{$seq_idx} = {
					'idx'		=> $seq_idx,
					'gene_id'	=> [{ $db => $gene_id }],
					'gene_name'	=> [{ $db => $gene_name }],
					'transc_id'	=> [{ $db => $transc_id }],
					'ccds_id'	=> [{ $db => $ccds_id }],
					'seq'		=> $seq_s
				}
			}
			else {
				if ( $seq_s eq $$ref_seqrep->{$seq_idx}->{'seq'} ) {
					push(@{$$ref_seqrep->{$seq_idx}->{'gene_id'}},   { $db => $gene_id });
					push(@{$$ref_seqrep->{$seq_idx}->{'gene_name'}}, { $db => $gene_name });
					push(@{$$ref_seqrep->{$seq_idx}->{'transc_id'}}, { $db => $transc_id });
					push(@{$$ref_seqrep->{$seq_idx}->{'ccds_id'}},   { $db => $ccds_id });
				}
				else { die "Diff protein sequence with equal MD5 IDX" }
			}
		}
	}
}

sub extract_seqfasta($$$)
{
	my ($seqreport_idx, $seqreport_id, $seqreport) = @_;
	my ($output) = '';	
	
	while (my ($seq_idx, $seqrep) = each(%{$seqreport}) )
	{
		#Êcreate list of identifiers from source db
		my ($seqrep_transc_id) = '';
		my ($seqrep_gene_id) = '';
		my ($seqrep_gene_name) = '';
		my ($seqrep_ccds_id) = '';
		foreach my $source_db (@{$SOURCE_DB_LIST}) {
			my ($srep_tid) = '';
			my ($srep_gid) = '';
			my ($srep_gn) = '';
			my ($srep_cid) = '';
			foreach my $grep (@{$seqrep->{'transc_id'}}) {
				my (@k) = keys(%{$grep});
				if ( $k[0] eq $source_db ) { $srep_tid .= $grep->{$source_db}.',' }				
			}
			foreach my $grep (@{$seqrep->{'gene_id'}}) {
				my (@k) = keys(%{$grep});
				if ( $k[0] eq $source_db ) { $srep_gid .= $grep->{$source_db}.',' }				
			}
			foreach my $grep (@{$seqrep->{'gene_name'}}) {
				my (@k) = keys(%{$grep});
				if ( $k[0] eq $source_db ) { $srep_gn .= $grep->{$source_db}.',' }				
			}
			foreach my $grep (@{$seqrep->{'ccds_id'}}) {
				my (@k) = keys(%{$grep});
				if ( $k[0] eq $source_db ) { $srep_cid .= $grep->{$source_db}.',' }				
			}
			$srep_tid =~ s/,$//g; $srep_gid =~ s/,$//g; $srep_gn =~ s/,$//g; $srep_cid =~ s/,$//g;
			if ( $srep_tid ne '' ) { $seqrep_transc_id .= $srep_tid.'/' } else { $seqrep_transc_id .= '-'.'/' }
			if ( $srep_gid ne '' ) { $seqrep_gene_id   .= $srep_gid.'/' } else { $seqrep_gene_id   .= '-'.'/' }
			if ( $srep_gn  ne '' ) { $seqrep_gene_name .= $srep_gn.'/'  } else { $seqrep_gene_name  .= '-'.'/' }
			if ( $srep_cid ne '' ) { $seqrep_ccds_id   .= $srep_cid.'/' } else { $seqrep_ccds_id   .= '-'.'/' }
		}
		$seqrep_transc_id =~ s/\/$//g; $seqrep_gene_id =~ s/\/$//g; $seqrep_gene_name =~ s/\/$//g; $seqrep_ccds_id =~ s/\/$//g;
		$output .= '>appris'.'|'.
					$seq_idx.' '.
					'xref_genes:'.$seqreport_idx.':'.$seqreport_id.' '.
					'transc:'.$seqrep_transc_id.' '.
					'genes:'.$seqrep_gene_id.' '.
					'gene_names:'.$seqrep_gene_name.' '.
					'ccds:'.$seqrep_ccds_id."\n".
					$seqrep->{'seq'}."\n";
	}
	
	
	return $output;
}

main();


1;

__END__

=head1 NAME

create_APPRISseqs_fromxref

=head1 DESCRIPTION

Create APPRIS fasta file from the Cross Reference Ensembl vs UniProt vs RefSeq 

=head1 SYNOPSIS

create_APPRISseqs_fromxref

=head2 Required arguments:

	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl create_APPRISseqs_fromxref.pl \

	-x  appris/workspaces/cmpEnsUniRef/xref/xref.ens84_rs107_UPentry-gene.txt \
	
	-e  features/homo_sapiens/e84_g24/Homo_sapiens.GRCh38.pep.all.extra.fa \
	
	-r  features/homo_sapiens/rs107/protein.extra.fa \
	
	-u  features/homo_sapiens/up201606/uniprot-proteome.extra.fasta \
	
	-o  features/appris_seqs.e84_rs107_up201606.transl.fa \
	
	--loglevel=debug --logfile=createAPPRISseqsfromXREF.log		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
