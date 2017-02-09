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

#use lib $ENV{APPRIS_SCRIPTS_DIR}."/lib";
#use appris;

###################
# Global variable #
###################
use vars qw(
	$SOURCE_LIST
);

$SOURCE_LIST = ['ensembl','refseq','uniprot'];

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
	my ($outfasta, $outmeta) = 	('','');
	
	$logger->info("-- create seqdata from ENSEMBL -------\n");
	my ($ens_seqdata) = create_seqdata($ensembl_seqfile);
	$logger->debug("ENS_SEQDATA:\n".Dumper($ens_seqdata)."\n");
	
	$logger->info("-- create seqdata from REFSEQ -------\n");
	my ($ref_seqdata) = create_seqdata($refseq_seqfile);
	$logger->debug("REF_SEQDATA:\n".Dumper($ref_seqdata)."\n");

	$logger->info("-- create seqdata from UNIPROT -------\n");
	my ($uni_seqdata) = create_seqdata($uniprot_seqfile);
	$logger->debug("UNI_SEQDATA:\n".Dumper($uni_seqdata)."\n");
	
	#Ensemb_ID	RefSeq_ID	UniProt_IDs
	$logger->info("-- create xrefdata -------\n");
	my ($flines) = getTotalStringFromFile($xref_file);
	for ( my $i=1; $i < scalar(@{$flines}); $i++ ) { # first line is the header
		my ($fline) = $flines->[$i];
		chomp($fline);
		my (@cols) = split('\t', $fline);
		my ($hgnc) = $cols[0];
		my ($ensembl_id) = ( defined $cols[1] ) ? $cols[1] : ''; 
		my ($entrez_id) = ( defined $cols[2] ) ? $cols[2] : '';
		my ($uniprot_list_ids) = join(',', @cols[3 .. $#cols]);
		my ($uniprot_id) = '-';
		$ensembl_id =~ s/\s*//g; $ensembl_id =~ s/\.[0-9]$//g;
		$entrez_id =~ s/\s*//g;
		$uniprot_list_ids =~ s/\s*//g;
		$uniprot_list_ids =~ s/^\,//g; $uniprot_list_ids =~ s/\,$//g;
		my ($xref_seqid) = join('/', @cols);
		
		#ÊCreate report with all the protein sequences from XREF gene
		my ($seqreport);
		if ( defined $ensembl_id and exists $ens_seqdata->{$ensembl_id} ) {
			create_seqrep('ensembl', $ens_seqdata->{$ensembl_id}, \$seqreport);
		}
		if ( defined $entrez_id and exists $ref_seqdata->{$entrez_id} ) {
			create_seqrep('refseq', $ref_seqdata->{$entrez_id}, \$seqreport);
		}
		foreach my $uniprot_ids (split(';', $uniprot_list_ids) ) {
			my ($up_ids) = ( $uniprot_ids =~ /([^\>]*)\>([^\$]*)/ ) ? $2 : $uniprot_ids;
			foreach my $up_id (split(',', $up_ids) ) {
				if ( exists $uni_seqdata->{$up_id} ) {
					create_seqrep('uniprot', $uni_seqdata->{$up_id}, \$seqreport);
				}
			}
		}
		
		# Create report id
		my ($seqreport_id) = create_seqrep_ids($seqreport);

		# Extract the Fasta sequences
		my ($outfa, $outmet) = extract_seqfasta($xref_seqid, $seqreport_id, $seqreport);
		$outfasta .= $outfa if ($outfa ne '' );
		$outmeta .= $outmet if ($outmet ne '' );
	}
	
	my ($p) = printStringIntoFile($outfasta, $outfile);
	my ($p2) = printStringIntoFile($outmeta, $outfile.'.meta');
	
	$logger->finish_log();
	
	exit 0;
}

sub create_seqdata($)
{
	my ($file) = @_;
	my ($data);

	if (-e $file and (-s $file > 0) ) {
		my ($in) = Bio::SeqIO->new(
							-file => $file,
							-format => 'Fasta'
		);
		while ( my $seq = $in->next_seq() ) {
			my ($s_id) = $seq->id;
			my ($s_desc) = $seq->desc;
			my ($s_seq) = $seq->seq;
			my ($isof_id);
			my ($transl_id);
			my ($gene_id);
			my ($gene_name);
			my ($ccds_id);
			my ($seq_length);
			my ($a_gene_ids);
			my ($a_transc_ids);
			
			# At the moment, only for UniProt/neXtProt cases
			if ( $s_id =~ /^(sp|tr)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences
				$isof_id = $2;
				my (@desc) = split('GN=', $s_desc);
				if ( scalar(@desc) >= 2 ) {
					if ( $desc[1] =~ /([^\s]*)/ ) { $gene_name = $1 }					
				}				
			}
			elsif ( $s_id =~ /^(sp_a|tr_a)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^\$]*)$/ ) { # UniProt sequences with extra values
				$isof_id = $2;
				my ($name) = $3;
				$gene_id = $isof_id;
				$gene_id =~ s/\-[0-9]*$//g; 
				$gene_name = $5;
				$ccds_id = $6;
				$seq_length = $7;
			}
			elsif ( $s_id =~ /^ge_a\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^\$]*)$/ ) { # GENCODE sequences with extra values
				my ($pep_id) = $1;
				$isof_id = $2;
				$gene_id = $3;
				$gene_name = $4;
				$ccds_id = $5;
				$seq_length = $6;
			}
			elsif ( $s_id =~ /^en_a\|([^\s]*)/ ) { # ENSEMBL sequences with extra values
				my ($pep_id) = $1;
				if ( $s_desc =~ /gene:([^\s]+).*transcript:([^\s]+).*gene_symbol:([^\s]+).*ccds:([^\s]+)/ ) {
					$gene_id = $1;
					$isof_id = $2;
					$gene_name = $3;
					$ccds_id = $4;
				}
			}
			elsif ( $s_id =~ /^gi_a\|[^|]*\|[^|]*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^\s]*)\s*([^\$]*)$/ ) { # RefSeq sequences with extra values
				my ($pep_id) = $1;
				$isof_id = $2;
				$gene_id = $3;
				$gene_name = $4;
				$ccds_id = $5;
				$seq_length = $6;
			}
			elsif ( $s_id =~ /^nxp:([^\s]*)/ ) { # neXtProt sequences
				$isof_id = $1;
				my (@desc) = split('Gname=', $s_desc);
				if ( scalar(@desc) >= 2 ) {
					if ( $desc[1] =~ /([^\s]*)/ ) { $gene_name = $1 }					
				}				
			}
			elsif ( $s_id =~ /^appris\|([^|]*)\|([^\s]*)/ ) { # APPRIS sequences FASTA file
				$isof_id = $1;
				$gene_id = $2;
				if ( $s_desc =~ /gene_ids\>([^\s]+)/ ) { $a_gene_ids = $1 }
				if ( $s_desc =~ /gene_names\>([^\s]+)/ ) { $gene_name = $1 }
				if ( $s_desc =~ /transc_ids\>([^\s]+)/ ) { $a_transc_ids = $1 }
				if ( $s_desc =~ /ccds_ids\>([^\s]+)/ ) { $ccds_id = $1 }
			}
			if ( defined $isof_id and defined $gene_id ) {
				$gene_id =~ s/\.[0-9]+$//g; $isof_id =~ s/\.[0-9]+$//g; # delete version suffix
				unless ( defined $gene_id ) {
					$gene_id = ( $isof_id =~ /([^\-]*)/ ) ? $1 : $isof_id;					
				}
				unless ( exists $data->{$gene_id} ) {
					$data->{$gene_id} = {
						'id'		=> $gene_id,
						'varsplic'	=> {}
					};													
					if ( defined $gene_name and $gene_name ne '' and $gene_name ne '-' ) {
						$data->{$gene_id}->{'name'} = $gene_name;
					}
				}
				$data->{$gene_id}->{'varsplic'}->{$isof_id} = {
					'desc'		=> $s_desc,
					'seq' 		=> $s_seq
				};
				if ( defined $ccds_id and $ccds_id ne '' and $ccds_id ne '-' and $ccds_id ne '?' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'ccds'} = $ccds_id;
				}
				if ( defined $a_gene_ids and $a_gene_ids ne '' and $a_gene_ids ne '-' and $a_gene_ids ne '?' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'gene_ids'} = $a_gene_ids;
				}
				if ( defined $a_transc_ids and $a_transc_ids ne '' and $a_transc_ids ne '-' and $a_transc_ids ne '?' ) {
					$data->{$gene_id}->{'varsplic'}->{$isof_id}->{'transc_ids'} = $a_transc_ids;
				}
			}
		}		
	}
	return $data;
	
} # End create_seqdata

sub create_seqrep($$\$)
{
	my ($db, $g_report, $ref_seqrep) = @_;
	my ($gene_id) = $g_report->{'id'};
	my ($gene_name) = ( exists $g_report->{'name'} ) ? $g_report->{'name'} : undef;
	
	while (my ($transc_id, $t_report) = each(%{$g_report->{'varsplic'}}) ) {
		if ( exists $t_report->{'seq'} ) {
			# create cache sequence idx
			my ($seq_s) = $t_report->{'seq'};
			my ($cache) = APPRIS::Utils::CacheMD5->new( -dat => $seq_s );		
			my ($seq_idx) = $cache->idx;
			my ($ccds_id) = ( exists $t_report->{'ccds'} ) ? $t_report->{'ccds'} : undef;
			unless ( exists $$ref_seqrep->{$seq_idx} ) {
				$$ref_seqrep->{$seq_idx} = {
					'idx'		=> $seq_idx,
					'seq'		=> $seq_s,
					'iden'		=> {
									$db	=> {
										$gene_id => {
											'gene_id'	=> $gene_id,
											'transc'	=> {
												$transc_id => {
													'transc_id'	=> $transc_id
												}
											}
										}	
									}
					}
				};
				if ( defined $gene_name ) {
					$$ref_seqrep->{$seq_idx}->{'iden'}->{$db}->{$gene_id}->{'gene_name'} = $gene_name;
				}
				if ( defined $ccds_id ) {
					$$ref_seqrep->{$seq_idx}->{'iden'}->{$db}->{$gene_id}->{'transc'}->{$transc_id}->{'ccds_id'} = $ccds_id;
				}
			}
			else {
				if ( $seq_s eq $$ref_seqrep->{$seq_idx}->{'seq'} ) {
					$$ref_seqrep->{$seq_idx}->{'iden'}->{$db}->{$gene_id}->{'gene_id'} = $gene_id;
					$$ref_seqrep->{$seq_idx}->{'iden'}->{$db}->{$gene_id}->{'transc'}->{$transc_id}->{'transc_id'} = $transc_id;
					if ( defined $gene_name ) {
						$$ref_seqrep->{$seq_idx}->{'iden'}->{$db}->{$gene_id}->{'gene_name'} = $gene_name;
					}
					if ( defined $ccds_id ) {
						$$ref_seqrep->{$seq_idx}->{'iden'}->{$db}->{$gene_id}->{'transc'}->{$transc_id}->{'ccds_id'} = $ccds_id;
					}
				}
				else { die "Diff protein sequence with equal MD5 IDX" }
			}
		}
	}
}

sub create_seqrep_ids($)
{
	my ($seqreport) = @_;
	my ($seqreport_ids);
	my ($seqreport_gid,$seqreport_gn) = ('','');
	my ($sep_rep);
		
	while (my ($seq_idx, $seqrep) = each(%{$seqreport}) )
	{
		if ( exists $seqrep->{'iden'} ) {
			my ($seprep_iden) = $seqrep->{'iden'};
			foreach my $source (@{$SOURCE_LIST}) {
				if ( exists $seprep_iden->{$source} ) {
					foreach my $gene_id ( keys(%{$seprep_iden->{$source}}) ) {
						$sep_rep->{$source}->{'gene_id'}->{$gene_id} = 1;
						if ( exists $seprep_iden->{$source}->{$gene_id}->{'gene_name'} ) {
							my ($gn) = $seprep_iden->{$source}->{$gene_id}->{'gene_name'};
							$sep_rep->{'gene_name'}->{$gn} = 1;
						}
					}
				}
			}
		}
	}
	foreach my $source (@{$SOURCE_LIST}) {
		my ($sr_i) = '';
		if ( exists $sep_rep->{$source} ) {
			foreach my $g ( sort {$a cmp $b} keys(%{$sep_rep->{$source}->{'gene_id'}}) ) { $sr_i .= $g.',' }
		}
		else { $sr_i = '?'  }
		$sr_i =~ s/\,$//g;
		$seqreport_gid .= $source.':'.$sr_i.'+';
	}
	$seqreport_gid =~ s/\+$//g;
	if ( exists $sep_rep->{'gene_name'} ) {
		 $seqreport_gn = join(',', sort {$a cmp $b} keys(%{$sep_rep->{'gene_name'}}) );
	}
	my ($seqreport_gidx) = '';
	eval {
		my ($ctx) = Digest::MD5->new;
		$ctx->add($seqreport_gid);
		($seqreport_gidx) = $ctx->hexdigest;
	};
	throw('Creating md5') if ($@);	
	$seqreport_ids->{'gene_idx'} = $seqreport_gidx if ( $seqreport_gidx ne '' );
	$seqreport_ids->{'gene_id'} = $seqreport_gid;
	$seqreport_ids->{'gene_name'} = $seqreport_gn; 

	return $seqreport_ids;
}

sub extract_seqfasta($$$)
{
	my ($xref_seqid, $seqreport_ids, $seqreport) = @_;
	my ($outfasta, $outmeta) = ('','');	
	
	#Êcreate FASTA sequence for seqreport
	while (my ($seq_idx, $seqrep) = each(%{$seqreport}) )
	{
		if ( exists $seqrep->{'iden'} ) {
			my ($seprep_iden) = $seqrep->{'iden'};
			#Êcreate list of transc identifiers from source db
			my ($seqrep_transc_id) = '';
			my ($seqrep_gene_name) = '';
			my ($seqrep_ccds_id) = '';
			my ($seqrep_ccds);
			foreach my $source (@{$SOURCE_LIST}) {
				my ($sr_tid) = '';
				if ( exists $seprep_iden->{$source} ) {
					my ($seprep_sc) = $seprep_iden->{$source};
					foreach my $gene_id ( sort {$a cmp $b} keys(%{$seprep_sc}) ) {
						foreach my $transc_id ( sort {$a cmp $b} keys(%{$seprep_sc->{$gene_id}->{'transc'}}) ) {
							my ($seprep_tr) = $seprep_sc->{$gene_id}->{'transc'}->{$transc_id};
							$sr_tid .= $transc_id.',';
							if ( exists $seprep_tr->{'ccds_id'} ) {
								my ($c) = $seprep_tr->{'ccds_id'};
								$c =~ s/\.[0-9]*$//g;
								$seqrep_ccds->{$c} = 1;
							}						
						}
					}
				}
				else { $sr_tid = '?' }
				$sr_tid =~ s/\,$//g;		
				$seqrep_transc_id .= $source.':'.$sr_tid.'+';
			}
			$seqrep_transc_id =~ s/\+$//g;
			$seqrep_ccds_id .= join(',', sort { $a cmp $b} keys(%{$seqrep_ccds}) ) if ( defined $seqrep_ccds );
			$outfasta .= '>appris'.'|'.
						$seq_idx.'|'.
						$seqreport_ids->{'gene_idx'}.' '.
						'gene_ids>'.$seqreport_ids->{'gene_id'}.' '.
						'gene_names>'.$seqreport_ids->{'gene_name'}.' '.
						'transc_ids>'.$seqrep_transc_id.' '.
						'ccds_ids>'.$seqrep_ccds_id."\n".
						$seqrep->{'seq'}."\n";
			$outmeta .= $xref_seqid."\t".
						$seq_idx."\t".
						$seqreport_ids->{'gene_idx'}."\t".
						$seqreport_ids->{'gene_id'}."\t".
						$seqreport_ids->{'gene_name'}."\t".
						$seqrep_transc_id."\t".
						$seqrep_ccds_id."\n";			
		}
		
		#Êcreate list of identifiers from source db
		my ($seqrep_transc_id) = '';
		my ($seqrep_gene_name) = '';
		my ($seqrep_ccds_id) = '';
	}
	
	return ($outfasta, $outmeta);
}

main();


1;

__END__

=head1 NAME

geneseqs

=head1 DESCRIPTION

Create APPRIS fasta file from the Cross Reference Ensembl vs UniProt vs RefSeq 

=head1 SYNOPSIS

geneseqs

=head2 Required arguments:

	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl geneseqs.pl \

	-x  appris/workspaces/cmpEnsUniRef/xref/xref.ens84_rs107_UPentry-gene.txt \
	
	-e  features/homo_sapiens/e84_g24/Homo_sapiens.GRCh38.pep.all.extra.fa \
	
	-r  features/homo_sapiens/rs107/protein.extra.fa \
	
	-u  features/homo_sapiens/up201606/uniprot-proteome.extra.fasta \
	
	-o  features/appris_seqs.e84_rs107_up201606.transl.fa \
	
	--loglevel=debug --logfile=createAPPRISseqsfromXREF.log		

or

perl geneseqs.pl \

	-x  features/mus_musculus/a1/xref_biomart.tsv \
	
	-e  features/mus_musculus/e87/mus_musculus.transl.extra.fa \
	
	-r  features/mus_musculus/rs106/protein.extra.fa \
	
	-u  features/mus_musculus/up201610/uniprot-proteome.extra.fasta \
	
	-o  features/mus_musculus/a1/appris_seqs.transl.fa \
	
	--loglevel=debug --logfile=createAPPRISseqsfromXREF.log		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
