#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Bio::SeqIO;
use LWP::UserAgent;
use Data::Dumper;

use APPRIS::Utils::File qw( prepare_workspace rm_dir getTotalStringFromFile getStringFromFile printStringIntoFile );
use APPRIS::Utils::Logger;
use APPRIS::Utils::Exception qw( info throw warning deprecate );

###################
# Global variable #
###################
use vars qw(
	$TMPDIR
	$TAXID
	$FTP_NAME_UP_IDMAP
);
$TMPDIR = '/tmp/appris_xref';
$TAXID = {
	'hsapiens'	=> 9606,
	'mmusculus'	=> 10090,
	'drerio'	=> 7955 
};
$FTP_NAME_UP_IDMAP = {
	'hsapiens'	=> 'HUMAN_9606',
	'mmusculus'	=> 'MOUSE_10090',
	'drerio'	=> 'DANRE_7955' 
};

# Input parameters
my ($str_params) = join "\n", @ARGV;
my ($species) = undef;
my ($refseq_hgnc_file) = undef;
my ($refseq_ens_file) = undef;
my ($refseq_uni_file) = undef;
my ($uniprot_file) = undef;
my ($outfile) = undef;
my ($logfile) = undef;
my ($logpath) = undef;
my ($logappend) = undef;
my ($loglevel) = undef;
&GetOptions(
	'species|s=s'		=> \$species,
	'refg|rg=s'			=> \$refseq_hgnc_file,
	'refe|re=s'			=> \$refseq_ens_file,
	'refu|ru=s'			=> \$refseq_uni_file,
	'uni|u=s'			=> \$uniprot_file,
	'outfile|o=s'		=> \$outfile,
	'loglevel=s'		=> \$loglevel,
	'logfile=s'			=> \$logfile,
	'logpath=s'			=> \$logpath,
	'logappend'			=> \$logappend,
);

# Required arguments
unless (
	defined $species and 
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
	# variables
	my ($gnXref_report);
	my ($up_report);
	
	# prepare workspace
	$logger->info("-- prepare workspace\n");
	rm_dir($TMPDIR);
	prepare_workspace($TMPDIR);	
	
	# Create xref report
	$logger->info("-- create cross reference data based on Ensembl \n");
	create_xref_ensembl(\$gnXref_report, \$up_report);
#	$logger->debug("gnXREF_ENS_REPORT:\n".Dumper($gnXref_report)."\n");
	
	$logger->info("-- create cross reference data based on RefSeq \n");
	create_xref_refseq(\$gnXref_report, \$up_report);
#	$logger->debug("gnXREF_REF_REPORT:\n".Dumper($gnXref_report)."\n");
	$logger->debug("UNI_REPORT:\n".Dumper($up_report)."\n");
	
	$logger->info("-- create cross reference data based on UniProt \n");
	create_xref_uniprot(\$gnXref_report, $up_report);
#	$logger->debug("gnXREF_UNI_REPORT:\n".Dumper($gnXref_report)."\n");

	# Create text report	
	$logger->info("-- create text report \n");
	my ($outcontent) = create_txt($gnXref_report);
	my ($printing_file_log) = printStringIntoFile($outcontent, $outfile);

	$logger->finish_log();
	
	exit 0;
}

sub export_xref_biomart($)
{
	my ($species) = @_;
	my ($report) = '';
	
	my $xml = '<?xml version="1.0" encoding="UTF-8"?>
	<!DOCTYPE Query>
	<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1">
	
	        <Dataset name = "'.$species.'_gene_ensembl" interface = "default" >
	                <Attribute name = "external_gene_name" />
	                <Attribute name = "ensembl_gene_id" />
	                <Attribute name = "entrezgene" />
	                <Attribute name = "uniprot_swissprot" />
	                <Attribute name = "uniprot_sptrembl" />
	        </Dataset>
	</Query>';
	
	my $path="http://www.ensembl.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
	my $ua = LWP::UserAgent->new;
	
	my $response;
	$ua->request($request, 
		     sub{   
			 my($data, $response) = @_;
			 if ($response->is_success) {
			 	$report .= $data;
			 }
			 else {
			     warn ("Problems with the web server: ".$response->status_line);
			 }
		     },1000);
	
	if ( $report =~ /\[success\]\n*$/ ) {
		$report =~ s/\[success\]\n*$//g;
	}
	else {
		warn ("Problems retrieving the data: It is not complete");
	}
	return $report;
}

sub create_xref_ensembl(\$\$)
{
	my ($ref_report, $ref_unirepot) = @_;
	my ($report);
	my ($file) = $TMPDIR.'/'.'xref.en.txt';
	
	#ÊDownlod xref from biomart
	eval {
		my ($cmd) = "perl $ENV{APPRIS_SCRIPTS_DIR}/create_xref_biomart.pl mmusculus > $file";
		system($cmd);
	};
	throw('Exporting xref.en from biomart') if ($@);
	
	# Associated Gene Name	Gene ID	EntrezGene ID	UniProt/SwissProt Accession	UniProt/TrEMBL Accession
	my ($flines) = getTotalStringFromFile($file);
	foreach my $fline ( @{$flines} ) {
		my (@cols) = split('\t', $fline);
		my ($gene_name) = $cols[0];
		my ($ensembl_id) = $cols[1];
		my ($entrez_id) = $cols[2];
		my ($uniprot_list_ids) = join(',', @cols[3 .. $#cols]);
		$gene_name =~ s/\s*//g;
		$ensembl_id =~ s/\s*//g; $ensembl_id =~ s/\.[0-9]$//g;
		$entrez_id =~ s/\s*//g;
		$uniprot_list_ids =~ s/\s*//g;
		$uniprot_list_ids =~ s/^\,//g; $uniprot_list_ids =~ s/\,$//g;
		
		if ( defined $gene_name and $gene_name ne '' ) {
			if ( !defined $ensembl_id or $ensembl_id eq '' ) {
				$ensembl_id = '-';
			}
			if ( !defined $entrez_id or $entrez_id eq '' ) {
				$entrez_id = '-';
			}			
			if ( !defined $uniprot_list_ids or $uniprot_list_ids eq '' ) {
				$uniprot_list_ids = '-';
			}
			foreach my $uniprot_id ( split(',', $uniprot_list_ids) ) {
				if ( $uniprot_id ne '' and $uniprot_id ne '-' ) {
					my ($label) = 'ENSEMBL';
					push( @{$$ref_report->{$gene_name}->{$ensembl_id}->{$entrez_id}->{$uniprot_id}}, $label );
					if ( $ensembl_id ne '' and $ensembl_id ne '-' ) {
						$$ref_unirepot->{$uniprot_id}->{'ensembl'}->{$ensembl_id} = $label;
					}
					if ( $entrez_id ne '' and $entrez_id ne '-' ) {
						$$ref_unirepot->{$uniprot_id}->{'refseq'}->{$entrez_id} = $label;
					}
				}
			}
		}
	}
}
sub create_xref_refseq(\$\$)
{
	my ($ref_report, $ref_unirepot) = @_;
	my ($report);
#	my ($refseq_hgnc_file) = $TMPDIR.'/'.'xref.rs.hgnc.txt';
	my ($refseq_ref_file) = $TMPDIR.'/'.'xref.rs.refseq.txt';
	my ($refseq_ens_file) = $TMPDIR.'/'.'xref.rs.ensembl.txt';
	my ($refseq_uni_file) = $TMPDIR.'/'.'xref.rs.uniprot.txt';
	my ($flines);
	my ($taxid) = $TAXID->{$species};
	
#	# Get HGNC
#	eval {
#		my ($cmd) = "wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz -P $TMPDIR && gzip -d $TMPDIR/gene_info.gz && grep '^$taxid' $TMPDIR/gene_info | cut -f 1,2,3,6 > $refseq_hgnc_file";
##		system($cmd);
#	};
#	throw('Exporting xref.rs for hgnc') if ($@);

	# Extract the Cross-Reference with RefSeq
	eval {
		my ($cmd) = "wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz -P $TMPDIR && gzip -d $TMPDIR/gene2refseq.gz && grep '^$taxid' $TMPDIR/gene2refseq | cut -f 1,2,3,4,6,16 > $refseq_ref_file";
		system($cmd);
	};
	throw('Exporting xref.rs for refseq') if ($@);
 
	# Extract the Cross-Reference with Ensembl
	eval {
		my ($cmd) = "wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz -P $TMPDIR && gzip -d $TMPDIR/gene2ensembl.gz && grep '^$taxid' $TMPDIR/gene2ensembl > $refseq_ens_file";
		system($cmd);
	};
	throw('Exporting xref.rs for ensembl') if ($@);

	# Extract RefSeq proteins (NP_ & YP_ & YP_)
	eval {
		my ($cmd) = "wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz -P $TMPDIR && gzip -d $TMPDIR/gene_refseq_uniprotkb_collab.gz && grep -e '^NP_\\|^XP_\\|^YP_' $TMPDIR/gene_refseq_uniprotkb_collab  > $refseq_uni_file";
		system($cmd);
	};
	throw('Exporting xref.rs for uniprot') if ($@);
	

	# Get Gene Symbols
#	$flines = undef; $flines = getTotalStringFromFile($refseq_hgnc_file);
#	foreach my $fline ( @{$flines} ) {
#		my (@cols) = split('\t', $fline);
#		my ($tax_id) = $cols[0];
#		my ($entrez_id) = $cols[1];
#		my ($gene_name) = $cols[2];
#		my ($mult_ids) = $cols[3];
#		my ($ensembl_id) = ( $mult_ids =~ /Ensembl:([^\|]*)/ ) ? $1 : ''; 
#		$entrez_id =~ s/\s*//g;
#		$gene_name =~ s/\s*//g;
#		$ensembl_id =~ s/\s*//g;
#		
#		if ( defined $gene_name and $gene_name ne '' ) {
#			if ( defined $entrez_id and $entrez_id ne '' ) {
#				if ( !exists $report->{$entrez_id}->{'hgnc'}->{$gene_name} ) {
#					$report->{$entrez_id}->{'hgnc'}->{$gene_name} = 1;					
#				}
#			}
#		}		
#	}
	# Get RefSeq Xref
	my ($for_uniprot);	
	$flines = undef; $flines = getTotalStringFromFile($refseq_ref_file);
	#10090	11287	VALIDATED	NM_007376.4	NP_031402.3	Pzp
	#10090	11298	MODEL	XM_017314223.1	XP_017169712.1	Aanat
	#10090	11298	MODEL	XM_017314224.1	XP_017169713.1	Aanat
	for (my $i=0; $i < scalar(@{$flines}); $i++ ) {
		my (@cols) = split('\t', $flines->[$i]);
		my ($tax_id) = $cols[0];
		my ($entrez_id) = $cols[1];
		my ($status) = $cols[2];
		my ($ref_transc_id) = $cols[3];
		my ($ref_transl_id) = $cols[4];
		my ($gene_name) = $cols[5];
		$entrez_id =~ s/\s*//g;
		$status =~ s/\s*//g;
		$ref_transc_id =~ s/\s*//g;
		$ref_transl_id =~ s/\s*//g;
		$ref_transl_id =~ s/\.[0-9]*$//g;
		$gene_name =~ s/\s*//g;

		if ( defined $entrez_id and $entrez_id ne '' ) {
			if ( defined $gene_name and $gene_name ne '' ) {
				if ( !exists $report->{$entrez_id}->{'hgnc'}->{$gene_name} ) {
					$report->{$entrez_id}->{'hgnc'}->{$gene_name} = 1;					
				}
			}
			if ( defined $ref_transl_id and $ref_transl_id ne '' and $ref_transl_id ne '-' ) {
				$for_uniprot->{$ref_transl_id} = $entrez_id;
			}			
		}
	}
	# Get Ensembl Xref
	$flines = undef; $flines = getTotalStringFromFile($refseq_ens_file);
	#9606	2	ENSG00000175899	NM_000014.4	ENST00000318602	NP_000005.2	ENSP00000323929
	#9606	3	ENSG00000256069	NR_040112.1	ENST00000543404	-	-
	#9606	9	ENSG00000171428	NM_000662.7	ENST00000307719	NP_000653.3	ENSP00000307218
	for (my $i=0; $i < scalar(@{$flines}); $i++ ) {
		my (@cols) = split('\t', $flines->[$i]);
		my ($tax_id) = $cols[0];
		my ($entrez_id) = $cols[1];
		my ($gene_id) = $cols[2];
		my ($ref_transc_id) = $cols[3];
		my ($transc_id) = $cols[4];
		my ($ref_transl_id) = $cols[5];
		$entrez_id =~ s/\s*//g;
		$gene_id =~ s/\s*//g;
		$ref_transc_id =~ s/\s*//g;
		$transc_id =~ s/\s*//g;
		$ref_transl_id =~ s/\s*//g;
		$ref_transl_id =~ s/\.[0-9]*$//g;
		
		if ( defined $entrez_id and $entrez_id ne '' ) {
			if ( defined $gene_id and $gene_id ne '' ) {
				if ( !exists $report->{$entrez_id}->{'ens'}->{$gene_id} ) {
					$report->{$entrez_id}->{'ens'}->{$gene_id} = 1;					
				}
			}
		}
	}	
	# Get UniProt Xref
	$flines = undef; $flines = getTotalStringFromFile($refseq_uni_file);
	#NP_000005       P01023
	#NP_000006       A4Z6T7
	#NP_000006       P11245
	for (my $i=0; $i < scalar(@{$flines}); $i++ ) {
		my (@cols) = split('\t', $flines->[$i]);
		my ($ref_transl_id) = $cols[0];
		my ($uniprot_acc) = $cols[1];
		$ref_transl_id =~ s/\s*//g;
		$ref_transl_id =~ s/\.[0-9]*$//g;
		$uniprot_acc =~ s/\s*//g;
		
		if ( defined $ref_transl_id and $ref_transl_id ne '' ) {
			if ( defined $uniprot_acc and $uniprot_acc ne '' ) {
				if ( exists $for_uniprot->{$ref_transl_id} ) {
					my ($entrez_id) = $for_uniprot->{$ref_transl_id};
					$report->{$entrez_id}->{'uni'}->{$uniprot_acc} = 1;					
				}
			}
		}		
	}	
	# Add to global report
	while ( my ($entrez_id,$rep) = each(%{$report}) ) {
		my (@hgnc);
		if ( exists $rep->{'hgnc'} ) { @hgnc = keys(%{$rep->{'hgnc'}}) } # required
		my (@ens);
		if ( exists $rep->{'ens'} ) { @ens = keys(%{$rep->{'ens'}}) } # optional
		else { push(@ens, '-') }
		my (@uni);
		if ( exists $rep->{'uni'} ) { @uni = keys(%{$rep->{'uni'}}) } # optional 
		else { push(@uni, '-') }		
		my (@keys);
		foreach my $gene_name (@hgnc) {			
			foreach my $ensembl_id (@ens) {
				foreach my $uniprot_id (@uni) {
					if ( $uniprot_id ne '' and $uniprot_id ne '-' ) {
						my ($label) = 'REFSEQ';
						push( @{$$ref_report->{$gene_name}->{$ensembl_id}->{$entrez_id}->{$uniprot_id}}, $label );
						if ( $ensembl_id ne '' and $ensembl_id ne '-' ) {
							$$ref_unirepot->{$uniprot_id}->{'ensembl'}->{$ensembl_id} = $label;
						}
						if ( $entrez_id ne '' and $entrez_id ne '-' ) {
							$$ref_unirepot->{$uniprot_id}->{'refseq'}->{$entrez_id} = $label;
						}
					}
				}					
			}				
		}
	}
}
sub create_xref_uniprot(\$$)
{
	my ($ref_report, $unirepot) = @_;
	my ($report);
	my ($for_uniprot);
	my ($file) = $TMPDIR.'/'.'xref.up.txt';
	my ($filename_idmapping) = $FTP_NAME_UP_IDMAP->{$species};
	
	# Extract Xref from UniProt	
	eval {
		my ($cmd) = "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/$filename_idmapping\_idmapping.dat.gz -P $TMPDIR && gzip -d $TMPDIR/$filename_idmapping\_idmapping.dat.gz ".
					"grep -e 'Gene_Name'       $TMPDIR/$filename_idmapping\_idmapping.dat >  $file && ".
					"grep -e 'Ensembl'         $TMPDIR/$filename_idmapping\_idmapping.dat >> $file && ".
					"grep -e 'GeneID\\|RefSeq' $TMPDIR/$filename_idmapping\_idmapping.dat >> $file && ".
					"grep -e 'CCDS'            $TMPDIR/$filename_idmapping\_idmapping.dat >>  $file";
		system($cmd);
	};
	throw('Exporting xref.rs for refseq') if ($@);	

	# Get Gene Symbol,Ensembl,RefSeq
	my ($flines) = getTotalStringFromFile($file);
	foreach my $line (@{$flines}) {
		my (@cols) = split('\t', $line);
		my ($uniprot_acc) = $cols[0];
		my ($ref_type) = $cols[1];
		my ($ref_id) = $cols[2];
		$uniprot_acc =~ s/\s*//g;
		$ref_type =~ s/\s*//g;		
		$ref_id =~ s/\s*//g;
		$ref_id =~ s/\.[0-9]*$//g;
		if ( defined $uniprot_acc and defined $ref_type and defined $ref_id and $ref_id ne '') {
			if ( $ref_type eq 'Gene_Name' ) {
				$for_uniprot->{$uniprot_acc} = $ref_id;
				$report->{$uniprot_acc}->{'hgnc'}->{$ref_id} = 1;
			}
			elsif ( $ref_type eq 'Ensembl' ) {				
				$report->{$uniprot_acc}->{'ens'}->{$ref_id} = 1;				
			}
			elsif ( $ref_type eq 'GeneID' ) {				
				$report->{$uniprot_acc}->{'ref'}->{$ref_id} = 1;
			}	
		}
	}
	# Add to global report
	while ( my ($uniprot_id,$rep) = each(%{$report}) ) {
		if ( exists $for_uniprot->{$uniprot_id} ) {
			my (@hgnc);
			if ( exists $rep->{'hgnc'} ) { @hgnc = keys(%{$rep->{'hgnc'}}) } # required
			my (@ens);
			if ( exists $rep->{'ens'} ) { @ens = keys(%{$rep->{'ens'}}) } # optional
			else { push(@ens, '-') }
			my (@ref);
			if ( exists $rep->{'ref'} ) { @ref = keys(%{$rep->{'ref'}}) }
			else { push(@ref, '-') }
			my (@keys);
			foreach my $gene_name (@hgnc) {
				foreach my $ensembl_id (@ens) {
					foreach my $entrez_id (@ref) {
						my ($label) = 'UNIPROT';
						if ( $ensembl_id eq '-' ) {
							foreach my $ensembl_id ( keys(%{$unirepot->{$uniprot_id}->{'ensembl'}}) ) {
								if ( $entrez_id eq '-' ) {
									foreach my $entrez_id ( keys(%{$unirepot->{$uniprot_id}->{'refseq'}}) ) {
										push( @{$$ref_report->{$gene_name}->{$ensembl_id}->{$entrez_id}->{$uniprot_id}}, $label );
									}
								}
								else {
									push( @{$$ref_report->{$gene_name}->{$ensembl_id}->{$entrez_id}->{$uniprot_id}}, $label );
								}								
							}
						}
						else {
							if ( $entrez_id eq '-' ) {
								foreach my $entrez_id ( keys(%{$unirepot->{$uniprot_id}->{'refseq'}}) ) {
									push( @{$$ref_report->{$gene_name}->{$ensembl_id}->{$entrez_id}->{$uniprot_id}}, $label );
								}
							}
							else {
								push( @{$$ref_report->{$gene_name}->{$ensembl_id}->{$entrez_id}->{$uniprot_id}}, $label );
							}							
						}
					}
				}				
			}
		}
	}
}
sub create_txt($)
{
	my ($report) = @_;
	my ($output) = 'HGNC'."\t".
					'ENsemb_ID'."\t".
					'RefSeq_ID'."\t".
					'UniProt_IDs'."\n";
	my (@genes) = keys( %{$report} );
	foreach my $gene_name ( sort {$a cmp $b} @genes ) {
		foreach my $ensembl_id ( sort {$a cmp $b} keys(%{$report->{$gene_name}}) ) {
			foreach my $entrez_id ( sort {$a cmp $b} keys(%{$report->{$gene_name}->{$ensembl_id}}) ) {
				my (@uniprot_list_ids) = sort {$a cmp $b} keys(%{$report->{$gene_name}->{$ensembl_id}->{$entrez_id}});
				my ($uniprot_ids) = join(',', @uniprot_list_ids);
				$output .= $gene_name."\t".$ensembl_id."\t".$entrez_id."\t".$uniprot_ids."\n";				
			}
		}
	}
	return $output;	
}


main();


1;

__END__

=head1 NAME

geneset

=head1 DESCRIPTION

Create Cross Reference HGNC vs Ensembl vs UniProt vs RefSeq 

=head1 SYNOPSIS

cmpEnsUniRef

=head2 Required arguments:

	-s,   --species= <Species name>

	-e,   --ens= <Ensembl - Biomart - >
	
	-r,   --ref= <RefSeq xref file>

	-u,   --uni= <UniProt xref file>
		
	-o,   --outfile= <Output file>	
	
=head2 Optional arguments:

	--loglevel=LEVEL <define log level (default: NONE)>	

	--logfile=FILE <Log to FILE (default: *STDOUT)>
	
	--logpath=PATH <Write logfile to PATH (default: .)>
	
	--logappend <Append to logfile (default: truncate)>


=head1 EXAMPLE

perl geneset.pl

	-s    mmusculus

	-e    xref/xref.ens84_refseq.biomart.PC.txt
	
	-r    xref/gene_RefSeqGene
	
	-u    xref/xref.uni201602_ens_refseq.txt
		
	-o    xref.ens84_rs107_uni201602.txt
	
	--loglevel=debug --logfile=cmpEnsUniRef.log
		

=head1 AUTHOR

Jose Manuel Rodriguez Carrasco -jmrodriguez@cnio.es- (INB-GN2,CNIO)

=cut
