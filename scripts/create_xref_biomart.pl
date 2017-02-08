#!/usr/bin/perl -w

# an example script demonstrating the use of BioMart webservice
use strict;
use LWP::UserAgent;


#ÊGive the query by input
#open (FH,$ARGV[0]) || die ("\nUsage: perl webExample.pl Query.xml\n\n");
#
#my $xml;
#while (<FH>){
#    $xml .= $_;
#}
#close(FH);

########################################################
# USAGE
#
my $USAGE =<<USAGE;

Usage:

  perl create_xref_biomart.pl mmusculus [-help]
    or
  perl create_xref_biomart.pl drerio [-help]
         
Help: Prints out this helpful message

USAGE
#
######################################################

if ( !defined $ARGV[0] || $ARGV[0] eq '-h' || $ARGV[0] eq '-help' ) {
    print "$USAGE\n";
	exit 0;
}

my $species = $ARGV[0];

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
my $all_data = '';
$ua->request($request, 
	     sub{   
		 my($data, $response) = @_;
		 if ($response->is_success) {
		 	$all_data .= $data;
		 	#print "$data";
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	     },1000);

if ( $all_data =~ /\[success\]\n*$/ ) {
	$all_data =~ s/\[success\]\n*$//g;
	print "$all_data";
}
else {
	warn ("Problems retrieving the data: It is not complete");
}



