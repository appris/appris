#!/usr/bin/perl -W

use strict;
use warnings;

use APPRIS::Registry;

sub main()
{
	# APPRIS Registry
	# The Registry system allows to tell your programs where to find the APPRIS databases and how to connect to them.
	# The following call will load all the local APPRIS database that was downloaded previously
	my($registry)=APPRIS::Registry->new();
	
	# By means of db configuration
	$registry->load_registry_from_db
	(
		-dbhost	=> 'jabba.cnio.es',
		-dbname	=> 'homo_sapiens_encode_3c',
		-dbuser	=> 'ensembl',
		-dbpass  => ''		
	);
	
	# Retrieve CDS positions
	my ($chromosome) = $registry->fetch_basic_by_region('M');
	foreach my $gene (@{$chromosome}) {
		foreach my $transcript (@{$gene->transcripts}) {			
				
			if ($transcript->translate and $transcript->exons) {
				my ($exons) = $transcript->exons;
				my ($translate) = $transcript->translate;
				for(my $icds=0;$icds<scalar(@{$translate->cds});$icds++) {
					my ($cds) = $translate->cds->[$icds];
					my ($exon) = $exons->[$icds];

					print 	"gene_id: ". $gene->stable_id. "transcript_id: ". $transcript->stable_id ."exon_id: ". $exon->stable_id .
							" cds_start: ".$cds->start." cds_end: ".$cds->end." cds_strand: ".$cds->strand."\n";

				}
			}
		}
	}
}

main();


1;
