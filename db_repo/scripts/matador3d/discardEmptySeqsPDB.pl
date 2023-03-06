#!/usr/bin/perl -W

use strict;
use Bio::SeqIO;

my ($fasta_file) = $ARGV[0];

my $fasta_object = Bio::SeqIO->new(
                        -file => $fasta_file,
                        -format => 'Fasta'
);
while ( my $seq = $fasta_object->next_seq() )
{
	my ($sequence_id) = $seq->id;
	my ($sequence_desc) = $seq->desc;
	my ($sequence) = $seq->seq;
	
	if ( length($sequence) != 0 ) {
		print STDOUT ">$sequence_id $sequence_desc\n$sequence\n";		
	}
}

exit 0;