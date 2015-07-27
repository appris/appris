#!/usr/bin/perl -w
#
# use this to remove stop codons from an alignment
# typically, this would be done to calculate dN/dS in HYPHY
# 
# Usage: perl ReplaceStopWithGaps.pl -pep 104D5_pep.fasta -nuc 104D5.fasta -output 104D5_nostop.fasta
# use this to replace stop codons from the nucleotide alignment
# the nucleotide and the peptide alignments are necessary 


use strict;
use Getopt::Long; 
use Bio::SeqIO;
use Data::Dumper;

my ($inpep,$innuc,$output, $i, %stop);
&GetOptions(
	    'pep:s'      => \$inpep,#
	    'nuc:s'      => \$innuc,
	    'out:s'   => \$output,#file without gaps
           );


my $pep = Bio::SeqIO->new(-file => "$inpep" , '-format' => 'fasta');
my $nuc  = Bio::SeqIO->new(-file => "$innuc" , '-format' => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$output" , '-format' => 'fasta');

while ( my $pepseq = $pep->next_seq() ) {
    my $pep_str=uc($pepseq->seq);
    if ($pep_str=~/\*/){
      my $pep_id=$pepseq->id();
      my @aa=split(//,uc($pepseq->seq));
      for ($i=0; $i<scalar(@aa); $i++){
        if ($aa[$i]=~/\*/){
      		$stop{$pep_id}{$i}++;
      		print STDERR "$pep_id peptide sequence has a stop $aa[$i] at ".($i+1)."\n";
      	}
      }
    }
}

while (my $nucseq = $nuc->next_seq()){
  my $nuc_id=$nucseq->id();
  my $nuc_str=uc($nucseq->seq);
  foreach my $pid (keys %stop){

    if ("$nuc_id" eq "$pid"){
      foreach my $site (keys %{$stop{$pid}}){
		  print STDERR "match $nuc_id and $pid\n";
		  print STDERR "The sequence for $nuc_id is \n$nuc_str\n";
		  my $nucpos=$site*3;
		  my $codon =  substr $nuc_str, $nucpos, 3;
		  print STDERR "$codon ";
		  if ($codon =~ /(((U|T)A(A|G|R))|((T|U)GA))/i){
			substr($nuc_str, $nucpos, 3) = '---';
			print "=> Match to a stop codon at nucleotide position ".($nucpos+1)."\nNew sequence for $nuc_id\n$nuc_str\n";
		  }else{
			print "Doesn't seem to match a stop codon at nucleotide position ".($nucpos+1)." in $nuc_id\n";
		  }
      }
    }
  }
  my $newseq = Bio::Seq->new(-seq => "$nuc_str",                           
                         -display_id => $nuc_id);
  $out->write_seq($newseq); 
}

