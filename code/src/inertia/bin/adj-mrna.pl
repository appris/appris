#!/usr/bin/perl
#get gff coordinates relative to cdna (i.e no introns, intron length 0)
#only give one set of gff (of 1 transcript) at the time
use strict;

my (@gff)=();
my ($line)="";
my ($i)=0;
my ($count)=0;
my ($diff)=0;
my ($j)=0;


open (IN, $ARGV[0]) or die;
while(<IN>){
    chomp;
    $line=$_;
    push(@gff,[split(/\s+/,$line)]);
    $count++;
}
close(IN)or die;

$diff=$gff[0][4]-$gff[0][3]+1;
$gff[0][3]=1;
$gff[0][4]=$diff;
for($i=1;$i<$count;$i++){
    $diff=$gff[$i][4]-$gff[$i][3]+1;
    $gff[$i][3]=$gff[$i-1][4]+1;
    $gff[$i][4]=$gff[$i][3]+$diff-1;
}

for($i=0;$i<$count;$i++){
    for($j=0;$j<9;$j++){
	print $gff[$i][$j]."\t";
    }
    print $gff[$i][9]."\n";
}
    
