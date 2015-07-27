#!/usr/bin/perl -w
use strict;
my ($line)="";
my (@gff)=();
my ($gff_count)=0;
my ($start)=0;
my ($stop)=0;
my ($help)=0;
my ($strand)="";
my ($id)="";
my ($i)=0;
my ($j)=0;
my (@neg)=();
my ($negs)=0;
my (@mrna)=();
my ($mrna_count)=0;
my ($k)=0;

open(GFF, "< $ARGV[0]") or die "cannot open $ARGV[0]: $!";
while(<GFF>){
    chomp;
    $line=$_;
    next if $line =~ /^\#/;
    push(@gff,[split(/\s+/,$line)]);
    $gff_count++;
}
close(GFF) or die "cannot close $ARGV[0]: $!";

open(MRNA, "$ARGV[1]") or die;
while(<MRNA>){
    chomp;
    $line=$_;
    next if $line =~ /^\#/;
    push(@mrna,[split(/\s+/,$line)]);
    $mrna_count++;

}


for($k=0;$k<$mrna_count;$k++){
	$id=$mrna[$k][8];
	$start= $mrna[$k][3];
	$stop= $mrna[$k][4];
	$strand=$mrna[$k][6];
	for($i=0;$i<$gff_count;$i++){
	  if(($gff[$i][8] eq $id)&&($gff[$i][6] eq "-")){#changed 26-11-07, was $gff[$i][0]
	    $gff[$i][3]-=$stop;
	    $gff[$i][3]=abs($gff[$i][3]);
	    $gff[$i][3]+=1;
	    $gff[$i][4]-=$stop;
	    $gff[$i][4]=abs($gff[$i][4]);
	    $gff[$i][4]+=1;
	    $help=$gff[$i][3];
	    $gff[$i][3]=$gff[$i][4];
	    $gff[$i][4]=$help;
	    $gff[$i][0]=$mrna[$k][0];#only 1-12
	    push(@neg,[$gff[$i][0],$gff[$i][1],$gff[$i][2],$gff[$i][3],$gff[$i][4],$gff[$i][5],$gff[$i][6],$gff[$i][7],$gff[$i][8]]);
	    $negs++;
	  }
	}
}
	for($j=0;$j<$#gff;$j++){
	    if($gff[$j][6] eq "+"){
		for($i=0;$i<8;$i++){
		    print STDOUT "$gff[$j][$i]\t";
		}
		print STDOUT "$gff[$j][8]\n";
	     }

	}
#neg strand records are printed with + because all coordinates are adjusted to + strand 

	for($i=0;$i<$negs;$i++){
	    $neg[$i][6]="+";
	}
	for($j=$negs-1;$j>=0;$j--){

	    for($i=0;$i<8;$i++){
		    print STDOUT "$neg[$j][$i]\t";
		}
		print STDOUT "$neg[$j][8]\n";
	}
