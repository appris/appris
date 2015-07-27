#!/usr/bin/perl -w
#to get "mrna" coordinates of the encode loci, get the lowest and  highest coordinate
#usage: perl get_locus_mrna.pl <gff_file>
use strict;

my(@loci)=();
my($i)=0;
my(@gff)=();
my($readline)="";
my($start)=0;
my($stop)=0;
my($grepline)= "";

if(scalar(@ARGV) != 1){
    print STDERR "usage: perl get_locus_mrna.pl <gff_file>","\n";
    exit(1);
}

open(GFF,$ARGV[0]) or die;
while($readline=<GFF>){
	chomp($readline);
	push(@gff,[split(/\t/,$readline)]);
}
#if ( $gff[0][6] eq '-') {
#	$start=$gff[$#gff][3];
#	$stop=$gff[0][4];
#} else {
#	$start=$gff[0][3];
#	$stop=$gff[$#gff][4];    	
#}
$start=$gff[0][3];
$stop=$gff[$#gff][4];
print "$gff[0][0]\t$gff[0][1]\ttranscript_coords\t$start\t$stop\t$gff[0][5]\t$gff[0][6]\t$gff[0][7]\t$gff[0][8]\n";