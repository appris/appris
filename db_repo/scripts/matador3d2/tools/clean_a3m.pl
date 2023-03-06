#!/usr/bin/perl -w
use strict;
use File::Basename;

my $input_file 	= $ARGV[0];
my $output_file = $ARGV[1];


if(!$input_file){
	die "ERROR: clean_a3m.pl script: An input file has to be specified."
}

if(!$input_file || !-e $input_file){
	die "ERROR: clean_a3m.pl script: The input file should exist."
}


open(FH, "$input_file") || die "Error while opening file $input_file";
open(FH_OUT, ">$output_file") || die "Error while opening file $output_file for writing";

my $pdb 	= basename($input_file);
my $print 	= 0;
while(my $line = <FH>){

	if($line =~ /^[A-Za-z-]/ && $print){
		print FH_OUT $line;
		$print = 0;
	}

	if($line =~ /^>ss_dssp/){
		$print = 0;
		next;
	}
	
	if($line =~ /^>ss/){
		$print = 0;
		next;
	}	
	
	if($line =~ /^>(\d(\w)+_\w)/){
		$print = 1;
		print FH_OUT $line;
		next;
	}
	
	if($line =~ /^>(sp|tr)/){
		print FH_OUT $line;
		$print = 1;
		next;
	}
}

close FH;
close FH_OUT;



exit;

