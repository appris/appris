#!/usr/bin/perl

open(LS,"listaguia");
while(<LS>){
	chomp $_;
	open(PN,">lista.pn");
	open(SN,">lista.sn");
	print PN "$_.chk";
	print SN "$_.chd";
	close PN;close SN;
	`makemat -P lista`
	}
close LS;
