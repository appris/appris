#!/usr/bin/perl

use strict;
use FindBin;
my $cwd=$FindBin::Bin;

use Config::IniFiles;
my $variables=Config::IniFiles->new(-file => "$cwd/../CONFIG_fire_var.ini");

my $home=$variables->val('PATHS','home');
my $faatmp=$home."/tmp/faatmp";

my %templates;
my %sequences;

unless(-e "$faatmp/FAA_LOG.txt"){die "Acá no está el fichero !!!!!\n";}
open(F,"$faatmp/FAA_LOG.txt");
while(<F>){
	chomp($_);
	my @lista=split(/\t/,$_);
	$templates{$lista[0]}=$lista[1];
}
close(F);

deleter();

foreach my $i(keys%templates){
	unless (-e "$faatmp/$i.faa" and -e "$faatmp/$i\_10.psi" and -s "$faatmp/$i\_10.psi"){delete $templates{$i};}
}

deleter();

foreach my $i(keys%templates){
	my $seq=$templates{$i};
        if (exists $sequences{$seq}){
		delete $templates{$i};
	}
	else {$sequences{$seq}=5;}
}

deleter();

unlink("$faatmp/FAA_LOG.txt");
open(F,'>',"$faatmp/FAA_LOG.txt");
foreach my $i(keys%templates){
	print F "$i\t$templates{$i}\n";
}

sub deleter{
	opendir(DIR,$faatmp);
	while(readdir(DIR)){
		if($_=~/(\d+)(\_\d+)?\..+/){
			my $temporal=$1;
			unless(exists$templates{$temporal}){unlink("$faatmp/$_");}
		}
	}
	closedir(DIR);
}
