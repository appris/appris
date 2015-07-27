#!/usr/unic/bin/perl5

## Version for TargetP

## Usage: getfasta.pl indexfile fastafile

$indexfile=$ARGV[0];
$fastafile=$ARGV[1];

open (INDEXFILE,$indexfile);
while (<INDEXFILE>) {
	chomp;
	($name_temp)=split;
	my $name= substr($name_temp,0,30);
	$isthere{$name}=1;
}
close INDEXFILE;

open (FASTAFILE,$fastafile);
while (<FASTAFILE>) {
#  chomp;
  if (substr($_,0,1) eq ">") {
    chomp;
    ($field1)=split;
    $oldname=$name;
    $name = substr($field1,1,30);
    if ($isthere{$oldname}) {
	print STDOUT ">",$oldname,"\n",  substr($seqoneline{$oldname},0,180),"\n";
    }
  }
  else {
    chomp;
    $seqoneline{$name}.= $_ if $isthere{$name};
  }
}
close FASTAFILE;

# last entry:
if ($isthere{$name}) {
        print STDOUT ">",$name,"\n",  substr($seqoneline{$name},0,180),"\n";
}

