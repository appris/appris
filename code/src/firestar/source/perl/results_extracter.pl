#!/usr/bin/perl

use strict;
use FindBin;
my $cwd;
BEGIN{
	$cwd=$FindBin::Bin;
}
use lib "$cwd/lib";
use Getopt::Long;
use results_parser;

my $info=results_parser->new();

my ($dir_target)=undef;
my ($reliability)=undef;
my ($tag)=undef;
my ($output)=undef;
my ($csa)=undef;
my ($binding)=undef;


&GetOptions(
        'r=i'           => \$reliability,
        'd=s'           => \$dir_target,
        't=s'           => \$tag,
	'c=s'		=> \$csa,
	'b=s'		=> \$binding
);
unless(defined $reliability){$reliability=0;}
unless(defined $tag){$tag="ALL";}
unless(defined $output){$output="text";}
unless(defined $csa){$csa="YES";}
unless(defined $binding){$binding="YES";}
uc($tag);uc($csa);uc($binding);

unless (defined $dir_target and -d $dir_target){help();print STDERR "\n\nERROR: I can't find the target directory\n\n";exit;}
unless ($reliability >=0 and $reliability<=100){help();print STDERR "\n\nERROR: The reliability cut-off has to be between 0 and 100\n\n";exit;}
unless ($tag eq "COGNATE" or $tag eq "POSSIBLE" or $tag eq "NON" or $tag eq "ALL"){
	help();
	print STDERR "\n\nERROR: in FireDB the available TAGS are COGNATE-POSSIBLE-NON (or ALL if you don't want to filter)\n\n";
	exit;
}
unless ($output eq "text" or $output eq "CAFA"){help();print STDERR "\n\nERROR: your output option is not valid\n\n";exit;}

$info->{dir_target}=$dir_target;

my $fail_counter=0;
opendir (D,$dir_target);
my @list_file=readdir(D);
foreach my $i(@list_file){
	if ($i=~/^\.+/){next;}
	if ($i=~/(.+)\.res/){
		my $name=$1;
		$info->parse($name);
		my %already_seen;
		foreach my $i(sort { $a <=> $b } keys%{$info->{$name}{SITE}}){
			if ($info->{$name}{SITE}{$i}{reli}>=60 and $info->{$name}{SITE}{$i}{tag} eq "COGNATE"){
				if ($info->{$name}{SITE}{$i}{GO_terms} ne ''){
					my @new_gos=split(' ',$info->{$name}{SITE}{$i}{GO_terms});
					foreach my $k(@new_gos){
						$k="GO:$k";
						unless(exists $already_seen{$k}){print "$name\t$k\tHAS\t$i\t$info->{$name}{SITE}{$i}{reli}\n";$already_seen{$k}=5;}
					}
				}
			}
		}
		delete $info->{$name};
	}
	else{$fail_counter++;}
#	print "Done $i NEXT ...";
#	<STDIN>
}
if (scalar@list_file == $fail_counter){print STDERR "The program can't find any result file of firestar (.res)\n\n";exit;}

sub help{
	my $usage ="This program parses and outputs firestar results. Usage:\n\n\t\tresults_extracter.pl -d=<directory_where_the_results_are_stored> <options>\n\n";
	$usage.="\t-d\t\tthe path of the directory where the results are;\n";
	$usage.="\t-r\t\treliability score cut-off [0-100]; all the results below the cut-off will not be shown (default=0)\n";
	$usage.="\t-t\t\tcompound biological relevance [COGNATE-POSSIBLE-NON-ALL]. Only the sites cointaining this type of ligands will be shown (default=ALL)\n";
	$usage.="\t-c\t\tcsa information [YES-NO]. Information from CSA predictions will be shown (default=YES)\n";
	$usage.="\t-b\t\tbinding sites information [YES-NO]. Information from binding sites predictions will be shown (default=YES)\n";
	print STDERR $usage."\n";
}


