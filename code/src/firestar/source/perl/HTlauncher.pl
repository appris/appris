#! /usr/bin/perl
# $Id: fire_launcher.pl
# Developed by: Paolo Maietta -pmaietta@cnio.es-
# Created:      25-Oct-2012

# This program is a simple connection between the user and firePredText.pl
# In order to pass all the parameters in the correct order to the program
# and to obtain the predictions. It can read a multi fasta file and returns
# single-file predictions
# _________________________________________________________________

use strict;
use Getopt::Long;
use FindBin;
my $cwd=$FindBin::Bin;
use Config::IniFiles;
my $variables=Config::IniFiles->new(-file => "$cwd/../CONFIG_fire_var.ini");


####################
# Input parameters #
####################

my ($file)=undef;
my ($out_dir)=undef;
my ($CSA)=undef;
my ($cut_off)=undef;
my ($cog)=undef;
my ($pro_opt)=undef;

&GetOptions(
        'f=s'           => \$file,
        'o=s'           => \$out_dir,
        'CSA=s'		=> \$CSA,
        'cut=s'		=> \$cut_off,
        'cog=s'		=> \$cog,
	'pro=s'		=> \$pro_opt
	
);

#########################################################
###############		VARIABLES	#################
#########################################################

my $home=$variables->val('PATHS','home');
my $temporal="$home/tmp/faatmp";
my $path_to_prof;

########################################################################
###############		INPUT PARAMETERS CONTROL	################
########################################################################

unless (defined $file && defined $out_dir){
	print "\n\n\tSorry, but some parameters are missing !!\n";
	print "Usage: fire_launcher.pl -f=[the name of input file] -o=[your output DIRECTORY]\n\n";
	die;
}

if (defined $cog){
	if ($cog!~/[YES|yes|Yes|No|NO|no]/){
		die "Sorry, but the possibilities for cog option are Yes or No; by default is set at Yes.\n\n";
	}
}
else{$cog="YES";}
if (defined $CSA){
	if ($CSA!~/[YES|yes|Yes|No|NO|no|ONLY|Only|only]/){
		die "Sorry, but the possibilities for CSA option are Yes, No or Only; by default is set at Yes.\n\n";
	}
}
else{$CSA="YES";}
if (defined $cut_off){
	if ($cut_off!~/\d+/ or $cut_off >100){
		die "Sorry, but the limits for cut option are 0-100; by default is set at NO.\n\n";
	}
}
else{$cut_off=0;}
unless (-d $out_dir){
	print "Output directory doesn't exist; Do you want to create it ? ";
	my $answer=<STDIN>;
	chomp($answer);
	if ($answer=~/[Y|Yes|YES|yes]/){mkdir $out_dir;}
	else{die;}
}
if (defined $pro_opt){
	if ($pro_opt !~ /[Y|YES|yes|Yes]/){
		die "Sorry the selected option for the profile is wrong; the only valid option is 'yes'.\n\n";
	}
	my @parking=split('/',$file);
	pop@parking;
	$path_to_prof=join('/',@parking);
}

############################################################################

my $name=undef;
my @cargador;
my $sequence='';
my $evalue=10;
my $tmpfname;
my @names;

open (F,$file);
while (<F>){
	if ($_=~/^>(.+)/){
		if (defined $name && defined $pro_opt){
			push(@cargador,"$name $sequence $evalue $out_dir/$name.res $cut_off $CSA $cog $path_to_prof/$name");
			$sequence='';
		}
		elsif (defined $name){
			push(@cargador,"-q=$name -s=$sequence -e=$evalue -o=$out_dir/$name.res -cut=$cut_off -csa=$CSA -cog=$cog");
			push(@names,$name);
			$sequence='';
		}
		my @chapalo=split(/\s+/,$1);
		$name=$chapalo[0];
		if (-e "$out_dir/$name.res"){$name=undef;}
	}
	else{
		if (defined $name){
			chomp($_);
			$sequence.=$_;
		}
	}
}

if (defined $name && defined $pro_opt){
	push(@cargador,"$name $sequence $evalue $out_dir/$name.res $cut_off $CSA $cog $path_to_prof/$name");
}
elsif (defined $name) {
	push(@cargador,"-q=$name -s=$sequence -e=$evalue -o=$out_dir/$name.res -cut=$cut_off -csa=$CSA -cog=$cog");
	push(@names,$name);
}

my $conteur=0;
foreach my $i (@cargador){
	my $flag="CLOSE";
	while ($flag eq "CLOSE"){
		my @processes=`ps aux | grep firestar | grep -v grep`;
		if ((scalar@processes/2) < 5){$flag="OPEN";}
		else {sleep(5);}
	}
	if (defined$pro_opt){`perl $cwd/firePredText_profile.pl $i >/dev/null 2>&1 &`;}
	else {`perl $cwd/firestar.pl $i > /dev/null 2> /tmp/$names[$conteur].err &`;}
	$conteur++;
}

