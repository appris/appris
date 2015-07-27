#! /usr/bin/perl -W

use strict;
use Getopt::Long;

####################
# Input parameters #
####################
my ($db_file) = undef;
my ($input_sequence) = undef;
my ($input_psiblast) = undef;
my ($output) = undef;

&GetOptions(
	'db=s'			=> \$db_file,
	'in-seq=s'		=> \$input_sequence,
	'in-psiblast=s'	=> \$input_psiblast,
	'output=s'		=> \$output,
);

###################
# Internal Values #
###################

my @fichero;
my @listado;
my @length_ali;
my @proteinas;
my @temporary;
my $aminoacids;
my $control;
my $line;
my $length;
my $length_seq;
my $perc;
my $perc_ali;
my $prot_name;
my $quantity=0;
my $reference;
my $score;
my $space_counter;
my $switch1="off";
my $i;
my $w;
my $z;
my $k;

##########
# PARSER #
##########

if ( -e $output and (-s $output > 0) ) {
	exit 0;
}


open (QUERY,"$input_sequence");
while (<QUERY>){
	if ($_=~/^>.+/){
		push(@fichero,$_);
	}
	elsif ($_=~/^([\s|\t]+)?[A-Z]+$/i) {
		$aminoacids=$aminoacids.$_;
	}
}
@temporary=split('',$aminoacids);
$aminoacids=undef;
foreach $i (@temporary){
	if ($i=~/[A-Z]/i){$aminoacids=$aminoacids.$i;$length_seq++;}
}
chomp($fichero[$#fichero]);
$fichero[$#fichero]=$fichero[$#fichero]."\t$length_seq\n";
push(@fichero,$aminoacids."\n");
close (QUERY);




open(FICH,"$input_psiblast");
@listado=<FICH>;
close(FICH);
for ($k=0;$k<=scalar@listado;$k++){
	if ($quantity > 28){last;} 
	if ($listado[$k]=~/^Results from round 2/){$switch1="on";}
	if ($listado[$k]=~/^>(.+)/ && $switch1 eq "on"){
		$score=undef;
		$aminoacids=undef;
		$perc=undef;
		@length_ali=undef;
		@temporary=undef;
		@temporary=split(/ /,$1);
		$prot_name=$temporary[0];	
		$space_counter=0;
		$z=$k+1;
		while ($z!=$k){
			if ($listado[$z]=~/^([\s|\t]+)?Length =\s+(\d+)/){
				$length=$2;$perc=($length/$length_seq)*100;
				$space_counter=0;
			}
			elsif ($listado[$z]=~/\s+Score =.+, Expect =\s+(.+),.+/){
				$space_counter=0;
				$score=$1;
				if ($score!~/[0\.0|\d?e\-\d+]/){$k=$z;last;}
			}
			elsif($listado[$z]=~/\s+Identities = (\d+)\/(\d+)\s+\(.+/){
				$space_counter=0;
				if ($1 == $2 && $1 == $length_seq){$k=$z;last;}
			}
			elsif ($listado[$z]=~/Sbjct:\s+\d+\s+([A-Z-]+)\s+\d+/i){
				$space_counter=0;
				@temporary=undef;
				@temporary=split('',$1);
				foreach $w(@temporary){
					if ($w=~/[A-Z]/i){$aminoacids=$aminoacids.$w;}
				}
			}
			elsif ($listado[$z]=~/^\n$/){$space_counter++;}
			else {$space_counter=0;}
			if (defined $aminoacids && $space_counter==2){
				@length_ali=split('',$aminoacids);
				$perc_ali=(scalar@length_ali/$length_seq)*100;
				if (($perc_ali>79 && $perc_ali<121 && $length_seq<401) || ($length_seq>400 && $length_seq<901 && $perc_ali>86 && $perc_ali<114) || ($length_seq>900 && $perc_ali>92 && $perc_ali<108)){
					push(@fichero,">".$prot_name."\t(aligned)\t",scalar@length_ali,"\n");
					push(@fichero,$aminoacids."\n");
					$k=$z-1;
					$quantity++;
					last;
				}
				elsif (defined $score && (($perc>79 && $perc<121 && $length_seq<401) || ($length_seq>400 && $length_seq<901 && $perc>86 && $perc<114) || ($length_seq>900 && $perc>92 && $perc<108))){
					push(@proteinas,$prot_name);
					$quantity++;
					$k=$z;
					last;
				}
				else {last;}
			}
			$z++;
			if ($listado[$z]=~/^>.+/){$k=$z-1;last;}
		}
	}
}

# con el listado obtenido, busco en el database y me quedo con un fichero que reune todas las secuencias para enviar al programa de alineamientos múltiples
# no cargo el database en un array sino lo leo linea por linea y me quedo sólo con las secuencias q me interesan
$w=0;
$aminoacids=undef;
open(DATA,$db_file);

if (scalar@proteinas!=0){
	$control='f';
	while ($line=<DATA>){
		if ($w==scalar@proteinas && $control eq 'f'){
			last;
		}
		if ($line=~/^>.+/){
			$control='f';
		}
		if ($control eq 't'){
			$aminoacids=$aminoacids.$line;
		}	
		if ($line=~/^>(.+)/){
			@temporary=undef;
			@temporary=split(/ /,$1);
			$reference=$temporary[0];	
			foreach $k (@proteinas){
				if ($reference eq $k && $w>0){
					@temporary=undef;
					@temporary=split('',$aminoacids);
					$aminoacids=undef;
					$length=undef;
					foreach $i (@temporary){
						if ($i=~/[A-Z]/i){$aminoacids=$aminoacids.$i;$length++}
					}
					chomp($fichero[$#fichero]);
					$fichero[$#fichero]=$fichero[$#fichero]."\t$length\n";
					push(@fichero,$aminoacids."\n");
					$aminoacids=undef;				
					$w++;
					push(@fichero,$line);
					$control='t';
					last;
				}
				elsif ($reference eq $k){
					$w++;
					push(@fichero,$line);
					$control='t';
					last;
				}
			}
		}
	}			
}
close(DATA);


@temporary="";
@temporary=split('',$aminoacids);
$aminoacids=undef;
$length=undef;
foreach $i (@temporary){
	if ($i=~/[A-Z]/i){$aminoacids=$aminoacids.$i;$length++;}
}
chomp($fichero[$#fichero]);
$fichero[$#fichero]=$fichero[$#fichero]."\t$length\n";
push(@fichero,$aminoacids."\n");	



open(RESU, ">$output");
foreach $i (@fichero) {
	print RESU $i;
}
close(RESU);


