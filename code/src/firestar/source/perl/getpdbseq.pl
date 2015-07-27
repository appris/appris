#!/usr/bin/perl

use strict;
use FindBin;
my $cwd=$FindBin::Bin;
use Config::IniFiles;
my $variables=Config::IniFiles->new(-file => "$cwd/../CONFIG_fire_var.ini");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $pdb_path=$variables->val('PATHS','PDB');$pdb_path.="/wwPDB";

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my @pname=split(/\//,$0);
my $help="#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Extraxt from pdb file the sequence in one letter aminoacid code
and print them in fasta format.
>
$pname[$#pname] -z 1tco 	-->	To get file from a pdb gziped DB
$pname[$#pname] -p /path/to/pdb -->	Path were pdb files are stored
$pname[$#pname] -f 1tco.pdb	-->	To open a text file
$pname[$#pname] -c		-->	You may want to select one chain or not
$pname[$#pname] -n		-->	Print line with PDB numbers
$pname[$#pname] -s		-->	Print line with PDB sequence
$pname[$#pname] -l length	-->	Min length of sequence
$pname[$#pname] -x		-->	Do not parse for multiple models (faster)
$pname[$#pname] -t 		-->	Print translation hash (3 letter to one)
$pname[$#pname] -/-h/-help	-->	This help
>
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";

if($ARGV[0]=~/-h/ or $ARGV[0] eq "-"){die "$help";}

#--------> global vars and lib arrays
my @PDB;
my $nums=0;
my $seqline=0;
my %allAA;
my $nmrParser=0;
my $pdb;
my %ARGS;
my $seq_length=2;
#--------<

#--------> pipeline
read_libs();
arguments();
arg_hash();
get_seq_nums();
#--------<

###-----------------------------------------------------------------------
sub arguments{

my $args=join("",@ARGV);
chomp $args;
my @args=split(/-/,$args);
for(@args){
	my $off=(length$_)-1;
	my $key=substr($_,0,1);
	my $val="";
	if($off > 0){
		$val=substr($_,1,$off);
		}
	$ARGS{$key}=$val;
	}
}


###----------------------------------------------------------------------
sub arg_hash{

if(exists$ARGS{t}){
	foreach my $key(keys%allAA){
		print "'$key'=>'$allAA{$key}'\n";
		}
	die "uncomplete?\n";
	}	
if(exists$ARGS{z}){
	if(exists$ARGS{p}){
		$pdb_path=$ARGS{p};
		}
	@PDB=`gunzip -c $pdb_path/pdb$ARGS{z}.ent.gz`;
	if($PDB[0]!~/\w/){print STDERR "no such pdb in /drives/databases/pdb/ !!!\n";}
	$pdb=$ARGS{z};
	}
if(exists$ARGS{f}){
	open(FILE,"$ARGS{f}") || die "$ARGS{f} no such file !!!\n";
	@PDB=<FILE>;	close FILE;
	my @pdb=split(/\//,$ARGS{f});
	$pdb=$pdb[$#pdb];$pdb=~s/\.pdb//;
	}
if(exists$ARGS{n}){$nums=1;}
if(exists$ARGS{s}){$seqline=1;}
if(exists$ARGS{l}){$seq_length=$ARGS{l};}
if(exists$ARGS{x}){$nmrParser=1;}
else{(@PDB);}
}

###-----------------------------------------------------------------------

sub get_seq_nums{
my @PDBCA=grep{$_=~/ATOM\s*\d+\s+CA\s*[A-Z]?[A-Z][A-Z][A-Z] / or $_=~/HETATM\s*\d+\s*CA\s+[A-Z]?[A-Z][A-Z][A-Z] /}@PDB;
my %CHAINS;
my $resnum;
my $resnumold;
my %NUMS;

for(@PDBCA){
	my $aa3=substr($_,17,3);
	my $chain=substr($_,21,1);
	$resnum=substr($_,22,5);
	$resnum=~s/ //g;
	my $aa1="";
	if(exists$allAA{$aa3}){
		$aa1=$allAA{$aa3};
		if($nums==0 and $resnum ne $resnumold){	
			$CHAINS{$chain}="$CHAINS{$chain}$aa1";
			}
		if($nums==1 and $resnum ne $resnumold){	
			$CHAINS{$chain}="$CHAINS{$chain}$aa1";
			$NUMS{$chain}="$NUMS{$chain} $resnum";
			}
		}
	$resnumold=$resnum;
	}
	
foreach my$key(keys%CHAINS){
	if(length$CHAINS{$key} > $seq_length){
		my @splNum=();
		if($nums==1){
			$NUMS{$key}=~s/^ +//;
			@splNum=split(/ /,$NUMS{$key});
			}
		my @splSeq=split(//,$CHAINS{$key});
		while($splSeq[$#splSeq] eq "X"){
			pop @splSeq;
			if($nums==1){pop @splNum;}
			}
		if($nums==1){
		if(scalar@splSeq != scalar@splNum){print STDERR "Warning!! $pdb: sequnce elements = $#splSeq; numeric elements = $#splNum\n";}
			}
		my $sequence=join("",@splSeq);
		my $numbers=join(" ",@splNum);
		$key=~s/\s+//;
		if(exists$ARGS{c}){
			if($key eq $ARGS{c}){
				if($nums==0 and $seqline==1 and $sequence=~/\w/){
					print ">$pdb\_$key\n$sequence\n";
					}
				elsif($nums==1 and $seqline==1 and $sequence=~/\w/){
					print ">$pdb\_$key\n$sequence\n$numbers\n";
					}
				elsif($nums==1 and $seqline==0 and $sequence=~/\w/){
					print ">$pdb\_$key\n$numbers\n";
					}
				}
			}
		else{	if($nums==0 and $seqline==1 and $sequence=~/\w/){
				print ">$pdb\_$key\n$sequence\n";
				}
			elsif($nums==1 and $seqline==1 and $sequence=~/\w/){
				print ">$pdb\_$key\n$sequence\n$numbers\n";
				}
			elsif($nums==1 and $seqline==0 and $sequence=~/\w/){
				print ">$pdb\_$key\n$numbers\n";
				}
			}
		}
	}
}

###-----------------------------------------------------------------------

sub nmr{
my $label=0;
my @PDBTMP=();	my @EXPDTA;
@EXPDTA=grep{$_=~/^MODEL/}@PDB;
if($EXPDTA[0]=~/^MODEL/){
	for(@PDB){
		if($_=~/^ENDMDL/ and $label==1){
			$label=2;
			last;
			}
		if($_=~/^MODEL/ and $label==0){
			$label=1;		
			}		
		if($label==1){	
			push(@PDBTMP,$_);
			}
		}
	@PDB=@PDBTMP;
	}
}

###------------------------------------------------###
# Translator hash from 3 to 1 aminoacid letter codes #
#  Source: the components.cif dictionary             #
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

sub read_libs{
%allAA=(
'ALA'=>'A','VAL'=>'V','LEU'=>'L','ILE'=>'I','GLY'=>'G',
'GLU'=>'E','GLN'=>'Q','ASP'=>'D','ASN'=>'N','ARG'=>'R',
'PRO'=>'P','TYR'=>'Y','LYS'=>'K','TRP'=>'W','SER'=>'S',
'THR'=>'T','PHE'=>'F','MET'=>'M','CYS'=>'C','HIS'=>'H',
'0A0'=>'D','0A1'=>'Y','0A2'=>'K','0A5'=>'N','0A8'=>'C','0A9'=>'F','0AA'=>'V','0AB'=>'V',
'0AC'=>'G','0AF'=>'W','0AG'=>'L','0AH'=>'S','0AK'=>'D','0AY'=>'K','0AZ'=>'P','0CS'=>'A',
'1LU'=>'L','1PA'=>'F','1TQ'=>'W','1TY'=>'Y','200'=>'F','23F'=>'F','2AG'=>'G','2AS'=>'D',
'2FM'=>'M','2LU'=>'L','2ML'=>'L','2MR'=>'R','2MT'=>'P','2TY'=>'Y','2VA'=>'V','3AH'=>'H',
'3MD'=>'D','4DP'=>'W','4FB'=>'P','4FW'=>'W','4HT'=>'W','4PH'=>'F','5CS'=>'C','5HP'=>'E',
'6CL'=>'K','6CW'=>'W','AA3'=>'A','AA4'=>'A','AAR'=>'R','ABA'=>'A','ACB'=>'D','ACL'=>'R',
'AEI'=>'D','AFA'=>'N','AGM'=>'R','AHB'=>'N','AHO'=>'A','AHP'=>'A','AIB'=>'A','AKL'=>'D',
'ALC'=>'A','ALG'=>'R','ALM'=>'A','ALN'=>'A','ALO'=>'T','ALS'=>'A','ALY'=>'K','APH'=>'A',
'API'=>'K','APK'=>'K','AR4'=>'E','ARM'=>'R','AS2'=>'D','ASA'=>'D','ASB'=>'D','ASI'=>'D',
'ASK'=>'D','ASL'=>'D','ASQ'=>'D','AYA'=>'A','AZK'=>'K','AZS'=>'S','AZY'=>'Y','B1F'=>'F',
'B2A'=>'A','B2F'=>'F','B2I'=>'I','B2V'=>'V','B3A'=>'A','B3D'=>'D','B3E'=>'E','B3K'=>'K',
'B3S'=>'S','B3X'=>'N','B3Y'=>'Y','BAL'=>'A','BBC'=>'C','BCS'=>'C','BCX'=>'C','BFD'=>'D',
'BHD'=>'D','BIF'=>'F','BLE'=>'L','BLY'=>'K','BMT'=>'T','BNN'=>'A','BOR'=>'R','BPE'=>'C',
'BSE'=>'S','BTA'=>'L','BTC'=>'C','BUC'=>'C','BUG'=>'L','C1X'=>'K','C3Y'=>'C','C5C'=>'C',
'C6C'=>'C','CAB'=>'A','CAF'=>'C','CAS'=>'C','CCL'=>'K','CCS'=>'C','CEA'=>'C','CGA'=>'E',
'CGU'=>'E','CHG'=>'A','CHP'=>'G','CIR'=>'R','CLB'=>'A','CLD'=>'A','CLE'=>'L','CLG'=>'K',
'CLH'=>'K','CME'=>'C','CMH'=>'C','CML'=>'C','CMT'=>'C','CR5'=>'G','CRU'=>'E','CS1'=>'C',
'CS3'=>'C','CS4'=>'C','CSA'=>'C','CSB'=>'C','CSD'=>'C','CSE'=>'C','CSI'=>'G','CSO'=>'C',
'CSP'=>'C','CSR'=>'C','CSS'=>'C','CSU'=>'C','CSW'=>'C','CSX'=>'C','CSZ'=>'C','CTH'=>'T',
'CWR'=>'S','CXM'=>'M','CY0'=>'C','CY1'=>'C','CY3'=>'C','CY4'=>'C','CYA'=>'C','CYD'=>'C',
'CYF'=>'C','CYG'=>'C','CYM'=>'C','CYQ'=>'C','CYR'=>'C','CZ2'=>'C','CZZ'=>'C','DAB'=>'A',
'DAH'=>'F','DBS'=>'S','DBU'=>'A','DBY'=>'Y','DBZ'=>'A','DDE'=>'H','DHA'=>'A','DHN'=>'V',
'DIR'=>'R','DLS'=>'K','DM0'=>'K','DMH'=>'N','DMK'=>'D','DNL'=>'K','DNP'=>'A','DNS'=>'K',
'DOH'=>'D','DON'=>'L','DPL'=>'P','DPN'=>'F','DPP'=>'A','DPQ'=>'Y','EFC'=>'C','FCL'=>'F',
'FGL'=>'G','FGP'=>'S','FLA'=>'A','FLE'=>'L','FME'=>'M','FOE'=>'C','FPA'=>'F','FT6'=>'W',
'FTR'=>'W','FTY'=>'Y','GAU'=>'E','GGL'=>'E','GHG'=>'Q','GHP'=>'G','GLH'=>'Q','GLQ'=>'E',
'GLZ'=>'G','GMA'=>'E','GPL'=>'K','GSC'=>'G','GSU'=>'E','GT9'=>'C','H5M'=>'P','HAC'=>'A',
'HAR'=>'R','HBN'=>'H','HIA'=>'H','HIC'=>'H','HIP'=>'H','HIQ'=>'H','HLU'=>'L','HMF'=>'A',
'HMR'=>'R','HPC'=>'F','HPE'=>'F','HPQ'=>'F','HRG'=>'R','HRP'=>'W','HSE'=>'S','HSL'=>'S',
'HSO'=>'H','HTI'=>'C','HTR'=>'W','HV5'=>'A','HY3'=>'P','HYP'=>'P','I58'=>'K','IAM'=>'A',
'IAS'=>'D','IGL'=>'G','IIL'=>'I','ILG'=>'E','ILX'=>'I','IML'=>'I','IT1'=>'K','IYR'=>'Y',
'K1R'=>'C','KCX'=>'K','KGC'=>'K','KOR'=>'M','KST'=>'K','KYN'=>'A','LA2'=>'K','LAL'=>'A',
'LCK'=>'K','LCX'=>'K','LDH'=>'K','LED'=>'L','LEF'=>'L','LLP'=>'K','LLY'=>'K','LME'=>'E',
'LPD'=>'P','LPG'=>'G','LPS'=>'S','LTR'=>'W','LVG'=>'G','LYM'=>'K','LYN'=>'K','LYR'=>'K',
'LYX'=>'K','LYZ'=>'K','M0H'=>'C','M3L'=>'K','MAA'=>'A','MAI'=>'R','MBQ'=>'Y','MC1'=>'S',
'MCL'=>'K','MCS'=>'C','MEA'=>'F','MEG'=>'E','MEN'=>'N','MEQ'=>'Q','MEU'=>'G','MGG'=>'R',
'MGN'=>'Q','MGY'=>'G','MHL'=>'L','MHO'=>'M','MHS'=>'H','MIS'=>'S','MLE'=>'L','MLL'=>'L',
'MLY'=>'K','MLZ'=>'K','MME'=>'M','MNL'=>'L','MNV'=>'V','MPQ'=>'G','MSE'=>'M','MSL'=>'M',
'MSO'=>'M','MVA'=>'V','N10'=>'S','N7P'=>'P','NAL'=>'A','NAM'=>'A','NBQ'=>'Y','NC1'=>'S',
'NCB'=>'A','NEM'=>'H','NEP'=>'H','NFA'=>'F','NIY'=>'Y','NLE'=>'L','NLN'=>'L','NLO'=>'L',
'NLP'=>'L','NLQ'=>'Q','NMC'=>'G','NMM'=>'R','NNH'=>'R','NPH'=>'C','NTY'=>'Y','NVA'=>'V',
'NZH'=>'H','OAS'=>'S','OCS'=>'C','OCY'=>'C','OHI'=>'H','OHS'=>'D','OMT'=>'M','OPR'=>'R',
'ORN'=>'A','ORQ'=>'R','OSE'=>'S','OTY'=>'Y','OXX'=>'D','P1L'=>'C','P2Y'=>'P','PAQ'=>'Y',
'PAS'=>'D','PAT'=>'W','PAU'=>'A','PBB'=>'C','PBF'=>'F','PCA'=>'E','PCC'=>'P','PCS'=>'F',
'PEC'=>'C','PF5'=>'F','PFF'=>'F','PG1'=>'S','PGY'=>'G','PHA'=>'F','PHD'=>'D','PHI'=>'F',
'PHL'=>'F','PHM'=>'F','PLE'=>'L','PM3'=>'F','POM'=>'P','PPH'=>'L','PPN'=>'F','PR3'=>'C',
'PRR'=>'A','PRS'=>'P','PSA'=>'F','PSH'=>'H','PTH'=>'Y','PTM'=>'Y','PTR'=>'Y','PVH'=>'H',
'PYA'=>'A','PYX'=>'C','R1A'=>'C','R1B'=>'C','R1F'=>'C','R7A'=>'C','RCY'=>'C','SAC'=>'S',
'SAH'=>'C','SAR'=>'G','SBD'=>'S','SBL'=>'S','SCH'=>'C','SCS'=>'C','SCY'=>'C','SDP'=>'S',
'SEB'=>'S','SEC'=>'A','SEG'=>'A','SEL'=>'S','SEP'=>'S','SET'=>'S','SHC'=>'C','SHP'=>'G',
'SHR'=>'K','SIB'=>'C','SLZ'=>'K','SMC'=>'C','SME'=>'M','SMF'=>'F','SNC'=>'C','SOC'=>'C',
'SOY'=>'S','STY'=>'Y','SVA'=>'S','T11'=>'F','TAV'=>'D','TBG'=>'G','TBM'=>'T','TFQ'=>'F',
'THC'=>'T','TIH'=>'A','TMB'=>'T','TMD'=>'T','TNB'=>'C','TNR'=>'S','TOX'=>'W','TPL'=>'W',
'TPO'=>'T','TPQ'=>'Y','TQQ'=>'W','TRG'=>'K','TRN'=>'W','TRO'=>'W','TRQ'=>'W','TRW'=>'W',
'TRX'=>'W','TTQ'=>'W','TTS'=>'Y','TYB'=>'Y','TYI'=>'Y','TYN'=>'Y','TYO'=>'Y','TYQ'=>'Y',
'TYS'=>'Y','TYT'=>'Y','TYY'=>'Y','UMA'=>'A','VAD'=>'V','VAF'=>'V','XX1'=>'K','YCM'=>'C',
);


}
