#! /usr/bin/perl
# This is a simple script which will carry out all of the basic steps
# required to make a MEMSAT V3 prediction. Note that it assumes that the
# following programs are in the appropriate directories:

###################
# Programs needed #
###################

# blastpgp - PSIBLAST executable (from NCBI toolkit)
# makemat - IMPALA utility (from NCBI toolkit)
# globmem - MEMSAT V3 program
# mem_pred - MEMSAT V3 program
# nnsat - MEMSAT V3 program

use strict;
use FindBin;
use Getopt::Long;

my ($db_file) = undef;
my ($name) = undef;
my ($input) = undef;
my ($out_memsat) = undef;
my ($out_chk) = undef;
my ($out_blast) = undef;
my ($tmp_dir) = undef;
my ($cache_dir) = undef;

&GetOptions(
	'db=s'			=> \$db_file,
	'name=s'		=> \$name,	
	'input=s'		=> \$input,
	'out-memsat=s'	=> \$out_memsat,
	'out-chk=s'		=> \$out_chk,
	'out-blast=s'	=> \$out_blast,
	'tmp-dir=s'		=> \$tmp_dir,
	'cache-dir=s'	=> \$cache_dir,
);

###################
# Internal values #
###################
my ($pwd) = "$FindBin::Bin";

# Where the PSIPRED V2 programs have been installed
my ($execdir) = $pwd."/bin";
my ($weights_dir) = $pwd."/data";
my ($glob_weight_file) = $weights_dir."/glob_weights.dat";
my ($weight_file) = $weights_dir."/weights.dat";

my ($pn_file) = $tmp_dir."/$name.pn";
my ($sn_file) = $tmp_dir."/$name.sn";
my ($mtx_file) = $cache_dir."/seq.mtx";
my ($globmem_file) = $tmp_dir."/$name.globmem";
my ($nn_file) = $tmp_dir."/$name.nn";
my ($memsat_file) = $tmp_dir."/$name.memsat";


unless (
	(-e $out_blast and (-s $out_blast > 0)) and
	(-e $out_chk and (-s $out_chk > 0))
) {
	eval {
		print "Running PSI-BLAST with sequence $input ...\n";
		print "blastpgp -b 50 -v 50 -j 2 -h 1e-4 -e 1e-3 -a 4 -d $db_file -i $input -C $out_chk -o $out_blast\n";
		system ("blastpgp -b 50 -v 50 -j 2 -h 1e-4 -e 1e-3 -a 4 -d $db_file -i $input -C $out_chk -o $out_blast");		
	};
	exit 1 if($@);
}

print "Predicting transmembrane topology... \n\n";

unless (
	(-e $pn_file and (-s $pn_file > 0)) and
	(-e $sn_file and (-s $sn_file > 0))
) {
	eval {
		open (FILE1,">$pn_file");
		print FILE1 "$out_chk\n";
		close (FILE1);	
	};
	exit 1 if($@);
	
	eval {
		open (FILE2,">$sn_file");
		print FILE2 "$input\n";
		close (FILE2);
	};
	exit 1 if($@);
}

unless ( -e $mtx_file and (-s $mtx_file > 0) ) {
	eval {
		print "cd $tmp_dir && makemat -P $name\n";
		system("cd $tmp_dir && makemat -P $name");
	};
	exit 1 if($@);
}

unless ( -e $globmem_file and (-s $globmem_file > 0) ) {
	eval {
		print "$execdir/globmem $glob_weight_file $mtx_file > $globmem_file\n";
		system("$execdir/globmem $glob_weight_file $mtx_file > $globmem_file");
	};
	exit 1 if($@);
}

unless ( -e $nn_file and (-s $nn_file > 0) ) {
	eval {
		print "$execdir/mem_pred $weight_file $mtx_file > $nn_file\n";
		system("$execdir/mem_pred $weight_file $mtx_file > $nn_file");
	};
	exit 1 if($@);
}

unless ( -e $out_memsat and (-s $out_memsat > 0) ) {
	eval {
		print "$execdir/nnsat $nn_file > $out_memsat\n";
		system("$execdir/nnsat $nn_file > $out_memsat");
	};
	exit 1 if($@);
}
