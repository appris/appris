#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $input_file 	= "";
my $output_file = "";
my $num_cpus 	= 10;
my $tmp_dir 	= "/tmp/matador3D2";
my $reset 		= 0;

# Needed binaries
my $hmmbuild_bin = "/home/jrodriguezr/apps/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild";
my $hmmpress_bin = "/home/jrodriguezr/apps/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmpress";
my $reformat_bin = "/home/jrodriguezr/apps/hhsuite-2.0.15-linux-x86_64/lib/hh/scripts/reformat.pl";


GetOptions (
	"i=s" 	=> \$input_file,
	"t=s" 	=> \$tmp_dir,
	"o=s" 	=> \$output_file,
	"n=i"   => \$num_cpus,
	"r=i"   => \$reset,
	);


check_input_file($input_file);

# Start from 0
if($reset){
	`rm -R "$tmp_dir/a3m" "$tmp_dir/a3m_clean" "$tmp_dir/a2m" "$tmp_dir/pdb_seqs" "$tmp_dir/hmms" $output_file*`;
}


# Untar and unzip input file:
my $a3m_dir = "$tmp_dir/a3m";
`mkdir -p $a3m_dir`;
system("tar -zxvf $input_file -C $a3m_dir");

# Clean a3m alis
my $a3m_clean_dir = "$tmp_dir/a3m_clean";
`mkdir -p $a3m_clean_dir`;
system("tools/clean_a3m.bash $a3m_dir $a3m_clean_dir $tmp_dir $num_cpus");

# From a3m to a2m removing gaps that are not in the target
my $a2m_dir = "$tmp_dir/a2m";
`mkdir -p $a2m_dir`;
system("tools/reformat.bash $a3m_clean_dir $a2m_dir $tmp_dir $num_cpus $reformat_bin");

# Extract the PDB sequence, i.e., extract the first two lines from the a2m alignment (PDB header + PDB sequence).
# This is done to generate a profile just based in each PDB sequence, to avoid that obvious
# hits are not taken into account due uninformative regions in profiles
my $pdb_seqs_dir = "$tmp_dir/pdb_seqs";
`mkdir -p $pdb_seqs_dir`;
my @files = `ls $a2m_dir`;
chomp @files;
foreach my $file (@files){
	my $input_file 	= "$a2m_dir/$file";
	my $root 		= substr($file, 0, -4);
	my $output_file = "$pdb_seqs_dir/$root\_pdb_seq.a2m";
	if(!-e $output_file){
		# Extract the first two lines (PDB header + PDB sequence), change the protein id in the header to template_pdb_seq to distinguish between pdb sequence only and the alignment with other sequences
		system("tools/extract_seq.bash $input_file $output_file")
	}
}




# Build HMM for each aligment
my $hmm_dir = "$tmp_dir/hmms";
`mkdir -p $hmm_dir`;

@files = `ls $a2m_dir`;
chomp @files;
foreach my $file (@files){
	if(!-e "$hmm_dir/$file"){ # Check profile does not exist
		system("$hmmbuild_bin --amino --cpu $num_cpus $hmm_dir/$file $a2m_dir/$file");
	}
}

# Add profile with only pdb sequence
@files = `ls $pdb_seqs_dir`;
chomp @files;
foreach my $file (@files){
	my $input_file 	= "$pdb_seqs_dir/$file";
	my $root 		= substr($file, 0, -4);
	my $output_file = "$hmm_dir/$file";
	if(!-e $output_file){ # Check profile does not exist
		system("$hmmbuild_bin --amino --cpu $num_cpus $output_file $input_file");
	}
}


# Create output directory if it does not exist
my $output_dir = dirname($output_file);
`mkdir -p $output_dir`;

# Concatenate alignments and compress database
`cat $hmm_dir/* > $output_file`;
system("$hmmpress_bin -f $output_file");




exit;



sub usage{
	
	print "\nUsage: generate_db_matador3D2.pl -i input_file -o output_file [-t tmp_dir] [-n num_cpus] [-r 0]\n";
	print "Description:\n";
	print "Given a tar.gz with the a3m HHPred alignments generate a profiles database.\n\n";
	print "-i input_file  -> tar.gz input file containing a3m alignments from HHPred.\n";
	print "-o output_file -> Output profile database that can be searched with HMMER tools.\n";
	print "-t tmp_dir     -> Optional, default /tmp. Temporal directorty where intermediate file and directories will be stored.\n";
	print "-n num_cpus    -> Optional, default 10. Number of processes that will be run.\n";
	print "-r reset       -> Optional, default 0. If different from 0, it will remove any intermediate file that has been previously generated. Repeat all the process from the scracth.\n";
	exit;
}



sub check_input_file{
	
	my $input_file = shift;
	
	if(-e $input_file){
		if(-z $input_file){
			print STDERR "Error: The input file $input_file is empty\n";
			usage();
		}
	}
	else{
		print STDERR "Error: The input file $input_file does not exist\n";
		usage();
	}
}
