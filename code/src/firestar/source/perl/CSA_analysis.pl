#!/usr/bin/perl

use strict;

use FindBin;
my $cwd;
BEGIN{
	$cwd=$FindBin::Bin;
}
use Config::IniFiles;

use lib "$cwd/lib";
use fire_libraries;
use fire_CSA_study;

#~~~~ Here we create the main object for the firestar analysis	~~~~~~~~

my $firestar=fire_CSA_study->new();

#~~~~~~~~~~~~~	 checking the input parameters	 ~~~~~~~~~~~~~~~~

$firestar->{cwd}=$cwd;
my $check=$firestar->input_parameters_eva(\@ARGV);

unless ($check eq "OK"){
	<STDIN>;
	system("perldoc $cwd/firestar.pl");
	exit;
}

#~~~~~~~~~~	LIBRARIES	~~~~~~~~~~~~~

my $fire_libs=fire_libraries->new(-file=>$firestar->{config_file});

$fire_libs->ec2go();
$fire_libs->tag(-type=>'COGNATE');			# aqui se guarda la informacion acerca de los cognate ligands
$fire_libs->tag(-type=>'POSSIBLE');		# aqui se guarda la informacion acerca de los possibles cognate ligands
$fire_libs->metal(-type=>'MET');			# aqui se guarda la informacion acerca de los metales
$fire_libs->cif2go();
$fire_libs->csa2ec();

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		BEGINNING of the ANALYSIS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~	HHSEARCH and PSIBLAST	    ~~~~~~~~~~~~~

$firestar->homologs_searching();

#~~~~~~~~~~	Extraction of the alignments from the output files	~~~~~~~~

$firestar->parse_psiblast();

if (length($firestar->{sequence})<2001){		$firestar->parse_hhsearch();}


my @output=$firestar->site_info_extractor($fire_libs);
my $OUTPUT;
$OUTPUT = \*STDOUT  unless(open($OUTPUT,'>',$firestar->{outfile}));

foreach my $i(@output){
	print $OUTPUT "$i\n";
}
close($OUTPUT);



$fire_libs->DESTROY();
$firestar->DESTROY();

END;

=head1 NAME

firestar.pl

=head1 DESCRIPTION
	
=head2 Required arguments:


	you can choose to run the script using an input file or introducing the information separately or introducing only the PDB code of your protein of interest


	-f= <Input file> input file

	-type= <str|fas> define the type of input file. Fasta and PDB format are supported


	or

	
	-q= <query name> the name or ID for your protein query

	-s= <sequence> the aminoacidic sequence of you protein
	

	or
	
	
	-pdb= <PDB code> the PDB accession code for your protein query
	

	or 


	-uni= <Uniprot accession code> the Uniprot accession for your protein of interest



=head2 Optional arguments :


	-e=<e-value>			you can choose the e-value used for the PSIBLAST results'extraction fltering. 			default = 10

	-o=<output_file>		the output file for the firestar results.							default = ./protein_name.res

	-cut=<reliability cut-off>	the minimum reliability score needed for a site to be showed.					deafult = 0

	-csa=<YES|NO|ONLY>		with this option you can control the output of the predictions based on catalytic site		default = YES
					atlas (CSA,<http://www.ebi.ac.uk/thornton-srv/databases/CSA/>).

	-cog=<YES|NO>			with this option set as "NO" you can exclude from the results printed out the predicted		default = YES
					sites tagged as non biological (NON_COGNATE).
	
	-opt=<output format>		with this option you can choose between different output formats. 'appris','text','SIAM'	default = text
					and 'web' are the formats available.

	-chain=<chain letter>		if you are running firestar using a PDB accession code or a PDB formatted file, you can		default = A
					specify the chain you're interested in.

	-conf=<path_to_Config file>	firestar gets some parameters (e.g paths, MySQL connection info, etc.) from a configuration	default = CONFIG_fire_var.ini
					file called CONFIG_fire_var.ini located in the main directory. if you want to use
					a different one, you can use this option


=head1 EXAMPLES

	perl firestar.pl -q=EXAMPLE -s=MKMASTRCKLARYLEDLEDVDLKKFKMHLEDYPPQKGCIPLPRGQTEKADHVDLATLMIDFNGEEKAWAMAVWIFAAINRRDLYEKAKRDE -csa=ONLY -e=0.01 -cut=60

	perl firestar.pl -f=../../firestar_analysis/example.fasta -type=fas -csa=NO -cog=NO
	
	perl firestar.pl -pdb=3rsv -chain=A -e=1O -opt=SIAM

=head1 AUTHOR
	
	
=head2 Paolo Maietta -paolo.maietta@gmail.com- (CNIO)


=cut
