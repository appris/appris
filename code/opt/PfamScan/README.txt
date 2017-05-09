This readme should help you get started with "pfam_scan.pl", which is for use
with the HMMER3 version of HMMER.

--------------------------------------------------------------------------------
- Setting up -------------------------------------------------------------------
--------------------------------------------------------------------------------


Unpack the script
=================

  shell% tar zxvf PfamScan.tar.gz
  ...
  shell% cd PfamScan


Install HMMER3
==============

1) Get the current tarball of the latest HMMER3 release from the HMMER site at
Janelia farm. For example:

  shell% wget ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-3.1b1.tar.gz

2) Uncompress the tarball

  shell% tar zxvf hmmer-3.1b1.tar.gz

3) Compile the HMMER3 source code (optional - skip to the section about adding
HMMER3 binaries to your path). You can also download precompiled binaries, but
it is recommmended that you compile your source code with your system compiler.
HMMER usually compiles in three easy steps.  For detailed instructions, see the
HMMER3 INSTALL guide.

i) Change the the HMMER source directory:

  shell% cd hmmer-3.1b1

ii) Run "configure". For best performance, and depending on the architecture of
the machine that will be running HMMER, it's recommended that you use the Intel
compiler, "icc" to build the HMMER binaries. You can also use "gcc" but note
that older versions (<3.3) will not correctly compile the code.

You can tell "configure" which compiler to use by adding the "CC" flag:

  shell% ./configure CC=icc LDFLAGS="-static" --prefix=/path/to/install/hmmer3

iii) Run make:

  shell% make

iv) Check that everything is compiled and working, by running:

  shell% make check

v) Install HMMER:

  shell% make install


Adding HMMER3 binaries to your path
===================================

You need to make sure that the HMMER binaries are found on your shell's
executable search path.

If you are using bash:

  bash% export PATH=/path/to/install/hmmer3:$PATH

or, if you are using csh (or tcsh):

  csh% setenv PATH /path/to/install/hmmer3:$PATH

Note: if you are using the pre-compiled binaries, just point to those, e.g.
/path/to/uncompressedTar/hmmer-3.0b3/binaries. 


Non-standard Perl dependencies
==============================

The PfamScan.pm module depends on several modules that don't come as part of a
standard Perl distribution, notably the Moose framework. Everything that you
need can be installed using the "cpan" tool. You'll need to make sure that
you've already configured your "cpan" environment, and then:

  shell% cpan Moose

Moose itself has quite a few dependencies, so don't worry if it looks like
you're installing half of CPAN !

PfamScan.pm also requires bioperl. We have currently only tested against
bioperl 1.4, and we believe it works with 1.6. You can install bioperl via
CPAN, or you can download bioperl 1.4 from here:

  http://bioperl.org/DIST/bioperl-1.4.tar.gz


Adding Pfam Modules to your PERL5LIB
====================================

If you are using bash:

  bash% export PERL5LIB=/path/to/pfam_scanDir:$PERL5LIB

or for C-shells:

  csh% setenv PERL5LIB /path/to/pfam_scanDir:$PERL5LIB


Change the path of Perl in pfam_scan.pl
=======================================

Open PfamScan/pfam_scan.pl in a text editor, and change the first line
of the code to point to your Perl.

You should be good to go now!


--------------------------------------------------------------------------------
- Running searches using "pfam_scan.pl" ----------------------------------------
--------------------------------------------------------------------------------


Download Pfam data files
========================

You will need to download the following files from the Pfam ftp site 
(ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/):

Pfam-A.hmm
Pfam-A.hmm.dat  
Pfam-B.hmm      
Pfam-B.hmm.dat  
active_site.dat

You will need to generate binary files for Pfam-A.hmm and Pfam-B.hmm by running the following commands:
 
hmmpress Pfam-A.hmm
hmmpress Pfam-B.hmm


Using pfam_scan.pl
==================

"pfam_scan.pl" is a program that searches a FASTA file against a library of
Pfam HMMs.


REQUIREMENTS
============

To recap, "pfam_scan.pl" requires:

  - Several modules written by Pfam, all of which are available as part of this
    tarball
  - Bioperl and HMMER3 installed (hmmscan and hmmalign should be in your path)
  - Pfam-A.hmm (and binaries): a data file that contains the Pfam-A library of HMMs
  - Pfam-B.hmm (and binaries): a data file that contains the Pfam-B library of HMMs
  - Pfam-A.hmm.dat: a data file that contains information about each Pfam-A
    family
  - Pfam-B.hmm.dat: a data file that contains information about each Pfam-B
    family
  - active_site.dat: a data file needed for the -as option, which contains active
    site information about each family
  - a FASTA-format file containing your query sequence(s)


Usage
=====

pfam_scan.pl -fasta <fasta_file> -dir <directory location of Pfam files>

Additonal options:

  -h              : show this help
  -o <file>       : output file, otherwise send to STDOUT
  -clan_overlap   : show overlapping hits within clan member families (applies to Pfam-A families only)
  -align          : show the HMM-sequence alignment for each match
  -e_seq <n>      : specify hmmscan evalue sequence cutoff for Pfam-A searches (default Pfam defined)
  -e_dom <n>      : specify hmmscan evalue domain cutoff for Pfam-A searches (default Pfam defined)
  -b_seq <n>      : specify hmmscan bit score sequence cutoff for Pfam-A searches (default Pfam defined)
  -b_dom <n>      : specify hmmscan bit score domain cutoff for Pfam-A searches (default Pfam defined)
  -pfamB          : search against Pfam-B HMMs (uses E-value sequence and domain cutoff 0.001),
                    in addition to searching Pfam-A HMMs
  -only_pfamB     : search against Pfam-B HMMs only (uses E-value sequence and domain cutoff 0.001)
  -as             : predict active site residues for Pfam-A matches
  -json [pretty]  : write results in JSON format. If the optional value "pretty" is given,
                    the JSON output will be formatted using the "pretty" option in the JSON
                    module


For more help, check the perldoc:

  shell% perldoc pfam_scan.pl


Output format
=============

The output should be familiar to anyone who's used the old "pfam_scan.pl"
script. Each line contains the following information:

<seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value><significance> <clan> <predicted_active_site_residues>


Example output (with -as option):

  Q5NEL3.1      2    224      2    227 PB01348.1    Pfam-B_13481      Pfam-B     1   184   226    358.5   1.4e-107  NA NA
  O65039.1     38     93     38     93 PF08246.5    Inhibitor_I29     Domain     1    58    58     45.9   2.8e-12    1 No_clan
  O65039.1    126    342    126    342 PF00112.16   Peptidase_C1      Domain     1   216   216    296.0   1.1e-88    1 CL0125   predicted_active_site[150,285,307]


FAQ
===

Q:  Why are the results from running hmmscan different from those I get when I run pfam_scan.pl?
A:  We group together families which we believe to have a common evolutionary ancestor in clans.  Where there 
    are overlapping matches within a clan, pfam_scan.pl will only show the most significant (the lowest E-value)
    match within the clan.  We perform the same clan filtering step on the Pfam website.  If you do want the
    script to report all the overlapping clan matches, you can use the -clan_overlap option.
