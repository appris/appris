					
					STEPS YOU HAVE TO DO TO RUN APPRIS PIPELINE

Using appris.pl
==================

"appris.pl" executes APPRIS pipeline for a given gene.


RELEASES
============
...

* v1:

* v2:
		
* v3: 5-Apr-2013
	- Take into account the specie of given input.
	
* v4: 19-Apr-2013
	- New output of appris.

* v5: 29-May-2013
	- Take into account old results.
	
* v6: 30-May-2013
	- Add INERTIA.
	
	
REQUIREMENTS
============


Log::Handler

PBS Platform
Sun Grid Engine
				

Prepare Methods

* SLR:

	1. Uncompress source files
			tar -zxvf slr_source.tgz
			
	2. Install program (for MacOSX system)
		cd slr/src
		cp -p Makefile.osx Makefile
		make

* prank:

	1. Create directory
			mkdir prank
			
	2. Uncompress source files
			tar -zxvf prank.src.091016.tgz -C prank
			
	3. Install program
			cd prank/src
			make
			
	4. Copy binary into prank/bin
			cd ..
			mkdir bin
			cp -p src/prank bin/.
			
* PfamScan:

	1. Uncompress source files
			tar -zxvf PfamScan.tar.gz
	
	2. Install HMMER3

	3. Download Pfam data files
	
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
	
	For more details:
		ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/README
		
* CExonic:


TODO: Check requirements
Bio::Seq
Bio::Das::Lite
Bio::Tools::BPlite
Bio::Graphics
GD::SVG
exonerate
muscle
compile indexfasta and getfasta files




TODO: Fill the documentation
		Copy files from "blast_segments" and "fastaDB" directory????
		
* CORSAIR:

TODO: Check requirements. At least, 'blastpgp' program

* SignalP

	1. Uncompress source files
		tar -zxvf signalp-3.0.Linux.tar.Z
		
	2. Change $SIGNALP variable from 'signalp' file
	
	For more details, read "signalp-3.0.readme" file
	
* TargetP

	1. Uncompress source files
		tar -zxvf targetp-1.1b.Linux.tar.Z
		
	2. Change $TARGETP variable from 'targetp' file
	
	For more details, read "targetp-1.1.readme" file



* Get the HTML documentation by means pdoc script:
	softs/pdoc-1.1> perl scripts/perlmod2www.pl -source /home/jmrodriguez/projects/Encode/release_7/lib/Perl/APPRIS/ -target /home/jmrodriguez/projects/Encode/release_7/lib/Perl/pdocs/ -css data/perl.css




#### EXTERNAL SOFTWARES THAT ARE USED BY APPRIS METHODS

> FIRESTAR:
	* blastpgp:
		1. sprot_clean_trembl_clean_70
		2. fdbTptDB_$release_date (specific db)
	* hhblits:
		1. hhblits_27Jan2012 (specific db)

> MATADOR3D:
	* blastpgp:
		1. PDB

> CORSAIR:
	* blastpgp:
		1. vertebrate database (specific db)
		
> SPADE:
	* pfam_scan:
		1. Pfam-A
		2. Pfam-B

> THUMP:
	* Memsat3
	* Phobius
	* Prodiv
	* PSIBLAST

> CRASH:
	* signalp
	* targetp
	
> INERTIA:
	* kalign:
	* prank:
	* slr:
	
> CEXONIC:
	* fastasubseq:		
	* indexfasta:
	* tblastx (local):
	* getfasta:
	* exonerate:
	* formatdb	
	* pressdb:
	* databases:
		1. H.sapiens/golden_path (Formatdb for mouse genome AND Pressdb for mouse genome)
		2. M.musculus/golden_path

