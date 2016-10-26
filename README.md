Welcome to APPRIS - A system for annotating alternative splice isoforms
=======================================================================
[APPRIS] (http://appris.bioinfo.cnio.es) [1] is a system that deploys a range of computational methods to provide value to the annotations of the human genome. APPRIS also selects one of the CDS for each gene as the principal isoform.

APPRIS defines the principal variant by combining protein structural and functional information and information from the conservation of related species.

The server has been used in the context of the scale up of the [GENCODE] (http://www.gencodegenes.org/), a sub-project of the ENCODE project to annotate the Human genome but APPRIS is being used for other species:
  * Human
  * Mouse
  * Rat
  * Zebrafish

System overview
===============
The goals of the APPRIS system are to annotate alternative variants with reliable biological data and to select the primary variant for each gene. APPRIS is based on a range of complementary computational methods.

The methods in APPRIS are the following:
  * Functionally important residues, firestar
  * Protein structural information, Matador3D
  * Presence of whole protein domains, SPADE
  * Conservation against vertebrates, CORSAIR
  * Presence of whole trans-membrane helices, THUMP
  * Prediction of signal peptide and sub-cellular location, CRASH
  * Selective pressure evolution of exons, INERTIA [Note: Currently, this method is not affecting in the selection of the principal isoform]

Data access
===========
The APPRIS websites offer an export option, suitable for small amounts of data. This is the ideal option if you want a protein sequence as a FASTA file, or a JSON/GTF file of a few features of a gene or transcript. Furthermore, you can get annotation tracks of gene/transcripts in the BED format. Simply click on one of the "Export" links in the right menu, and select the output options. If you wish to extract multiple features, we recommend the following alternatives.

  * Web Services, APPRIS data can be returned remotely using web services.
  http://apprisws.bioinfo.cnio.es/

  * Downloads, APPRIS data text files.
  http://appris.bioinfo.cnio.es/#/downloads

References
==========
[1] Rodriguez JM, Maietta P, Ezkurdia I, Pietrelli A, Wesselink JJ, Lopez G, Valencia A, Tress ML.
APPRIS: annotation of principal and alternative splice isoforms. 
Nucleic Acids Res. 2013 Jan;41(Database issue):D110-7.

[2] Rodriguez JM, Carro A, Valencia A, Tress ML. APPRIS WebServer and WebServices.
Nucleic Acids Res. 2015 Jul 1;43(W1):W455-9. doi: 10.1093/nar/gkv512.

Contact
=======
This APPRIS website is powered by the Structural Computational Biology Group at
	Centro Nacional de Investigaciones Oncologicas, [CNIO] (http://www.cnio.es)
		and
	Instituto Nacional de Bioinformatica, [INB] (http://www.inab.org)

If you have questions or comments, please write to:
	Jose Manuel Rodríguez, jmrodriguez@cnio.es
	Michael Tress, mtress@cnio.es.

Installing APPRIS
=================

Steps you have to do to acquire APPRIS system

1. Clone APPRIS code:

	git clone https://github.com/appris/appris.git
	
2. Download databases for APPRIS code:

TODO!! Think where locate 70Gb of data.

3. Setting up the environment variables:

  3.1. Add in your bash profile the following lines:
  
		export APPRIS_HOME="APPRIS HOME"
		source ${APPRIS_HOME}/conf/apprisrc
		source ${APPRIS_HOME}/conf/apprisrc.WS

4. Setting up environment vars for "firestar" ("conf/code/fire_var.ini"):

	1. Change the env vars:
		[PATHS]
			home
			DB
			tmp
			AFM
		[CLUSTER_PATHS]
			home
			root
			dir
			DB
			
	2. Add FireDB database:
		database: FireDB
		user: firedb
		pwd:
					
Setting up variables of Ensembl database:

	1. Modify variables of config file that sets up Ensembl database, "conf/code/ensembldb.ini":
		ENSEMBL_CORE_REGISTRY, ENSEMBL_COMPARA_REGISTRY
	
Setting up environment vars of "appris" for each specie you will execute:
	
	1. Modify APPRIS_HOME enviroment variable that is saved in the configuration file "conf/apprisrc.*":
		APPRIS_HOME -> directory of appris
		APPRIS_WORKSPACE -> workspace of current project
		
Take into account the temporal files coming from "code/opt" programs
	
	1. SignalP and TargetP: signalp-3.0/tmp and targetp-1.1/tmp
		chmod -R +t tmp
		chmod -R 777 tmp
		

Requirements
============

- Perl requirements (recommendation, use [CPAN] https://www.perl.org/cpan.html)):

  * for global scripts of appris:
	- FindBin
	- Getopt::Long
	- Config::IniFiles
	- Bio::SeqIO
	- Bio::SearchIO
	- File::Temp
	- File::Basename
	- Data::Dumper

  * for firestar scripts:
	- DBI
	- POSIX

  * for inertia scripts:
	- Statistics::Descriptive

  * for spade scripts:
	- Moose
	- Class::Load
	- Data::OptList
	- Module::Implementation
	- Class::Load::XS
	- MRO::Compat
	- Data::Printer
	- IPC::Run
	
  * for web services:
	- CGI
	- HTTP::Status
	- Email::Valid
	- MIME::Lite
	
  * for ensembl scripts:
	- Parse::RecDescent

- [BioPerl] (http://bioperl.org/) (at least, 1.2.3)

- [Ensembl API] (http://www.ensembl.org/info/docs/api/api_installation.html)

- Softwares:

  * for crash scripts:
	- gawk
		
