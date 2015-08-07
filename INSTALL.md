Steps you have to do to acquire APPRIS system
=============================================

Add APPRIS subversion code, http://sourceforge.net/projects/appris

	1. Export the last release of APPRIS subversion:
		svn co https://svn.code.sf.net/p/appris/code/trunk appris
	
Setting up environment vars of "appris":

	1. Add into your bash profile:
		export APPRIS_HOME="APPRIS HOME"
		source ${APPRIS_HOME}/conf/apprisrc
		source ${APPRIS_HOME}/conf/apprisrc.WS

Setting up databases for "appris":

	1. Add databases into "code/db". For more information: Read "appris/code/db/README.txt"

Setting up environment vars for "firestar" ("conf/code/fire_var.ini"):

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

Perl requirments:

	>> for global scripts of appris:
		- FindBin
		- Getopt::Long
		- Config::IniFiles
		- Bio::SeqIO
		- Bio::SearchIO
		- File::Temp
		- File::Basename
		- Data::Dumper

	>> for firestar scripts:
		- DBI
		- POSIX

	>> for inertia scripts:
		- Statistics::Descriptive

	>> for spade scripts:
		- Moose
		- Class::Load
		- Data::OptList
		- Module::Implementation
		- Class::Load::XS
		- MRO::Compat
		- Data::Printer
		- IPC::Run
		
	>> for web services:
		- CGI
		- HTTP::Status
		- Email::Valid
		- MIME::Lite
		
	>> for ensembl scripts:
		- Parse::RecDescent

BioPerl API (at least, 1.2.3)

Ensembl API

	- http://www.ensembl.org/info/docs/api/api_installation.html

Softwares:

	>> for crash scripts:
		- gawk
		
