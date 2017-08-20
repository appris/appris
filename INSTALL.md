Installation from the Docker repository
=======================================
We recommend to use our Docker repository to avoid conflicts.

Install Docker program
----------------------

The _getting started_ guide on Docker has detailed instructions for setting up Docker on [Mac](https://docs.docker.com/docker-for-mac/), [Linux](https://docs.docker.com/engine/installation/linux/) and [Windows](https://docs.docker.com/docker-for-windows/).

Running APPRIS-Docker container
-------------------------------
You are going to run an __APPRIS__ container (with an Ubuntun linux) on your system and get a taste of the `docker run` command.

To get started, let's run the following in our terminal:
```
$ docker pull appris/appris
```

The `pull` command fetches the appris **image** from the **Docker registry** and saves it in our system. You can use the `docker images` command to see a list of all images on your system.
```
$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
appris/appris       latest              84a2cef92a68        45 minutes ago      3.224 GB
appris/appris       2017_08.v24         12c4188c29c7        50 minutes ago      3.224 GB
```

Great! Let's now run a APPRIS-Docker **container** based on this image. To do that you are going to use the `docker
run` command.

> __Note:__ For more information about how to run *__APPRIS-Docker container__*,
read the [DOCKER readme file](DOCKER.md)


Installation from Scratch
=========================

Requirements
------------

- Ubuntu Linux x86_64 machine

- Softwares:
  	```
  	- gawk
  	- git
  	- npm
  	- nodejs
  	```

- Perl requirements (recommendation, use [CPAN] (https://www.perl.org/cpan.html)):

	```
	- FindBin
	- Getopt::Long
	- Config::IniFiles
	- Bio::SeqIO
	- Bio::SearchIO
	- File::Temp
	- File::Basename
	- Data::Dumper
	- JSON
	- DBI
	- POSIX
	- Statistics::Descriptive
	- Moose
	- Class::Load
	- Data::OptList
	- Module::Implementation
	- Class::Load::XS
	- MRO::Compat
	- Data::Printer
	- IPC::Run
	- CGI
	- HTTP::Status
	- Email::Valid
	- MIME::Lite
	- Parse::RecDescent
	```
	Note: See the lib/appris_perllib/Makefile.PL file

- [MySQL Client] (http://dev.mysql.com/doc/refman/5.7/en/linux-installation.html)

- [BioPerl] (http://bioperl.org/) (at least, 1.2.3)

- [Ensembl API] (http://www.ensembl.org/info/docs/api/api_installation.html)


Installation
------------

Steps you have to do to acquire APPRIS system

1. Clone APPRIS code:
	```
	git clone https://github.com/appris/appris.git
	```

2. Setting up the environment variables:

	1. Add in your bash profile the following lines:
		```
		export APPRIS_HOME="APPRIS HOME"
		source ${APPRIS_HOME}/conf/apprisrc
		source ${APPRIS_HOME}/conf/apprisrc.WS
		```

	2. Setting up the environment variables for each execution in the sample configuration file "${APPRIS_HOME}/conf/scripts/apprisrc.*"

3. Download databases for APPRIS code:
	```
	wget
	cd ${APPRIS_HOME}
	tar -xf appris_local_db.${date_version}.tar.gz
	```

4. Setting up 'firestar' method:

	1. Create FireDB database:
		```
		mysql> create user 'firedb'@'localhost';
		mysql> grant all on FireDB.* to 'firedb'@'localhost';
		mysql> create database FireDB;
		```

	2. Import FireDB database:
		```
		mysql FireDB -h localhost -u firedb < ${APPRIS_HOME}/db/FireDB_*.sql
		```

	3. Setting up the environment variables in the file "${APPRIS_HOME}/conf/code/fire_var.ini":

		1. Include database variables:
			```
			database: FireDB
			user: firedb
			pwd:
			```

		2. Change the env vars:
			```
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
			```

5. Setting up the variables for Ensembl database in the "${APPRIS_HOME}/conf/code/ensembldb.ini":
	```
	ENSEMBL_CORE_REGISTRY, ENSEMBL_COMPARA_REGISTRY
	```

6. Take into account the temporal files coming for "${APPRIS_HOME}/code/opt" programs SignalP and TargetP:
	```
	cd signalp-3.0/tmp and targetp-1.1/tmp
	chmod -R +t tmp
	chmod -R 777 tmp
	```


