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
$ docker pull appris/core
```

The `pull` command fetches the appris **image** from the **Docker registry** and saves it in our system. You can use the `docker images` command to see a list of all images on your system.
```
$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
appris/core       latest              84a2cef92a68        45 minutes ago      3.224 GB
appris/core       2017_08.v24         12c4188c29c7        50 minutes ago      3.224 GB
```

Great! Let's now run a APPRIS-Docker **container** based on this image. To do that you are going to use the `docker
run` command.

Create a working directory
--------------------------
Create a working directory where the input file for APPRIS, and the results from APPRIS will be locate (with
read/write permissions)
```
$ mkdir ${working_dir} && chmod -R og+w ${working_dir}

```

Download databases for APPRIS code
----------------------------------
Download the database files that the code needs into the **working directory**

> __WARNING__: The database file is quite big (around 35Gb). Take care!

```
$ cd ${working_dir} && \
  wget http://apprisws.bioinfo.cnio.es/archives/db/appris_db_archives.${appris_version}.tar.gz && \
  tar -xf appris_db_archives.${appris_version}.tar.gz
```

In windows, the symbolic links don't work. For that reason, you have to change them:
- Delete the symbolic links:
```
$ cd db
$ rm firestar
$ rm pdb_{version}
```

- Rename the files
```
$ mv firestar_{version} firestar
$ mv pdb_{version}.emptySeqs.fa pdb_{version}
```

> __WARNING__: Docker in Windows is limited to c:\Users folder. If you want to create a volume using a directory outside of c:\Users try this ([ref](https://stackoverflow.com/questions/33126271/how-to-use-volume-option-with-docker-toolbox-on-windows/42435077#42435077)):

>    In windows 7, I used docker toolbox. It used Virtual Box.
    1. Open virtual box
    2. Select the machine (in my case default).
    3. Right clicked and select settings option
    4. Go to Shared Folders
    5. Include a new machine folder.
    For example, in my case I have included:
        **Name**: e:\{working_dir}
        **Path**: e/{working_dir}
    6. Click and close
    7. Open "Docker Quickstart Terminal" and restart the docker machine.
    Use this command:
        $ docker-machine restart
        ...
        $ docker-machine env

Create/Download the **features** input files
--------------------------------------------
Create (or download) the input data files.

> __Note:__ For more information, read the [FEATURES readme file](FEATURES.md)

Running APPRIS-Docker container
-------------------------------
The next step is to run the *__appris/core__* image mounting the database directory and the _working directory_

```
$ docker run -itd \
    -v ${database_dir}:/opt/appris/db \
    -v ${working_dir}:/home/appris/ws \
    appris/code
```

Now, we run APPRIS commands in a running container with _appris_ user.

```
$ docker exec --user appris -it {CONTAINER_ID} \
    bash -c "source /opt/appris/conf/apprisrc.docker && \
            appris_run_appris \
                -c /home/appris/ws/test.ws.env \
                -m firestar,matador3d,matador3d2,spade,corsair,thump,crash,proteo,appris \
                -l info"
```

> __Note:__ For more detail, read the related section in the [DOCKER readme file](DOCKER.md#running-appris-docker-container)


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
	- DBD::mysql
	- POSIX
	- IO::String
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
	- XML::DOM
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
    wget http://apprisws.bioinfo.cnio.es/archives/db/appris_db_archives.${appris_version}.tar.gz && \
    tar -xf appris_db_archives.${appris_version}.tar.gz

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

	3. Setting up the environment variables in the file "${APPRIS_HOME}/conf/code/firestar.ini":

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

5. Take into account the temporal files coming for "${APPRIS_HOME}/code/opt" programs SignalP and TargetP
(__DEPRECATED__):
	```
	cd signalp-3.0/tmp and targetp-1.1/tmp
	chmod -R +t tmp
	chmod -R 777 tmp
	```

5. Setting up the variables for Ensembl database in the "${APPRIS_HOME}/conf/code/ensembldb.ini" (__DEPRECATED__):
	```
	ENSEMBL_CORE_REGISTRY, ENSEMBL_COMPARA_REGISTRY
	```


Build the images of APPRIS
==========================


Build the image of APPRIS/core
------------------------------
Now that you have your Dockerfile, you can build your image. The docker build command does the heavy-lifting of creating a docker image from a Dockerfile.

```
$ docker build -t appris/core -f build/appris_core.dockerfile .
```
The following files are required:

+ build/appris_core.dockerfile
+ build/setup.sh
+ build/entrypoint.sh
+ build/FireDB_{data_version}.sql.gz


Build the image of APPRIS/server
------------------------------
Now that you have your Dockerfile, you can build your image. The docker build command does the heavy-lifting of creating a docker image from a Dockerfile.

```
$ docker build -t appris/server -f build/appris_server.dockerfile .
```

### Download data files for APPRIS server that server needs into the {data_directory}
We have to ways to obtain the data files:

1. Download the data files from external repository
```
$ wget http://apprisws.bioinfo.cnio.es/archives/data/appris_data_archives.${appris_version}.tar.gz && \
  tar -xf appris_data_archives.${appris_version}.tar.gz
```
> __WARNING__: The data file is quite big (around 7Gb). Take care!

2. Where *__appristools__* saves the annotations


### Import databases of APPRIS annotations
1. Run the *__appris/server__* image mounting the data directory
```
$ docker run -itd \
    -v ${data_dir}:/opt/appris/data \
    appris/server
```

2. Import the databases files within the Docker image

> __WARNING__: You will need around 50Gb in the system to import database files. Think on that!

    2.1 Run a *_bash_* in the running *__appris/server__* container
    ```
    $ docker exec --user appris -it {CONTAINER_ID} bash
    ```

    2.1 Import databases
    ```
    $ docker exec --user appris -it {CONTAINER_ID} \
        bash -c "source /opt/appris/conf/apprisrc.docker && \
                perl /opt/appris/docker/scripts/import_appris_dbs.pl \
                    -c /opt/appris/conf/config_{appris_version}.json \
                    -d /opt/appris/data \
                    -l info"
    ```