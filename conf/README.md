Conf files
==========
Config files used in the APPRIS system.

This directory contains the following subdirectories:

```
   conf
	|
	|___ apprisrc
	|
	|___ apprisrc.WS
	|
	|___ config.json
	|
	|___ methods.json

            /* directories */
	|
	|___ code
          |
          |__ /config files for the source code of APPRIS pipeline/
	|
	|___ db
          |
          |__ /config files for the database of APPRIS pipeline/
	|
	|___ ws
          |
          |__ /config files for the APPRIS web-server/
	|
	|___ scripts
            |
            |__ /config files for the scripts that execute the APPRIS pipeline/
	|
	|___ ucsc
           |
           |__ /config files for the trackHubs of UCSC Genome Browser (bed files)/
```

Config files
============
The following files are the scripts to operate with the code of APPRIS pipeline.

+ *__apprisrc__*,
    Bash file with the environment variables used by the APPRIS pipeline/database.
    > Note: Exists several files for the configuration in the production server (pro), development server (dev),
    and computational machine (loc), and docker machine (docker).

+ *__apprisrc.WS__*,
    Bash file with the environment variables used by the APPRIS server.

+ *__config.json__*,
    File that sets up the APPRIS server and the involved databases.

+ *__methods.json__*,
    File that describes the methods in the APPRIS server.

### code/

+ *__pipeline.ini__*,
    File with the parameters that will be applied into the methods in the APPRIS pipeline.

+ *__firestar.ini__*,
    Fle with the *_internal_* parameters that will be applied to *__firestar__* method in the APPRIS pipeline.

+ *__corsair_species.json__*,
    File with the *_internal_* parameters that will be applied to *__CORSAIR__* method in the APPRIS pipeline.

+ *__ensembl.ini__*,
    File with Ensembl configuration for the *__ensembl/compara__* scripts.
    >Note: *__Obsoleted__*

### db/

+ *__apprisdb.sql__*,
    SQL structure of APPRIS databases.

+ *__report.sql__*, *__annotation.sql__*, *__g_annotation.sql__*,
    SQL queries to retrieve information from APPRIS databases.
    >Note: *__Obsoleted__*

### ws/

+ *__apprisdb.ini__*,
    File that sets up the databases within the APPRIS server.

+ *__config.json__*,
    File that sets up the APPRIS server and the involved databases. By default,
    it is a symbolic link to the main *_config.json_*.

+ *__methods.json__*,
    File that describes the methods in the APPRIS server. By default, it is a symbolic link to the main *_methods.json_*.

+ *__cluster.ini__*,
    File that controls the cluster configuration of APPRIS pipeline/system (server side).
    >Note: *__Obsoleted__*

### scripts/

+ *__apprisdb.ini__*,
    File that sets up the databases within the APPRIS pipeline.

+ *__apprisrc.{species}.{gene_set}__*,
    Set-up the script (*_appristools_*) that executes the APPRIS pipeline depending on the input species and gene sets.

+ *__email.ini__*,
    File that sets up the Email configuration.

+ *__ensembldb.ini__*,
    File with Ensembl configuration for the *__ensembl/compara__* scripts.
    >Note: *__Obsoleted__*

+ *__cluster.ini__*,
    File that controls the cluster configuration of APPRIS pipeline/system (client side).
    >Note: *__Obsoleted__*

+ *__apprisdiff.ini__*,
    >Note: *__Deprecated__*

### ucsc/

+ *__bedAs.{methods}__*,
    Config files to the create of track Hubs (in bed12 format) for the UCSC Genome Browser.

