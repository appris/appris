Description of binary files
===========================

The following binary files help you to execute the multiple operations evolving in the APPRIS system.

 * [Files that execute the APPRIS pipeline (appristools*)](#files-that-execute-the-appris-pipeline-appristools-)

 * [Files that operate with the APPRIS database (appris\_db\_*)](#files-that-operate-with-the-appris-database-appris-_db-_-)

 * [Files that retrieves information from APPRIS annotations (appris\_retrieve\_*)](#files-that-retrieves-information-from-appris-annotations-appris-_retrieve-_-)

 * [Files that creates the APPRIS gene sets (appris\_gs\_*)](#files-that-creates-the-appris-gene-sets-appris-_gs-_-)



Files that execute the APPRIS pipeline (appristools*)
-----------------------------------------------------

* __appristools__, executes all steps in the APPRIS pipeline from a config file
    1. Executes APPRIS pipeline
    2. Retrieves the list of principal isoforms and Compare stats between versions
    3. Inserts the annotations into database
    4. Retrieves the data files of methods

```
appristools -p 1234 -c ws/config.json -m fmsctrpa -e 'example@appris-tools.org' -f gtf

appristools -p 1234 -d conf/scripts/apprisrc.Hsap -m fmsctrpa -e 'example@appris-tools.org' -f 'gtf,bed,bed12'

```

* __appristools_srv__, executes the steps that operate with the APPRIS server
    1. Download gene data files
	2. Copy gene data to workspace space
	3. Upload annotation files to server
	4. Upload gene data files to server

```
appristools_srv -p 3 -r 2016_06.v17  -n changelog.md -c ws/config.json
```


Files that retrieves information from APPRIS annotations (appris\_retrieve\_*)
------------------------------------------------------------------------------

* __appris_retrieve_main_data__, retrieves statistics and data from APPRIS methods

* __appris_retrieve_method_data__, retrieves annotations from APPRIS methods. These files will be downloaded from the
 website.

> __Note:__ Theses files are executed within __appristools__ bin file.


Files that retrieves information from APPRIS annotations (appris\_retrieve\_*)
------------------------------------------------------------------------------

* __appris_retrieve_main_data__, retrieves statistics and data from APPRIS methods

* __appris_retrieve_method_data__, retrieves annotations from APPRIS methods. These files will be downloaded from the
 website.

> __Note:__ Theses files are executed within __appristools__ bin file.


Files that operate with the APPRIS database (appris\_db\_*)
-----------------------------------------------------------

* __appris_db_create__, creates the APPRIS database

```
appris_db_create -h
```

* __appris_db_delete__, removes APPRIS databases

```
appris_db_delete -h
```

* __appris_db_import__, imports APPRIS database from mysql dump file

```
appris_db_import -d appris_homo_sapiens_gencode_19_dev -h localhost -u appris -i appris_db_dump.mysql.gz
```

* __appris_db_backup__, creates mysql dump file for an APPRIS database

```
appris_db_backup -d appris_homo_sapiens_gencode_19_dev -h localhost -u appris -o appris_db_dump.mysql.gz
```


Files that creates the APPRIS gene sets (appris\_gs\_*)
-------------------------------------------------------

* __appris_gs_create__, creates the gene set for APPRIS from the cross-reference of Ensembl, RefSeq and UniProt
```
 appris_gs_create \
 -s danio_rerio \
 -ie /home/jmrodriguez/projects/APPRIS/features/danio_rerio/e87 \
 -ir /home/jmrodriguez/projects/APPRIS/features/danio_rerio/rs105 \
 -iu /home/jmrodriguez/projects/APPRIS/features/danio_rerio/up201610 \
 -o /home/jmrodriguez/projects/APPRIS/features/danio_rerio/a1
 ```
> __Note:__ the gene sets of human have been obtained manually. For that reaso,
the script is the following __appris_gs_create_human__


