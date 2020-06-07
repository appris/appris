Code
====
Source code of APPRIS pipeline

This directory contains the following subdirectories:

```
   opt /* External programs */
	|
	|___ ...


   src
	|
	|___ appris.pl
	|
	|___ cappris.pl /* obsoleted */
	|
	|___ iappris.pl

            /* directories */
	|
	|___ appris
	|
	|___ corsair
	|
	|___ crash
	|
	|___ ensembl /* obsoleted */
	|
	|___ firestar
	|
	|___ geneset
	|
	|___ inertia /* obsoleted */
	|
	|___ matador3d
	|
	|___ matador3d2
	|
	|___ proteo
	|
	|___ spade
	|
	|___ thump
	|
	|___ trifid
	|
	|___ ucsc /* obsoleted */
	|
	|___ uniprot
```

Source code files
=================
The following files are the scripts to operate with the code of APPRIS pipeline.

+ *__appris.pl__*,
    It executes APPRIS pipeline for a given gene. For more detail, execute the script (*_perl appris.pl_*) to see
    the needed parameters.

+ *__iappris.pl__*,
    It inserts the annotations of a gene from APPRIS pipeline. For more detail, execute the script (*_perl iappris.pl_*) to see the needed parameters.

+ *__cappris.pl__* (*__Obsolete__*),
    It checks the annotations from APPRIS pipeline.

