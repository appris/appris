Features files
==============
Files with the inputs for the APPRIS pipeline.

There are three (four) different gene sets:
 1. Ensembl/GENCODE
 2. RefSeq Gene
 3. UniProt
 4. Intersection

Files you need for each gene sets:

1.


2.


3.

In a general way, you can download the *_gene data files_* using the *__appristools_srv__* script.
```
appristools_srv -p 1 -c ws/config.json
```


4. In the creation of the intersection gene set, you have to execute the *__appris_gs_create__* script.
```
 appris_gs_create \
 -s danio_rerio \
 -ie /home/jmrodriguez/projects/APPRIS/features/danio_rerio/e87 \
 -ir /home/jmrodriguez/projects/APPRIS/features/danio_rerio/rs105 \
 -iu /home/jmrodriguez/projects/APPRIS/features/danio_rerio/up201610 \
 -o /home/jmrodriguez/projects/APPRIS/features/danio_rerio/a1
 ```
 > Note: There is specific sccript for the manual cross-references in human *__appris_gs_create_human__*