Features files
==============
Files with the inputs for the APPRIS pipeline.

There are three (four) different gene sets:
 1. Ensembl/GENCODE
 2. RefSeq Gene
 3. UniProt
 4. Intersection

In a general way, you can download the *_gene data files_* using the *__appristools_srv__* script.
```
appristools_srv -p 1 -c ws/config.json
```

In the creation of the intersection gene set, you have to execute the *__appris_gs_create__* script.
```
 appris_gs_create \
 -s danio_rerio \
 -ie /home/jmrodriguez/projects/APPRIS/features/danio_rerio/e87 \
 -ir /home/jmrodriguez/projects/APPRIS/features/danio_rerio/rs105 \
 -iu /home/jmrodriguez/projects/APPRIS/features/danio_rerio/up201610 \
 -o /home/jmrodriguez/projects/APPRIS/features/danio_rerio/a1
 ```
 > Note: There is specific sccript for the manual cross-references in human *__appris_gs_create_human__*


Download the features file in manual way
========================================

Files you need for each gene sets:

1. __Ensembl/GENCODE__
2. __RefSeq Gene__
3. __UniProt__

    1. Download the following files using the web site of UniProt (from the PROTEOME section!!!)
        1. uniprot-proteome.fasta (FASTA canonical & isoforms)
        2. uniprot-proteome.txt   (Text format)
        3. uniprot-proteome.tab   (Tab-separated)

    2. Gunzip files

    3. Create links:
        1. uniprot-proteome.fasta
        2. uniprot-proteome.txt
        3. uniprot-proteome.tab

    4. Run the script
    ```perl
    perl scripts/add_extraVals_into_UPseq.pl \
    -f uniprot-proteome.fasta \
    -t uniprot-proteome.tab \
    -d uniprot-proteome.txt \
    -o uniprot-proteome.extra.fasta \
    --loglevel=info \
    --logfile=add_extraVals_into_UPseq.log
    ```

4. __Intersection__



