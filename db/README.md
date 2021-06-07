Database files of Code's APPRIS
===============================

### For firestar: [FMI](#firestar)
  + BLOSUM62
  + firestar_{version}/fdbTptDB_{version}
  + firestar_{version}/chads_{version}
  + firestar_{version}/hhblits_{version}_a3m_db
  + firestar_{version}/hhblits_{version}.cs219
  + firestar_{version}/hhblits_{version}_hhm_db
  + firestar_{version}/sprot_clean_trembl_clean_70*
  + local Mysql __FireDB__ database

> Note: Deprecated files:
    firestar_{version}/nr20_12Aug11_a3m_db*
    firestar_{version}/nr20_12Aug11.cs219*
    firestar_{version}/nr20_12Aug11_hhm_db*

### For Matador3D: [FMI](#matador3d)
  + components.cif
  + pdb_{version}.fasta
  + pdb_{version}.emptySeqs.fa
  + pdb_{version}
  + pdb_{version}.phr
  + pdb_{version}.psq
  + pdb_{version}.pin

### For Matador3D2: [FMI](#matador3d2)
  + pdb_70.with_pdb_seq
  + pdb_70.with_pdb_seq.h3m
  + pdb_70.with_pdb_seq.h3p
  + pdb_70.with_pdb_seq.h3f
  + pdb_70.with_pdb_seq.h3i

### For SPADE: [FMI](#spade)
  + pfam_{version}/active_site.dat
  + pfam_{version}/Pfam.version
  + pfam_{version}/Pfam-A.hmm.dat
  + pfam_{version}/Pfam-A.hmm
  + pfam_{version}/Pfam-A.hmm.h3p
  + pfam_{version}/Pfam-A.hmm.h3m
  + pfam_{version}/Pfam-A.hmm.h3i
  + pfam_{version}/Pfam-A.hmm.h3f

### For CORSAIR: [FMI](#corsair)
  + refseq_vert
  + refseq_vert.pin
  + refseq_vert.phr
  + refseq_vert.psq
  + refseq_invert
  + refseq_invert.pin
  + refseq_invert.phr
  + refseq_invert.psq

### For THUMP: [FMI](#thump)
  + sprot_clean_trembl_clean_90
  + sprot_clean_trembl_clean_90.*.phr
  + sprot_clean_trembl_clean_90.*.pin
  + sprot_clean_trembl_clean_90.*.psq
  + sprot_clean_trembl_clean_90.pal

### For PROTEO: [FMI](#proteo)
  + proteo_20210604.csv

### For TRIFID: [FMI](#trifid)
  There are typically multiple TRIFID data files, with each one containing data
  for a specific species, assembly, feature set, and TRIFID prediction set. For
  example, the TRIFID file created on 12th May 2021 for GENCODE v37 annotations
  of human assembly `GRCh38` would be named `trifid_GRCh38_g37_20210512.tsv.gz`.


Creation of the Code's Databases
================================

firestar
--------
1. Download and create the databases for firestar (*__Paolo Maietta__* developer,
http://firedb.bioinfo.cnio.es/Php/Help.php)
```
    - sprot_clean_trembl_clean_70*
    - fdbTptDB_*
    - chads_*
    - hhblits_*
    - nr20_*
```
2. Install FireDB database:
	mysql FireDB -h localhost -u firedb < FireDB_*.sql


Matador3D
---------

1. Generate Matador3D BLAST database

To generate a new Matador3D BLAST database,
execute a command such as the following:
```shell
python appris/scripts/prep_pdbsum_blastdb.py --output-dir appris/db
```
This downloads PDB sequence data from [PDBsum](https://www.ebi.ac.uk/pdbsum/),
processes it, and creates a BLAST database in the specified output directory.

2. Configure Matador3D BLAST database

In section `[MATADOR3D_VARS]` of the config file `appris/conf/code/pipeline.ini`,
set `db` to the name of the new database. For example, given a BLAST database in
`appris/db` whose filenames start with `pdbsum_20210427`, the Matador3D config
section might look like this:
```ini
[MATADOR3D_VARS]
  name=matador3d
  program=blastpgp
  db=pdbsum_20210427
  evalue=0.01
  cache=yes
```

Matador3D2
----------
1. Create the following PDB files. Developed by *__Juan Rodríguez__* (read **actualizar_db_matador3D2.odt** file)
```
    - pdb_70.with_pdb_seq*
```


SPADE
-----
Extract the databases from Pfam.

The URL from the PfamScan code is in ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/

1. Create directory with date version. Eg. Pfam_201706
```
	$ mkdir Pfam_201706 && cd Pfam_201706
```
2. Download release note file:
```
	$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
	$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
```
3. Get the Pfam database in a raw dir (and unzip them):
```
	$ mkdir raw && cd raw && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam24.0/Pfam-A.hmm.* && gunzip Pfam-A.hmm.*.gz
```
4. 	To use the active site option you will also need to	download the active site alignments (and unzip them):
```
	$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz && gunzip active_site.dat.gz
```
5. Grab and install HMMER, NCBI BLAST and Bioperl, and make sure your paths etc are set up properly.
	TO REMAIND YOU, the last version of HMMER has to be: 3.1b2
6. Create Pfam database in HMM format with a given list of domains:
```
	$ perl scripts/spade/selectGivenDomains.pl \
		-d scripts/spade/AllPfamDomainData.20170621.ID.txt \
		-i1 pfam_201706/raw/Pfam-A.hmm \
		-i2 pfam_201706/raw/Pfam-A.hmm.dat \
		-o1 pfam_201706/Pfam-A.hmm \
		-o2 pfam_201706/Pfam-A.hmm.dat
```
7. You will need to generate binary file for Pfam-A.hmm by running the following command:
```
	$ cd pfam_201706 && hmmpress Pfam-A.hmm
```



CORSAIR
-------
__RefSeq__ databases ["vertebrate_mamalian"+"vertebrate_other", "invertebrate"]  (ftp://ftp.ncbi.nlm.nih.gov/refseq/release/)

1. Prepare workspace for the data files
```
  $ mkdir refseq_{dataversion}
  $ cd refseq_{dataversion}
  $ mkdir raw
  $ cd raw
```
2. Get RefSeq database for 'vertebrate' and 'invertebrate' from
```
  $ wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.*.protein.faa.gz
  $ wget ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other.*.protein.faa.gz
  $ wget ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate.*.protein.faa.gz
```
3. Unzip them
```
  $ gzip -d vertebrate_*
  $ gzip -d invertebrate.*
```
4. Concatenate them
```
  $ cat vertebrate_mammalian.* vertebrate_other.* >> ../refseq_vert
  $ cat invertebrate.* >> ../refseq_invert
  $ cat vertebrate_mammalian.* vertebrate_other.* invertebrate.* >> ../refseq
```
5. Index database
```
  $ cd ..
  $ formatdb -i refseq_vert -p T
  $ formatdb -i refseq_invert -p T
  $ formatdb -i refseq -p T
```

Another way to download the databases:

1. Prepare workspace for the data files
```
  $ mkdir refseq_{dataversion}
  $ cd refseq_{dataversion}
  $ mkdir raw
  $ cd raw
```
2. Get databases
```
  $ perl appris/db/scripts/corsair/download_refseq.pl -c appris/conf/code/corsair_alt.diverge_time.human.json -o db/refseq_{dataversion}/raw --loglevel=debug
```
3. Unzip them
```
  $ gzip -d *.gz
```
4. Concatenate them
```
  $ cat *.fasta >> ../refseq
```
5. Index database
```
  $ cd ..
  $ formatdb -i refseq -p T
```

> Note: There are some species that have been downloaded manually.


__UniProt__ databases ["vertebrate", "invertebrate"] (ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/):

1. Prepare workspace for the data files
```
  $ mkdir uniprot_{dataversion}
  $ cd uniprot_{dataversion}
  $ mkdir raw
  $ cd raw
```
2. Get UniProt databases
```
  $ perl appris/db/scripts/corsair/download_uniprot.pl -c appris/conf/code/corsair_alt.diverge_time.human.json -o db/uniprot_{dataversion}/raw --loglevel=debug
```
3. Unzip them
```
  $ gzip -d *.gz
```
4. Concatenate them
```
  $ cat *.fasta >> ../uniprot
```
5. Index database
```
  $ cd ..
  $ formatdb -i uniprot -p T
```


THUMP
-----
The **sprot_clean_trembl_clean_70** database (fasta and 'psq', 'pin', and 'phr') formatted for Blast search has been
developed by *__José María Fernández__*
> Note: sprot_clean_trembl_clean_90* files (fasta and 'psq', 'pin', and 'phr'). You can find it in the /drives/databases/BlastDB folder


PROTEO
------
Developed by *__Iakes Ezcurdia Garmendia__* and *__Michael Tress__*.

1. Generate PROTEO data file

Given a file `peptide_results.tsv` of data output from a proteomics peptide
search workflow, the following command will generate a PROTEO data file:
```shell
python prep_proteo_data.py -i peptide_results.tsv -o proteo_20210604.csv
```

2. Configure PROTEO data file

In section `[PROTEO_VARS]` of the config file `appris/conf/code/pipeline.ini`,
set `db` to the name of the new PROTEO data file. For example, with PROTEO data
file `appris/db/proteo_20210604.csv`, the PROTEO config section might look like
this:
```ini
[PROTEO_VARS]
  name=proteo
  db=proteo_20210604.csv
```

TRIFID
------

In typical use, TRIFID predictions are downloaded from the
[TRIFID GitLab repository](https://gitlab.com/bu_cnio/trifid).
When preparing TRIFID data in this way, it is only necessary to
provide the key details for the TRIFID prediction set (i.e. its
genome assembly, and the source name and version number of its
genome annotation).

First, it is necessary to download the relevant annotation data
and translation sequence files. For example, the following
commands would download these files for GENCODE v35:
```shell
curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz \
  --output gencode.v35.annotation.gtf.gz

curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.pc_translations.fa.gz \
   --output gencode.v35.pc_translations.fa.gz
```

Then, the Python script `prep_trifid_data.py` may be
used to prepare a TRIFID data file for input to APPRIS.
Again, using GENCODE v35 as an example:
```shell
python prep_trifid_data.py \
  --assembly GRCh38 --source-name gencode --source-version 35 \
  --transl-file gencode.v35.pc_translations.fa.gz \
  --annot-file gencode.v35.annotation.gtf.gz \
  --output-dir appris/db/trifid/homo_sapiens/GRCh38/g35
```

A TRIFID data file must be configured for each dataset. This may
be done, for example, by specifying an individual TRIFID file with
the `-trifid` parameter of `appris_run_appris`, or by setting the
TRIFID release for a given dataset in the global JSON config file
used by `appristools`.


Create database for APPRIS-Docker
=================================
To execute a *__Doker__* container, APPRIS needs the database files. Here we show an example how to copy the files for the version **2017_08.v24**
```
$ tar -cf - \
  db/sprot_clean_trembl_clean_90* \
  db/refseq_vert* \
  db/refseq_invert* \
  db/APD8_a2_sortu.2.csv \
  db/firestar_22Aug2013 \
  db/firestar \
  db/pdb_70.with_pdb_seq* \
  db/components.cif \
  db/BLOSUM62 \
  db/pdb_20170420.* \
  db/pdb_20170420 \
  db/pfam_201706 \
  db/README.md \
  db/scripts | gzip -9c > ../appris_db_archives/appris_db_archives.2017_08.v24.tar.gz
```
