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
  + proteo_20200604.csv



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
1. Download wwwPDB. Check the web: https://www.wwpdb.org/download/downloads. As of today, the way to download the files is:
```
    $ rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./wwwPDB
```
2. Download the 'monomers components' from wwwPDB:
```
	$ wget ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
```
3. Execute the script (the script needs Java 1.8) (created by *__José María Fernández__*):
```
	$ touch ficherovacio.txt && touch ficherovacio2.txt
	
	$ java -jar scripts/matador3d/GOPHERPrepare.jar \
	    -s ficherovacio.txt \
	    ficherovacio2.txt \
	    wwwPDB pdb_{DATE_VERSION}.fasta \
	    components.cif \
	    /tmp/matador3d_update_wwwPDB >& /tmp/matador3d_update_wwwPDB.log
```
4. Then, discard empty sequence
```
	$ perl scripts/matador3d/discardEmptySeqsPDB.pl pdb_{DATE_VERSION}.fasta 1> pdb_{DATE_VERSION}.emptySeqs.fa
	
	$ ln -s pdb_{DATE_VERSION}.emptySeqs.fa pdb_{DATE_VERSION}
```
5. Index database
```
	$ formatdb -i pdb_{DATE_VERSION} -p T
```
6. Delete tmp files
```
	$ rm -rf /tmp/matador3d_update_wwwPDB && rm /tmp/matador3d_update_wwwPDB.log && rm ficherovacio.txt && rm ficherovacio2.txt && rm formatdb.log
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
Developed by *__Iakes Ezcurdia Garmendia__* and *__Michael Tress__*




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
