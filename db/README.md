#####################
# LIST OF DATABASES #
#####################

	sprot_clean_trembl_clean_90
	active_site.dat
	
	# for pfam method
	Pfam-A.hmm
	Pfam-B.hmm
		
	# for matador method
	pdb.*
	
	# for corsair method
	refseq_vert.*
	refseq_vert_files
	
	# for firestar method
	BLOSUM62
	firestar_{version}/fdbTptDB_{version}
	firestar_{version}/chads_{version}
	firestar_{version}/hhblits_{version}_a3m_db
	firestar_{version}/hhblits_{version}.cs219
	firestar_{version}/hhblits_{version}_hhm_db
	firestar_{version}/sprot_clean_trembl_clean_70
	
	firestar_{version}/nr20_12Aug11_a3m_db # DEPRECATED
	firestar_{version}/nr20_12Aug11.cs219  # DEPRECATED
	firestar_{version}/nr20_12Aug11_hhm_db # DEPRECATED
	
	# for proteo method
 	APD8_a2_sortu.2.csv
		
############################
# DESCRIPTION OF DATABASES #
############################


## PDB created by Jose Maria:

1. Download wwwPDB. Check the web: https://www.wwpdb.org/download/downloads. As of today, the way to download the files is:
	
	$ rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./wwwPDB	
	 
2. Download the 'monomers components' from wwwPDB:

	$ wget ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
		
3. Execute the script (the script needs Java 1.8):

	$ touch ficherovacio.txt && touch ficherovacio2.txt
	
	$ java -jar scripts/matador3d/GOPHERPrepare.jar -s ficherovacio.txt ficherovacio2.txt wwwPDB pdb_{DATE_VERSION}.fasta components.cif /tmp/matador3d_update_wwwPDB >& /tmp/matador3d_update_wwwPDB.log
	
4. Then, discard empty sequence
	
	$ perl scripts/matador3d/discardEmptySeqsPDB.pl pdb_{DATE_VERSION}.fasta 1> pdb_{DATE_VERSION}.emptySeqs.fa
	
	$ ln -s pdb_{DATE_VERSION}.emptySeqs.fa pdb_{DATE_VERSION}
	
5. Index database
	
	$ formatdb -i pdb_{DATE_VERSION} -p T
	
6. Delete tmp files

	$ rm -rf /tmp/matador3d_update_wwwPDB && rm /tmp/matador3d_update_wwwPDB.log && rm ficherovacio.txt && rm ficherovacio2.txt && rm formatdb.log



## RefSeq Vertebrate database comes from "vertebrate_mamalian" and "vertebrate_other" (ftp://ftp.ncbi.nlm.nih.gov/refseq/release/)

1. Get RefSeq database for 'vertebrate' and 'invertebrate' from

	$ wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.*.protein.faa.gz
	$ wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_other.*.protein.faa.gz

	$ wget ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate.*.protein.faa.gz
	
2. Unzip them
	
	$ gzip -d vertebrate_*
	
	$ gzip -d invertebrate.*
	
3. Concatenate them
	
	$ cat vertebrate_mammalian.* vertebrate_other.* >> refseq_vert
	
	$ cat invertebrate.* >> refseq_invert

4. Index database

	$ formatdb -i refseq_vert -p T
	
	$ formatdb -i refseq_invert -p T



## PfamScan, ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/

1. Create directory with date version. Eg. Pfam_201706

	$ mkdir Pfam_201706 && cd Pfam_201706 

2. Download release note file:

	$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
	$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
		
3. Get the Pfam database in a raw dir (and unzip them):

	$ mkdir raw && cd raw && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam24.0/Pfam-A.hmm.* && gunzip Pfam-A.hmm.*.gz

4. 	To use the active site option you will also need to	download the active site alignments (and unzip them):

	$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz && gunzip active_site.dat.gz

5. Grab and install HMMER, NCBI BLAST and Bioperl, and make sure your paths etc are set up properly.
	TO REMAIND YOU, the last version of HMMER has to be: 3.1b2
	
6. Create Pfam database in HMM format with a given list of domains:

	$ perl scripts/spade/selectGivenDomains.pl \
		-d AllPfamDomainData.20170529.ID.txt \
		-i1 ../../pfam_201706/raw/Pfam-A.hmm \
		-i2 ../../pfam_201706/raw/Pfam-A.hmm.dat \
		-o1 ../../pfam_201706/Pfam-A.hmm \
		-o2 ../../pfam_201706/Pfam-A.hmm.dat
	
7. You will need to generate binary file for Pfam-A.hmm by running the following command:
	
	$ hmmpress Pfam-A.hmm
	

    
## sprot_clean_trembl_clean_70, fdbTptDB_*, chads_*, hhblits_*, nr20_*, FireDB for firestar

1. Install FireDB database:
	mysql FireDB -h localhost -u firedb < FireDB_*.sql


	
## PDB for Matador3D2


TODO!!! Check the documentation of Juan Rodriguez in 'scripts/matador3d2':


	
		
