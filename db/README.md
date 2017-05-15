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
		
>> PDB created by Jose Maria:

1. Download wwwPDB. Check the web: https://www.wwpdb.org/download/downloads. As of today, the way to download the files is:
	
	rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./wwwPDB	
	 
2. Download the 'monomers components' from wwwPDB:

	wget ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
		
3. Execute the script (the script needs Java 1.8):

	touch ficherovacio.txt && touch ficherovacio2.txt
	
	java -jar scripts/matador3d/GOPHERPrepare.jar -s ficherovacio.txt ficherovacio2.txt wwwPDB pdb_{DATE_VERSION}.fasta components.cif /tmp/matador3d_update_wwwPDB >& /tmp/matador3d_update_wwwPDB.log
	
4. Then, discard empty sequence
	
	perl scripts/matador3d/discardEmptySeqsPDB.pl pdb_{DATE_VERSION}.fasta 1> pdb_{DATE_VERSION}.emptySeqs.fa
	
	ln -s pdb_{DATE_VERSION}.emptySeqs.fa pdb_{DATE_VERSION}
	
5. Index database
	
	formatdb -i pdb_{DATE_VERSION} -p T
	
6. Delete tmp files

	rm -rf /tmp/matador3d_update_wwwPDB && rm /tmp/matador3d_update_wwwPDB.log && rm ficherovacio.txt && rm ficherovacio2.txt && rm formatdb.log

>> RefSeq Vertebrate database comes from "vertebrate_mamalian" and "vertebrate_other" (ftp://ftp.ncbi.nlm.nih.gov/refseq/release/)

1. Get RefSeq database for 'vertebrate' and 'invertebrate' from
	wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.*.protein.faa.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_other.*.protein.faa.gz

	wget ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate.*.protein.faa.gz
	
2. Unzip them	
	gzip -d vertebrate_*
	
	gzip -d invertebrate.*
	
3. Concatenate them	
	cat vertebrate_mammalian.* vertebrate_other.* >> refseq_vert
	
	cat invertebrate.* >> refseq_invert

4. Index database
	formatdb -i refseq_vert -p T
	
	formatdb -i refseq_invert -p T

>> Pfam, ftp://ftp.ebi.ac.uk/pub/databases/Pfam/

ftp://ftp.sanger.ac.uk/pub/databases/Pfam/ (DEPRECATED)


1. Get the Pfam database from:
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam*

	In particular you need the files 
	Pfam-A.fasta, Pfam_ls, Pfam_fs, Pfam_ls.bin, Pfam_fs.bin,
	Pfam_ls.ssi Pfam_fs.bin.ssi, and Pfam-A.seed, and optionally
	Pfam-C.
	To use the active site option you will also need to
	download the active site alignments which are available as a
	tarball (active_site.tgz).

2. Unzip them if necessary
    $ gunzip Pfam*.gz

3. Grab and install HMMER, NCBI BLAST and Bioperl, and make sure your
   paths etc are set up properly.

4. Index Pfam-A.fasta for BLAST searches
    $ formatdb -i Pfam-A.fasta -p T

5. Index the Pfam_ls and Pfam_fs libraries for HMM fetching
    $ hmmindex Pfam_ls
    $ hmmindex Pfam_fs
    
>> sprot_clean_trembl_clean_70, fdbTptDB_*, chads_*, hhblits_*, nr20_*, FireDB for firestar

6. Install FireDB database:
	mysql FireDB -h localhost -u firedb < FireDB_*.sql

>> TSL annotations, which will be used in appris methods, is downloaded from BioMart
	
7. TSL.annots.e80_g22_gM5.txt #DEPRACATED#!!!
	$ BioMart > Download TSL annotations for human filtering the gene by: "Only" Transcript Support Level (TSL) and "Only" APPRIS annotations
	$ mv mart_export.txt TSL.annots.e80_g22.txt
	
	$ BioMart > Download TSL annotations for mouse filtering the gene by: "Only" Transcript Support Level (TSL) and "Only" APPRIS annotations
	$ mv mart_export.txt TSL.annots.e80_gM5.txt
	
	concatenate files together except first line
	$ awk 'FNR>1' TSL.annots.e80_g22.txt TSL.annots.e80_gM5.txt > TSL.annots.e80_g22_gM5.txt
	
	
>> PDB for Matador3D2


TODO!!! Check the documentation of Juan Rodriguez in 'scripts/matador3d2':


	
		
