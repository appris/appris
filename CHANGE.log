#####################################
# APPRIS RELEASES:					
#									
# A-RELEASE: Annotation release		
#		{DATE}.v{ANNOT_VERSION} 	
#									
# C-RELEASE: Code release			
# 		{WS}.{C}					
#	WS -> Web server release		
# 	C  -> Code release 				
#									
# I-DATE: Date when input data files have obtained
#									
#####################################


#############################
# A-RELEASE: 	2016_11.v19	#
# C-RELEASE: 	4.4.2.11	#
# I-DATE: 		01Nov2016	#
#############################

* HIGHLIGHTS:
	
* APPRIS-CODE, v2.11

	- CORSAIR:
		- A variant with "start/stop codons not found" can not to win.
		 
* APPRIS-DATA, v17
		- Human
			> Ensembl86/Gencode25 (GRCh38)		ens86.v19
			> Ensembl84/Gencode24 (GRCh38)		ens84.v19
			> Ensembl81/Gencode23 (GRCh38)		ens81.v19 ???
			> Ensembl79/Gencode22 (GRCh38)		ens79.v19 ???
			> RefSeq107							rs107.v19
			> Ensembl74/Gencode19 (GRCh37)		gen19.v19
			> RefSeq105							rs105.v19
			> UniProt 2016_06					up201606.v19
		- Mouse
			> Ensembl86/GencodeM11 (GRCm38)		ens86.v19
			> Ensembl84/GencodeM9  (GRCm38)		ens84.v19
		- Zebra-fish
			> Ensembl86 (GRCz10)				ens86.v19
			> Ensembl77 (Zv9)					ens77.v19
		- Rat
			> Ensembl86 (Rn6)					ens86.v19
			> Ensembl77 (Rn5)					ens77.v19
		- Pig
			> Ensembl86 (Sscrofa10.2)			ens86.v19
		- Chimp
			> Ensembl86 (CHIMP2.1.4)			ens86.v19
		- Fruitfly
			> Ensembl86 (BDGP6)					ens86.v19		
		- C-elegans
			> Ensembl86 (WBcel235)				ens86.v19
			

#############################
# A-RELEASE: 	2016_07.v18	#
# C-RELEASE: 	4.4.2.10	#
# I-DATE: 		01Nov2016	#
#############################

* HIGHLIGHTS:
	
* APPRIS-CODE, v2.10

	- Matador3D v2: Has been reported to UniProt.
	
	

#############################
# A-RELEASE: 	2016_06.v17	#
# C-RELEASE: 	4.4.2.9		#
# I-DATE: 		01Jun2016	#
#############################

* HIGHLIGHTS:
	- Freeze the code. 4.4.2.9


#############################
# A-RELEASE: 	2016_03.v15	#
# C-RELEASE: 	4.4.2.9		#
# I-DATE: 		01Mar2016	#
#############################

* HIGHLIGHTS:
	- New name for each release... Forget the date
	- New algorithm for SPADE: SPADE2
	
* APPRIS-CODE, v2.9

	- APPRIS
		- Fragments (CDS start/stop not found) changes (eg, ADGRD2, NEDD8-MDP1),
			if all variants has fragments, we don't reject them (and keep going with pipeline).
		- We keep discarding the NMDs
		- In the step_length: (Eg, GNAO1),
			the sequences are different but with the same length we choose the winner based on the appris-score.
		- Now appris scripts contains the decision of all methods.
		- New columns for the outputs including the TSL annotation
		- Increase the weight of SPADE
		- For the moment, discard THUMP method
	
	- Matador3D
		- Fixed bug with the number of residues. 

	- SPADE2
		- Decision with the bitscore.
		- First of all, discard overlapping domains coming from pfamscan
		- Herence only domains with bigger bitscore. Eg, ENSG00000164692

		 
* APPRIS-DATA, v15
		- Human
			> Ensembl84/Gencode24 (GRCh38),	ens84.v15
			> Ensembl81/Gencode23 (GRCh38),	ens81.v15
			> Ensembl79/Gencode22 (GRCh38),	ens79.v15
			> Ensembl74/Gencode19 (GRCh37), gen19.v15
		- Mouse
			> Ensembl84/GencodeM9 (GRCm38),	ens84.v15
		- Zebra-fish
			> Ensembl84 (GRCz10), 			ens84.v15
		- Rat
			> Ensembl84 (Rn6),    			ens84.v15
		- Pig
			> Ensembl84 (Sscrofa10.2),		ens84.v15
		- Chimp
			> Ensembl84 (CHIMP2.1.4), 		ens84.v15
		- Fruitfly
			> Ensembl84 (BDGP6), 			ens84.v15		
		- C-elegans
			> Ensembl84 (WBcel235), 		ens84.v15
			
	
#############################
# A-RELEASE: 	2016_01.v14	#
# C-RELEASE: 	4.3.2.8		#
# I-DATE: 		29Jan2016	#
#############################

* HIGHLIGHTS:

* APPRIS-CODE, v2.8
	- APPRIS
		- Discard the filter of "start/stop" codon not found.
		- Add new cutoff for CORSAIR method.
		If the aa length is longer than len-cutoff but the max score is bigger than score-cutoff, then continue; Otherwise, discard the scores of all variants. Eg, FRYL
		
* APPRIS-DATA, v14
		- Human
			> Ensembl84/Gencode24 (GRCh38),	ens84.v14.29Jan2016
			> Ensembl74/Gencode19 (GRCh37), gen19.v14.29Jan2016
		- Mouse
			> Ensembl84/GencodeM9 (GRCm38),	ens84.v14.29Jan2016
		- Zebra-fish
			> Ensembl84 (GRCz10), 			ens84.v14.29Jan2016


		
#############################
# A-RELEASE: 	2015_12.v13	#
# C-RELEASE: 	4.3.2.7		#
# I-DATE: 		14Dec2015	#
#############################

* HIGHLIGHTS:

* APPRIS-SERVER, v4.3
	- Add TSL annots in the report web page.

* APPRIS-CODE, v2.7
	- Matador3D
		- Fixed bug: Take into account the -1 and -0.5 score. 
	
	- APPRIS
		- Fragments (CDS start/stop not found) changes (eg, ADGRD2, NEDD8-MDP1):
			- If all variants has fragments, we don't reject them.
			- If all variants are NM, we don't reject them.
			- If all variants are NM or fragments, we don't reject them (and keep going with pipeline).
		- If fragments have CCDS, we don't reject them (eg, EIF4G2)
					
* APPRIS-DATA, v13
		- Human
			> Ensembl83/Gencode24 (GRCh38),	ens83.v13.14Dec2015
			> Ensembl74/Gencode19 (GRCh37), gen19.v13.14Dec2015
		- Mouse
			> Ensembl83/GencodeM8 (GRCm38),	ens83.v13.14Dec2015
		- Rat
			> Ensembl83 (Rn6),    			ens83.v13.14Dec2015
		- Zebra-fish
			> Ensembl82 (GRCz10), 			ens83.v13.14Dec2015
		- Pig
			> Ensembl82 (Sscrofa10.2),		ens83.v13.14Dec2015
		- Chimp
			> Ensembl82 (CHIMP2.1.4), 		ens83.v13.14Dec2015
		- Fruitfly
			> Ensembl82 (BDGP6), 			ens83.v13.14Dec2015		
		- C-elegans
			> Ensembl82 (WBcel235), 		ens83.v13.14Dec2015


#############################
# A-RELEASE: 	2015_11.v12	#
# C-RELEASE: 	4.3.2.6		#
# I-DATE: 		26Oct2015	#
#############################

* HIGHLIGHTS:
	- Hand-over for e83: Human, Mouse and Rat.
		
* APPRIS-DATA, v12
		- Human
			> Ensembl83/Gencode24 (GRCh38),	ens83.v12.26Oct2015
		- Mouse
			> Ensembl83/GencodeM8 (GRCm38),	ens83.v12.26Oct2015
		- Rat
			> Ensembl83 (Rn6),    			ens83.v12.26Oct2015
		- Zebra-fish
			> Ensembl82 (GRCz10), 			ens82.v11.01Sep2015
		- Pig
			> Ensembl82 (Sscrofa10.2),		ens82.v11.01Sep2015
		- Chimp
			> Ensembl82 (CHIMP2.1.4), 		ens82.v11.01Sep2015
		- Fruitfly
			> Ensembl82 (BDGP6), 			ens82.v11.01Sep2015		
		- C-elegans
			> Ensembl82 (WBcel235), 		ens82.v11.01Sep2015
	
	- New gene data set: RefSeq!!! (version r107)
	
#############################
# A-RELEASE: 	2015_10.v11	#
# C-RELEASE: 	4.3.2.6		#
# I-DATE: 		01Sep2015	#
#############################

* HIGHLIGHTS:
	
* APPRIS-CODE, v2.6
	- APPRIS
		- The filter of 'codon not found' has been deleted. For cases like GTF3A, or CG43172 gene for Drosophila. 
	
	- Scripts than download the features files for species and Ensembl version.
	
* APPRIS-DATA, v11
	- Annotations for all species using the code version v2.6.
		- Human
			> Ensembl82/Gencode23 (GRCh38),	ens82.v11.01Sep2015
			> Ensembl74/Gencode19 (GRCh37), gen19.v9.17Jul2015
		- Mouse
			> Ensembl82/GencodeM7 (GRCm38),	ens82.v11.01Sep2015
		- Zebra-fish
			> Ensembl82 (GRCz10), 			ens82.v11.01Sep2015
		- Rat
			> Ensembl82 (Rn6),    			ens82.v11.01Sep2015
		- Pig
			> Ensembl82 (Sscrofa10.2),		ens82.v11.01Sep2015
		- Chimp
			> Ensembl82 (CHIMP2.1.4), 		ens82.v11.01Sep2015
		- Fruitfly
			> Ensembl82 (BDGP6), 			ens82.v11.01Sep2015		
		- C-elegans
			> Ensembl82 (WBcel235), 		ens82.v11.01Sep2015
	
#############################
# A-RELEASE: 	2015_09.v10	#
# C-RELEASE: 	4.3.2.5		#
# I-DATE: 		21Aug2015	#
#############################

* HIGHLIGHTS:
	- Hand-over for e82: Mouse.
	- NEW!!! Creation of UCSC 'trackHub' for APPRIS species:
		human(hg38,hg19), mouse(mm10), zebra-fish(danRer10), rat(rn6), pig(susScr3), chimp(panTro4), fruitfly(dm6), C.elegans(ce10)   
	
* APPRIS-DATA, v10
	New gene data-set for Mouse:
	- Mouse
		> Ensembl82/GencodeM7 (GRCm38),	ens82.v10.21Aug2015
	The same annotations for the rest of species 	
	
	NEW SPECIES:
	- Fruitfly
		> Ensembl81 (BDGP6), 			ens81.v10.21Aug2015
	- C-elegans
		> Ensembl81 (WBcel235), 		ens81.v10.21Aug2015
	
#############################
# A-RELEASE: 	2015_08.v10 #
# C-RELEASE: 	4.3.2.5		#
# I-DATE: 		21Aug2015	#
#############################

* APPRIS-CODE, v2.5
	- APPRIS 
		- Horrible bug fixed. The decision in length was bad :-(
		- Bug fixed in the APPRIS decision:
        	Eg. for CEP152 gene, we decide as ALTernative transcript that the biotype is not protein_coding and codons not found
        	But we check also if the gene has unique transcripts, eg. CTD-3222D19.2
	- Change the name of databases. Add suffix for version

* APPRIS-DATA, v10
	- Human
		> Ensembl81/Gencode23 (GRCh38),	ens81.v10.21Aug2015
		> Ensembl76/Gencode20 (GRCh38),	ens76.v10.21Aug2015
		> Ensembl74/Gencode19 (GRCh37), gen19.v10.21Aug2015
	- Mouse
		> Ensembl81/GencodeM6 (GRCm38),	ens81.v10.21Aug2015
		> Ensembl74/GencodeM2 (GRCm38),	ens74.v10.21Aug2015
	- Zebra-fish
		> Ensembl81 (GRCz10), 			ens81.v10.21Aug2015
	- Rat
		> Ensembl81 (Rn6),    			ens81.v10.21Aug2015
	- Pig
		> Ensembl81 (Sscrofa10.2),		ens81.v10.21Aug2015
	- Chimp
		> Ensembl81 (CHIMP2.1.4), 		ens81.v10.21Aug2015
		
#############################
# A-RELEASE: 	2015_08.v9  #
# C-RELEASE: 	4.3.2.4 	#
# I-DATE: 		17Jul2015	#
#############################

* HIGHLIGHTS:
	- Hand-over for e81: Human, and Mouse.
	- APPRIS has included invertebrate species: Drosophila melanogaster and Caenorhabditis elegans!!!	
	- TSL method (http://www.ensembl.org/Help/Glossary?id=492) is included into APPRIS decision.
	
* APPRIS-SERVER, v4.3
	- Now Sequence Browser retrieves annotations in detail for residue.
	- Improve/fixing the annotations within sequence Browser
	- Change of 'viewer' web service that prints the UCSC image to 'sequencer'
	- Include the highlight of CDS regions within Sequence Browser
	- Row panel with the following frames: GitHub, Google Group news-email, and Twitter.
	- Add the new species into frontpage.
	- Include version label.
	- Include 'Changelogs' section.
	- Maintain templates 'app_mnt'
	
* APPRIS-CODE, v2.4
	- MATADOR3D
		- The score of total gaps was not working. There were not cases with -1,-0.5 score. E.g. PAX6 gene
	- SPADE
  		- It does not take into account "insertions" within HMM sequences. Fixed. E.g. PAX6 gene
		- Fixing the range of residues when a transcripts receives exported alignments
  	- APPRIS
		- TSL annotation (tsl1) is included into APPRIS decision		
	- APPRIS pipeline has included invertebrate species: Drosophila melanogaster and Caenorhabditis elegans.
	- CORSAIR
		- Include new RefSeq database for "Invertebrates".
	- Discard the Ensembl version only in the cases where the identifier starts with 'ENS'.
	- In the download directory, we have classified data files into directories based on format type.
	- APPRIS data file for method have changed their prefix name from 'appris_data.{method_name}' to 'appris_method.{method_name}'
	- Print the alignments of Pfam separately in BED format
	- 'retrieve_method_data', Methods as input parameter
	- Using fork to spread load to multiple cores	

* APPRIS-DATA, v9
	- Human
		> Ensembl81/Gencode23 (GRCh38),	ens81.v9.17Jul2015
		> Ensembl76/Gencode20 (GRCh38),	ens76.v9.17Jul2015 <- FREEZE!!
		> Ensembl74/Gencode19 (GRCh37), gen19.v9.17Jul2015 <- FREEZE!!
	- Mouse
		> Ensembl81/GencodeM6 (GRCm38),	ens81.v9.17Jul2015
		> Ensembl74/GencodeM2 (GRCm38),	ens74.v9.17Jul2015 <- FREEZE!!
	- Zebra-fish
		> Ensembl81 (GRCz10), 			ens81.v9.17Jul2015
	- Rat
		> Ensembl81 (Rn6),    			ens81.v9.17Jul2015
	- Pig
		> Ensembl81 (Sscrofa10.2),		ens81.v9.17Jul2015

	NEW SPECIES:
	- Chimp
		> Ensembl81 (CHIMP2.1.4), 		ens81.v9.17Jul2015
	- Fruitfly
		> Ensembl81 (BDGP6), 			ens81.v9.17Jul2015
		
#############################
# A-RELEASE: 	2015_05.v8  #
# C-RELEASE: 	4.2.2.3 	#
# I-DATE: 		16Apr2015	#
#############################

* HIGHLIGHTS:
	- Hand-over for e80: Mouse, Zebra-fish and Rat.
	- New assemblies for Zebra-fish (GRCz10) and Rat (Rnor_6.0).

* APPRIS-SERVER, v4.2
	- Fix a bug when UCSC image is shown.	
	- Include AWStats is a free powerful and featureful tool that generates advanced web, streaming, ftp or mail server statistics, graphically.
	http://apprisws.bioinfo.cnio.es/awstats/awstats.pl?month=all&year=2015&output=main&config=appris&framename=index
		
* APPRIS-CODE, v2.3
	- APPRIS
		- Fixing a bug where we detect "start/stop" codons
		
* APPRIS-DATA, v8
	- Human
		> Ensembl80/Gencode22 (GRCh38),	ens80.v8.16Apr2015 = ens79.v8.2Apr2015 
		> Ensembl76/Gencode20 (GRCh38),	ens76.v8.16Apr2015 <- FREEZE!!
	- Mouse
		> Ensembl80/GencodeM5 (GRCm38),	ens80.v8.16Apr2015
		> Ensembl74/GencodeM2 (GRCm38),	ens74.v8.16Apr2015 <- FREEZE!!
	- Zebra-fish
		> Ensembl80 (GRCz10),			ens80.v8.16Apr2015
	- Rat
		> Ensembl80 (Rn6),				ens80.v8.16Apr2015
	- Pig
		> Ensembl80 (Sscrofa10.2),		ens80.v8.16Apr2015
		
#############################
# A-RELEASE: 	2015_03.v7  #
# C-RELEASE: 	4.1.2.2 	#
# I-DATE: 		9Feb2015	#
#############################

* HIGHLIGHTS:
	- Fixing a bug in the CORSAIR filter.
	- We have increase the number of the version for the release. The reason is the same number.
	- The annotations for GRCh37 of human have been frozen: gen19.v7.9Feb2015

* APPRIS-SERVER, v4.1

* APPRIS-CODE, v2.2
	- APPRIS 
		- Fixing a bug in the CORSAIR filter
		- Fixing a bug in the filter of length sequences.
		
* APPRIS-DATA, v7
	- Human
		> Ensembl77/Gencode21 (GRCh38), ens77.v7.9Feb2015
		> Ensembl76/Gencode20 (GRCh38), ens76.v7.9Feb2015
		> Ensembl74/Gencode19 (GRCh37), gen19.v7.9Feb2015 <- FREEZE!!
	- Mouse
		> Ensembl78/GencodeM4 (GRCm38), ens78.v7.9Feb2015
	- Zebra-fish
		> Ensembl77 (Zv9),    ens77.v7.9Feb2015 <- FREEZE!!
	- Rat
		> Ensembl77 (Rn5),    ens77.v7.9Feb2015 <- FREEZE!!
	- Pig
		> Ensembl77 (Sscrofa10.2), ens77.v7.9Feb2015
		
#############################
# A-RELEASE: 	2014_02.v1  #
# C-RELEASE: 	4.0.2.1 	#
# I-DATE: 		27Jan2014	#
#############################

* HIGHLIGHTS:
	- NEW APPRIS SERVER
	- New APPRIS algorithm using oldest CCDS ids as principal isoforms:

* APPRIS-SERVER, v4.0

* APPRIS-CODE, v2.1
	- New tags of PRINCIPAL ISOFORMS
	- New script to retrieve the list of principal isforms for Ensembl
	- APPRIS
		- New APPRIS algorithm using oldest CCDS ids as principal isoforms:
			- With the 'appris_candiates' we check:
				1. unique CCDS
				2. the oldest CCDS (smaller than 10numbers comparing with the next CCDS number)
				3. the longest CCDS
				4. the longes protein sequence
				
* APPRIS-DATA, v1
	- Human
		> Ensembl79 (GRCh38), ens79.v1.26Jan2015 # TEMPORAL
		> Ensembl77/Gencode21 (GRCh38), ens77.v5.27Jan2015
		> Ensembl74/Gencode19 (GRCh37), gen19.v6.27Jan2015
	- Mouse
		> Ensembl78/GencodeM4 (GRCm38), ens78.v3.27Jan2015
	- Zebra-fish
		> Ensembl77 (Zv9),    ens77.v3.27Jan2015
	- Rat
		> Ensembl77 (Rn5),    ens77.v3.27Jan2015
	- Pig
		> Ensembl77 (Sscrofa10.2), ens77.v2.27Jan2015
		
###############################################
# DESCRIPTION OF APPRIS RELEASES: WS.C1.D	  #
# WS -> Release of web server				  #
# C1 -> Release of code (the second number)	  #
# D  -> Release of data files			 	  #
###############################################
		 
#######################
# DATE: 	1Nov2014  #
# RELEASE:  3.4.4     #
#######################

* HIGHLIGHTS:
	- APPRIS improvements included in for Dominant Isoform paper.
	- Frozen version (tag release-1.4 in svn repository).

* APPRIS-DATA
	- Human
		> Ensembl77/Gencode21 (GRCh38), ens77.v4.31Oct2014
		> Gencode19/Ensembl74 (GRCh37), gen19.v5.31Oct2014
	- Mouse
		> Ensembl78 (GRCm38), ens78.v2.31Oct2014
		> GencodeM3/Ensembl77 (GRCm38), genM3.v4.31Oct2014
	- Zebra-fish
		> Ensembl77 (Zv9),    ens77.v2.31Oct2014
	- Rat
		> Ensembl77 (Rn5),    ens77.v2.31Oct2014
	- Pig
		> Ensembl77 (Sscrofa10.2), ens77.v1.13Jan2015
		
* APPRIS-CODE, v1.4 (-> v.2.0)
	- APPRIS
		- New APPRIS algorithm:
			- new weights for methods.
			
#######################
# DATE: 	31Oct2014 #
# RELEASE:  3.4.3     #
#######################

* HIGHLIGHTS:
	- New filter in APPRIS in the case of CORSAIR decision.
	- Release for Ensembl.

* APPRIS-DATA
	- Human
		> Ensembl77/Gencode21 (GRCh38), ens77.v3.28Oct2014
	- Mouse
		> Ensembl78 (GRCm38), ens78.v1.28Oct2014
	- Zebra-fish
		> Ensembl77 (Zv9),    ens77.v1.28Oct2014
	- Rat
		> Ensembl77 (Rn5),    ens77.v1.30Oct2014
		
		
* APPRIS-CODE, v1.4
	- APPRIS
		- New APPRIS algorithm:
			- filter by aa length in CORSAIR
			- new weights for CORSAIR.
			- new normalization for methods. 
	
#######################
# DATE: 	09Oct2014 #
# RELEASE:  3.4.3     #
#######################

* HIGHLIGHTS:
	- Include CCDS correctly for Ensembl77/Gencode21.

* APPRIS-DATA
	- Human
		> Gencode21/Ensembl77 (GRCh38), ens77.v2.08Oct2014
		
#######################
# DATE: 	08Oct2014 #
# RELEASE:  3.4.2     #
#######################

* HIGHLIGHTS:
	- Include CCDS correctly for Gencode20 and GencodeM3.

* APPRIS-DATA
	- Human
		> Gencode21/Ensembl77 (GRCh38), ens77.v1.29Aug2014
		> Gencode20 (GRCh38), gen20.v2.13Aug2014
	- Mouse
		> GencodeM3 (GRCm38), genM3.v2.13Aug2014
	- Zebra-fish
		> Ensembl76 (Zv9),    ens76.v1.03Oct2014
		
* APPRIS-CODE, v1.4
	- GENCODE has changed the comment line of FASTA translation file.
	- Include method which parser the input files coming from Ensembl (gtf, pep.fa, cdna.fa).
	- Cached files for gene. 
						
#######################
# DATE: 	16Jul2014 #
# RELEASE:  3.4.1     #
#######################

* APPRIS-DATA
	- Human
		> Gencode20 (GRCh38), gen20.v1.16Jul2014
		> Gencode19, gen19.v4.16Jul2014
		> Gencode18, gen18.v1.16Jul2014
	- Mouse
		> GencodeM3 (GRCm38), genM3.v1.16Jul2014
		
* APPRIS-CODE, v1.4
	- Print parameters within log file
	- Gencode methods do not need connect to Ensembl (Parser.pm)
	- Include RESTful web services of APPRIS Runner.
	- New Export format (TXT.pm)	
	- apprisall.pl executes in a cluster when is located with itself host.
	- apprisall.pl and appris.pl accepts new parameter of alignments. Now, 'compara' alignments needs ensembl version.
	- Skeleton of TEST battery
	
	- Update scripts of stats. CExonic output of appris rst has been deleted.
	
	- SPADE (NOT USED YET: IN TEST)
		- Normalized lengths for each DOMAIN, plus template gap penalties.
		
	- APPRIS
		- Score output prints for the following methods:
			-- spade: "(Num. Whole Domains+Num. Possible Domains) - (Num. Wrong Domains+Num. Damaged Domains)"
			-- thump: "Num. TMH - Num. Damaged TMH" 

	- PROTEO
		- New script that saves peptides by gene
		
	- APPRIS
		- New algorithm to select the main variant. We normalize the scores of methods.
		- Delete CExonic output.
		- We take into account the CLASS of transcript for decisition but we don't use "codon start/end not found" label yet.
		
	- CEXONIC
		- Has been deleted
				
	- Export PROTEO features by means RESTful services.
	
* APPRIS-WEBSITE, v1.3
	- CExonic is DEPRECATED
	- Add PROTEO annotations
	- Improve Web Site features.
	- Include PROTEO RESTful service into web site.
	
	- New Archives: Gencode18, and Gencode19... (Finish as soon as possible the new web site because it is a nightmare!!!)	
					
#######################
# SKIP      --------- #
# DATE: 	--------- #
# RELEASE:  3.3.1     # 
#######################

* APPRIS-DATA
	- Human
		> Gencode19, gen19.v3
	- Mouse
		> GencodeM2, genM2.v2
		
* APPRIS-CODE, v1.3
	- Print parameters within log file
	- Gencode methods do not need connect to Ensembl (Parser.pm)
	- Include RESTful web services of APPRIS Runner.
	- New Export format (TXT.pm)
	
	- INERTIA
		- INERTIA is executed only for the transcripts with CDS (protein-coding), and CDS start-end have been founded.

	- CORSAIR
		- Exon data is not required.
		- We have added the conservation score by exons by exon.

	- MATADOR3D
		- New threshold of identity/gaps when we parse the alignment of "Blast".
		
	- FIRESTAR
		- New version of firestar.
		- New FireDB.
		- Now, the score of SQUARE is a median of site.
		- appris-firestar takes into account all residues coming from firestar's core.
		- appris-firestar: I have reduced the difference between the number of residues when a variant wins. From 6 to 2.
		
	- THUMP
		- Fix the count of helixes.
		
	- PROTEO
		- Add new method into APPRIS database (but REMEMBER it is NOT within appris decision)
		
#######################
# DATE: 	20Dec2013 #
# RELEASE:  3.2.3     #
#######################

* MAIN DESCRIPTION
	- New release of annotations of Gencode19. Preliminar gencode data did not contain external name of genes.

* APPRIS-DATA
	- Human
		> Gencode19, gen19.v2.18Dec2013

#######################
# DATE: 	 2Dec2013 #
# RELEASE:  3.2.2     #
#######################

* MAIN DESCRIPTION

* APPRIS-DATA
	- Change the format of list of principal isoforms.
	For more information:
		http://appris.bioinfo.cnio.es/download/README.txt 
	- Human
		> Gencode19, gen19.v1.25Nov2013
	- Mouse
		> GencodeM2, genM2.v1.2Dec2013
	- Rat
		> Ensembl70, e70.v4.29Oct2013
	- Zebra-fish
		> Ensembl70, e70.v1.25Nov2013		


#######################
# DATE: 	25Nov2013 #
# RELEASE:  3.2.1     #
#######################

* MAIN DESCRIPTION

* APPRIS-DATA
	- Include "chromosome" column within "conserved_exons" data.
	- Human
		> Gencode15, g15.v4.16Oct2013		
	- Mouse
		> Ensembl70, e70.v4.29Oct2013
	- Rat
		> Ensembl70, e70.v4.29Oct2013
	- Lynx (private data)
		> p23A.v3.14Oct2013 

* APPRIS-CODE, v1.2
	- Put together README files
	- Include all code of APPRIS in GitHub (https://github.com/inab/appris)
	
	- (main) APPRIS
		- Modify the way to obtain the inputs
			- Type of inputs: 
				- directory data (cached)
				- gencode data
				- ensembl data
				- sequence (protein) data
			- Type of alignment:
				- compara
		- Delete runner files				

		- Add INERTIA script
		- Accepts new parameters. 't-align' that describes the type of alignment.
			- Type of alignment:
				- compara		

	- CORSAIR
		- New scores for the Lynx pardinus
		- New cutoff of 1.5
		- Include the list of species for 'Danio rerio'
	
	- MATADOR3D
		- The gene ENSG00000160404 is rejecting ENST00000373284 with an negative alignment -0.20125 because 
		the percentage of gaps is bigger than 25.
		We have modify the threshold of percentages of gaps to:
		
				elsif ($gaps > 25)
					{ $totalgaps = 0 }
				elsif ($gaps > 33)
					{ $totalgaps = -0.5 }
				elsif ($gaps > 40)
				
		- Take the correct transcript		
		- Change the threshold of Matador3D to 0.65

	- APPRIS
		- Reliability scores:
			If transcript is UNKOWN, its reliability score is NR (not rejected).
			And if the transcript has the longest protein sequences, its reliability score is NR*
			
			Only the transcripts with PRICIPAL have reliability scores.
		- If we don't have information of codons, we continue saving the appris score.

	- THUMP		
		- Changes the regular expresion to parse transcript for any kind of identifier. 
			

#########################
# DATE: 	16Aug2013   #
# RELEASE:  3.1.1       #
#########################

* APPRIS-WEBSITE, v1.3

- Links of archives redirect to more specific website.

- Rename the suffix of data files of 'rattus_norvegicus'.

- Change the annotations of SPADE and THUMP. Now, we print their "scores" in the website:
	- SPADE shows an float number where the absolute number is the pfam domains that are detected. Then, the number after the decimal indicates the numbers of partial domains.
	- THUMP shows the number of TMH detected. Again, the number after the decimal indicates the partial TMH.
 
- Delete the phrase "on Human (GRCh37/hg19) Assembly" from the 'report' web page.

- Fix the URL that retrieves proteomics data as BED format.

- Update the proteomic data (gencode12):
	- Peptides were assembled from five previously available proteomics data sets. Three of the peptide data sets came from published large-scale experiments (Ezkurdia, Nagaraj, Geiger) 
	and the others were from two large spectra libraries, PeptideAtlas (Farrah) and NIST (http://peptide.nist.gov/).
	These tracks show where peptide evidence exists (0 or 1) , but there is no means of scoring these tracks.  


#########################
# DATE: 	25Jul2013   #
# RELEASE:  3.1         #
#########################


* MAIN DESCRIPTION:


* APPRIS-CODE, v1.1

We have decided to establish new nomenclature of version of appris's code.

	- SPADE does not migrate domains of variants with the same sequence where pfamscan does not provide results. 
		For example, ENSG00000021300.7 (Gencode15)
		
	- FIRESTAR prints results whose cutoff is smaller than 65%. This is bad!!! 
	For example, ENSG00000184361 (Gencode15)

	- CORSAIR, has new scores of the rest of species (including Lynx pardinus):
		'Rest of species'
			'Homo sapiens'				=> 1,
			'Pan troglodytes'			=> 1,
			'Mus musculus'				=> 1.1,
			'Rattus norvegicus'			=> 1.1,
			'Bos taurus'				=> 1.2,
			'Canis lupus familiaris'	=> 1.2,
			'Sus scrofa'				=> 1.2,
			'Monodelphis domestica'		=> 1.3,
			'Gallus gallus'				=> 2,
			'Taeniopygia guttata'		=> 2,
			'Anolis carolinensis'		=> 2,
			'Xenopus tropicalis'		=> 2.2,
			'Tetraodon nigroviridis'	=> 2,
			'Danio rerio'				=> 2.5,

* APPRIS-DATA:

- Human
	> Gencode12, rel12.27Jun2013.v5
		- SPADE did not migrate variant domains with the same sequence.
		- INERTIA new results.

	> Gencode15, g15.v3.15Jul2013
		- SPADE did not migrate variant domains with the same sequence.		
		- FIRESTAR printed results whose cutoff was smaller than 65%.

* NEW GENOMES:

- Mouse
	> Ensembl70, e70.v3.15Jul2013
		- SPADE did not migrate variant domains with the same sequence.		
		- FIRESTAR printed results whose cutoff was smaller than 65%.

- Rat
	> Ensembl70, e70.v3.10Jul2013
		- SPADE did not migrate variant domains with the same sequence.		
		- FIRESTAR printed results whose cutoff was smaller than 65%.


* APPRIS-WEBSITE, v1.3
	
	- Add script for the rest of species
	
	- Include new data files into website:
		- appris data exon by exon
		
	- Modify "download" web page.




