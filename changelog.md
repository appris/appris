___
## 2017_01.v20
```
SERVER-RELEASE:   2017_01.v20
CODE-RELEASE:     4.6.2.13
```

### Highlights
- New release of Ensembl (e87) for all species.
- APPRIS annotations for UniProt (2016_10 version) in Mouse, Zebra-fish, and Rat proteomes.

### Code (2.13)
- APPRIS:
	Changes in the weight of CORSAIR.
	Changes in the weight of firestar (e.g. bad result for ANAPC15 - ENSG00000110200).
- Matador3D v2: Has been included in APPRIS code.


___
## 2016_11.v19
```
SERVER-RELEASE:   2016_11.v19
CODE-RELEASE:     4.6.2.12
```

### Highlights
- New release of Ensembl (e86 version).
- APPRIS annotations for UniProt (2016_06 version).

### Server (4.6)
- Include Advanced Search panel.
- CNIO proteomic evidences (PROTEO) have been included into the SequenceBrowser

### Code (2.12)
- CORSAIR:
	A variant with "start/stop codons not found" can not to win.


___
## 2016_07.v18
```
SERVER-RELEASE:   2016_07.v18
CODE-RELEASE:     4.5.2.11b
```

### Highlights
- New release of Ensembl (Ensembl 85) for human, mouse, and rat.
	
### Code (2.11b)
- Matador3D v2: Has been reported to UniProt.



___
## 2016_07.v17
```
SERVER-RELEASE:   2016_07.v17
CODE-RELEASE:     4.5.2.11
```

### Highlights
- New release of Ensembl (Ensembl 85) for human, mouse, and rat.
- Freeze the code. [4.5.2.11] (https://github.com/appris/appris/tree/freeze_4.5.2.11)



___
## 2016_06.v17
```
SERVER-RELEASE:   2016_06.v17
CODE-RELEASE:     4.5.2.11
```

### Highlights
- RefSeq versions 105 and 107 for the assemblies GRCh37 and GRCh38, respectively. 
	
### Server (4.5)
- Change in the web services to accepts multiple datasets sources.

- Now the search web services retrieves the results for a specific query.

- Add TSL annots in the report web page:

### Code (2.11)
- The file names of method results have changed.

- Add TSL and TAG features in "entity" table.

- Change the "namespace" of "datasource" table.

- Delete HAVANA datasource

- APPRIS
	- Fixed bug. No decision when the sequence length is equal, with the same appris score. The list of genes was:
		ENSG00000070831 ENSG00000103248 ENSG00000132639 ENSG00000163632 ENSG00000265590 ENSG00000267127 ENSG00000268434
		
- Create "server.json" config file for server.

- Create script "appristools_srv" to update the APPRIS release


___
## 2016_04.v16
```
SERVER-RELEASE:   2016_04.v16
CODE-RELEASE:     4.4.2.10
```

### Highlights
- RefSeq versions 105 and 107 for the assemblies GRCh37 and GRCh38, respectively.
	
### Code (2.10)
- Add "tag" field into database and return its values like "readthought_transcripts"

- Script that add extra values into dataset file

- APPRIS
	- Use the tag features to identify the "readthought_transcripts"
	
	- Fixed bug: It is wrong the assignation of Principal in the step of eldest CCDS and TSL. Eg. BAIAP2 TODO????? It's working correctly
	
	- Check tsl order in ensemble and gencode. Eg. BBS9
	
	- Now TSL annotations is obtained from gene annotation file. The function that extracts the TSL is deprecated.
		

___
## 2016_03.v15
```
SERVER-RELEASE:   2016_03.v15
CODE-RELEASE:     4.4.2.9
```

### Highlights
- New name for each release... Forget the date
- New algorithm for SPADE: SPADE2
	
### Code (2.9)
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

	

___
## 2016_01.v14
```
A-RELEASE: 	2016_01.v14
C-RELEASE: 	4.3.2.8
I-DATE: 		29Jan2016
```

### Highlights

### Code (2.8)
- APPRIS
	- Discard the filter of "start/stop" codon not found.
	- Add new cutoff for CORSAIR method.
	If the aa length is longer than len-cutoff but the max score is bigger than score-cutoff, then continue; Otherwise, discard the scores of all variants. Eg, FRYL
		
### Annotations (v14)
- Human
	* Ensembl84/Gencode24 (GRCh38),	ens84.v14.29Jan2016
	* Ensembl74/Gencode19 (GRCh37), gen19.v14.29Jan2016
- Mouse
	* Ensembl84/GencodeM9 (GRCm38),	ens84.v14.29Jan2016
- Zebra-fish
	* Ensembl84 (GRCz10), 			ens84.v14.29Jan2016


___
## 2015_12.v13
```
A-RELEASE: 	2015_12.v13
C-RELEASE: 	4.3.2.7
```

### Highlights

### APPRIS-SERVER, v4.3
	- Add TSL annots in the report web page.

### Code (2.7)
- Matador3D
	- Fixed bug: Take into account the -1 and -0.5 score. 

- APPRIS
	- Fragments (CDS start/stop not found) changes (eg, ADGRD2, NEDD8-MDP1):
		- If all variants has fragments, we don't reject them.
		- If all variants are NM, we don't reject them.
		- If all variants are NM or fragments, we don't reject them (and keep going with pipeline).
	- If fragments have CCDS, we don't reject them (eg, EIF4G2)
					
### Annotations (v13)
- Human
	* Ensembl83/Gencode24 (GRCh38),	ens83.v13.14Dec2015
	* Ensembl74/Gencode19 (GRCh37), gen19.v13.14Dec2015
- Mouse
	* Ensembl83/GencodeM8 (GRCm38),	ens83.v13.14Dec2015
- Rat
	* Ensembl83 (Rn6),    			ens83.v13.14Dec2015
- Zebra-fish
	* Ensembl82 (GRCz10), 			ens83.v13.14Dec2015
- Pig
	* Ensembl82 (Sscrofa10.2),		ens83.v13.14Dec2015
- Chimp
	* Ensembl82 (CHIMP2.1.4), 		ens83.v13.14Dec2015
- Fruitfly
	* Ensembl82 (BDGP6), 			ens83.v13.14Dec2015		
- C-elegans
	* Ensembl82 (WBcel235), 		ens83.v13.14Dec2015

___
## 2015_11.v12
```
A-RELEASE: 	2015_11.v12
C-RELEASE: 	4.3.2.6
```

### Highlights
- Hand-over for e83: Human, Mouse and Rat.
		
### Annotations (v12)
- Human
	* Ensembl83/Gencode24 (GRCh38),	ens83.v12.26Oct2015
- Mouse
	* Ensembl83/GencodeM8 (GRCm38),	ens83.v12.26Oct2015
- Rat
	* Ensembl83 (Rn6),    			ens83.v12.26Oct2015
- Zebra-fish
	* Ensembl82 (GRCz10), 			ens82.v11.01Sep2015
- Pig
	* Ensembl82 (Sscrofa10.2),		ens82.v11.01Sep2015
- Chimp
	* Ensembl82 (CHIMP2.1.4), 		ens82.v11.01Sep2015
- Fruitfly
	* Ensembl82 (BDGP6), 			ens82.v11.01Sep2015		
- C-elegans
	* Ensembl82 (WBcel235), 		ens82.v11.01Sep2015
	
- New gene data set: RefSeq!!! (version r107)

___
## 2015_10.v11
```
A-RELEASE: 	2015_10.v11
C-RELEASE: 	4.3.2.6
I-DATE: 	01Sep2015
```

### Highlights
	
### Code (2.6)
- APPRIS
	- The filter of 'codon not found' has been deleted. For cases like GTF3A, or CG43172 gene for Drosophila. 

- Scripts than download the features files for species and Ensembl version.
	
### Annotations (v11)
- Annotations for all species using the code version v2.6.
	- Human
		* Ensembl82/Gencode23 (GRCh38),	ens82.v11.01Sep2015
		* Ensembl74/Gencode19 (GRCh37), gen19.v9.17Jul2015
	- Mouse
		* Ensembl82/GencodeM7 (GRCm38),	ens82.v11.01Sep2015
	- Zebra-fish
		* Ensembl82 (GRCz10), 			ens82.v11.01Sep2015
	- Rat
		* Ensembl82 (Rn6),    			ens82.v11.01Sep2015
	- Pig
		* Ensembl82 (Sscrofa10.2),		ens82.v11.01Sep2015
	- Chimp
		* Ensembl82 (CHIMP2.1.4), 		ens82.v11.01Sep2015
	- Fruitfly
		* Ensembl82 (BDGP6), 			ens82.v11.01Sep2015		
	- C-elegans
		* Ensembl82 (WBcel235), 		ens82.v11.01Sep2015

___
## 2015_09.v10
```
A-RELEASE: 	2015_09.v10
C-RELEASE: 	4.3.2.5
I-DATE: 	21Aug2015
```

### Highlights
- Hand-over for e82: Mouse.
- NEW!!! Creation of UCSC 'trackHub' for APPRIS species:
	human(hg38,hg19), mouse(mm10), zebra-fish(danRer10), rat(rn6), pig(susScr3), chimp(panTro4), fruitfly(dm6), C.elegans(ce10)   
	
### Annotations (v10)
New gene data-set for Mouse:
- Mouse
	* Ensembl82/GencodeM7 (GRCm38),	ens82.v10.21Aug2015
The same annotations for the rest of species 	

NEW SPECIES:
- Fruitfly
	* Ensembl81 (BDGP6), 			ens81.v10.21Aug2015
- C-elegans
	* Ensembl81 (WBcel235), 		ens81.v10.21Aug2015
	
___
## 2015_08.v10
```
A-RELEASE: 	2015_08.v10
C-RELEASE: 	4.3.2.5
I-DATE: 	21Aug2015
```

### Code (2.5)
- APPRIS 
	- Horrible bug fixed. The decision in length was bad :-(
	- Bug fixed in the APPRIS decision:
    	Eg. for CEP152 gene, we decide as ALTernative transcript that the biotype is not protein_coding and codons not found
    	But we check also if the gene has unique transcripts, eg. CTD-3222D19.2
- Change the name of databases. Add suffix for version

### Annotations (v10)
- Human
	* Ensembl81/Gencode23 (GRCh38),	ens81.v10.21Aug2015
	* Ensembl76/Gencode20 (GRCh38),	ens76.v10.21Aug2015
	* Ensembl74/Gencode19 (GRCh37), gen19.v10.21Aug2015
- Mouse
	* Ensembl81/GencodeM6 (GRCm38),	ens81.v10.21Aug2015
	* Ensembl74/GencodeM2 (GRCm38),	ens74.v10.21Aug2015
- Zebra-fish
	* Ensembl81 (GRCz10), 			ens81.v10.21Aug2015
- Rat
	* Ensembl81 (Rn6),    			ens81.v10.21Aug2015
- Pig
	* Ensembl81 (Sscrofa10.2),		ens81.v10.21Aug2015
- Chimp
	* Ensembl81 (CHIMP2.1.4), 		ens81.v10.21Aug2015
		
___
## 2015_08.v9
```
A-RELEASE: 	2015_08.v9
C-RELEASE: 	4.3.2.4
I-DATE: 		17Jul2015
```

### Highlights
- Hand-over for e81: Human, and Mouse.
- APPRIS has included invertebrate species: Drosophila melanogaster and Caenorhabditis elegans!!!	
- TSL method (http://www.ensembl.org/Help/Glossary?id=492) is included into APPRIS decision.
	
### APPRIS-SERVER, v4.3
- Now Sequence Browser retrieves annotations in detail for residue.
- Improve/fixing the annotations within sequence Browser
- Change of 'viewer' web service that prints the UCSC image to 'sequencer'
- Include the highlight of CDS regions within Sequence Browser
- Row panel with the following frames: GitHub, Google Group news-email, and Twitter.
- Add the new species into frontpage.
- Include version label.
- Include 'Changelogs' section.
- Maintain templates 'app_mnt'
	
### Code (2.4)
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

### Annotations (v9)
- Human
	* Ensembl81/Gencode23 (GRCh38),	ens81.v9.17Jul2015
	* Ensembl76/Gencode20 (GRCh38),	ens76.v9.17Jul2015 <- FREEZE!!
	* Ensembl74/Gencode19 (GRCh37), gen19.v9.17Jul2015 <- FREEZE!!
- Mouse
	* Ensembl81/GencodeM6 (GRCm38),	ens81.v9.17Jul2015
	* Ensembl74/GencodeM2 (GRCm38),	ens74.v9.17Jul2015 <- FREEZE!!
- Zebra-fish
	* Ensembl81 (GRCz10), 			ens81.v9.17Jul2015
- Rat
	* Ensembl81 (Rn6),    			ens81.v9.17Jul2015
- Pig
	* Ensembl81 (Sscrofa10.2),		ens81.v9.17Jul2015

NEW SPECIES:
- Chimp
	* Ensembl81 (CHIMP2.1.4), 		ens81.v9.17Jul2015
- Fruitfly
	* Ensembl81 (BDGP6), 			ens81.v9.17Jul2015


___
## 2015_05.v8
```
A-RELEASE: 	2015_05.v8
C-RELEASE: 	4.2.2.3
I-DATE: 	16Apr2015
```

### Highlights
- Hand-over for e80: Mouse, Zebra-fish and Rat.
- New assemblies for Zebra-fish (GRCz10) and Rat (Rnor_6.0).

### APPRIS-SERVER, v4.2
- Fix a bug when UCSC image is shown.	
- Include AWStats is a free powerful and featureful tool that generates advanced web, streaming, ftp or mail server statistics, graphically.
http://apprisws.bioinfo.cnio.es/awstats/awstats.pl?month=all&year=2015&output=main&config=appris&framename=index
		
### Code (2.3)
- APPRIS
	- Fixing a bug where we detect "start/stop" codons
		
### Annotations (v8)
- Human
	* Ensembl80/Gencode22 (GRCh38),	ens80.v8.16Apr2015 = ens79.v8.2Apr2015 
	* Ensembl76/Gencode20 (GRCh38),	ens76.v8.16Apr2015 <- FREEZE!!
- Mouse
	* Ensembl80/GencodeM5 (GRCm38),	ens80.v8.16Apr2015
	* Ensembl74/GencodeM2 (GRCm38),	ens74.v8.16Apr2015 <- FREEZE!!
- Zebra-fish
	* Ensembl80 (GRCz10),			ens80.v8.16Apr2015
- Rat
	* Ensembl80 (Rn6),				ens80.v8.16Apr2015
- Pig
	* Ensembl80 (Sscrofa10.2),		ens80.v8.16Apr2015


___
## 2015_03.v7
```
A-RELEASE: 	2015_03.v7
C-RELEASE: 	4.1.2.2
I-DATE: 		9Feb2015
```

### Highlights
- Fixing a bug in the CORSAIR filter.
- We have increase the number of the version for the release. The reason is the same number.
- The annotations for GRCh37 of human have been frozen: gen19.v7.9Feb2015

### APPRIS-SERVER, v4.1

### Code (2.2)
- APPRIS 
	- Fixing a bug in the CORSAIR filter
	- Fixing a bug in the filter of length sequences.
		
### Annotations (v7)
- Human
	* Ensembl77/Gencode21 (GRCh38), ens77.v7.9Feb2015
	* Ensembl76/Gencode20 (GRCh38), ens76.v7.9Feb2015
	* Ensembl74/Gencode19 (GRCh37), gen19.v7.9Feb2015 <- FREEZE!!
- Mouse
	* Ensembl78/GencodeM4 (GRCm38), ens78.v7.9Feb2015
- Zebra-fish
	* Ensembl77 (Zv9),    ens77.v7.9Feb2015 <- FREEZE!!
- Rat
	* Ensembl77 (Rn5),    ens77.v7.9Feb2015 <- FREEZE!!
- Pig
	* Ensembl77 (Sscrofa10.2), ens77.v7.9Feb2015


___
## 2014_02.v1
```
A-RELEASE: 	2014_02.v1
C-RELEASE: 	4.0.2.1
I-DATE: 		27Jan2014
```

### Highlights
- NEW APPRIS SERVER
- New APPRIS algorithm using oldest CCDS ids as principal isoforms:

### APPRIS-SERVER, v4.0

### Code (2.1)
- New tags of PRINCIPAL ISOFORMS
- New script to retrieve the list of principal isforms for Ensembl
- APPRIS
	- New APPRIS algorithm using oldest CCDS ids as principal isoforms:
		- With the 'appris_candiates' we check:
			1. unique CCDS
			2. the oldest CCDS (smaller than 10numbers comparing with the next CCDS number)
			3. the longest CCDS
			4. the longest protein sequence
				
### Annotations (v1)
- Human
	* Ensembl79 (GRCh38), ens79.v1.26Jan2015 # TEMPORAL
	* Ensembl77/Gencode21 (GRCh38), ens77.v5.27Jan2015
	* Ensembl74/Gencode19 (GRCh37), gen19.v6.27Jan2015
- Mouse
	* Ensembl78/GencodeM4 (GRCm38), ens78.v3.27Jan2015
- Zebra-fish
	* Ensembl77 (Zv9),    ens77.v3.27Jan2015
- Rat
	* Ensembl77 (Rn5),    ens77.v3.27Jan2015
- Pig
	* Ensembl77 (Sscrofa10.2), ens77.v2.27Jan2015
		
