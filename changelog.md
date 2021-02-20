___
## 2021_02.v39
```
SERVER-RELEASE:   2021_02.v39
CODE-RELEASE:     4.11.2.25
```

### Highlights

+ New release of Ensembl (e102 version), for the species:
	- Homo sapiens

### Server (4.11)

+ Added an interim news section to the web front page.

### Code (2.25)

+ APPRIS:
	- Added checksum calculation step to the appristools pipeline. This will
	  enable APPRIS releases to include a checksum file for each new dataset.
	- Fixed counting of initial partial codons when the partial codon is not
	  an ambiguous 'X' amino acid.
	- Fixed output of revised APPRIS score in cases where APPRIS advances to
	  the second scoring phase.

+ CORSAIR:
	- Fixed issue with comparison of BLAST query lengths that contain commas.

+ FIRESTAR:
	- Integrated Firestar software updates.
	- Fixed Firestar BED export.

+ SPADE:
	- Fixed issue with identification of domain overlaps.

___
## 2020_11.v37
```
SERVER-RELEASE:   2020_11.v37
CODE-RELEASE:     4.10.2.24
```

### Highlights

+ New release of GENCODE human (version 35).
+ New release of Ensembl (e101 version), for the species:
	- Caenorhabditis elegans
+ Rereleases of GENCODE datasets for human version 34,
  human version 27 and mouse version vM25.

### Server (4.10)

+ Addition of a table of important functional isoforms.
+ Addition of TRIFID scores to gene report in human,
  highlighting of high-scoring isoforms in bold text.
+ Update of the PDB chain ID regex to reflect current PDB chain IDs,
  also allowing for the modified PDB chain IDs used by Matador3D2.
+ Removal of the embedded Google Group from APPRIS website front page.
+ Removal of hyperlink from the Cat_Site_Atl PDB ligand
  in the Firestar annotations section of the gene report.
+ Addition of APPRIS WebServices landing page.
+ Update of documentation and website content,
  primarily to include information on TRIFID.

### Code (2.24)

+ APPRIS:
	- Integration of Alt-CORSAIR method for human and mouse.
	- Recalibration of APPRIS method scoring and weighting.
	- Fixed data concatenation issue affecting chromosomes with many genes.
	- Updated chromosome sizes for rat and Drosophila melanogaster.
	- Matador3D results are used in preference to
	  Matador3D2 results only if they exist.
	- Harmonisation of APPRIS method parameter codes.
	- Fixed error whereby redundant BED attributes
	  were being output for some features.
	- Fixed error whereby some transcripts with multiple values for a given
	  attribute were being output with some attribute values missing.
	- Changed the name of the APPRIS bigBed autoSql config
	  file to reflect the actual number of fields.
	- Fixed bug whereby conversion of residue positions to genomic
	  coordinates was not taking account of the initial and final
	  reading frames of a coding sequence.
	- Changed the code converting residue positions to genomic
	  coordinates, so that microexons shorter than one codon and with
	  an end phase of zero are assigned the residue at that position.
	  Previously this residue would have been assigned to the preceding
	  exon, resulting in a slight misalignment of subsequent exons.

+ CORSAIR:
	- Integration of Alt-CORSAIR method for human and mouse.
	- Improved error-handling for BLAST search.
	- Update of vertebrate database to refseq_201811,
	  invertebrate database to refseq_202010.
	- Cached BLAST result files are specific
	  to the configured BLAST database.
	- Configured E-value is passed to BLAST.

+ FIRESTAR:
	- Fixed error whereby functional residues with multiple ligands
	were being exported to GTF with only one of those ligands.

+ MATADOR3D:
	- Added a minimum length constraint to mini-exons, so that
	  every mini-exon has at least one full codon, except in
	  cases where the corresponding exon is itself shorter
	  than the minimum length constraint.
	- Fixed error whereby Matador3D exons were being exported
	  to BED output alongside the mini-exon alignments, resulting
	  in overlapping BED blocks.
	- For Matador3D and Matador3D2, added a bed3+10 format option
	  for bigBed files, to allow for longer fields where needed.

+ PROTEO:
	- Peptides only contribute to the PROTEO score if they map
	  to the amino-acid sequence of the given isoform.

+ SPADE:
	- When transferring domains between transcripts, preference
	  is given to the external transcript whose transferred domains
	  confer the greatest score increase on the receiving transcript.
	- PfamScan results are parsed using a regex tailored to Pfam-A output.
	- Isoforms with internal stop sites are skipped.
	- Fixed BED export error whereby records were being duplicated.

+ TRIFID:
	- Streamlined TRIFID integration, with the option to choose
	  the 'TRIFID' or 'classic' transcript selection process.
	- Independent configuration of TRIFID data for each dataset.
	- Within APPRIS, the TRIFID module's primary input
	  is taken from the translation sequences file.
	- The identity of each isoform is confirmed by matching the translation
	  sequence against the sequence in the TRIFID prediction file.

___
## 2020_06.v32
```
SERVER-RELEASE:   2020_06.v32
CODE-RELEASE:     4.9.2.23
```

### Highlights

+ New release of Ensembl (e100 version), for the species:
	- Drosophila melanogaster

### Server (4.9)

+ Fixed display of start/stop codon flags.

### Code (2.23)

+ APPRIS:
	- Updated chromosome lists in rat and fruitfly.
	- Fixed bug in handling of start/stop codon flags where transcripts
	  with both start and stop missing were not being flagged.
	- Fixed extraction of transcript annotation from GTF file.
	- Fixed transcript ID filter in gene report.

___
## 2020_06.v31
```
SERVER-RELEASE:   2020_06.v31
CODE-RELEASE:     4.8.2.22
```

### Highlights

+ New release of Ensembl (e100 version), for the species:
	- Homo sapiens
	- Mus musculus

+ New APPRIS annotations using the gene dataset of RefSeq for the species:
	- Homo sapiens (version 109.20200228)

### Server (4.8)

+ Updated to display only queryable datasets in the principal APPRIS web interface.

### Code (2.22)

+ TRIFID:
	- TRIFID predictions are used in the APPRIS decider steps as follows: a transcript
	  is classified as PRINCIPAL:2 if it is assigned the best TRIFID score of the given
	  gene among those transcripts that passed stage 1 and its TRIFID score is at least
	  0.25 greater than any other transcript that passed stage 1. TRIFID predictions are
	  further used to select PRINCIPAL:4 transcripts, such that the transcript with the
	  best TRIFID score among those transcripts that passed stage 3 is selected as the
	  PRINCIPAL:4 transcript.

+ APPRIS:
	- Added queryable dataset option to denote datasets that are accessible in the
	  principal APPRIS web interface.
	- Updated RefSeq data download code to reflect changes to RefSeq FTP site layout.
	- Added optional bed3+9 format for APPRIS bigBed file to accommodate longer merged
	  transcript IDs.
	- Updated proteomics data for GENCODE v33, added upper limit to PROTEO score in BED
	  output.
	- Added filter to exclude PAR_Y genes from GENCODE input.
	- Added compression of APPRIS main GTF and FASTA files.
	- Fixed bugs in Matador3D mini-CDS alignment, tab-delimited file parsing
	  and APPRIS::Gene module.
	- Enabled independent configuration of SignalP-3.0 temp directory.
	- Updated dependency information and documentation.

___
## 2020_02.v30
```
SERVER-RELEASE:   2020_02.v30
CODE-RELEASE:     4.7.2.21
```

### Highlights
+ New release of Ensembl (e99 version), for the species:
	- Homo sapiens
	- Mus musculus

### Server (4.7)
+ No changes

### Code (2.21)

+ APPRIS:
	- Calculation of APPRIS score in 2 phases: first using Spade bitscore,
	  then using Spade domain integrity score if necessary to break ties.

+ FIRESTAR:
	- Effective range used to penalise low-scoring transcripts.

+ MATADOR3D:
	- Effective range used to penalise low-scoring transcripts.
	- Recalibrated alignment identity and gap scores to improve principal transcript identification.
	- Retention of only the first HSP for each PDB BLAST hit to avoid spurious secondary matches
	  near repeat regions.
	- Revised calculation of transcript sequence region position and reading frame to improve
	  transfer of Matador3D annotations between transcripts.

+ PROTEO:
	- Proteomics data incorporates results from an analysis of the data originally produced by:
	  Kim, M., Pinto, S., Getnet, D. et al. (2014) A draft map of the human proteome.
	  *Nature.* 509:575-581. [doi:10.1038/nature13302](https://doi.org/10.1038/nature13302)
	- Leucine and isoleucine are treated as interchangeable during
	  mapping of peptide fragments onto the transcript sequence.

+ SPADE:
	- Effective range of Spade bitscore used to penalise low-scoring transcripts.
	- Spade domain integrity score (effective range 2.5) reintroduced to break ties.

___
## 2019_11.v29
```
SERVER-RELEASE:   2019_11.v29
CODE-RELEASE:     4.7.2.20
```

### Highlights
+ New release of Ensembl (e99 version), for the species:
	- Homo sapiens
	- Mus musculus
	- Drosophila melanogaster

### Server (4.7)
+ No changes

### Code (2.20)
+ No changes

___
## 2019_09.v29
```
SERVER-RELEASE:   2019_09.v29
CODE-RELEASE:     4.7.2.20
```

### Highlights
+ New release of Ensembl (e98 version), for the species:
	- Homo sapiens
	- Mus musculus
	- Sus scrofa

### Server (4.7)
+ No changes

### Code (2.20)
+ No changes

___
## 2019_07.v29
```
SERVER-RELEASE:   2019_07.v29
CODE-RELEASE:     4.7.2.20
```

### Highlights
+ New release of Ensembl (e97 version), for the species:
	- Homo sapiens
	- Mus musculus
	- Caenorhabditis elegans

### Server (4.7)
+ No changes

### Code (2.20)
+ No changes

___
## 2019_02.v29
```
SERVER-RELEASE:   2019_02.v29
CODE-RELEASE:     4.7.2.20
```

### Highlights
+ New release of Ensembl (e96 version), for the species:
	- Homo sapiens
	- Mus musculus
	- Bos taurus
	- Drosophila melanogaster

### Server (4.7)
+ No changes

### Code (2.20)
+ Matador3D: Fixing some bugs:
	- Changing the ranges of decision in the gaps
	- Take into account the frame for the comparision between alignments coming from different isoforms.
	- Bug fixed getting the CDS for Celegans and fruitfly, (eg. FBtr0305654)


___
## 2018_12.v28
```
SERVER-RELEASE:   2018_12.v28
CODE-RELEASE:     4.7.2.19
```

### Highlights
+ RefSeq gene set has been updated for Human. Release 109.

### Server (4.7)
+ No changes

### Code (2.19)
+ No changes
___
## 2018_02.v26
```
SERVER-RELEASE:   2018_02.v26
CODE-RELEASE:     4.7.2.19
```

### Highlights
+ Ensembl 91!! New annotations for the following species:
	Human
	Mouse
	Chimp		(new assembly Pan_tro_3.0)

### Server (4.7)
+ No changes

### Code (2.19)
+ CORSAIR:
	- New release of RefSeq database (201709)

___
## 2017_10.v25
```
SERVER-RELEASE:   2017_10.v25
CODE-RELEASE:     4.7.2.18
```

### Highlights
+ New annotations for the following species and gene sets
	Human: e90v25, rs108v25, up201703v25, a1v25
	Mouse: e90v25, rs106v25, up201610v25, a1v25
	Zebra-fish: e90v25, rs105v25, up201610v25, a1v25
	Rat: e88v25, rs106v25, up201610v25, a1v25
	Pig: e90v25, rs105v25, up201702v25, a1v25
	Chimp: rs104v25, e88v25, up201702v25, a1v25
	Fruitfly: e88v25
	C.elegans: e88v25

### Server (4.7)

### Code (2.18)

+ CORSAIR:
	- Bug fixed: Taking into account the "gaps" for the all length of subject sequence.
	- Now, it does not take into account the exon coordinates.
	- New release of RefSeq database (201709)
+ APPRIS:
	- Now, it does not take into account the exon coordinates.

___
## 2017_08.v24
```
SERVER-RELEASE:   2017_08.v24
CODE-RELEASE:     4.7.2.17
```

### Highlights

+ Ensembl 90. New annotations for the following species:
	Human
	Mouse
	Zebra-fish
	Pig		(new assembly Sscrofa11.1)


### Server (4.7)

### Code (2.17)

+ APPRIS:
	- Now the decision of which Matador3D to use is automatic.
	By default, we use Matador3D2 but If we have genome coordinates, then we use Matador3D.
	- Bug fixed: although the variant has start/stop codon not found, the correct score of CORSAIR will be printed.
	But the final annotation of APPRIS will be the same.

+ SPADE
	- New database with the curated domains of Pfam-A (pfam_201706)

___
## 2017_06.v23
```
SERVER-RELEASE:   2017_06.v23
CODE-RELEASE:     4.7.2.16
```

### Highlights

+ Creation of APPRIS gene-set (version a1) for all vertebrate species located: human, mouse, zebra-fish, rat, pig, chimp.
The APPRIS gene-set has been developed from the join of the gene sets of Ensembl and RefSeq, and enrich with the isoforms of UniProt.

+ New APPRIS annotations using the gene dataset of RefSeq for the following species:
	Human       (version 108)
	Mouse       (version 106)
	Zebra-fish  (version 105)
	Rat         (version 106)
	Pig         (version 105)
	Chimp       (version 104)

+ New APPRIS annotations using the gene names annotated in UniProt for the following species:
	Human       (version 2017_03)
	Pig         (version 2017_02)
	Chimp       (version 2017_02)


### Server (4.7)

+ Improvements in the website to include the annotations for the APPRIS gene data-set

+ PROTEO
	- Bug fixed printing the number of peptides. Eg. PDE3B

+ THUMP
	- Bug Fixed: printing the TMH
	Eg, http://appris-dev.bioinfo.cnio.es/#/database/id/homo_sapiens/ensembl:ENSG00000156869+refseq:391059?as=hg38&sc=appris&ds=a1v23		

+ Bug fixed: printing annotations in the method pop-ups of Sequence Browser when the transcript id is too long.
Eg. "ensembl:ENSG00000180509+refseq:3753"

+ Print the whole identifiers in the "Principal Isoforms" panel of the website report.
Eg. http://appris-dev.bioinfo.cnio.es/#/database/id/homo_sapiens/ensembl:ENSG00000224712+refseq:642778?as=hg38&sc=appris&ds=a1v23

+ Matador3D: print the Matador3D and Matador3D2 residues when the dataset comes from Ensembl and UniProt/APPRIS, respectively.
Eg. from ensembl dataset http://appris-dev.bioinfo.cnio.es/#/database/id/homo_sapiens/ENSG00000196323?as=hg38&sc=ensembl&ds=e88v22
Eg. from appris dataset http://appris-dev.bioinfo.cnio.es/#/database/id/homo_sapiens/ensembl:ENSG00000196323+refseq:29068?as=hg38&sc=appris&ds=a1v23
Now, in this examples, is printing Matador3D2 in both cases. Not Matador3D and Matador3D2


### Code (2.16)

+ SPADE
	- Bugs fixed: parsing the Pfam outputs for APPRIS databases.

+ Bug fixed: inserting genes as "ensembl:ENSG00000180509+refseq:3753" from the appris gene dataset (MSG: No a correct stable id).

+ Bug fixed: inserting firestar residues into database. for genes in appris dataset.

+ Matador3D2
	- New DB (PDB_70_with_pdb_seq)

+ APPRIS
	- Matador3D2 is used for appris gene set and UniProt. Matador3D is used for the rest of the data-sets (Ensembl, RefSeq).

+ PROTEO
	- Bug fixed: saving the protemic experiments into appris results
	Eg. http://appris-dev.bioinfo.cnio.es/#/database/id/homo_sapiens/ENSG00000049540?as=hg38&sc=ensembl&ds=e88v22

+ CORSAIR
	- Tuning the conservation score when there are not exon information
	Eg. NM_001012720 in ensembl:ENSG00000148604+refseq:5995
	Eg. NM_001122842 in ensembl:ENSG00000111912+refseq:135112
	Eg. XM_006722069 in ensembl:ENSG00000184922+refseq:752

+ THUMP
	- Bug Fixed: getting on consensus section.
	Eg, http://appris-dev.bioinfo.cnio.es/#/database/id/homo_sapiens/ensembl:ENSG00000156869+refseq:391059?as=hg38&sc=appris&ds=a1v23



___
## 2017_05.v22
```
SERVER-RELEASE:   2017_05.v22
CODE-RELEASE:     4.6.2.15
```

### Highlights

+ New release of APPRIS for Ensembl (e88) -all local species-.

+ New APPRIS annotations using the gene dataset of RefSeq for the following species:
	Mouse 		(version 106)
	Zebra-fish	(version 105)
	Rat 		(version 106)
	Pig			(version 105)
	Chimp		(version 104)

+ New APPRIS annotations using the gene names annotated in UniProt for the following species:
	Human 		(version 2017_03)
	Pig			(version 2017_02)
	Chimp		(version 2017_02)

### Code (2.15)

* The pipeline is running based on the config.json file, from now.

* New cache files based on date.

* Matador3D:
	- Update of PDB database

* THUMP
	- Bug Fix on consensus section. Eg, PDE3B, ENSMUSG00000098306.

* APPRIS
	- Includes both Matador3D and Matador3D2 methods. Then, appris chooses from the 'pipeline.ini' file.
	- New algorithm. Take into account the transcripts that have passed the filters of methods previously.
	- Decrease the weight of Matador3D.



___
## 2017_04.v21
```
SERVER-RELEASE:   2017_04.v21
CODE-RELEASE:     4.6.2.14
```

### Highlights

- Only Annotations of human have been updated for Ensembl87, Gencode19, RefSeq 105 and RefSeq 107.


### Code (2.14)

- Matador3D (version 1) has been reconnected in the APPRIS pipeline

- Include the two methods of Matador3D but Matador3D (v1) is used in the pipeline.

- Matador3D2: Fixing a bug when a domain sequences was not defined, it creates very high scores.

- SPADE: Fixing a bug with the overlapping of exported alignments.

- APPRIS: Improving the algorithm
In the cases when all sequences are PRINCIPAL:5, APPRIS firstly selects the principal isoform by the following decision:
	The transcript validated manually (NM_ in RefSeq) and the transcript has not to be automatic transcript (with -2XX transcript name in Ensembl).
In addition, APPRIS makes a decision when the all final transcripts have the same length of sequences and the same APPRIS score based on the smaller id.
Eg, '162966' for human RefSeq 107

- CRASH
	- Fix the bug why sometimes it executed and others not. (Eg, ENSG00000066654/ENST00000381337)

- Fixed bug reporting the APPRIS annotations for UniProt

- Update the code for the creation of the APPRIS gene data-sets.

- Now APPRIS accepts the few shorted sequences.


___
## 2017_01.v20
```
SERVER-RELEASE:   2017_01.v20
CODE-RELEASE:     4.6.2.13
```

### Highlights
- New release of Ensembl (e87) for all species.
- New annotations for RefSeq human gene data-set (versions 105 and 107 for the assemblies GRCh37 and GRCh38, respectively).
- APPRIS annotations for UniProt human, mouse, zebra-fish, and rat. Proteome versions:
	Human (2016_06)
	Mouse (2016_10)
	Zebra-fish (2016_10)
	Rat (2016_10)
- New version of Matado3D method (v2).

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

	- Fixed bug: It is wrong the assignation of Principal in the step of eldest CCDS and TSL. Eg. BAIAP2

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
