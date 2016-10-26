## 2016_11.v19
```
#################################
# SERVER-RELEASE:   2016_11.v19 #
# CODE-RELEASE:     4.5.2.12    #
#################################
```

#### Highlights
	
#### Code (2.11)
- CORSAIR:
A variant with "start/stop codons not found" can not to win.

__

## 2016_07.v18
```
SERVER-RELEASE:   2016_07.v18
CODE-RELEASE:     4.4.2.10
```

#### Highlights
- New release of Ensembl (Ensembl 85) for human, mouse, and rat.
	
#### Code (2.11)
- Matador3D v2: Has been reported to UniProt.
	
__

## 2016_07.v17
```
SERVER-RELEASE:   2016_07.v17
CODE-RELEASE:     4.5.2.11
```

#### Highlights
- New release of Ensembl (Ensembl 85) for human, mouse, and rat.
- Freeze the code. [4.4.2.9] (https://github.com/appris/appris/tree/freeze_4.4.2.9)

__

## 2016_06.v17
```
SERVER-RELEASE:   2016_06.v17
CODE-RELEASE:     4.5.2.11
```

#### Highlights
- RefSeq versions 105 and 107 for the assemblies GRCh37 and GRCh38, respectively. 
	
#### APPRIS-SERVER, v4.5
- Change in the web services to accepts multiple datasets sources.

- Now the search web services retrieves the results for a specific query.

- Add TSL annots in the report web page:

#### APPRIS-CODE, v2.11
- The file names of method results have changed.

- Add TSL and TAG features in "entity" table.

- Change the "namespace" of "datasource" table.

- Delete HAVANA datasource

- APPRIS
	- Fixed bug. No decision when the sequence length is equal, with the same appris score. The list of genes was:
		ENSG00000070831 ENSG00000103248 ENSG00000132639 ENSG00000163632 ENSG00000265590 ENSG00000267127 ENSG00000268434
		
- Create "server.json" config file for server.

- Create script "appristools_srv" to update the APPRIS release

==APPRIS_RELEASE:RELNOTES:2016_04.v16==
```
A-RELEASE: 	2016_04.v16
C-RELEASE: 	4.4.2.10
```

#### Highlights
- RefSeq versions 105 and 107 for the assemblies GRCh37 and GRCh38, respectively.
	
#### APPRIS-CODE, v2.10

- Add "tag" field into database and return its values like "readthought_transcripts"

- Script that add extra values into dataset file

- APPRIS
	- Use the tag features to identify the "readthought_transcripts"
	
	- Fixed bug: It is wrong the assignation of Principal in the step of eldest CCDS and TSL. Eg. BAIAP2 TODO????? It's working correctly
	
	- Check tsl order in ensemble and gencode. Eg. BBS9
	
	- Now TSL annotations is obtained from gene annotation file. The function that extracts the TSL is deprecated.
		
==APPRIS_RELEASE:DATAFILES:2016_04.v16==
homo_sapiens				GRCh38			rs107.v16				appris_homo_sapiens_rs107v16
homo_sapiens				GRCh37			rs105.v16				appris_homo_sapiens_rs105v16
homo_sapiens				GRCh38			ens84.v16				appris_homo_sapiens_e84v16
homo_sapiens				GRCh37			gen19.v16				appris_homo_sapiens_g19v16
homo_sapiens				GRCh38			uni201509.v15			appris_homo_sapiens_uni201509v15
homo_sapiens				GRCh38			ens81.v15				appris_homo_sapiens_e81v15
homo_sapiens				GRCh38			ens79.v15				appris_homo_sapiens_e79v15
mus_musculus				GRCm38			ens84.v15 				appris_mus_musculus_e84v15
danio_rerio					GRCz10			ens84.v15 				appris_danio_rerio_e84v15
danio_rerio					Zv9				ens77.v7.9Feb2015 		appris_danio_rerio_e77
rattus_norvegicus			Rnor_6.0		ens84.v15				appris_rattus_norvegicus_e84v15
rattus_norvegicus			Rnor_5.0		ens77.v7.9Feb2015		appris_rattus_norvegicus_e77
sus_scrofa					Sscrofa10.2		ens84.v15				appris_sus_scrofa_e84v15
pan_troglodytes				CHIMP2.1.4		ens84.v15				appris_pan_troglodytes_e84v15
drosophila_melanogaster		BDGP6			ens84.v15				appris_drosophila_melanogaster_e84v15
caenorhabditis_elegans		WBcel235		ens84.v15				appris_caenorhabditis_elegans_e84v15

==APPRIS_RELEASE:RELNOTES:2016_03.v15==
```
A-RELEASE: 	2016_03.v15
C-RELEASE: 	4.4.2.9
```

#### Highlights
- New name for each release... Forget the date
- New algorithm for SPADE: SPADE2
	
#### APPRIS-CODE, v2.9

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
		
==APPRIS_RELEASE:DATAFILES:2016_03.v15==
homo_sapiens				ens84.v15	GRCh38		appris_homo_sapiens_e84v15
homo_sapiens				gen19.v15	GRCh37		appris_homo_sapiens_g19v15
homo_sapiens				ens81.v15	-		appris_homo_sapiens_e81v15
homo_sapiens				ens79.v15	-		appris_homo_sapiens_e79v15
mus_musculus				ens84.v15 	GRCm38		appris_mus_musculus_e84v15
danio_rerio					ens84.v15 	GRCz10		appris_danio_rerio_e84v15
danio_rerio					ens77.v7.9Feb2015 	Zv9			appris_danio_rerio_e77
rattus_norvegicus			ens84.v15	Rnor_6.0	appris_rattus_norvegicus_e84v15
rattus_norvegicus			ens77.v7.9Feb2015	Rnor_5.0	appris_rattus_norvegicus_e77
sus_scrofa				ens84.v15	Sscrofa10.2	appris_sus_scrofa_e84v15
pan_troglodytes				ens84.v15	CHIMP2.1.4	appris_pan_troglodytes_e84v15
drosophila_melanogaster			ens84.v15	BDGP6		appris_drosophila_melanogaster_e84v15
caenorhabditis_elegans			ens84.v15	WBcel235	appris_caenorhabditis_elegans_e84v15
