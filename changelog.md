___
## 2016_11.v19
```
SERVER-RELEASE:   2016_11.v19
CODE-RELEASE:     4.5.2.12
```

### Highlights
	
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
