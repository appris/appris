Directory structure
===================
This directory contains the following subdirectories:

current_release
	|	
	|___ ensembl_datafiles		* APPRIS annotations from Ensembl/GENCODE gene datasets
			|
			|__ species			* annotations per species
	|
	|___ releases				* Dir with all releases
			|
			|__ old				* Old releases. The structure of dirs will not maintain
	

data (species-release directories) # DEPRECATED
			|
			|__ homo_sapiens
					|
					|__ {version directory}	({GENCODE/Ensembl version}.{annotation version}.{date version})
			|
			|__ mus_musculus
					|
					|__ {version directory}	({GENCODE/Ensembl version}.{annotation version}.{date version})
			|
			|__ danio_rerio
					|
					|__ {version directory}	({GENCODE/Ensembl version}.{annotation version}.{date version})
			|
			|__ rat_norvegicus
					|
					|__ {version directory}	({GENCODE/Ensembl version}.{annotation version}.{date version})
			|
			|__ sus_scrofa
					|
					|__ {version directory}	({GENCODE/Ensembl version}.{annotation version}.{date version})

'data' directory is deprecated. The directory with the last annotations are in: 'download/current_release'.
If you want the directory with all releases are in: 'download/releases'
