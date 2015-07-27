Scripts/Binaries description
============================

* BIN

	>> apprisall, global script that executes APPRIS

	>> appris_retrieve_stats, retrieves statistics from APPRIS methods

	>> appris_insert_appris, inserts annotations into APPRIS database
	
	>> appris_retrieve_annots, retrieves annotations from APPRIS methods. Files that they will be downloaded from the website

	>> appris_check_appris, checks and gets runtimes from annotations.
	
	>> appris_export_exon_data, exports exon list whose annotations are the following:
		- PRINCIPAL/ALTER: exon that belongs to principal isoform/alternative isoform.
		- CONSERVE/NO_CONSERVE: exon that belongs isoform with (without) vertebrate conservation.
		- OVERLAP/NO_OVERLAP: exon that overlaps with (without) the exons of all isoforms.
	
	>> appris_export_method_annots, exports annotations for the following method:
		- Matador3D: list of pdb's, and coordinate alignments (at the moment, it's reported cds coordinates).	
		- SPADE: list of domains, and alignment coordinates.
		- INERTIA: list of exons and the best score of slr (omega), the smaller score.
	
	>> appris_gencode_report, retrieve statistics from gencode report
	
	>> appris_ls_annots, count the result files of appris
	
	
* SCRIPTS

	>> get_diff_genes.pl: retrives the list of genes that are differents among versions of GENCODE.

