Welcome to APPRIS - A system for annotating alternative splice isoforms
=======================================================================
APPRIS [1] is a system that deploys a range of computational methods to provide value to the annotations of the human genome. APPRIS also selects one of the CDS for each gene as the principal isoform.

APPRIS defines the principal variant by combining protein structural and functional information and information from the conservation of related species.



Directory structure
===================
download
	|
	|___ data (species-release directories)
			|
			|__ homo_sapiens
					|
					|__ ens77.v4.31Oct2014	({GENCODE/Ensembl version}.{annotation version}.{date version})
					|
					|__ gen19.v5.31Oct2014	({GENCODE/Ensembl version}.{annotation version}.{date version})
			|
			|__ mus_musculus
					|
					|__ genM3.v4.31Oct2014	({GENCODE/Ensembl version}.{annotation version}.{date version})
			|
			|__ danio_rerio
					|
					|__ ens77.v2.31Oct2014	({Ensembl version}.{annotation version}.{date version})
			|
			|__ rat_norvegicus
					|
					|__ ens77.v2.31Oct2014	({Ensembl version}.{annotation version}.{date version})
			|
			|__ sus_scrofa
					|
					|__ ens77.v1.13Jan2015	({Ensembl version}.{annotation version}.{date version})

Note: the name of release directories are samples that ilustrate the name structure ({Ensembl version}.{annotation version}.{date version})


Data files
==========
Inside this directory there are some directories coming from different releases of several species (human, mouse, rat, etc).

For each "specie-release" directory, there are the following annotation files:

	== appris_data.principal.txt
		It prints a list of the principal isoforms selected by APPRIS based on a range of protein features.
		
	== appris_data.appris.txt
		APPRIS detects principal isoforms based on a range of methods whose scores are described.

	== appris_data.firestar.gtf.gz
		GTF annotation file with Functional Residues information (see 'GTF files' section for more information)
		
	== appris_data.matador3d.gtf.gz
		GTF annotation file with Tertiary Structure information (see 'GTF files' section for more information)
		
	== appris_data.corsair.gtf.gz
		GTF annotation file with Vertebrate Conservation information (see 'GTF files' section for more information)
		
	== appris_data.spade.gtf.gz
		GTF annotation file with Domain information (see 'GTF files' section for more information)
			
	== appris_data.thump.gtf.gz
		GTF annotation file with Transmembrane Helices information (see 'GTF files' section for more information)
			
	== appris_data.proteo.gtf.gz (only for human)
		GTF annotation file with Proteomic evidences (see 'GTF files' section for more information)
	


List of Principal Isoforms (appris_data.principal.txt)
==========================
APPRIS (Nucleic Acids Res. 2013 41:D110-7) is a system that deploys a range of computational methods to provide value to the annotations of the human genome.
APPRIS also defines the principal variant by combining protein structural and functional information and information from the conservation of related species.

APPRIS has selected a single CDS variant for each gene as the 'PRINCIPAL' isoform based on the range of protein features. 
Principal isoforms are tagged with the numbers 1 to 5, with 1 being the most reliable. The definition of the flags are as follows:

 * "PRINCIPAL:1", (This flag corresponds to the older flag "appris_principal")
 Transcript(s) expected to code for the main functional isoform based solely on the core modules in the APPRIS database. 
 The APPRIS core modules map protein structural and functional information and cross-species conservation to the annotated variants.
 
 * "PRINCIPAL:2", (This flag corresponds to the older flag "appris_candidate_ccds")
 Where the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes), 
 the database chooses two or more of the CDS variants as "candidates" to be the principal variant. 
 
 If one (but no more than one) of these candidates has a distinct CCDS identifier it is selected as the principal variant for that gene. 
 A CCDS identifier shows that there is consensus between RefSeq and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support. 
 
 * "PRINCIPAL:3", (This is a new flag)
 Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants have distinct CCDS identifiers, 
 APPRIS selects the variant with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the earlier it was annotated. 
 
 Consensus CDS annotated earlier are likely to have more cDNA evidence. 
 Consecutive CCDS identifiers are not included in this flag, since they will have been annotated in the same release of CCDS. These are distinguished with the next flag.
 
 * "PRINCIPAL:4", (This flag corresponds to the older flag "appris_candidate_longest_ccds")
 Where the APPRIS core modules are unable to choose a clear principal CDS and there is more than one variant with a distinct (but consecutive) CCDS identifiers, 
 APPRIS selects the longest CCDS isoform as the principal variant. 
 
 * "PRINCIPAL:5", (This flag corresponds to the older flag "appris_candidate_longest_seq")
 Where the APPRIS core modules are unable to choose a clear principal variant and none of the candidate variants are annotated by CCDS, 
 APPRIS selects the longest of the candidate isoforms as the principal variant. 

For genes in which the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes) 
the "candidate" variants not chosen as principal are labeled in the following way:

 * "ALTERNATIVE:1", (This is a new flag)
 Candidate transcript(s) models that are conserved in at least three tested non-primate species.
 
 * "ALTERNATIVE:2", (This is a new flag)
 Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.

Non-candidate transcripts are not flagged and are considered as "Minor" transcripts.


Scores files (appris_data.appris.txt)
============
Tabular file that prints the scores of APPRIS methods. The description of the columns are the following:

* Gene identifier, Ensembl id.

* Transcript identifier, Ensembl id.

* Protein coding label,
	** TRANSLATION, transcript translates to protein.
	** NO_TRANSLATION, transcript does not translate to protein.

* Transcript Status:
	** known transcript is 100% Identical to RefSeq NP or Swiss-Prot entry.
	** A novel transcript shares >60% length with known coding sequence from RefSeq or Swiss-Prot or has cross-species/family support or domain evidence.
	** A putative shares <60% length with known coding sequence from RefSeq or Swiss-Prot, or has an alternative first or last coding exon.
	** A unknown transcript comes from the Ensembl automatic annotation pipeline.

* Transcript Class (or Biotype):
	** A protein coding transcript is a spliced mRNA that leads to a protein product.
	** A processed transcript is a noncoding transcript that does not contain an open reading frame (ORF). This type of transcript is annotated by the VEGA/Havana manual curation project.
	** Nonsense-mediated decay indicates that the transcript undergoes nonsense mediated decay, a process which detects nonsense mutations and prevents the expression of truncated or erroneous proteins.
	** Transcribed pseudogenes and other non-coding transcripts do not result in a protein product.

* Start/Stop codons do not found.

* Consensus CDS identifier of transcript.

* The absolute numbers of functional residues detected (firestar).

* Score related to the number of exon that map to protein structure (Matador3D):
	Each exon that maps to a structure gets a score of 1 if it maps with high identity and no gaps.
	The score is lower if the mapping includes gaps/has low identity.
	Insertions into structures are penalised with -1.

* The number of vertebrate species that have an isoform that aligns to the human isoform over the whole sequence and without gaps (CORSAIR).
	Alignments with the same species scores just 0.5
	We generate multiple alignments against orthologues from a vertebrate protein database.
	We only align a limited number of vertebrate species, chimp, mouse, rat, dog, cow etc.

* The absolute numbers of pfam domains that are detected (SPADE).
	The numbers after the '-' indicate the numbers of partial domains: 'Whole Pfam domains'-'Partial Pfam domains'
	By partial we could mean "broken" or not whole. Some domains will be broken by a splicing event, but many domains are not whole because the pfam domain boundaries do not always describe the domain well.
	
	The reason these domains come after the '-' is to allow us to compare. For example, an alternative isoform may have an insertion that breaks the pfam domain from the main isoform.
	Here it would be possible for pfamscan to detect a single pfam domain for the main variant and two partial domains in the alternative isoform. But in this case a single domain is clearly better than two domains.

* The number of TMH detected (THUMP).
	The numbers after the '-' indicate the numbers of partial TMH: 'Whole TMH'-'Partial TMH (damaged)'	
	By partial we could mean "broken" or not whole. Some TMH will be broken by a splicing event, but many TMH are not whole because the TMH domain boundaries do not always describe the domain well.

* Reliability score for signal peptides (CRASH-SP).
	We use a score of 3 or above as a reliable signal peptide.

* Reliability score for mitochondrial signal sequences (CRASH-TP).
	We use a score of 3 or above as a reliable mitochondrial signal sequences.

* APPRIS score (deprecated).
	
* Reliability score for APPRIS (deprecated).


GTF data files (*.gtf.gz)
==============

Description of the method scores:

	* firestar are the absolute numbers of functional residues detected.
	* Matador3D is a score related to the number of exon that map to structure.
	* CORSAIR shows the number of vertebrate species that have an isoform that aligns to the human isoform over the whole sequence and without gaps (human scores just 0.5).
	* SPADE shows the absolute numbers of pfam domains that are detected. The numbers after the decimal indicate the numbers of partial domains.
	* THUMP shows the number of TMH detected. Again numbers after the decimal indicate partial TMH.
	* CRASH gives a reliability score for signal peptides (we use a score of 3 or above as a reliable signal peptide). Crash-M does the same thing for mitochondrial signal sequences.

These GTF files are text/plain files with tabular format whose columns are the following:

	== FIRESTAR, No. Functional Residues
		- Chromosome
		- Method name
		- Type of annotation: 'functional_residue'
		- Start position of functional residue
		- End position of functional residue
		- Score: for these annotations it is always '0'
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon. Note: peptide position of functional residue

	== MATADOR3D, Tertiary Structure Score
		- Chromosome
		- Method name: 'MATADOR3D'
		- Type of annotation: 'homologous_structure'
		- Start position of tertiary structure
		- End position of tertiary structure
		- Score: based on the number of regions that can be mapped to structural homologues
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon. Note: pdb id.

	== CORSAIR, Conservation Score
		- Chromosome
		- Method name: 'CORSAIR'
		- Type of annotation: 'no_conservation', 'doubtful_conservation', and 'conservation'
		- Start position of transcript
		- End position of transcript
		- Score: that is approximately the number of vertebrate species that can be aligned without introducing gaps
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon

	== SPADE, Whole Domains
		- Chromosome
		- Method name: 'SPADE'
		- Type of annotation: 'domain', 'domain_possibly_damaged', 'domain_damaged', and 'domain_wrong'
		- Start position of domain
		- End position of domain
		- Score: a local pfam domain integrity score which decides whether a domain is damaged or not
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon

	== THUMP, No. Transmembrane Helices
		- Chromosome
		- Method name: 'THUMP'
		- Type of annotation: 'tmh_signal', and 'damaged_tmh_signal'
		- Start position of transmembrane helix
		- End position of transmembrane helix
		- Score: for these annotations it is always '0'
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon


References
==========
[1] Rodriguez JM, Maietta P, Ezkurdia I, Pietrelli A, Wesselink JJ, Lopez G, Valencia A, Tress ML.
APPRIS: annotation of principal and alternative splice isoforms. 
Nucleic Acids Res. 2013 Jan;41(Database issue):D110-7.



Contact
=======
This APPRIS website is powered by the Structural Computational Biology Group at
	Centro Nacional de Investigaciones Oncologicas, (CNIO, http://www.cnio.es)
		and
	Instituto Nacional de Bioinformatica, (INB, http://www.inab.org)

If you have questions or comments, please write to:
	Jose Rodriguez, jmrodriguez@cnio.es
	Michael Tress, mtress@cnio.es.



