Data files
==========
In this directory we would save the *__data files__* from APPRIS annotations.


Directory structure
===================
This directory contains the following subdirectories:

```
current_release /* Release in production */
	|
	|___ datafiles
			|
            |___ {species_name}
                    |
                    |__ {datasource_version}
                                |
                                |__ /data_files/

releases        /* Dir with all releases */
    |
    |___ {APPRIS_version}
            |
            |___ datafiles
                    |
                    |___ {species_name}
                                |
                                |__ {datasource_version}
                                            |
                                            |__ /data_files/

trackHub        /* Dir with all releases */
    |
    |___ {genome_assembly_versions_for_species}
                        |
                        |__ /trackhub_data_files/


trackHubPROTEO  /* Dir with all releases */
    |
    |___ {human_genome_assemblies}
                    |
                    |__ /trackhub_data_files/

```

> __Note__: Internally, when APPRIS pipeline is running, the data directory would be different. See below under section [Directory structure running the pipeline](#directory-structure-running-the-pipeline)

APPRIS data files
=================
The following files are the data files retrieved by APPRIS pipeline.

+ *__appris_data.principal.txt__*,
    It prints a list of the principal isoforms selected by APPRIS based on a range of protein features.

+ *__appris_data.appris.txt__*,
    APPRIS detects principal isoforms based on a range of methods whose scores are described.

+ *__appris_data.firestar.gtf.gz__*,
    GTF annotation file with Functional Residues information (see 'GTF data files' section for more information)

+ *__appris_data.matador3d.gtf.gz__*,
    GTF annotation file with Tertiary Structure information (FMI: see 'GTF data files' section)

+ *__appris_data.corsair.gtf.gz__*,
    GTF annotation file with Vertebrate Conservation information (FMI: see 'GTF data files' section)

+ *__appris_data.spade.gtf.gz__*,
    GTF annotation file with Domain information (FMI: see 'GTF data files' section)

+ *__appris_data.thump.gtf.gz__*,
    GTF annotation file with Transmembrane Helices information (FMI: see 'GTF data files' section)

+ *__appris_data.proteo.gtf.gz__* (only for human GENCODE/Ensembl),
    GTF annotation file with Proteomic evidences (FMI: see 'GTF data files' section)


List of Principal Isoforms (appris_data.principal.txt)
==========================
APPRIS (Nucleic Acids Res. 2013 41:D110-7) is a system that deploys a range of computational methods to provide value to the annotations of the human genome.
APPRIS also defines the principal variant by combining protein structural and functional information and information from the conservation of related species.

For more information, see 'Principal Isoforms Flags' section.


Principal Isoforms Flags
========================
APPRIS has selected a single CDS variant for each gene as the 'PRINCIPAL' isoform based on the range of protein features.
Principal isoforms are tagged with the numbers 1 to 5, with 1 being the most reliable. The definition of the flags are as follows:

 * "PRINCIPAL:1",
 Transcript(s) expected to code for the main functional isoform based solely on the core modules in the APPRIS database.
 The APPRIS core modules map protein structural and functional information and cross-species conservation to the annotated variants.

 * "PRINCIPAL:2",
 Where the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes),
 the database chooses two or more of the CDS variants as "candidates" to be the principal variant.

 If one (but no more than one) of these candidates has a distinct CCDS identifier it is selected as the principal variant for that gene.
 A CCDS identifier shows that there is consensus between RefSeq and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support.

 * "PRINCIPAL:3",
 Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants have distinct CCDS identifiers,
 APPRIS selects the variant with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the earlier it was annotated.

 Consensus CDS annotated earlier are likely to have more cDNA evidence.
 Consecutive CCDS identifiers are not included in this flag, since they will have been annotated in the same release of CCDS. These are distinguished with the next flag.

 * "PRINCIPAL:4",
 Where the APPRIS core modules are unable to choose a clear principal CDS and there is more than one variant with a distinct (but consecutive) CCDS identifiers,
 APPRIS selects the longest CCDS isoform as the principal variant.

 * "PRINCIPAL:5",
 Where the APPRIS core modules are unable to choose a clear principal variant and none of the candidate variants are annotated by CCDS,
 APPRIS selects the longest of the candidate isoforms as the principal variant.

For genes in which the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes)
the "candidate" variants not chosen as principal are labeled in the following way:

 * "ALTERNATIVE:1",
 Candidate transcript(s) models that are conserved in at least three tested non-primate species.

 * "ALTERNATIVE:2",
 Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.

Non-candidate transcripts are not flagged and are considered as "MINOR" transcripts.


Scores files (appris_data.appris.txt)
============
Tabular file that prints the scores of APPRIS methods. The description of the columns are the following:

* Gene identifier:
	Ensembl id, RefSeq id, or UniProt entry.

* Transcript identifier:
	Ensembl id, RefSeq id, or UniProt entry.

* Protein coding label,
	** TRANSLATION, transcript translates to protein.
	** NO_TRANSLATION, transcript does not translate to protein.

* Transcript Status (DEPRECATED):
	** known transcript is 100% Identical to RefSeq NP or Swiss-Prot entry.
	** A novel transcript shares >60% length with known coding sequence from RefSeq or Swiss-Prot or has cross-species/family support or domain evidence.
	** A putative shares <60% length with known coding sequence from RefSeq or Swiss-Prot, or has an alternative first or last coding exon.
	** A unknown transcript comes from the Ensembl automatic annotation pipeline.

* Transcript Class (or Biotype):
	** A protein coding transcript is a spliced mRNA that leads to a protein product.
	** A processed transcript is a noncoding transcript that does not contain an open reading frame (ORF). This type of transcript is annotated by the VEGA/Havana manual curation project.
	** Nonsense-mediated decay indicates that the transcript undergoes nonsense mediated decay, a process which detects nonsense mutations and prevents the expression of truncated or erroneous proteins.
	** Transcribed pseudogenes and other non-coding transcripts do not result in a protein product.

* Start/Stop codons do not found:
	** start means 'Start codon does not found'
	** stop means 'Start codon does not found'

* Consensus CDS identifier (CCDS):
 	The Consensus CDS (CCDS) project is a collaborative effort to identify a core set of human and mouse protein coding regions that are consistently annotated and of high quality.
 	The long term goal is to support convergence towards a standard set of gene annotations.

 	For more information:
	https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi

* Transcript Support Level (TSL):
	The method relies on the primary data that can support full-length transcript structure: mRNA and EST alignments supplied by UCSC and Ensembl.
	The following categories are assigned to each of the evaluated annotations:
		tsl1 � all splice junctions of the transcript are supported by at least one non-suspect mRNA
		tsl2 � the best supporting mRNA is flagged as suspect or the support is from multiple ESTs
		tsl3 � the only support is from a single EST
		tsl4 � the best supporting EST is flagged as suspect
		tsl5 � no single transcript supports the model structure
		tslNA � the transcript was not analysed for one of the following reasons:
			pseudogene annotation, including transcribed pseudogenes
			human leukocyte antigen (HLA) transcript
			immunoglobin gene transcript
			T-cell receptor transcript
			single-exon transcript (will be included in a future version)

	For more information:
	http://www.ensembl.org/Help/Glossary?id=492

* The absolute numbers of functional residues detected (firestar).

* Score related to the number of exon that map to protein structure (Matador3Dv2):
	Matador3Dv2 analyses protein structural information for each variant. The number represents the sum of bitscores in PDB alignment.

* The number of vertebrate species that have an isoform that aligns to the human isoform over the whole sequence and without gaps (CORSAIR).
	Alignments with the same species scores just 0.5
	We generate multiple alignments against orthologues from a vertebrate protein database.
	We only align a limited number of vertebrate species, chimp, mouse, rat, dog, cow etc.

* The absolute numbers of pfam domains that are detected (SPADE):
	SPADE identifies the functional domains present in a transcript and detects those domains that are damaged (not whole). The number represents the sum of bitscores in Pfam alignment.

* The number of TMH detected (THUMP).
	The numbers after the '-' indicate the numbers of partial TMH: 'Whole TMH'-'Partial TMH (damaged)'
	By partial we could mean "broken" or not whole. Some TMH will be broken by a splicing event, but many TMH are not whole because the TMH domain boundaries do not always describe the domain well.

* Reliability score for signal peptides and mitochondrial signal sequences (CRASH):
	We use a score of 3 or above as a reliable signal peptide, and mitochondrial signal sequences (separated by comma).

* The number of exons with unusual evolutionary rats (INERTIA): DEPRECATED!!
	INERTIA uses three alignment methods to generate cross-species alignments, from which SLR identifies exons with unusual evolutionary rates.

* APPRIS score:
	Reliability score for the variants based on the scores of methods and a weight for them.

* Reliability labels for APPRIS:
	See 'Principal Isoforms Flags' section.


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

	+ FIRESTAR, No. Functional Residues
		- Chromosome
		- Method name
		- Type of annotation: 'functional_residue'
		- Start position of functional residue
		- End position of functional residue
		- Score: for these annotations it is always '0'
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon. Note: peptide position of functional residue

	+ MATADOR3D, Tertiary Structure Score
		- Chromosome
		- Method name: 'MATADOR3D'
		- Type of annotation: 'homologous_structure'
		- Start position of tertiary structure
		- End position of tertiary structure
		- Score: based on the number of regions that can be mapped to structural homologues
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon. Note: pdb id.

	+ CORSAIR, Conservation Score
		- Chromosome
		- Method name: 'CORSAIR'
		- Type of annotation: 'no_conservation', 'doubtful_conservation', and 'conservation'
		- Start position of transcript
		- End position of transcript
		- Score: that is approximately the number of vertebrate species that can be aligned without introducing gaps
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon

	+ SPADE, Whole Domains
		- Chromosome
		- Method name: 'SPADE'
		- Type of annotation: 'domain', 'domain_possibly_damaged', 'domain_damaged', and 'domain_wrong'
		- Start position of domain
		- End position of domain
		- Score: a local pfam domain integrity score which decides whether a domain is damaged or not
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon

	+ THUMP, No. Transmembrane Helices
		- Chromosome
		- Method name: 'THUMP'
		- Type of annotation: 'tmh_signal', and 'damaged_tmh_signal'
		- Start position of transmembrane helix
		- End position of transmembrane helix
		- Score: for these annotations it is always '0'
		- Strand
		- Frame: for these annotations it is always '.'
		- Attributes end in a semicolon


Directory structure running the pipeline
========================================
Internally, when APPRIS pipeline is running, the data directory would be different:

```
{APPRIS_version}
        |
        |___ {species_name}
                    |
                    |__ {datasource_version}
                                |
                                |__ /data_files/
```
