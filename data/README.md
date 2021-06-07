Data files
==========
In this directory we would save the *__data files__* from APPRIS annotations.

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
    It prints a list of the principal isoforms selected by APPRIS based on a range of protein features (see [*__Principal Isoforms Flags__*](#principal-isoforms-flags) section, for more information)

+ *__appris_data.appris.txt__*,
    APPRIS detects principal isoforms based on a range of methods whose scores are described (see [*__Score files__*](#score-files) section, for more information).

+ *__appris_data.appris_label.txt__*,
    Labels (_YES,NO,UNKNOWN_) of methods that establish the decision for the principal isoforms.

+ *__appris_stats.ccds.txt__*,
    Statistics with the gene coverage of methods, comparison with the CCDS, and number with the principal isoform
    decisions.

+ *__appris_data.{method}.gtf.gz__*,
    Annotations of *_methods_* mapped into genome based on *_GTF_* format (see [*__GTF data files__*](#gtf-data-files) section for more information)

+ *__appris_data.{method}.bed.gz__*,
    Annotations of *_methods_* mapped into genome based on *_BED_* format (see [*__BED data files__*](#bed-data-files) section for more information)

+ *__appris_data.{method}.bed12.gz__*,
    Annotations of *_methods_* mapped into genome based on *_BED12_* format (see [*__BED12 data files__*](#bed12-data-files) section for more information)

+ *__appris_data.{method}.bb__*,
    Binary files with the annotations of *_methods_* mapped into genome based on *_bigBed_* format (see [*__bigBed data files__*](#bigbed-data-files) section for more information)


Principal Isoform Flags
=======================

APPRIS selects a single CDS variant for each gene as the 'PRINCIPAL' isoform based on the range of
protein features. Note that since multiple transcripts may have the same CDS, it is possible for
more than one transcript to be annotated as the principal variant.

Principal isoforms are tagged with the numbers 1 to 5, with 1 being the most reliable. Note that
under the classic transcript selection system (see below) we also regard PRINCIPAL:2 isoforms as
highly reliable, while under the TRIFID transcript selection system (see below) we regard
PRINCIPAL:2 and PRINCIPAL:3 isoforms as highly reliable.

+ __PRINCIPAL:1__
 Transcript(s) predicted to code for the main functional isoform based solely on the core modules
 in the APPRIS database. Core modules map protein structural and functional information and
 cross-species conservation to the annotated transcripts.

At this stage many transcripts are flagged as __"MINOR"__ transcripts. MINOR transcripts are not
considered in any subsequent steps.

When APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human
protein-coding genes), the database chooses two or more of the CDS variants as "candidates" to be
the principal variant. These are assigned tags from PRINCIPAL:2 to PRINCIPAL:5, or ALTERNATIVE:1 or
ALTERNATIVE:2, according to one of two selection processes The 'TRIFID' selection process is used
on datasets for which TRIFID scores are available, (e.g. Gencode34/Ensembl100 human annotation),
while the 'Classic' selection process is used for all other datasets.

### TRIFID transcript selection

+ __PRINCIPAL:2__
 Given two or more principal transcript candidates identified by the APPRIS core modules, APPRIS
 can tag a candidate transcript as PRINCIPAL:2 if it has the best raw TRIFID score and its raw TRIFID
 score is greater than any other candidate by a sufficiently wide margin. TRIFID is a machine learning
 classifier based on APPRIS and external database input which we have found to be a reliable
 indicator of isoform functionality. For more information about TRIFID, please consult the [TRIFID
 GitLab repository](https://gitlab.com/bu_cnio/trifid).

+ __PRINCIPAL:3__
 Where the APPRIS core modules are unable to choose a clear principal variant and no single variant
 has a TRIFID score exceeding that of the other candidates by a sufficient margin, APPRIS selects
 from among the remaining candidates the variant with the most supporting proteomics evidence.

+ __PRINCIPAL:4__
 Where the APPRIS core modules are unable to choose a clear principal variant and no single isoform
 from among the remaining candidates has strong proteomics evidence support or a dominant TRIFID
 score, APPRIS selects the remaining candidate with the best raw TRIFID score.

+ __PRINCIPAL:5__
 Where the APPRIS core modules are unable to choose a clear principal variant and none of the
 candidate variants have a sufficiently strong proteomics evidence base or TRIFID score, APPRIS
 selects the longest of the candidate isoforms as the principal variant.

 In the small number of cases where two or more candidates remain that cannot be distinguished even
 by isoform length, a simple string sort is applied to the transcript identifiers of remaining
 candidates, and the candidate with the first-sorting transcript ID is arbitrarily selected as the
 principal variant.

### Classic transcript selection

+ __PRINCIPAL:2__
 When APPRIS core modules are unable to choose a clear principal variant and one of the principal
 transcript candidates has a distinct CCDS identifier, it is selected as the principal variant for
 that gene. A CCDS identifier shows that there is consensus between RefSeq and GENCODE/Ensembl for
 that variant, guaranteeing that the variant has cDNA support.

+ __PRINCIPAL:3__
 Where the APPRIS core modules are unable to choose a clear principal variant and multiple
 candidates have distinct CCDS identifiers, APPRIS selects the variant with the lowest CCDS
 identifier as the principal variant. The lower the CCDS identifier, the earlier it was annotated.

 Consensus CDS annotated earlier are likely to have more cDNA evidence. Consecutive CCDS
 identifiers are not included in this flag, since they will have been annotated in the same
 release of CCDS. These are distinguished with the next flag.

 In addition, if there are multiple variants with distinct (but consecutive) CCDS identifiers,
 APPRIS chooses the variant for which all splice junctions are supported by at least one
 non-suspect mRNA. This information is reported by the method Transcript Support Level (TSL),
 which is a method to highlight well-supported and poorly-supported transcript models for users.
 The method relies on the primary data that can support full-length transcript structure: mRNA and
 EST alignments supplied by UCSC and Ensembl.

+ __PRINCIPAL:4__
 Where the APPRIS core modules are unable to choose a clear principal CDS, there are multiple
 variants with distinct (but consecutive) CCDS identifiers, and all the splice junctions are not
 well-supported, APPRIS selects the longest CCDS isoform as the principal variant.

+ __PRINCIPAL:5__
 Where the APPRIS core modules are unable to choose a clear principal variant and none of the
 candidate variants are annotated by CCDS, APPRIS selects the longest of the candidate isoforms
 as the principal variant.

 In the small number of cases where two or more candidates remain that cannot be distinguished even
 by isoform length, a simple string sort is applied to the transcript identifiers of remaining
 candidates, and the candidate with the first-sorting transcript ID is arbitrarily selected as
 the principal variant.

For genes with PRINCIPAL:2 to PRINCIPAL:5 isoforms "candidate" variants not chosen as principal
are labeled in the following way:

+ __ALTERNATIVE:1__
 Candidate transcript(s) models that are conserved in at least three tested non-primate species.

+ __ALTERNATIVE:2__
 Candidate transcript(s) models that appear to be conserved in fewer than three tested
 non-primate species.


Score files
============
Tabular file that prints the scores of APPRIS methods. The description of the columns are the following:

+ __Gene identifier__:
	Ensembl id, RefSeq id, or UniProt entry.

+ __Gene name__

+ __Transcript identifier__:
	Ensembl id, RefSeq id, or UniProt entry.

+ __Protein identifier__:
	Ensembl id, RefSeq id, or UniProt entry.

+ __Protein coding label__:
	- TRANSLATION, transcript translates to protein.
	- NO_TRANSLATION, transcript does not translate to protein.

+ __Transcript Class (or Biotype)__:
	- A protein coding transcript is a spliced mRNA that leads to a protein product.
	- A processed transcript is a noncoding transcript that does not contain an open reading frame (ORF). This type of transcript is annotated by the VEGA/Havana manual curation project.
	- Nonsense-mediated decay indicates that the transcript undergoes nonsense mediated decay, a process which detects nonsense mutations and prevents the expression of truncated or erroneous proteins.
	- Transcribed pseudogenes and other non-coding transcripts do not result in a protein product.

+ __Start/Stop codons do not found__:
	- start means 'Start codon does not found'
	- stop means 'Start codon does not found'

+ __Consensus CDS identifier (CCDS)__:
 	The Consensus CDS ([CCDS](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi)) project is a collaborative effort to identify a core set of human and mouse protein coding regions that are consistently annotated and of high quality.
 	The long term goal is to support convergence towards a standard set of gene annotations.

+ __Transcript Support Level (TSL)__:
	The method relies on the primary data that can support full-length transcript structure: mRNA and EST alignments supplied by UCSC and Ensembl.
	The following categories are assigned to each of the evaluated annotations:
	- _tsl1_,  all splice junctions of the transcript are supported by at least one non-suspect mRNA
    - _tsl2_,  the best supporting mRNA is flagged as suspect or the support is from multiple ESTs
    - _tsl3_,  the only support is from a single EST
    - _tsl4_,  the best supporting EST is flagged as suspect
    - _tsl5_,  no single transcript supports the model structure
    - _tslNA_, the transcript was not analysed for one of the following reasons:
		- pseudogene annotation, including transcribed pseudogenes
        - human leukocyte antigen (HLA) transcript
		- immunoglobin gene transcript
		- T-cell receptor transcript
		- single-exon transcript (will be included in a future version)

	For more information:
	https://www.ensembl.org/Help/Glossary?id=492

+ The absolute numbers of __functional residues__ detected (_firestar_)

+ Score related to the number of exon that map to __protein structure__. Whether we have genomic information or not, we use _Matador3D_ or _Matador3Dv2_.

    _Matador3D_, the score is based on the number of exons that can be mapped to structural homologues.

    _Matador3Dv2_, the number represents the sum of bit-scores in PDB alignment.

+ The __number of vertebrate species__ that have an isoform that aligns to the human isoform over the whole sequence
and without gaps (_CORSAIR_).

    Alignments with the same species scores just 0.5.

    We generate multiple alignments against orthologues from a vertebrate protein database.

	We only align a limited number of vertebrate species, chimp, mouse, rat, dog, cow etc.

+ The absolute numbers of __pfam domains__ that are detected (_SPADE_):

    SPADE identifies the functional domains present in a transcript and detects those domains that are damaged (not whole). The number represents the sum of bitscores in Pfam alignment.

+ The number of __TMH__ detected (_THUMP_).

    The numbers after the '-' indicate the numbers of partial TMH: 'Whole TMH'-'Partial TMH (damaged)'. By partial we could mean "broken" or not whole. Some TMH will be broken by a splicing event, but many TMH are not whole because the TMH domain boundaries do not always describe the domain well.

+ Reliability score for __signal peptides and mitochondrial signal__ sequences (_CRASH_).

    We use a score of 3 or above as a reliable signal peptide, and mitochondrial signal sequences (separated by comma).

+ The number of exons with unusual evolutionary rats (_INERTIA_) *__DEPRECATED!!__*

    INERTIA uses three alignment methods to generate cross-species alignments, from which SLR identifies exons with unusual evolutionary rates.

+ __APPRIS score__

    Reliability score for the variants based on the scores of methods and a weight for them.

+ __No. mapping peptides__

    Proteomic evidence *__only for the human genome (GENCODE gene set)__*.

    This proteomic evidence has been collected from various mass
    spectrometry (MS) sources, covering a range of tissues and cell types.

    Prior to GENCODE v33, peptides were collected from 8 sources
    ([Ezkurdia et al. 2015](https://doi.org/10.1021/pr501286b);
    [Farrah et al. 2013](https://doi.org/10.1021/pr301012j);
    [Ezkurdia et al. 2012](https://doi.org/10.1093/molbev/mss100);
    [Munoz et al. 2011](https://doi.org/10.1038/msb.2011.84);
    [Nagaraj et al. 2011](https://doi.org/10.1038/msb.2011.81);
    [Geiger et al. 2012](https://doi.org/10.1074/mcp.m111.014050);
    [Kim et al. 2014](https://doi.org/10.1038/nature13302);
    [Wilhelm et al. 2014](https://doi.org/10.1038/nature13319)).
    To improve reliability, peptides from each of these sources were filtered,
    eliminating non-tryptic and semi-tryptic peptides and peptides containing
    missed cleavages, and where possible only considering peptides identified
    by multiple search engines.

    APPRIS releases from GENCODE v33 onwards have
    used peptide evidence from 2 proteomics studies
    ([Kim et al. 2014](https://doi.org/10.1038/nature13302);
    [Wang et al. 2019](https://doi.org/10.15252/msb.20188503)),
    and from GENCODE v38, these have been supplemented by peptides from a further 5 studies
    ([Zhang et al. 2014](https://doi.org/10.1038/nature13438);
    [Bekker-Jensen et al. 2017](https://doi.org/10.1016/j.cels.2017.05.009);
    [Carlyle et al. 2017](https://doi.org/10.1038/s41593-017-0011-2);
    [Schiza et al. 2019](https://doi.org/10.1074/mcp.RA118.001170);
    [Jiang et al. 2020](https://doi.org/10.1016/j.cell.2020.08.036)).
    For these studies, the Comet ([Eng et al. 2012](https://doi.org/10.1002/pmic.201200439))
    search engine was used with default parameters, then post-processed with Percolator
    ([Käll et al. 2007](https://doi.org/10.1038/nmeth1113);
    [The et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27572102/)).
    Peptide-spectrum matches (PSMs) that had a Posterior Error Probability (PEP)
    of lower than 0.001 were allowed as long as they were fully tryptic peptides
    and had no more than 2 missed cleavages.

+ __Reliability labels__ of APPRIS

    See [*__Principal Isoforms Flags__*](#principal-isoforms-flags) section, for more information.


GTF data files
==============
Description of the method scores:

+ *__firestar__* are the absolute numbers of functional residues detected.
+ *__Matador3D__* is a score related to the number of exon that map to structure.
+ *__CORSAIR__* shows the number of vertebrate species that have an isoform that aligns to the human isoform over the whole sequence and without gaps (human scores just 0.5).
+ *__SPADE__* shows the absolute numbers of pfam domains that are detected. The numbers after the decimal indicate the
numbers of partial domains.
+ *__THUMP__* shows the number of TMH detected. Again numbers after the decimal indicate partial TMH.
+ *__CRASH__* gives a reliability score for signal peptides (we use a score of 3 or above as a reliable signal peptide). Crash-M does the same thing for mitochondrial signal sequences.

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


BED data files
==============
BED (Browser Extensible Data) format provides a flexible way to define the data lines that are displayed in an
annotation track (in principle for [UCSC Genome Browser](https://genome.ucsc.edu/)). There is one bed file for each APPRIS method, and the data represents the genome regions where the annotations of methods are located.

>For more information:
https://genome.ucsc.edu/FAQ/FAQformat.html#format1


BED12 data files
================
This is an extension of BED format. BED detail uses the first 4 to 12 columns of BED format, plus 2 additional fields that are used to enhance the track details pages. The first additional field is an ID, and the second additional field is a description of the item. For the APPRIS method, we have used the description field to include information as the PDB ligands (firestar), the Pfam domains (SPADE),
and the PDB structure (Matador3D), etc.

>For more information:
https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7


bigBed data files
=================
The bigBed format stores annotation items that can be either a simple or a linked collection of exons, much as BED files do. BigBed files are created from BED12 type files using the program bedToBigBed. The resulting bigBed files are in an indexed binary format. These files have been used for the creation of *__APPRIS Track Hub__* in the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgHubConnect).


>For more information:
https://genome.ucsc.edu/goldenPath/help/bigBed.html
https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html


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
