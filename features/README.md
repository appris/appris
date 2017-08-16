UNIPROT DATA SOURCE
===================

#Êinfo files

ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README

# Download data files

http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=proteome:UP000005640&force=no&preview=false&format=tab&columns=id,entry%20name,genes(PREFERRED),database(CCDS),length

http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=proteome:UP000005640&fil=&force=no&preview=false&format=fasta&include=yes


REFSEQ DATA SOURCE
==================


# info files

## release
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/README_CURRENT_RELEASE

##Êchromosome info
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/chr_accessions_GRCm38.p4

## Get RefSeq gene to GeneName Symbol
ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene


# Download data files

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/GFF/ref_GRCm38.p4_top_level.gff3.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/RNA/rna.fa.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/protein/protein.fa.gz


# Extract the Cross-Reference with Ensembl (for mouse)
wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz
zgrep "^10090" gene2ensembl > xref.refseq_ensembl.txt

# Add extra values (DEPRECATED)
#perl /local/jmrodriguez/appris/scripts/add_extraVals_into_GTF.pl \
#	-r \
#	--xref=xref.refseq_ensembl.txt \
#	--data=ref_GRCm38.p4_top_level.gff3 \
#	--seq-data=protein.fa \
#	--extra-data=../e86_gM11/gencode.vM11.annotation.gtf \
#	--outfile=ref_GRCm38.p4_top_level.extra.gff3 \
#	--loglevel=info \
#	--logfile=add_extraVals_into_GTF.log
#
#
## Get the ReadThrought tags
#perl /local/jmrodriguez/appris/scripts/add_RT_into_GFF3.pl --gdata=e81_g23/gencode.v23.annotation.gtf --rdata=rs105/ref_GRCh37.p13_top_level.gff3 --xref=rs105/xref.rs107_ensembl.txt --outfile=rs105/ref_
#GRCh37.p13_top_level.RT.gff3 --loglevel=info
#perl /local/jmrodriguez/appris/scripts/add_RT_into_GFF3.pl --gdata=e81_g23/gencode.v23.annotation.gtf --rdata=rs107/ref_GRCh38.p2_top_level.gff3 --xref=rs107/xref.rs107_ensembl.txt --outfile=rs107/ref_G
#RCh38.p2_top_level.RT.gff3 --loglevel=info



