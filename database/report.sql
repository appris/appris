-- 	
-- Count the number of genes
-- 	
SELECT
	count(*) AS num_genes
FROM
	entity e
WHERE
	datasource_id=5;

-- 	
-- Count the number of genes coming FROM HAVANA
-- 	
SELECT
	count(*) AS num_genes_havana
FROM
	entity e
WHERE
	datasource_id=5 AND
	source='HAVANA';

-- 	
-- Count the number of genes coming FROM ENSEMBL
-- 	
SELECT
	count(*) AS num_genes_ensembl
FROM
	entity e
WHERE
	datasource_id=5 AND 
	source='ENSEMBL';

-- 	
-- Count the number of transcripts
-- 	
SELECT
	count(*) AS num_transcripts
FROM
	entity e
WHERE
	datasource_id=6;

-- 	
-- Count the number of transcripts coming FROM HAVANA
-- 	
SELECT
	count(*) AS num_transcripts_havana
FROM
	entity e
WHERE
	datasource_id=6 AND
	source='HAVANA';

-- 	
-- Count the number of transcripts coming FROM ENSEMBL
-- 	
SELECT
	count(*) AS num_transcripts_ensembl
FROM
	entity e
WHERE
	datasource_id=6 AND
	source='ENSEMBL';


-- 	
-- Count transcripts that at least it has got sequence
-- 	
SELECT
	count(distinct(e1.identifier)) AS num_transcripts_with_sequence
FROM
	entity e1, sequence s
WHERE 	
	e1.datasource_id=6 AND 
	e1.entity_id=s.entity_id AND
	s.type_id=1;

-- 	
-- Count transcripts that at least it has got translation
-- 
SELECT
	count(distinct(e1.identifier)) AS num_transcripts_with_translation
FROM
	entity e1, sequence s
WHERE 	
	e1.datasource_id=6 AND 
	e1.entity_id=s.entity_id AND
	s.type_id=2;

-- 	
-- Count exons
-- 
SELECT
	count(entity_id) AS num_exons
FROM
	exon;

-- 	
-- Count CDS
-- 
SELECT
	count(entity_id) AS num_cds
FROM
	cds;
	
-- 
-- Number of genes by chromosome
-- 
SELECT
	c.chromosome AS chromosome,
	count(distinct(e1.identifier)) AS num_genes
FROM
	entity e1, coordinate c
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 
-- Number of transcripts by chromosome
-- 
SELECT
	c.chromosome AS chromosome,
	count(distinct(e1.identifier)) AS num_transcripts
FROM
	entity e1, coordinate c
WHERE
	e1.datasource_id=6 AND
	e1.entity_id=c.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 
-- Number of genes by chromosome that at least one transcript has sequence 
-- 
SELECT
	c.chromosome AS chromosome,
	count(distinct(e1.identifier)) AS num_genes_with_sequence
FROM
	entity e1, coordinate c, xref_identify x, entity e2, sequence s
WHERE
	e1.entity_id=x.entity_id AND
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	x.identifier=e2.identifier AND
	x.datasource_id=6 AND	
	e2.entity_id=s.entity_id AND
	s.type_id=1
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 
-- Number of transcripts by chromosome that have sequence 
-- 
SELECT
	c.chromosome AS chromosome,
	count(distinct(e1.identifier)) AS num_transcripts_with_sequence
FROM
	entity e1, coordinate c, sequence s
WHERE
	e1.datasource_id=6 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id AND
	s.type_id=1
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;


-- 
-- Number of genes by chromosome that at least one transcript has translation 
-- 
SELECT
	c.chromosome AS chromosome,
	count(distinct(e1.identifier)) AS num_genes_with_translation
FROM
	entity e1, coordinate c, xref_identify x, entity e2, sequence s
WHERE
	e1.entity_id=x.entity_id AND
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	x.identifier=e2.identifier AND
	x.datasource_id=6 AND	
	e2.entity_id=s.entity_id AND
	s.type_id=2
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 
-- Number of transcripts by chromosome that have translation
-- 
SELECT
	c.chromosome AS chromosome,
	count(distinct(e1.identifier)) AS num_transcripts_with_translation
FROM
	entity e1, coordinate c, sequence s
WHERE
	e1.datasource_id=6 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id AND
	s.type_id=2
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 
-- Number of genes by biotype
-- 
SELECT
	e1.biotype AS biotype,
	count(distinct(e1.identifier)) AS num_genes 
FROM
	entity e1
WHERE
	e1.datasource_id=5
GROUP BY e1.biotype;

-- 
-- Number of transcripts by biotype
-- 
SELECT
	e1.biotype AS biotype,
	count(distinct(e1.identifier)) AS num_transcripts 
FROM
	entity e1
WHERE
	e1.datasource_id=6
GROUP BY e1.biotype;

-- 
-- Number of transcripts by biotype that have sequence 
-- 
SELECT
	e1.biotype AS biotype,
	count(distinct(e1.identifier)) AS num_transcripts_with_sequence 
FROM
	entity e1, sequence s
WHERE
	e1.datasource_id=6 AND 
	e1.entity_id=s.entity_id AND
	s.type_id=1
GROUP BY e1.biotype;

-- 
-- Number of transcripts by biotype that have translation
-- 
SELECT
	e1.biotype AS biotype,
	count(distinct(e1.identifier)) AS num_transcripts_with_translation 
FROM
	entity e1, sequence s
WHERE
	e1.datasource_id=6 AND 
	e1.entity_id=s.entity_id AND
	s.type_id=2
GROUP BY e1.biotype;