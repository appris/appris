-- 
-- Number of genes by chromosome that have translation
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
-- Get the number of genes by chromosome that has FIRESTAR result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_firestar
FROM
	entity e1, coordinate c, firestar f
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=f.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

--
-- Get the number of genes by chromosome that has MATADOR3D result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_matador3d
FROM
	entity e1, coordinate c, matador3d m
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=m.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

--
-- Get the number of genes by chromosome that has SPADE result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_spade
FROM
	entity e1, coordinate c, spade s
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of genes by chromosome that has CORSAIR result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_corsair
FROM
	entity e1, coordinate c, corsair s
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of genes by chromosome that has THUMP result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_thump
FROM
	entity e1, coordinate c, thump t
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=t.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

--
-- Get the number of genes by chromosome that has CRASH result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_crash
FROM
	entity e1, coordinate c, crash s
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of genes by chromosome that has INERTIA result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_inertia
FROM
	entity e1, coordinate c, inertia i
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=i.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of genes by chromosome that has PROTEO result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_proteo
FROM
	entity e1, coordinate c, proteo e
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=e.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of genes by chromosome that has APPRIS result
-- 
SELECT
	c.chromosome AS chromosome,
	count(e1.identifier) AS num_genes_appris
FROM
	entity e1, coordinate c, appris a
WHERE
	e1.datasource_id=5 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=a.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;
