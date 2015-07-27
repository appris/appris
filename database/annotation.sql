-- 
-- Number of transcripts by chromosome that have translation
-- 
SELECT
	count(distinct(e1.identifier)) AS num_trans,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, sequence s
WHERE
	e1.datasource_id=6 AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id AND
	s.type_id=2
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of transcripts by chromosome that has FIRESTAR result
-- 
SELECT
	count(e1.identifier) AS num_trans_firestar,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, firestar f
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=f.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

--
-- Get the number of transcripts by chromosome that has MATADOR3D result
-- 
SELECT
	count(e1.identifier) AS num_trans_matador3d,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, matador3d m
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=m.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

--
-- Get the number of transcripts by chromosome that has SPADE result
-- 
SELECT
	count(e1.identifier) AS num_trans_spade,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, spade s
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of transcripts by chromosome that has CORSAIR result
-- 
SELECT
	count(e1.identifier) AS num_trans_corsair,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, corsair s
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of transcripts by chromosome that has THUMP result
-- 
SELECT
	count(e1.identifier) AS num_trans_thump,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, thump t
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=t.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

--
-- Get the number of transcripts by chromosome that has CRASH result
-- 
SELECT
	count(e1.identifier) AS num_trans_crash,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, crash s
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=s.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of transcripts by chromosome that has INERTIA result
-- 
SELECT
	count(e1.identifier) AS num_trans_inertia,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, inertia i
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=i.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of transcripts by chromosome that has PROTEO result
-- 
SELECT
	count(e1.identifier) AS num_trans_proteo,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, proteo e
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=e.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;

-- 	
-- Get the number of transcripts by chromosome that has APPRIS result
-- 
SELECT
	count(e1.identifier) AS num_trans_appris,
	c.chromosome AS chromosome
FROM
	entity e1, coordinate c, appris a
WHERE
	(e1.datasource_id=2 OR e1.datasource_id=6) AND
	e1.entity_id=c.entity_id AND
	e1.entity_id=a.entity_id
GROUP BY c.chromosome ORDER BY CAST(c.chromosome AS SIGNED) ASC;
