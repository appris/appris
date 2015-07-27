/*
 * Exporter RESTful web services
 *
 * This is a RESTful web services that export APPRIS annotations.
 */

/**
 * @api {get} /position/:specie/:name/:method Retrieves data from a genomic region filtering by methods.
 * @apiDescription Retrieves the APPRIS data from a given genomic region filtering by the given list of methods.
 * @apiName exporterRegionMethod
 * @apiGroup Exporter
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} position	Genomic region (chr:start-end).
 * @apiParam {String} methods 	Methods to execute (appris,firestar,matador3d,spade,corsair,crash,thump).
 * @apiParam {Integer} ens		Ensembl version.
 * @apiParam {String} format 	Format result.
 *
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 *
 * 		no: no head information
 * 		yes:ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} jobid A plain text document containing the current query.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/exporter/position/homo_sapiens/chr22:30773835-30821305/spade?format=json&ens=67
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *		[{"gene_id":"ENSG00000100003","source":"SPADE","frame":".","score":"1","note":"hmm_name:CRAL_TRIO_N,evalue:0.00000000042,pep_start:13,pep_end:66","end":"30803107","seqname":"22","strand":"+","type":"domain_wrong","transcript_id":"ENST00000405717","start":"30793142"},
 *		 {"gene_id":"ENSG00000100003","source":"SPADE","frame":".","score":"14","note":"hmm_name:CRAL_TRIO,evalue:2.5e-38,pep_start:111,pep_end:244","end":"30811815","seqname":"22","strand":"+","type":"domain_wrong","transcript_id":"ENST00000405717","start":"30803500"}]
 *
 * @apiError BadRequest The request is wrong.
 * @apiError JobNotFound The <code>id</code> of the Job was not found.
 * @apiError MethodNotAllowed The <code>method</code> is not allowed.
 *
 * @apiErrorExample Error-Response:
 *     HTTP/1.1 400 Bad Request The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
 *     HTTP/1.1 404 Not Found The server has not found anything matching the Request-URI.
 *     HTTP/1.1 405 Method Not Allowed: The parameter XXX is not allowed. The parameter must be: YYY.
 *
 */

/**
 * @api {get} /position/:specie/:position Retrieves data from a genomic region.
 * @apiDescription Retrieves the APPRIS data from a given genomic region.
 * @apiName exporterRegion
 * @apiGroup Exporter
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} position	Genomic region (chr:start-end).
 * @apiParam {Integer} ens		Ensembl version.
 * @apiParam {String} format 	Format result.
 *
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 *
 * 		no: no head information
 * 		yes:ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} jobid A plain text document containing the current query.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/exporter/position/homo_sapiens/chr22:30773835-30821305?format=gtf
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *		22	SPADE	domain_wrong	30793142	30803107	1	+	.	gene_id "ENSG00000100003"; transcript_id "ENST00000405717"; note "hmm_name:CRAL_TRIO_N,evalue:0.00000000042,pep_start:13,pep_end:66"
 *		22	SPADE	domain_wrong	30803500	30811815	14	+	.	gene_id "ENSG00000100003"; transcript_id "ENST00000405717"; note "hmm_name:CRAL_TRIO,evalue:2.5e-38,pep_start:111,pep_end:244"
 *		22	APPRIS	functional_domain	30793026	30812964	0.2	+	.	gene_id "ENSG00000100003"; transcript_name "SEC14L2-004"; transcript_id "ENST00000405717"; annotation "UNKNOWN"
 *
 * @apiError BadRequest The request is wrong.
 * @apiError JobNotFound The <code>id</code> of the Job was not found.
 * @apiError MethodNotAllowed The <code>method</code> is not allowed.
 *
 * @apiErrorExample Error-Response:
 *     HTTP/1.1 400 Bad Request The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
 *     HTTP/1.1 404 Not Found The server has not found anything matching the Request-URI.
 *     HTTP/1.1 405 Method Not Allowed: The parameter XXX is not allowed. The parameter must be: YYY.
 *
 */

/**
 * @api {get} /name/:specie/:name/:method Retrieves data from a given gene/transcript name filtering by methods.
 * @apiDescription Retrieves the APPRIS data from a given gene/transcript name filtering by the given list of methods.
 * @apiName exporterNameMethod
 * @apiGroup Exporter
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} name 		Gene identifier/name.
 * @apiParam {String} methods 	Methods to execute (appris,firestar,matador3d,spade,corsair,crash,thump).
 * @apiParam {Integer} ens		Ensembl version.
 * @apiParam {String} format 	Format result.
 *
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 *
 * 		no: no head information
 * 		yes:ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} jobid A plain text document containing the current query.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/exporter/name/homo_sapiens/RNF215/firestar,spade?format=bed&ens=70
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *		browser position chr22:30773835-30817760
 *		browser pix 800
 *		browser hide all
 *		track name=Known_functional_residues description='Known functional residues' visibility=2 color='210,145,35' group='0'
 *		chr22	30775715	30775727	ENST00000421022	0	-	30775715	30775727	0	2	3,3	0,9
 *		chr22	30775713	30776084	ENST00000215798	0	-	30775713	30776084	0	8	3,3,3,3,3,3,3,3	0,9,42,51,60,66,359,368
 *		chr22	30775715	30776086	ENST00000382363	0	-	30775715	30776086	0	8	3,3,3,3,3,3,3,3	0,9,42,51,60,66,359,368
 *		track name=Functional_domains description='Whole Pfam functional domains' visibility=2 color='118,156,2' group='0'
 *		chr22	30775713	30776084	ENST00000215798	0	-	30775713	30776084	0	2	89,34	0,337
 *		chr22	30775709	30776086	ENST00000382363	0	-	30775709	30776086	0	2	93,36	0,341
 *
 * @apiError BadRequest The request is wrong.
 * @apiError JobNotFound The <code>id</code> of the Job was not found.
 * @apiError MethodNotAllowed The <code>method</code> is not allowed.
 *
 * @apiErrorExample Error-Response:
 *     HTTP/1.1 400 Bad Request The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
 *     HTTP/1.1 404 Not Found The server has not found anything matching the Request-URI.
 *     HTTP/1.1 405 Method Not Allowed: The parameter XXX is not allowed. The parameter must be: YYY.
 *
 */

/**
 * @api {get} /name/:specie/:name Retrieves data from a given gene/transcript name.
 * @apiDescription Retrieves the APPRIS data from a given gene/transcript name.
 * @apiName exporterName
 * @apiGroup Exporter
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} name 		Gene identifier/name.
 * @apiParam {Integer} ens		Ensembl version.
 * @apiParam {String} format 	Format result.
 *
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 *
 * 		no: no head information
 * 		yes:ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} jobid A plain text document containing the current query.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/exporter/name/homo_sapiens/RNF215?format=tsv&ens=67
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     	== APPRIS > FIRESTAR
 *     	== >sequence_id
 *     	== residue	amino_acid	ligand	reliability_score (1-->6)
 *     	>ENST00000382363
 *     	325	C	ZN	6
 *     	328	C	ZN	6
 *     	343	C	ZN	6
 *     	345	H	ZN	6
 *     	348	H	ZN	4
 *     	351	C	ZN	4
 *     	362	C	ZN	6
 *     	365	C	ZN	6
 *     	== APPRIS > SPADE
 *     	== >sequence_id
 *     	== start	end	domain_name	best_e-value
 *     	>ENST00000382363
 *     	325	367	zf-C3HC4	0.00000025
 *
 * @apiError BadRequest The request is wrong.
 * @apiError JobNotFound The <code>id</code> of the Job was not found.
 * @apiError MethodNotAllowed The <code>method</code> is not allowed.
 *
 * @apiErrorExample Error-Response:
 *     HTTP/1.1 400 Bad Request The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
 *     HTTP/1.1 404 Not Found The server has not found anything matching the Request-URI.
 *     HTTP/1.1 405 Method Not Allowed: The parameter XXX is not allowed. The parameter must be: YYY.
 *
 */

/**
 * @api {get} /id/:specie/:id/:method Retrieves data from a given gene/transcript identifier filtering by methods.
 * @apiDescription Retrieves the APPRIS data from a given Ensembl gene/transcript identifier filtering by the given list of methods.
 * @apiName exporterIdMethod
 * @apiGroup Exporter
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} id 		Gene identifier/name.
 * @apiParam {String} methods 	Methods to execute (appris,firestar,matador3d,spade,corsair,crash,thump).
 * @apiParam {Integer} ens		Ensembl version.
 * @apiParam {String} format 	Format result.
 *
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 *
 * 		no: no head information
 * 		yes:ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} jobid A plain text document containing the current query.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/exporter/id/homo_sapiens/ENST00000382363/firestar,spade?format=tsv
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     	== APPRIS > FIRESTAR
 *     	== >sequence_id
 *     	== residue	amino_acid	ligand	reliability_score (1-->6)
 *     	>ENST00000382363
 *     	325	C	ZN	6
 *     	328	C	ZN	6
 *     	343	C	ZN	6
 *     	345	H	ZN	6
 *     	348	H	ZN	4
 *     	351	C	ZN	4
 *     	362	C	ZN	6
 *     	365	C	ZN	6
 *     	== APPRIS > SPADE
 *     	== >sequence_id
 *     	== start	end	domain_name	best_e-value
 *     	>ENST00000382363
 *     	325	367	zf-C3HC4	0.00000025
 *
 * @apiError BadRequest The request is wrong.
 * @apiError JobNotFound The <code>id</code> of the Job was not found.
 * @apiError MethodNotAllowed The <code>method</code> is not allowed.
 *
 * @apiErrorExample Error-Response:
 *     HTTP/1.1 400 Bad Request The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
 *     HTTP/1.1 404 Not Found The server has not found anything matching the Request-URI.
 *     HTTP/1.1 405 Method Not Allowed: The parameter XXX is not allowed. The parameter must be: YYY.
 *
 */

/**
 * @api {get} /id/:specie/:id Retrieves data from a given gene/transcript identifier.
 * @apiDescription Retrieves the APPRIS data from a given Ensembl gene/transcript identifier.
 * @apiName exporterId
 * @apiGroup Exporter
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} id 		Gene identifier/name.
 * @apiParam {Integer} ens		Ensembl version.
 * @apiParam {String} format 	Format result.
 *
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 *
 * 		no: no head information
 * 		yes:ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} jobid A plain text document containing the current query.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/exporter/id/homo_sapiens/ENST00000382363?format=tsv
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     	== APPRIS > FIRESTAR
 *     	== >sequence_id
 *     	== residue	amino_acid	ligand	reliability_score (1-->6)
 *     	>ENST00000382363
 *     	325	C	ZN	6
 *     	328	C	ZN	6
 *     	343	C	ZN	6
 *     	345	H	ZN	6
 *     	348	H	ZN	4
 *     	351	C	ZN	4
 *     	362	C	ZN	6
 *     	365	C	ZN	6
 *     	== APPRIS > SPADE
 *     	== >sequence_id
 *     	== start	end	domain_name	best_e-value
 *     	>ENST00000382363
 *     	325	367	zf-C3HC4	0.00000025
 *
 * @apiError BadRequest The request is wrong.
 * @apiError JobNotFound The <code>id</code> of the Job was not found.
 * @apiError MethodNotAllowed The <code>method</code> is not allowed.
 *
 * @apiErrorExample Error-Response:
 *     HTTP/1.1 400 Bad Request The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
 *     HTTP/1.1 404 Not Found The server has not found anything matching the Request-URI.
 *     HTTP/1.1 405 Method Not Allowed: The parameter XXX is not allowed. The parameter must be: YYY.
 *
 */