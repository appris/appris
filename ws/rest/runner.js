/*
 * Runner RESTful web services
 *
 * This is a RESTful web services that execute APPRIS pipeline.
 */

/**
 * @api {post} /run run > Submit a job with the specified parameters
 * @apiName run
 * @apiGroup Runner
 * @apiVersion 0.1.0
 *
 * @apiParam {String} species 	Species name (homo_sapiens, hsap, mus_musculus, mmus, danio_rerio, drer).
 * @apiParam {String} sequences	Fasta protein sequence.
 * @apiParam {String} id 		Gene identifier/name.
 * @apiParam {String} e_version	Ensembl version.
 * @apiParam {String} methods 	Methods to execute (appris,firestar,matador3d,spade,corsair,crash,thump).
 * @apiParam {String} outfile 	Output file name.
 * @apiParam {String} format 	Output format.
 * 
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format 
 *
 * @apiSuccess {String} jobid A plain text document containing the job identifier.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/runner/run
 *	POST
 *	{
 *	'species' => 'Homo sapiens',
 *	'id' => 'ENSG00000099999',
 *	'e_version' => '70',
 *	'methods' => 'firestar,matador3d,corsair,spade,thump,crash,appris',
  *	'format' => 'raw',
 *	}; 
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     c5fcac3d75a7e6d56b753c2265a73ad1
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
 * @api {get} /status/:jobid status > Get the status of a submitted job.
 * @apiName status
 * @apiGroup Runner
 * @apiVersion 0.1.0
 *
 * @apiParam {String} jobid 	Job identifier.
 *
 * @apiSuccess {String} status A plain text document containing the job status. The values for the status are:
 * 
 * RUNNING: the job is currently being processed.
 * FINISHED: job has finished, and the results can then be retrieved.
 * ERROR: an error occurred attempting to get the job status.
 * FAILURE: the job failed.
 * NOT_FOUND: the job cannot be found.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/runner/status/2eb1bb92f32a0384b1b417cb7c9a02ff
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     	RUNNING
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
 * @api {get} /resulttypes/:jobid resulttypes > Get available result types for a finished job.
 * @apiName resulttypes
 * @apiGroup Runner
 * @apiVersion 0.1.0
 *
 * @apiParam {String} jobid 	Job identifier.
 *
 * @apiSuccess {JSON} resulttypes 	An JSON detailing the available result types.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/runner/resulttypes/2eb1bb92f32a0384b1b417cb7c9a02ff
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     	{
 *     	'identifier' => 'appris',
 *     	'fileSuffix' => 'raw',
 *     	'mediaType' => 'text/plain',
 *     	'description' => ''
 *     	};
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
 * @api {get} /result/:jobid result > Get the job result of the specified type.
 * @apiName result
 * @apiGroup Runner
 * @apiVersion 0.1.0
 *
 * @apiParam {String} jobid 		Job identifier.
 * @apiParam {String} format 		Format result.
 * 
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 * 
 * 		no: no head information
 * 		yes: ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} result 	A document containing the result in the requested format. The MIME type of the returned document is set according to the format.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/runner/result/c5fcac3d75a7e6d56b753c2265a73ad1?format=tsv
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *     == APPRIS > FIRESTAR
 *     == >sequence_id
 *     == residue	amino_acid	ligand	reliability_score (1-->6)
 *     >ENST00000405930
 *     106	C	ZN	4
 *     109	C	ZN	4
 *     120	C	ZN,ZN	6,6
 *     123	C	ZN,ZN	6,4
 *     132	H	ZN	6
 *     134	C	ZN	6
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
 * @api {get} /result/:jobid/:method result > Get the job result of the specified type filtering by a given method.
 * @apiName result
 * @apiGroup Runner
 * @apiVersion 0.1.0
 *
 * @apiParam {String} jobid 		Job identifier.
 * @apiParam {String} methods 	Methods to execute (appris,firestar,matador3d,spade,corsair,crash,thump).
 * @apiParam {String} format 		Format result.
 * 
 * 		raw: result of every method in text plain format
 * 		json: completed result of every method in json format
 *		tsv: highlight results of every method in tabular
 *		gtf: genomic data of every method in gtf format
 * 		bed: genomic data of every method in bed format
 * @apiParam {String} headbed 		Add head info into bed file. The values are:
 * 
 * 		no: no head information
 * 		yes: ensGene,ccdsGene,burgeRnaSeqGemMapperAlign add head tracks to bed file
 *		only: print only the head of bed format
 *
 * @apiSuccess {String} result 	A document containing the result in the requested format. The MIME type of the returned document is set according to the format.
 *
 * @apiExample Example usage:
 * http://apprisws.bioinfo.cnio.es/ws/rest/runner/result/d2b90a9c5a7bd1b0272c930eccabc3f2/corsair?format=tsv
 *
 * @apiSuccessExample Success-Response:
 *     HTTP/1.1 200 OK
 *		== APPRIS > CORSAIR
 *		== >sequence_id
 *		== nearest_homologue	%ID
 *		>ENST00000405930
 *		Homo sapiens	100.00
 *		>ENST00000436518
 *		>ENST00000320602
 *		>ENST00000334554
 *		Homo sapiens	100.00
 *		Pan troglodytes	96.49
 *		Sus scrofa	94.12
 *		Canis lupus familiaris	93.86
 *		Rattus norvegicus	91.37
 *		Mus musculus	91.37
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