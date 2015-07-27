'use strict';

/* Services */

var apprisServices = angular.module('apprisServices', ['ngResource']);

/* Example query for Run Server */
apprisServices.factory('Examples', ['$resource', function($resource){
    return $resource('ws/examples.json');
}]);

/* CONFIG file of DOWNLOAD files */
apprisServices.factory('Downloads', ['$resource', 'consBaseUrlWS', function($resource, consBaseUrlWS){
    return $resource(consBaseUrlWS+'/downloads.json', {
        query: {
            method:'GET',
            //cache: false,
            isArray: false,
            params: { format:'json' }
        }
    });
}]);

/* CONFIG file of used SPECIES */
apprisServices.factory('Species', ['$resource', 'consBaseUrlWS', function($resource, consBaseUrlWS){
    return $resource(consBaseUrlWS+'/species.json', {
        query: {
            method:'GET',
            //cache: false,
            isArray: false,
            params: { format:'json' }
        }
    });
}]);

/* CONFIG file of used METHODS */
apprisServices.factory('Methods', ['$resource', 'consBaseUrlWS', function($resource, consBaseUrlWS){
    return $resource(consBaseUrlWS+'/methods.json', {
        query: {
            method:'GET',
            //cache: false,
            isArray: false,
            params: { format:'json' }
        }
    });
}]);

/* RUNNER services */
apprisServices.factory('Runner', ['$q', '$resource', 'consUrlRunnerRunWS', function ($q, $resource, consUrlRunnerRunWS) {
//    return $resource(consUrlRunnerRunWS, { species: '@species', methods: '@methods' }, {
//        save: {
//            method:'POST'
//        }
//    });

// TODO!!! PROBLEMS WHEN THE REQUEST-URI IS TOO LONG!!
    var runnerResource = $resource(consUrlRunnerRunWS);
    return {
        /* method: create */
        create: function (data) {
                    var delay = $q.defer();
                    runnerResource.get(data, function (resp) {
                        delay.resolve(resp);
                    }, function(err) {
                        delay.reject(err);
                    });
                    return delay.promise;
        }
    }
}]);

apprisServices.factory('Status', ['$resource', 'consUrlRunnerStatusWS', function($resource, consUrlRunnerStatusWS){
    return $resource(consUrlRunnerStatusWS+'/:jobid', { jobid:'@jobid' });
}]);

//apprisServices.factory('ResultTypes', ['$resource', 'consUrlRunnerResultTypesWS', function($resource, consUrlRunnerResultTypesWS){
//    return $resource(consUrlRunnerResultTypesWS+'/:jobid', { jobid:'@jobid' }, {
//        /* method: query */
//        query: {
//            method:'GET',
//            cache: false,
//            isArray: true,
//            params: { format:'json' }
//        }
//    });
//}]);

apprisServices.factory('ResultTypes', ['$http', '$q', '$filter', 'consUrlRunnerResultTypesWS', function($http, $q, $filter, consUrlRunnerResultTypesWS){
    return({
        query: getResultTypes
    });
    // PUBLIC METHODS
    function getResultTypes( iparams ) {
        var url = getURL(iparams);
        var params = getParams(iparams);
        var request = $http({
            method: "get",
            cache: false,
            url: url,
            params: params
        });
        return request.then(handleSuccess, handleError) ;
    }

    // PRIVATE METHODS
    // get path
    function getURL( iparams ) {
        var url = null;
        if ( angular.isDefined(iparams.jobid) && iparams.jobid !== null ) {
            url = consUrlRunnerResultTypesWS+'/'+iparams.jobid;
        }
        else {
            url = consUrlRunnerResultTypesWS;
        }
        return url;
    }
    // get params
    function getParams( iparams ) {
        var params = {};
        if ( angular.isDefined(iparams.species) ) { params.species = iparams.species }
        if ( angular.isDefined(iparams.format) ) { params.format = iparams.format }
        return params;
    }

    // I transform the error response, unwrapping the application dta from the API response payload.
    function handleError(response) {
        if ( !(angular.isDefined(response.data))  ) {
            return $q.reject( "An unknown error occurred." );
        }
        else {
            if ( response.status === 404 ) {
                return $q.reject( "The annotations were not found for your query. Try to add at least one sequence or one method." );
            }
            else if ( response.status === 405 ) {
                return $q.reject( "The methods were not selected. Try to add some methods." );
            }
        }
        return $q.reject(response.data);
    }

    // I transform the successful response, unwrapping the application data from the API response payload.
    function handleSuccess( response ) {
        if ( angular.isDefined(response.data) ) {
            var data = response.data;
            var isSeqRunner = false;
            var resultDetailAnnotHeads = [];
            var currentMethods = [];

            if ( $filter('isSeqRunner')(data) ) {
                isSeqRunner = true;
            }
            // Add seq. id and methods
            resultDetailAnnotHeads.push({
                id: "transcript_id",
                label: "Seq. id"
            },{
                id: "transcript_name",
                label: "Seq. name"
            },{
                id: "length_aa",
                label: "Length (aa)"
            });
            angular.forEach(data, function(item) {
                // discard methods for detailed panel!!!
                if ( item.id !== 'principal_isoform' ) {
                    resultDetailAnnotHeads.push({
                        id: item.id,
                        label: item.label
                    });
                }
                // discard methods for browser panel !!!
                currentMethods.push({
                    id: item.id,
                    name: item.name,
                    label: item.label,
                    desc: item.desc
                });
            });
            return {
                isSeqRunner: isSeqRunner,
                resultDetailAnnotHeads: resultDetailAnnotHeads,
                currentMethods: currentMethods
            }
        }
        else { return {} }

    }

}]);

//apprisServices.factory('Result', ['$resource', 'consUrlRunnerResultWS', function($resource, consUrlRunnerResultWS){
//    return $resource(consUrlRunnerResultWS+'/:jobid', { jobid:'@jobid', format: '@format'}, {
//        getJSON: {
//            method:'GET',
//            isArray: true,
//            params: { format:'json' }
//        },
//        getGTF: {
//            method:'GET',
//            isArray: false,
//            params: { format:'gtf' }
//        },
//        getTXT: {
//            method:'GET',
//            isArray: false,
//            params: { format:'tsv' },
//            transformResponse: transformTsvData
//        }
//        /* method: query */
////        query: {
////            method:'GET',
////            isArray: true,
////            params: { format:'json' }
////        }
//    });
//}]);

/* EXPORTER services */
//apprisServices.factory('Export', ['$resource', 'consUrlExporterWS', function($resource, consUrlExporterWS){
//    return $resource(consUrlExporterWS+'/:tid/:species/:id', { tid:'@tid', species:'@species', id:'@id', ens: '@ens', format: '@format', methods: '@methods' }, {
//        getJSON: {
//            method:'GET',
//            isArray: true,
//            params: { format:'json' }
//        },
//        getGTF: {
//            method:'GET',
//            isArray: false,
//            params: { format:'gtf' }
//        },
//        getTXT: {
//            method:'GET',
//            isArray: false,
//            params: { format:'tsv' },
//            transformResponse: transformTsvData
//        }
//        /* method: query */
////        query: {
////            method:'GET',
////            isArray: true,
////            params: { format:'json' }
////        }
//    });
//}]);


/* RUNNER/EXPORTER services */
apprisServices.factory('Retriever', ['$http', '$q', 'consUrlRunnerResultWS', 'consUrlExporterWS', function($http, $q, consUrlRunnerResultWS, consUrlExporterWS){
    return({
        getJSON: getJSON,
        getGTF: getTXT,
        getTXT: getTXT
    });

    // PUBLIC METHODS
    function getJSON( iparams ) {
        var url = getURL(iparams);
        var params = getParams(iparams);
        var request = $http({
            method: "get",
            cache: false,
            url: url,
            params: params
        });
//        return request.then(handleSuccess) ;
        return request.then(handleSuccess, handleError) ;
    }

    function getTXT( iparams ) {
        var url = getURL(iparams);
        var params = getParams(iparams);
        var request = $http({
            method: "get",
            cache: false,
            url: url,
            params: params
        });
//        return request.then(handleSuccess) ;
        return request.then(handleSuccess, handleError) ;
    }

    // PRIVATE METHODS
    // get path
    function getURL( iparams ) {
        var url = null;
        if ( angular.isDefined(iparams.jobid) ) {
            url = consUrlRunnerResultWS+'/'+iparams.jobid;
        }
        else {
            url = consUrlExporterWS+'/'+iparams.tid+'/'+iparams.species+'/'+iparams.id;
        }
        return url;
    }
    // get params
    function getParams( iparams ) {
        var params = {};
        if ( angular.isDefined(iparams.operation) ) { params.operation = iparams.operation }
        if ( angular.isDefined(iparams.methods) ) { params.methods = iparams.methods }
        if ( angular.isDefined(iparams.ids) ) { params.ids = iparams.ids }
        if ( angular.isDefined(iparams.ens) ) { params.ens = iparams.ens }
        if ( angular.isDefined(iparams.db) ) { params.db = iparams.db }
        if ( angular.isDefined(iparams.format) ) { params.format = iparams.format }
        return params;
    }
    // I transform the error response, unwrapping the application dta from the API response payload.
    function handleError(err) {
        var response = {
            "status": null,
            "data": null
        };
        if ( !(angular.isDefined(err.data))  ) {
            response.status = 500;
            response.data = "An unknown error occurred.";
        }
        else if ( err.status === 0 ) {
                response.status = 500;
                response.data = "An unknown error occurred.";
        }
        else if ( err.status === 404 ) {
            response.status = 404;
            response.data = "Your query is malformed. Please, rephrase your query.";
        }
        else if ( err.status === 405 ) {
                response.status = 405;
                response.data = "Your query is malformed. Please, rephrase your query.";
        }
        else if ( err.status === 500 ) {
            response.status = 500;
            response.data = "An unknown error occurred.";
        }
        else {
            response.status = 500;
            response.data = "An unknown error occurred.";
        }
        return $q.reject(response);
    }
    // I transform the successful response, unwrapping the application data from the API response payload.
    function handleSuccess( response ) {
        if ( angular.isDefined(response.data) ) {
            if ( (angular.isArray(response.data) || angular.isString(response.data) ) && (response.data.length > 0) ) {
                return response.data;
            }
            else { return "The annotations were not found. Try to add another methods." }
        }
        else { return "Empty response" }
    }

}]);

/* SEEKER services */
apprisServices.factory('Seeker', ['$resource', 'consUrlSeekerWS', function($resource, consUrlSeekerWS){
    return $resource(consUrlSeekerWS+'/:id', { id:'@id', methods:'@methods', format: '@format'}, {
        /* method: query */
        query: {
            method:'GET',
            cache: false,
            isArray: false,
            params: { format:'json' }
        }
    });
}]);

/* SEQUENCER services */
apprisServices.factory('Sequencer', ['$http', '$q', 'consUrlSequencerWS', function($http, $q, consUrlSequencerWS){
    return({
        get: getSequencer
    });
    function getSequencer( params ) {
        var urlParams = '';
        var qParams = {};
        if ( angular.isDefined(params.jobid) ) {
            urlParams = params.jobid;
        }
        else {
            urlParams = params.tid+'/'+params.species+'/'+params.id;
        }
        if ( angular.isDefined(params.operation) ) { qParams.operation = params.operation }
        if ( angular.isDefined(params.methods) && (params.methods != '') ) { qParams.methods = params.methods }
        if ( angular.isDefined(params.ids) ) { qParams.ids = params.ids }
        if ( angular.isDefined(params.ens) ) { qParams.ens = params.ens }
        if ( angular.isDefined(params.db) ) { qParams.db = params.db }
        if ( angular.isDefined(params.format) ) { qParams.format = params.format }
        if ( angular.isDefined(params.r) ) { qParams.r = params.r }
        var request = $http({
            method: "get",
            cache: false,
            url: consUrlSequencerWS+'/'+urlParams,
            params: qParams
        });
//        return request.then(handleSuccess) ;
        return request.then(handleSuccess, handleError) ;

        // PRIVATE METHODS

        // I transform the error response, unwrapping the application dta from the API response payload.
        function handleError(response) {
            if ( !(angular.isDefined(response.data))  ) {
                return $q.reject( "An unknown error occurred." );
            }
            else {
                if ( response.status === 404 ) {
                    return $q.reject( "The annotations were not found for your query. Try to add at least one sequence or one method." );
                }
                else if ( response.status === 405 ) {
                    return $q.reject( "Your query is malformed. Please, rephrase your query. Try to add at least one sequence or one method." );
                }
            }
            return $q.reject(response.data);
        }

        // I transform the successful response, unwrapping the application data from the API response payload.
        function handleSuccess( response ) {
            if ( angular.isDefined(response.data) ) {
                //if ( (angular.isArray(response.data) || angular.isString(response.data) || angular.isObject(response.data) ) && (response.data.length > 0) ) {
                if ( angular.isArray(response.data) || angular.isString(response.data) || angular.isObject(response.data) ) {
                    return response.data;
                }
                else { return "The annotations were not found. Try to add another methods." }
            }
            else { return "Empty response" }
        }
    }
}]);

/* VIEWER services */
apprisServices.factory('Viewer', ['$http', '$q', 'consUrlViewerWS', function($http, $q, consUrlViewerWS){
    return({
        get: getSequencer
    });
    var url = null;
    var params = {};

    function getSequencer( iparams ) {
        url = getURL(iparams);
        params = getParams(iparams);
        var request = $http({
            method: "get",
            cache: false,
            url: url,
            params: params
        });
//        return request.then(handleSuccess) ;
        return request.then(handleSuccess, handleError) ;
    }

    // PRIVATE METHODS
    // get path url
    function getURL( iparams ) {
        var urlParams = '';
        if ( angular.isDefined(iparams.jobid) ) {
            urlParams = consUrlViewerWS+'/'+iparams.jobid;
        }
        else {
            urlParams = consUrlViewerWS+'/'+iparams.tid+'/'+iparams.species+'/'+iparams.id;
        }
        return urlParams;
    }
    // get params
    function getParams( iparams ) {
        var qParams = {};
        if ( angular.isDefined(iparams.operation) ) { qParams.operation = iparams.operation }
        if ( angular.isDefined(iparams.methods) ) { qParams.methods = iparams.methods }
        if ( angular.isDefined(iparams.ids) ) { qParams.ids = iparams.ids }
        if ( angular.isDefined(iparams.ens) ) { qParams.ens = iparams.ens }
        if ( angular.isDefined(iparams.db) ) { qParams.db = iparams.db }
        if ( angular.isDefined(iparams.format) ) { qParams.format = iparams.format }
        return qParams;
    }
    // I transform the error response, unwrapping the application dta from the API response payload.
    function handleError(response) {
        if ( !(angular.isDefined(response.data))  ) {
            return $q.reject( "An unknown error occurred." );
        }
        else {
            if ( response.status === 404 ) {
                return $q.reject( "The annotations were not found for your query. Try to add at least one sequence or one method." );
            }
            else if ( response.status === 405 ) {
                return $q.reject( "Your query is malformed. Please, rephrase your query. Try to add at least one sequence or one method." );
            }
        }
        return $q.reject(response.data);
    }
    // I transform the successful response, unwrapping the application data from the API response payload.
    function handleSuccess( response ) {
        if ( angular.isDefined(response.data) ) {
            if ( (angular.isArray(response.data) || angular.isString(response.data) || angular.isObject(response.data) ) && (response.data.length > 0) ) {
                return response.data;
            }
            else { return "The annotations were not found. Try to add another methods." }
        }
        else { return "Empty response" }
    }

}]);

/**
 * COMMON functions
 **/
// Transform the response from Exporter services
function transformExportData(data) {
    if ( angular.isDefined(data) ) {
        var dataJson = angular.fromJson(data);
        if ( angular.isArray(dataJson) ) {
            return dataJson;
        }
    }
}
// Transform the TSV Data
function transformTsvData(data) {
    var filtered = [];
    var arrayData = data.split(/^>/mg);
    if ( angular.isArray(arrayData) && (arrayData.length > 0) ) {
        angular.forEach(arrayData, function(item) {
            if ( angular.isDefined(item) && (item != '') ) {
                var arrayIdMethods = item.split(/^([^\n]*)\n*/);
                if ( angular.isArray(arrayIdMethods) && (arrayIdMethods.length >= 2 ) ) {
                    var seqId = arrayIdMethods[1];
                    var methods = [];
                    var arrayMethods = arrayIdMethods[2];
                    var arrayMets = arrayMethods.split(/^==\s*/mg);
                    if ( angular.isArray(arrayMets) && (arrayMets.length > 0) ) {
                        angular.forEach(arrayMets, function(item) {
                            if ( angular.isDefined(item) && (item != '') ) {
                                var mets = item.split(/^\s*([^:]*):\s*/);
                                if ( angular.isArray(mets) && (mets.length >= 3 ) ) {
                                    var metId = mets[1];
                                    var metBody = mets[2];
                                    var residues = [];
                                    var resRows = metBody.split(/\n/);
                                    for (var i = 1; i < resRows.length; i++) {
                                        if ( angular.isDefined(resRows[i]) ) {
                                            var resTitle = resRows[0].split(/\t/);
                                            var resCols = resRows[i].split(/\t/);
                                            var residue = {};
                                            for (var j = 0; j < resTitle.length; j++) {
                                                residue[ resTitle[j] ] = resCols[j];
                                            }
                                            if ( angular.isDefined(residue) ) {
                                                residues.push(residue);
                                            }
                                        }
                                    }
                                    if ( residues.length > 0 ) {
                                        methods.push({
                                            'id':       metId,
                                            'body':     metBody,
                                            'title':    metBody,
                                            'residues': residues
                                        });
                                    }
                                }
                            }
                        });
                    }
                    if ( methods.length > 0 ) {
                        filtered.push({
                            'id': seqId,
                            'methods': methods
                        });
                    }
                }
            }
        });
    }
    return filtered;
}