/**
 *
 * Report - AngularJS module for visualization of gene report
 *
 */
var module = angular.module('appris.report', []);

/* CONTROLLERS */

module.controller('ReportController', ['consPageError', '$rootScope', '$scope', '$routeParams', '$location', '$filter', 'Methods', 'ResultTypes', 'Retriever', 'Sequencer', function (consPageError, $rootScope, $scope, $routeParams, $location, $filter, Methods, ResultTypes, Retriever, Sequencer) {

    // init vars
    $rootScope.isLoadingScreen = true;
    $scope.alert = {};

    $scope.runnerid = null;
    $scope.exporterid = null;
    $scope.jobid = null;
    $scope.isSeqRunner = false;
    $scope.assembly = null;
    $scope.assemblies = [];

    $scope.currentMethods = [];
    $scope.resultDetailAnnotHeads = [];
    $scope.resultMainAnnotHeads = [];
    $scope.resultMainResidues = [];

    // get route parameters from type of mode
    // EXPORTER mode
    if ( angular.isDefined($routeParams.tid) && angular.isDefined($routeParams.species) && angular.isDefined($routeParams.id) ) {
        $scope.exporterid = $routeParams.tid+'/'+$routeParams.species+'/'+$routeParams.id;
        $scope.jobid = $scope.exporterid;
    }
    // RUNNER mode (idRUNNER and seqRUNNER)
    else if ( angular.isDefined($routeParams.jobid) ) {
        $scope.runnerid = $routeParams.jobid;
        $scope.jobid = $scope.runnerid;
    }

    // APPRIS annots

    // load result types (method results). BY DEFAULT: retrieves all
    ResultTypes.query({ jobid: $scope.runnerid, species: $routeParams.species }).then( function(data) {
        $scope.isSeqRunner = data.isSeqRunner;
        $scope.resultDetailAnnotHeads = data.resultDetailAnnotHeads;
        $scope.currentMethods = data.currentMethods;
    });

    // create annotation table depending on mode (EXPORTER or RUNNER) at the moment we have already the resultHeadLbs
    $scope.$watchCollection('[resultDetailAnnotHeads, currentMethods]', function(newValues) { // For AngularJS 1.1.4
        if ( angular.isDefined(newValues[0]) && angular.isDefined(newValues[1]) ) {
            if ( angular.isArray(newValues[0]) && newValues[0].length > 0 && angular.isArray(newValues[1]) && newValues[1].length > 0 ) {
            var queryData = null;
            if ( $scope.exporterid ) {
                queryData = {
                    tid: $routeParams.tid,
                    species: $routeParams.species,
                    id: $routeParams.id,
                    as: $routeParams.as,
                    sc: $routeParams.sc,
                    methods: 'appris'
                };
            }
            else if ( $scope.runnerid ) {
                queryData = {
                    jobid: $scope.runnerid,
                    methods: 'appris'
                };
            }

            var rawResults = Retriever.getJSON(queryData);
            rawResults
                .then( function(data) {
                    if ( $scope.exporterid ) {
                        var resultData = $filter('convertTransScoreObj')(data);
                        $scope.resultMainAnnotHeads.push({
                            id: "transcript_id",
                            label: "Seq. id"
                        },{
                            id: "transcript_name",
                            label: "Seq. name"
                        },{
                            id: "length_aa",
                            label: "Length (aa)"
                        },{
                            id: "biotype",
                            label: "Biotype"
                        },{
                            id: "ccds_id",
                            label: "CCDS"
                        },{
                            id: "flags",
                            label: "Flags"
                        },{
                            id: "principal_isoform",
                            label: "Principal Isoform"
                        });
                    }
                    else if ( $scope.runnerid ) {
                        if ( $scope.isSeqRunner ) {
                            var resultData = $filter('convertSeqRunnerObj')(data);
                            $scope.resultMainAnnotHeads.push({
                                id: "transcript_id",
                                label: "Seq. id"
                            },{
                                id: "principal_isoform",
                                label: "Principal Isoform"
                            });
                        }
                        else {
                            var resultData = $filter('convertTransScoreObj')(data);
                            $scope.resultMainAnnotHeads.push({
                                id: "transcript_id",
                                label: "Seq. id"
                            },{
                                id: "transcript_name",
                                label: "Seq. name"
                                },{
                                    id: "length_aa",
                                    label: "Length (aa)"
                            },{
                                id: "principal_isoform",
                                label: "Principal Isoform"
                            });
                        }
                    }
                    $rootScope.isLoadingScreen = false;
                    if ( angular.isUndefined(resultData) || (angular.isArray(resultData) && resultData.length == 0 ) ) {
                        $scope.alert.enable = true;
                        $scope.alert.type = "Warning";
                        $scope.alert.message = "Your query were not found. Refreshing page helps to fix the problem (clearing the cache). If problem persistent, Please, contact the administrators";
                    }
                    else {
                        $scope.resultAnnots = resultData[1];
                        $scope.browserDataMethods = angular.copy($scope.currentMethods);
                        $scope.addBrowserDatalSeqs();
                    }
                }, function(error) {
                    $rootScope.isLoadingScreen = false;
                    if ( error.status >= 400 && error.status < 500 ) {
                        $scope.alert.enable = true;
                        $scope.alert.type = "Error";
                        $scope.alert.message = error.data;
                    }
                    else if ( error.status >= 500 ) {
                        $location.url(consPageError);
                    }
                });
        }
    }
    });

    $scope.addBrowserDatalSeqs = function () {
        if ( $scope.browserDataSeqs == null ) {
            $scope.browserDataSeqs = true;
        }
        else if ( $scope.browserDataSeqs ) {
            $scope.browserDataSeqs = false;
        }
        else {
            $scope.browserDataSeqs = true;
        }
        angular.forEach($scope.resultAnnots, function (item) {
            item.selected = $scope.browserDataSeqs;
        });
    };

    // Export data
    $scope.saveAs = [{
        type: "annotation",
        name: "Annotations as",
        format: "json",
        mimetype:   "application/json"
    },{
        type: "annotation",
        name: "Annotations as",
        format: "gtf",
        mimetype:   "text/plain"
    },{
        type: "annotation",
        name: "Annotations as",
        format: "tsv",
        mimetype:   "text/plain"
    },{
        type: "sequence",
        name: "Sequence annot. as",
        format: "json",
        mimetype:   "application/json"
    },{
        type: "sequence",
        name: "Protein seq. as",
        format: "fasta",
        mimetype:   "text/plain"
    },{
        type: "genome",
        name: "Genome annot. as",
        format: "bed",
        mimetype:   "text/plain"
    }];
    $scope.exportData = function(save) {
        // init vars
        var inputs = $scope.resultAnnots;
        var methods = $scope.browserDataMethods;
        var rawResults = null;
        var browserID = save.type;
        var format = save.format;
        // create data (path)
        var browserData = null;
        if ( $scope.exporterid ) {
            browserData = {
                tid: $routeParams.tid,
                species: $routeParams.species,
                id: $routeParams.id,
                as: $routeParams.as,
                sc: $routeParams.sc
            };
        }
        else if ( $scope.runnerid ) {
            browserData = {
                jobid: $scope.runnerid
            };
        }
        // create data (parameters)
        var ids = '';
        angular.forEach($filter('orderBy')(inputs, 'transcript_id'), function(item) {
            if ( item.selected ) {
                ids += item.transcript_id+',';
            }
        });
        ids = ids.replace(/,$/g,'');
        var mets = '';
        angular.forEach(methods, function(item) {
            mets += item.name+',';
        });
        mets = mets.replace(/,$/g,'');
        browserData.ids = ids;
        browserData.format = format;
        browserData.methods = mets;
        // call services
        if ( browserID == 'annotation' ) {
            if ( format === 'json' ) {
                rawResults = Retriever.getJSON(browserData);
            }
            else if ( format === 'gtf' ) {
                rawResults = Retriever.getGTF(browserData);
            }
            else if ( format === 'tsv' ) {
                rawResults = Retriever.getTXT(browserData);
            }
            $scope.saveName = $routeParams.id + '.' + browserID;
        }
        else if ( (browserID == 'sequence') ) {
            if ( format === 'json' ) {
                browserData.operation = 'residues';
                rawResults = Sequencer.get(browserData);
            }
            else if ( format === 'fasta' ) {
                browserData.operation = 'sequences';
                rawResults = Sequencer.get(browserData);
            }
            $scope.saveName = $routeParams.id + '.' + browserID;
        }
        else if ( (browserID == 'genome') ) {
            if ( format === 'bed' ) {
                rawResults = Retriever.getTXT(browserData);
            }
            $scope.saveName = $routeParams.id + '.' + browserID;
        }
        return rawResults;
    };
}
]);


/* FILTERS */

// Checks if the result is "SeqRUNNER mode"
apprisFilters.filter('isSeqRunner', function() {
    return function(input) {
        var isSeqResult = null;
        if ( angular.isDefined(input) && input.length > 0 ) {
            isSeqResult = false;
            var data = input[0];
            if ( data.mode == "sequence" ) {
                isSeqResult = true;
            }
        }
        return isSeqResult;
    };
});

// Creates report from results with genome input
apprisFilters.filter('convertTransScoreObj', function() {
    return function(input) {
        var filtered = {};
        var idList = [];
        var selected = [];
        angular.forEach(input, function(item) {
            if ( angular.isDefined(item.annotation) ) {
                var iTrans = item.transcript_id;
                var sLabel = item.type;
                var sScore = item.score;
                if ( !(filtered[iTrans]) ) {
                    filtered[iTrans] = {};
                    filtered[iTrans]['transcript_id'] = iTrans;
                    filtered[iTrans]['flags'] = '';
                    if ( angular.isDefined(item.transcript_name) ) {
                        filtered[iTrans]['transcript_name'] = item.transcript_name;
                    }
                    if ( angular.isDefined(item.ccds_id) ) {
                        filtered[iTrans]['ccds_id'] = item.ccds_id;
                    }
                    if ( angular.isDefined(item.length_aa) ) {
                        filtered[iTrans]['length_aa'] = item.length_aa;
                    }
                    if ( angular.isDefined(item.biotype) ) {
                        filtered[iTrans]['biotype'] = item.biotype;
                    }
                    if ( angular.isDefined(item.no_codons) ) {
                        filtered[iTrans]['no_codons'] = item.no_codons;
                    }
                    if ( angular.isDefined(item.tsl) ) {
                        filtered[iTrans]['flags'] += ', '+'TSL'+item.tsl;
                    }
                    if ( angular.isDefined(item.no_codons) ) {
                        filtered[iTrans]['flags'] += ', '+item.no_codons+'_codon_NF';
                    }
                    if ( angular.isDefined(item.tag) ) {
                        if ( item.tag.indexOf('readthrough_transcript') > -1 ) {
                            filtered[iTrans]['flags'] += ', '+'ReadThrough';
                        }
                    }
                    filtered[iTrans]['flags'] = filtered[iTrans]['flags'].replace(/^,\s*/g,'');
                    var id = {}
                    id['id'] = iTrans;
                    if ( angular.isDefined(item.length_aa) ) {
                        id['length_aa'] = item.length_aa;
                    }
                    idList.push(id);
                }
                if ( sLabel == "principal_isoform") {
                    var sAnnot = '-';
                    if ( angular.isDefined(item.reliability) ) {
                        sAnnot = item.reliability;
                    }
                    filtered[iTrans][sLabel] = sAnnot;
                }
                else if ( sLabel == "peptide_signal" || sLabel == "mitochondrial_signal" ) {
                    filtered[iTrans][sLabel] = sAnnot;
                }
                else {
                    filtered[iTrans][sLabel] = sScore;
                }
            }
        });

        Object.keys(filtered).forEach(function (key) {
            var val = filtered[key];
            selected.push(val);
        });
        return [idList, selected];
    };
});

// Creates report from results with seq input (SeqRUNNER mode)
apprisFilters.filter('convertSeqRunnerObj', function() {
    return function(input) {
        var idList = [];
        var selected = [];
        if ( angular.isDefined(input) && input.length > 0 ) {
            var data = input[0];
            if ( data.appris && angular.isObject(data.appris) ) {
                var dat = data.appris;
                angular.forEach(dat, function(item, key) {
                    if ( angular.isObject(item) ) {
                        var filtered = {};
                        var id = {}
                        id['id'] = key;
                        idList.push(id);
                        filtered['transcript_id'] = key;
                        if ( angular.isDefined(item.principal_isoform_signal) ) {
                            var sAnnot = '-';
                            if ( angular.isDefined(item.reliability) ) {
                                sAnnot = item.reliability;
                            }
                            filtered['principal_isoform'] = sAnnot;
                        }
                        if ( angular.isDefined(item.functional_residues_score) ) {
                            filtered['functional_residue'] = item.functional_residues_score;
                        }
                        if ( angular.isDefined(item.homologous_structure_score) ) {
                            filtered['homologous_structure'] = item.homologous_structure_score;
                        }
                        if ( angular.isDefined(item.vertebrate_conservation_score) ) {
                            filtered['vertebrate_conservation'] = item.vertebrate_conservation_score;
                        }
                        if ( angular.isDefined(item.domain_score) ) {
                            filtered['functional_domain'] = item.domain_score;
                        }
                        if ( angular.isDefined(item.transmembrane_helices_score) ) {
                            filtered['transmembrane_signal'] = item.transmembrane_helices_score;
                        }
                        if ( angular.isDefined(item.peptide_signal) ) {
                            filtered['signal_peptide_mitochondrial'] = item.peptide_signal;
                        }
                        if ( angular.isDefined(item.mitochondrial_signal) ) {
                            filtered['signal_peptide_mitochondrial'] += '/'+item.mitochondrial_signal;
                        }
                        selected.push(filtered);
                    }
                });
            }
        }
        return [idList, selected];
    };
});

// Replaces the appris labels for "css class"
apprisFilters.filter('activeAnnotClass', function(principal1, principal2, principal3, principal4, principal5, alternative1, alternative2) {
    return function(input, type){
        var filtered = '';
        if ( angular.isDefined(input) ) {
            if ( (input == principal1) ||  (input == principal2) || (input == principal3) || (input == principal4) || (input == principal5) ) {
                if ( angular.isDefined(type) && type == 'label' ) { filtered = "success" }
                else { filtered = "principal" }
            }
            else if ( (input == alternative1) || (input == alternative2) ) {
                if ( angular.isDefined(type) && type == 'label' ) { filtered = "warning" }
                else { filtered = "candidate" }
            }
        }
        return filtered;
    };
});