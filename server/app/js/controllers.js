'use strict';

/* Controllers */

var apprisControllers = angular.module('apprisControllers', []);

/* APPRIS (init) Controllers */
apprisControllers.controller('ApprisController', ['$rootScope', 'Species', 'Methods',
    function($rootScope, Species, Methods) {
//        $timeout(function() { alert.expired = true; }, 2000);

        // load list of species
        if ( angular.isUndefined($rootScope.species) ) {
            $rootScope.species = null;
            Species.get().$promise.then( function(data) {
                $rootScope.species = data;
            });
        }

        // load list of methods
        if ( angular.isUndefined($rootScope.methods) || (angular.isArray($rootScope.methods) && $rootScope.methods.length == 0 ) ) {
            $rootScope.methods = null;
            Methods.get().$promise.then( function(data) {
                $rootScope.methods = data;
            });
        }
    }
]);

/* The ROOT CONTROLLER */
apprisControllers.controller('RootController', ['$rootScope', 'Species', 'Methods',
    function($rootScope, Species, Methods) {
        // loading screen
        //$rootScope.isLoadingScreen = false;

        // load list of species
        Species.get().$promise.then( function(data) {
            $rootScope.species = data;
        });

        // load list of methods
        Methods.get().$promise.then( function(data) {
            $rootScope.methods = data;
        });
    }
]);

/* SPECIES Controllers */
apprisControllers.controller('SpeciesController', ['$rootScope', 'Species',
    function($rootScope, Species) {
        // load list of species
        if ( angular.isUndefined($rootScope.species) ) {
            $rootScope.species = null;
            Species.get().$promise.then( function(data) {
                $rootScope.species = data;
            });
        }
    }
]);

///* RUNNER Controllers */
//apprisControllers.controller('RunnerController', ['consPathServerStatus', '$rootScope', '$scope', '$location', '$filter', '$cacheFactory', 'Runner', 'Species', 'Methods', 'Examples',
//    function(consPathServerStatus, $rootScope, $scope, $location, $filter, $cacheFactory, Runner, Species, Methods, Examples) {
//
//        // init form var
//        $rootScope.isLoadingScreen = false;
//        $scope.runnerForm = {};
//        $scope.runnerForm.species = {};
//        $scope.runnerForm.methods = [];
//        $scope.runnerForm.sequences = null;
//        $scope.runnerForm.ensembl = null;
//        $scope.runnerForm.gene = null;
//
////        // load list of species
////        if ( angular.isUndefined($rootScope.species) ) {
////            $rootScope.species = null;
////            Species.get().$promise.then( function(data) {
////                $rootScope.species = data;
////            });
////        }
////
////        // load list of methods
////        if ( angular.isUndefined($rootScope.methods) || (angular.isArray($rootScope.methods) && $rootScope.methods.length == 0 ) ) {
////            $rootScope.methods = null;
////            Methods.get().$promise.then( function(data) {
////                $rootScope.methods = data;
////            });
////        }
//
//        // load examples
//        if ( angular.isUndefined($scope.examples) ) {
//            $scope.examples = Examples.query();
//        }
//
//        // submit runnerForm: run a Job
//        $scope.runJob = function (form) {
//
//            // Validate form
//
//            // init flag values
//            $scope.speInvalid = false;
//            $scope.genInvalid = false;
//            $scope.seqInvalid = false;
//            $scope.dupInvalid = false;
//            $scope.metInvalid = false;
//            $scope.formInvalid = false;
//
//            // empty specie input
//            if ( angular.isUndefined(form.species) ) {
//                $scope.speInvalid = true;
//            }
//            // empty data input
//            if ( angular.isUndefined(form.ensembl) && angular.isUndefined(form.gene) && angular.isUndefined(form.sequences) ) {
//                $scope.inpInvalid = true;
//            }
//            else if ( ( angular.isUndefined(form.ensembl) || angular.isUndefined(form.gene)) && angular.isUndefined(form.sequences) ) {
//                $scope.genInvalid = true;
//            }
//            else if ( angular.isUndefined(form.sequences) && angular.isUndefined(form.ensembl) && angular.isUndefined(form.gene) ) {
//                $scope.seqInvalid = true;
//            }
//            // duplicate data input
//            if ( angular.isObject(form.sequences) && angular.isObject(form.gene) ) {
//                $scope.dupInvalid = true;
//            }
//            // empty data met
//            if ( angular.isObject(form.methods) && !( $filter('isSelectedMethod')(form.methods) ) ) {
//                $scope.metInvalid = true;
//            }
//
//            // messages
//            if ( $scope.speInvalid || $scope.seqInvalid || $scope.genInvalid || $scope.dupInvalid || $scope.metInvalid ) {
//                $scope.formInvalid = true;
//                return;
//            }
//
//            // create postData
//            var postData = {};
//            if (    angular.isDefined(form.species) && angular.isObject(form.species) && form.species &&
//                    angular.isDefined(form.sequences) && form.sequences &&
//                    $filter('isSelectedMethod')(form.methods) )
//            {
//                postData.species = form.species.scientific;
//                postData.sequences = form.sequences;
//                //postData.methods = Object.keys(form.methods).join(',')+',appris' // TODO: Improve this 'hard-code shit'
//                postData.methods = $filter('selectedMethods')(form.methods).join(',');
//            }
//            else if (   angular.isDefined(form.species) && angular.isObject(form.species) && form.species &&
//                        angular.isDefined(form.ensembl) && form.ensembl &&
//                        angular.isDefined(form.gene) && form.gene &&
//                        $filter('isSelectedMethod')(form.methods) )
//            {
//                postData.species = form.species.scientific;
//                postData.id = form.gene;
//                postData.e_version = form.ensembl.dataset;
//                //postData.methods = Object.keys(form.methods).join(',')+',appris' // TODO: Improve this 'hard-code shit'
//                postData.methods = $filter('selectedMethods')(form.methods).join(',');
//            }
//            else {
//                form.$invalid = true;
//                return;
//            }
//
//            // Submit job
//            var newRunner = Runner.create(postData);
//
//            // Trigger validation flag.
//            $rootScope.isLoadingScreen = true;
//            $scope.submitted = true;
//
//            newRunner.then(function(data) {
//                if ( angular.isDefined(data.jobid) ) {
//                    var jobId = data.jobid;
//                    $rootScope.isLoadingScreen = false;
//                    $location.url(consPathServerStatus+jobId);
//                }
//                else {
//                    $location.url(consPathServerStatus);
//                }
//            }, function (err) {
//                //console.log('Request failed: '+JSON.stringify(err));
//                //console.log('Request failed: '+err);
//                $rootScope.isLoadingScreen = false;
//                $location.url(consPathServerStatus);
//            });
//        };
//
//        // fill with an example
//        $scope.exRunnerController = function (inEx) {
//            $scope.runnerForm.species = null;
//            $scope.runnerForm.sequences = null;
//            $scope.runnerForm.ensembl = null;
//            $scope.runnerForm.gene = null;
//            if ( inEx == 'seq' ) {
//                var ex = $scope.examples[0];
//                var spe_id = ex.species_id;
//                $scope.runnerForm.species = $scope.species[spe_id];
//                $scope.runnerForm.sequences = ex.sequences;
//            }
//            else if ( inEx == 'id' ) {
//                var ex = $scope.examples[1];
//                var spe_id = ex.species_id;
//                $scope.runnerForm.species = $scope.species[spe_id];
//                $scope.runnerForm.ensembl = $scope.species[spe_id].assemblies[0];
//                $scope.runnerForm.gene = ex.gene;
//            }
//        };
//
//    }
//]);

//apprisControllers.controller('StatusController', ['consPathServerResult', '$scope', '$location', '$routeParams', '$timeout', 'Status',
//    function(consPathServerResult, $scope, $location, $routeParams, $timeout, Status) {
//        // recursive method that retrieves status
//        $scope.control = function(jobid, data){
//            if ( !angular.isDefined(data) ) {
//                $scope.status = {
//                    type:   'danger',
//                    status: 'ERROR',
//                    header: "An error occurred attempting the job: "+jobid,
//                    log: "The status is not defined",
//                    progress: 0
//                };
//            }
//            else {
//                if ( data.status == 'NOT_FOUND' ){
//                    $scope.status = {
//                        type:   'danger',
//                        status: 'ERROR',
//                        header: "An error occurred attempting the job: "+jobid,
//                        log: data.log,
//                        progress: 0
//                    };
//                }
//                else if ( data.status == 'FINISHED' || data.status == 'FAILURE' || data.status == 'ERROR' ){
//                    $scope.status = {
//                        status: data.status,
//                        log:    data.log,
//                        progress: 100
//                    };
//                    if ( data.status == 'FINISHED' ) {
//                        $scope.status.type = 'success';
//                        $scope.status.header = "Your job has finished succesfully";
//                        $location.url(consPathServerResult+jobid);
//                    }
//                    else if ( data.status == 'FAILURE' ) {
//                        $scope.status.type = 'warning';
//                        $scope.status.header = "Your job has finished with some errors";
//                        $location.url(consPathServerResult+jobid);
//                    }
//                    else if ( data.status == 'ERROR' ) {
//                        $scope.status.type = 'danger';
//                        $scope.status.header = "Your job has failed";
//                    }
//                }
//                else if ( data.status == 'PENDING' || data.status == 'RUNNING' || data.status == 'FAILURE' ){
//                    $scope.status = {
//                        type:   'info',
//                        header: "Your job is running to be processed... Please be patient",
//                        status: data.status,
//                        log:    data.log,
//                        progress: (data.progress*100)
//                    };
//                    Status.get({ jobid: jobid },
//                        function(dat) {
//                            $timeout(function() { $scope.control(jobid, dat); }, 5000);
//                        },
//                        function(err) {
//                            $scope.status = {
//                                type: 'danger',
//                                status: 'ERROR',
//                                header:    "Failed to received: "+jobid,
//                                progress: 0
//                            };
//                        }
//                    );
//                }
//            }
//        };
//
//        // init status control
//        if ( angular.isUndefined($routeParams.jobid) ) {
//            $scope.status = {
//                type: 'danger',
//                status: "ERROR",
//                header: "An error occurred attempting the job",
//                log: "The job is not defined",
//                progress: 0
//            };
//        }
//        else {
//            var jobid = $routeParams.jobid;
//            $scope.status = {
//                type: 'info',
//                status: "PENDING",
//                header: "Your job is waiting to be processed... Please be patient",
//                log: "The job is waiting to be processed",
//                progress: 0
//            };
//            $scope.control(jobid, $scope.status);
//        }
//    }
//]);

/* EXPORTER/RUNNER result Controllers */
//apprisControllers.controller('ResultController', ['consPageError', '$rootScope', '$scope', '$routeParams', '$location', '$filter', 'Species', 'Methods', 'ResultTypes', 'Retriever', 'Sequencer', 'Viewer', function (consPageError, $rootScope, $scope, $routeParams, $location, $filter, Species, Methods, ResultTypes, Retriever, Sequencer, Viewer) {
//
//    // init vars
//    $rootScope.isLoadingScreen = true;
//    $scope.alert = {};
//
//    $scope.runnerid = null;
//    $scope.exporterid = null;
//    $scope.jobid = null;
//    $scope.isSeqRunner = false;
//    $scope.assembly = null;
//    $scope.assemblies = [];
//
//    $scope.currentMethods = [];
//    $scope.resultDetailAnnotHeads = [];
//    $scope.resultMainAnnotHeads = [];
//    $scope.resultMainResidues = [];
//    //$scope.browserDataMethods = [];
//    //$scope.browserDataAlign = false;
//    //$scope.browserDataSeqs = null;
//
//    // get route parameters from type of mode
//    // EXPORTER mode
//    if ( angular.isDefined($routeParams.tid) && angular.isDefined($routeParams.species) && angular.isDefined($routeParams.id) ) {
//        $scope.exporterid = $routeParams.tid+'/'+$routeParams.species+'/'+$routeParams.id;
//        $scope.jobid = $scope.exporterid;
//    }
//    // RUNNER mode (idRUNNER and seqRUNNER)
//    else if ( angular.isDefined($routeParams.jobid) ) {
//        $scope.runnerid = $routeParams.jobid;
//        $scope.jobid = $scope.runnerid;
//    }
//
////    if ( angular.isDefined($routeParams.species) ) {
////        $scope.assemblies = $rootScope.species[$routeParams.species].assemblies;
////        $scope.assembly = $rootScope.species[$routeParams.species].assembly;
////    }
//
////    // load list of species
////    if ( angular.isUndefined($rootScope.species) ) {
////        $rootScope.species = null;
////        Species.get().$promise.then( function(data) {
////            $rootScope.species = data;
////            if ( angular.isDefined($routeParams.species) ) {
////                $scope.assemblies = $rootScope.species[$routeParams.species].assemblies;
////                $scope.assembly = $rootScope.species[$routeParams.species].assembly;
////            }
////        });
////    }
////    else { // add assembly
////        if ( angular.isDefined($routeParams.species) ) {
////            $scope.assemblies = $rootScope.species[$routeParams.species].assemblies;
////            $scope.assembly = $rootScope.species[$routeParams.species].assembly;
////        }
////    }
//
////    // load list of methods
////    if ( angular.isUndefined($rootScope.methods) || (angular.isArray($rootScope.methods) && $rootScope.methods.length == 0 ) ) {
////        $rootScope.methods = null;
////        Methods.get().$promise.then( function(data) {
////            $rootScope.methods = data;
////        });
////    }
//
//    // APPRIS annots
//
//    // load result types (method results). BY DEFAULT: retrieves all
//    ResultTypes.query({ jobid: $scope.runnerid }).$promise.then( function(data) {
//        if ( $filter('isSeqRunner')(data) ) {
//            $scope.isSeqRunner = true;
//        }
//        // Add seq. id and methods
//        $scope.resultDetailAnnotHeads.push({
//            id: "transcript_id",
//            label: "Seq. id"
//        },{
//            id: "length_aa",
//            label: "Length (aa)"
//        });
//        angular.forEach(data, function(item) {
//            // discard methods for detailed panel!!!
//            if ( item.id !== 'principal_isoform' ) {
//                $scope.resultDetailAnnotHeads.push({
//                    id: item.id,
//                    label: item.label
//                });
//            }
//            // discard methods for browser panel !!!
//            $scope.currentMethods.push({
//                id: item.id,
//                name: item.name,
//                label: item.label,
//                desc: item.desc
//            });
//        });
//    });
//
//    // create annotation table depending on mode (EXPORTER or RUNNER) at the moment we have already the resultHeadLbs
//    $scope.$watch('resultDetailAnnotHeads', function (value) {
//        if ( angular.isDefined(value) && value !== null ) {
//            var queryData = null;
//            if ( $scope.exporterid ) {
//                queryData = {
//                    tid: $routeParams.tid,
//                    species: $routeParams.species,
//                    id: $routeParams.id,
//                    //ens: $routeParams.ens,
//                    db: $routeParams.db,
//                    methods: 'appris'
//                };
//            }
//            else if ( $scope.runnerid ) {
//                queryData = {
//                    jobid: $scope.runnerid,
//                    methods: 'appris'
//                };
//            }
//
//            var rawResults = Retriever.getJSON(queryData);
//            rawResults
//                .then( function(data) {
//                    if ( $scope.exporterid ) {
//                        var resultData = $filter('convertTransScoreObj')(data);
//                        //$scope.resultAnnots = resultData[1];
//                        $scope.resultMainAnnotHeads.push({
//                            id: "transcript_id",
//                            label: "Seq. id"
//                        },{
//                            id: "transcript_name",
//                            label: "Seq. name"
//                        },{
//                            id: "biotype",
//                            label: "Biotype"
////                            },{
////                                id: "length_aa",
////                                label: "Length (aa)"
//                        },{
//                            id: "no_codons",
//                            label: "Codons not found"
//                        },{
//                            id: "ccds_id",
//                            label: "CCDS"
//                        },{
//                            id: "principal_isoform",
//                            label: "Principal Isoform"
//                        });
//                    }
//                    else if ( $scope.runnerid ) {
//                        if ( $scope.isSeqRunner ) {
//                            var resultData = $filter('convertSeqRunnerObj')(data);
//                            $scope.resultMainAnnotHeads.push({
//                                id: "transcript_id",
//                                label: "Seq. id"
//                            },{
//                                id: "principal_isoform",
//                                label: "Principal Isoform"
//                            });
//                        }
//                        else {
//                            var resultData = $filter('convertTransScoreObj')(data);
//                            $scope.resultMainAnnotHeads.push({
//                                id: "transcript_id",
//                                label: "Seq. id"
//                            },{
//                                id: "transcript_name",
//                                label: "Seq. name"
////                                },{
////                                    id: "length_aa",
////                                    label: "Length (aa)"
//                            },{
//                                id: "principal_isoform",
//                                label: "Principal Isoform"
//                            });
//                        }
//                    }
//                    $rootScope.isLoadingScreen = false;
//                    if ( angular.isUndefined(resultData) || (angular.isArray(resultData) && resultData.length == 0 ) ) {
//                        $scope.alert.enable = true;
//                        $scope.alert.type = "Warning";
//                        $scope.alert.message = "Your query were not found. Refreshing page helps to fix the problem (clearing the cache). If problem persistent, Please, contact the administrators";
//                    }
//                    else {
//                        $scope.resultAnnots = resultData[1];
//                        $scope.addBrowserDatalSeqs();
//                        $scope.browserDataMethods = angular.copy($scope.currentMethods);
//                        //var inputs = $scope.resultAnnots;
//                        //var methods = 'appris';
//                        //if ( angular.isDefined($scope.browserDataMethods) ) { methods = $scope.browserDataMethods; }
//                        //var align = $scope.browserDataAlign;
//                        $scope.submitBrowserData($scope.resultAnnots, $scope.currentMethods, $scope.browserDataAlign);
//                    }
//                }, function(error) {
//                    $rootScope.isLoadingScreen = false;
//                    if ( error.status >= 400 && error.status < 500 ) {
//                        $scope.alert.enable = true;
//                        $scope.alert.type = "Error";
//                        $scope.alert.message = error.data;
//                    }
//                    else if ( error.status >= 500 ) {
//                        //$location.path(consPageError);
//                        $location.url(consPageError);
//                    }
//                });
//        }
//    });
//
//
////    // When is ready the result annotations
////    $scope.$watch('resultAnnots', function (inputs) {
////        if ( angular.isArray(inputs) && inputs.length > 0 ) {
////            $scope.addBrowserDatalSeqs();
////            $scope.browserDataMethods = angular.copy($scope.currentMethods);
////            var inputs = $scope.resultAnnots;
////            var methods = 'appris';
////            if ( angular.isDefined($scope.browserDataMethods) ) { methods = $scope.browserDataMethods; }
////            var align = $scope.browserDataAlign;
////            $scope.submitBrowserData(inputs, methods, align);
////        }
////    });
//    $scope.addBrowserDatalSeqs = function () {
//        if ( $scope.browserDataSeqs == null ) {
//            $scope.browserDataSeqs = true;
//        }
//        else if ( $scope.browserDataSeqs ) {
//            $scope.browserDataSeqs = false;
//        }
//        else {
//            $scope.browserDataSeqs = true;
//        }
//        angular.forEach($scope.resultAnnots, function (item) {
//            item.selected = $scope.browserDataSeqs;
//        });
//    };
//
//    // APPRIS Browser tabs
//    $scope.browserTabs = [{
//        id:         'annotation',
//        title:      'Annotation Browser',
//        refresh:    false,
//        info:       '#infoBrowserAnnot',
//        panel:      null
//    },{
//        id:         'sequence',
//        title:      'Sequence Browser',
//        refresh:    false,
//        info:       '#infoBrowserSeq',
//        panel:      null
//    },
//        {
//            id:         'genome',
//            title:      'Genome Browser',
//            refresh:    false,
//            info:       '#infoBrowserGen',
//            panel:      null
//        }];
//    $scope.browserActive = function() {
//        return $scope.browserTabs.filter(function(tab){
//            return tab.active;
//        })[0];
//    };
//
//    // submit data
//    $scope.submitBrowserData = function(inputs, methods, align) {
//        angular.forEach($scope.browserTabs, function(browser) {
//            // init vars
//            browser.refresh = true;
//            browser.frame = null;
//            // creaate data
//            var browserData = null;
//            if ( $scope.exporterid ) {
//                browserData = {
//                    tid: $routeParams.tid,
//                    species: $routeParams.species,
//                    id: $routeParams.id,
//                    //ens: $routeParams.ens
//                    db: $routeParams.db
//                };
//            }
//            else if ( $scope.runnerid ) {
//                browserData = {
//                    jobid: $scope.runnerid
//                };
//            }
//            var idList = '';
//            angular.forEach($filter('orderBy')(inputs, 'transcript_id'), function(item) {
//                if ( item.selected ) {
//                    idList += item.transcript_id+',';
//                }
//            });
//            idList = idList.replace(/,$/g,'');
//            browserData.ids = idList;
//            var methodList = 'appris,';
//            angular.forEach(methods, function(item) {
//                methodList += item.name+',';
//            });
//            methodList = methodList.replace(/,$/g,'');
//            browserData.methods = methodList;
//            // call services
//            if ( (browser.id == 'annotation') ) {
//                browserData.operation = 'residues';
//                var View = Sequencer.get(browserData);
//                View.then(function (result) {
//                    $scope.resultMainResidues = result;
//                    browser.frame = result;
//                    browser.refresh = false;
//                }, function(error) {
//                    browser.frame = error;
//                    browser.refresh = false;
//                });
//            }
//            else if ( (browser.id == 'sequence') ) {
////                browserData.operation = 'sequences';
////                if ( angular.isDefined(align) && align ) {
////                    browserData.operation = 'align';
////                }
////                var View = Viewer.get(browserData);
////                View.then(function (result) {
////                    browser.frame = result;
////                    browser.refresh = false;
////                }, function(error) {
////                    browser.frame = error;
////                    browser.refresh = false;
////                });
//                browserData.operation = 'sequences';
//                if ( angular.isDefined(align) && align ) {
//                    browserData.operation = 'align';
//                }
//                var View = Sequencer.get(browserData);
//                View.then(function (result) {
//                    $scope.resultMainResidues = result;
//                    browser.frame = result;
//                    browser.refresh = false;
//                }, function(error) {
//                    browser.frame = error;
//                    browser.refresh = false;
//                });
//            }
//            else if ( (browser.id == 'genome') ) {
//                browserData.operation = 'genome';
//                var View = Viewer.get(browserData);
//                View.then(function (result) {
//                    browser.frame = result;
//                    browser.refresh = false;
//                }, function(error) {
//                    browser.frame = error;
//                    browser.refresh = false;
//                });
//            }
//        });
//    };
//
//
//    // Export data
//    $scope.saveAs = [{
//        type: "annotation",
//        name: "Annotations as",
//        format: "json",
//        mimetype:   "application/json"
//    },{
//        type: "annotation",
//        name: "Annotations as",
//        format: "gtf",
//        mimetype:   "text/plain"
//    },{
//        type: "annotation",
//        name: "Annotations as",
//        format: "tsv",
//        mimetype:   "text/plain"
//    },{
//        type: "sequence",
//        name: "Sequence annot. as",
//        format: "json",
//        mimetype:   "application/json"
//    },{
//        type: "sequence",
//        name: "Protein seq. as",
//        format: "fasta",
//        mimetype:   "text/plain"
//    },{
//        type: "genome",
//        name: "Genome annot. as",
//        format: "bed",
//        mimetype:   "text/plain"
//    }];
//    $scope.exportData = function(save) {
//        // init vars
//        var inputs = $scope.resultAnnots;
//        var methods = $scope.browserDataMethods;
//        var rawResults = null;
//        var browserID = save.type;
//        var format = save.format;
//        // create data (path)
//        var browserData = null;
//        if ( $scope.exporterid ) {
//            browserData = {
//                tid: $routeParams.tid,
//                species: $routeParams.species,
//                id: $routeParams.id,
//                //ens: $routeParams.ens
//                db: $routeParams.db
//            };
//        }
//        else if ( $scope.runnerid ) {
//            browserData = {
//                jobid: $scope.runnerid
//            };
//        }
//        // create data (parameters)
//        var ids = '';
//        angular.forEach($filter('orderBy')(inputs, 'transcript_id'), function(item) {
//            if ( item.selected ) {
//                ids += item.transcript_id+',';
//            }
//        });
//        ids = ids.replace(/,$/g,'');
//        var mets = '';
//        angular.forEach(methods, function(item) {
//            mets += item.name+',';
//        });
//        mets = mets.replace(/,$/g,'');
//        //if ( ids !== '' ) browserData.ids = ids;
//        //if ( mets !== '' ) browserData.methods = mets;
//        //if ( format !== '' ) browserData.format = format;
//        browserData.ids = ids;
//        browserData.format = format;
//        browserData.methods = mets;
//        // call services
//        if ( browserID == 'annotation' ) {
//            if ( format === 'json' ) {
//                rawResults = Retriever.getJSON(browserData);
//            }
//            else if ( format === 'gtf' ) {
//                rawResults = Retriever.getGTF(browserData);
//            }
//            else if ( format === 'tsv' ) {
//                rawResults = Retriever.getTXT(browserData);
//            }
//            $scope.saveName = $routeParams.id + '.' + browserID;
//        }
//        else if ( (browserID == 'sequence') ) {
//            if ( format === 'json' ) {
//                browserData.operation = 'residues';
//                rawResults = Sequencer.get(browserData);
//            }
//            else if ( format === 'fasta' ) {
//                browserData.operation = 'sequences';
//                rawResults = Sequencer.get(browserData);
//            }
//            $scope.saveName = $routeParams.id + '.' + browserID;
//        }
//        else if ( (browserID == 'genome') ) {
//            if ( format === 'bed' ) {
//                rawResults = Retriever.getTXT(browserData);
//            }
//            $scope.saveName = $routeParams.id + '.' + browserID;
//        }
//        return rawResults;
//    };
//}
//]);

/* The Gene Info Controllers */
apprisControllers.controller('GeneResultController', ['consUrlEnsembl', '$rootScope', '$scope', '$routeParams', '$filter', 'Species', 'Seeker',
    function(consUrlEnsembl, $rootScope, $scope, $routeParams, $filter, Species, Seeker) {

        // init vars
        var rawResults = null;
        var ens = null;

        // get route parameters from type of mode
        // EXPORTER mode
        if ( angular.isDefined($routeParams.tid) && angular.isDefined($routeParams.species) && angular.isDefined($routeParams.id) ) {
            rawResults = Seeker.query({
                id: $routeParams.id
            });
            if ( angular.isDefined($routeParams.ens) ) {
                ens = $routeParams.ens;
            }
        }
        // RUNNER mode
        else if ( angular.isDefined($routeParams.jobid) ) {
            rawResults = Seeker.query({
                id: $routeParams.jobid
            });
        }

        // create info
        if ( angular.isDefined(rawResults) ) {
            rawResults.$promise.then( function(data) {
                if ( angular.isDefined(data) && angular.isArray(data.match) && data.match.length > 0 ) {
                    var item = data.match[0]; // get the first
                    var specie_id = item.species;
                    var rst = [];
                    if ( angular.isDefined($rootScope.species[specie_id]) ) {
                        var speciesInfo = $rootScope.species[specie_id];
                        if ( angular.isDefined(item.namespace) && item.namespace !== null ) {
                            rst.push({
                                "label": "Id",
                                "value": item.id,
                                "link_ensembl": consUrlEnsembl + '/' + specie_id + '/Gene/Summary?db=core;g=' + item.id
                            });
                        }
                        else {
                            rst.push({
                                "label": "Id",
                                "value": item.id
                            });
                        }
                        if ( angular.isDefined(speciesInfo.scientific) && speciesInfo.scientific !== null ) {
                            rst.push({
                                "label": "Species",
                                "value": $rootScope.species[specie_id].scientific
                            });
                        }
                        if ( angular.isDefined(item.dblink) && item.dblink !== null && angular.isDefined($filter('filter')(item.dblink, { "namespace": 'External_Id' })) && angular.isArray($filter('filter')(item.dblink, { "namespace": 'External_Id' })) && ($filter('filter')(item.dblink, { "namespace": 'External_Id' }).length > 0) ) {
                            rst.push({
                                "label": "Name",
                                "value": $filter('filter')(item.dblink, { "namespace": 'External_Id' })[0].id
                            });
                        }
                        if ( angular.isDefined(item.biotype) && item.biotype !== null ) {
                            rst.push({
                                "label": "Biotype",
                                "value": item.biotype
                            });
                        }
                        if ( angular.isDefined(speciesInfo.assemblies) && angular.isArray(speciesInfo.assemblies) && (speciesInfo.assemblies.length > 0)) {
                            var assembly = $filter('filter')(speciesInfo.assemblies, { "dataset": ens });
                            if ( angular.isDefined(assembly) && angular.isArray(assembly) && (assembly.length > 0) ) {
                                rst.push({
                                    "label": "Assembly",
                                    "value": assembly[0].name
                                });
                            }
                        }
                        if ( angular.isDefined(item.chr) && item.chr !== null && angular.isDefined(item.start) && item.start !== null && angular.isDefined(item.end) && item.end !== null ) {
                            rst.push({
                                "label": "Location",
                                "value": item.chr+':'+item.start+'-'+item.end
                            });
                        }
                    }
                    $scope.geneInfo = rst;
                }
            });
        }
    }
]);

/*
 *
 * SEEKER Controllers
 *
 */
apprisControllers.controller('SeekerController', ['consPathSeeker', '$scope', '$location',
    function(consPathSeeker, $scope, $location) {

        // init form var
        $scope.seekerForm = {};

        // submit search
        $scope.runJob = function (form) {

            // Validate
            // empty specie input
            if ( angular.isUndefined(form.inQuery) && (form.inQuery != '')) {
                $scope.queInvalid = true;
            }
            // messages
            if ( $scope.queInvalid ) {
                return;
            }
            // Redirect URL
            //$location.path(consPathSeeker+form.inQuery);
            $location.url(consPathSeeker+form.inQuery);
        };
    }
]);

apprisControllers.controller('SeekerResultController', ['consQueryNotMatch', '$rootScope', '$scope', '$routeParams', '$location', '$filter', 'Species', 'Seeker',
    function(consQueryNotMatch, $rootScope, $scope, $routeParams, $location, $filter, Species, Seeker) {

        // init vars
        $rootScope.isLoadingScreen = true;
        var id = $routeParams.id;

        // load specie list
        if ( angular.isUndefined($rootScope.species) ) {
            $rootScope.species = null;
            Species.get().$promise.then( function(data) {
                $rootScope.species = data;
            });
        }

        // load results

        // create search summary
        var specieResult = {};
        var rawResults = Seeker.query({
            id: id
        });
        rawResults.$promise.then( function(data) {
            angular.forEach(data.match, function(item) {
                if (  angular.isUndefined(specieResult[item.species]) ) {
                    specieResult[item.species] = [];
                }
                var assemblies = $rootScope.species[item.species].assemblies;
                angular.forEach(assemblies, function(assembly) {
                    if ( assembly.dataset === parseInt(item.dataset) ) {
                        var speRst = {
                            "species":      item.species,
                            "assembly": {
                                "id": $filter('filter')($rootScope.species[item.species].assemblies, { "name": item.assembly })[0].id,
                                "name": $filter('filter')($rootScope.species[item.species].assemblies, { "name": item.assembly })[0].name
                            },
                            "dataset":      item.dataset,
                            "label":        item.label,
                            "namespace":    item.namespace,
                            "id":           item.id,
                            "biotype":      item.biotype,
                            "dblink":       item.dblink,
                            "position":     'chr'+item.chr+':'+item.start+'-'+item.end
                        }
                        specieResult[item.species].push(speRst);
                    }
                });
            });

            // sort by specie
            $scope.seekerResult = [];
            angular.forEach($rootScope.species, function(item, key) {
                if (  angular.isDefined(specieResult[key]) ) {
                    var rst = {
                        "species":   item.scientific,
                        "match":    specieResult[key]
                    };
                    $scope.seekerResult.push(rst);
                    $rootScope.isLoadingScreen = false;
                }
            });

        }, function(error) {
            $rootScope.isLoadingScreen = false;
            if ( error.status >= 400 && error.status < 500 ) {
                //$location.path(consQueryNotMatch);
                $location.url(consQueryNotMatch);
            }
            else if ( error.status >= 500 ) {
                //$location.path(consQueryNotMatch);
                $location.url(consQueryNotMatch);
            }
        });

    }
]);


/* The REST of Controllers */
apprisControllers.controller('DownloadsController', ['consBaseUrlWS', '$scope', 'Downloads',
    function(consBaseUrlWS, $scope, Downloads) {
        // init pagination
        $scope.currentPage = 1;
        $scope.pageSize = 4;

        // load all data files
        $scope.urlREADME = consBaseUrlWS + '/download/README.txt';
        $scope.baseUrlDownload = consBaseUrlWS + '/download/data';
        $scope.headers = [];
        $scope.datafiles = [];
        $scope.species = [];
        $scope.assemblies = [];
        $scope.datasets = [];
        Downloads.get().$promise.then( function(data) {
            if (angular.isDefined(data.headers) ) {
                $scope.headers = data.headers;
            }
            if (angular.isDefined(data.datafiles) ) {
                $scope.datafiles = data.datafiles;
                angular.forEach(data.datafiles, function(val) {
                    var specie = val['specie'];
                    var assembly = val['assembly'];
                    var dataset = val['dataset'];
                    if ( $scope.species.indexOf(specie) == -1 ) {
                        $scope.species.push(specie);
                    }
                    if ( $scope.assemblies.indexOf(assembly) == -1 ) {
                        $scope.assemblies.push(assembly);
                    }
                    if ( $scope.datasets.indexOf(dataset) == -1 ) {
                        $scope.datasets.push(dataset);
                    }
                });
            }
        });
    }
]);

apprisControllers.controller('NavTopController', ['$scope',
    function($scope) {
        $scope.apiLink = 'http://apprisws.cnio.es/apidoc';
    }
]);

apprisControllers.controller('AboutController', ['$scope',
    function($scope) {
        $scope.logoCNIO = 'img/CNIO-logo_small.png';
        $scope.logoINB = 'img/INB-logo_small.png';
    }
]);

apprisControllers.controller('pagDownloadsController', ['$scope',
    function($scope) {
        $scope.currentPage = 1;
        $scope.pageSize = 4;
        $scope.pageChangeHandler = function(num) {
//            console.log('going to page ' + num);
        };
    }
]);

/* Alert Controllers */
apprisControllers.controller('AlertController', ['$rootScope',
    function($rootScope) {

//        $timeout(function() { alert.expired = true; }, 2000);

    }
]);

apprisControllers.controller('HelpController', ['$rootScope', '$scope', '$location', '$routeParams', 'Species',
    function($rootScope, $scope, $location, $routeParams, Species) {
        $scope.help = $routeParams.help;
        $scope.isActive = function (viewLocation) {
            return viewLocation === $location.path();
        };

        // load list of species
        if ( angular.isUndefined($rootScope.species) ) {
            $rootScope.species = null;
            Species.get().$promise.then( function(data) {
                $rootScope.species = data;
            });
        }
    }
]);
