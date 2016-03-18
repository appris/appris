'use strict';

/* Controllers */

var apprisControllers = angular.module('apprisControllers', []);

apprisControllers.run(function($rootScope, serverName, serverType, serverHost, serverHostWS, Species, Methods, Downloads) {
    $rootScope.serverConf = {
        "name":    serverName,
        "type":    serverType,
        "host":    serverHost,
        "hostWS":  serverHostWS
    };
    $rootScope.downloads = Downloads.get();
    $rootScope.species = Species.get();
    $rootScope.methods = Methods.get();
    $rootScope.isLoadingScreen = false;
})

/* The Gene Info Controllers */
apprisControllers.controller('GeneResultController', ['consUrlEnsembl', '$rootScope', '$scope', '$routeParams', '$filter', 'Seeker',
    function(consUrlEnsembl, $rootScope, $scope, $routeParams, $filter, Seeker) {

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

apprisControllers.controller('SeekerResultController', ['consQueryNotMatch', '$rootScope', '$scope', '$routeParams', '$location', '$filter', 'Seeker',
    function(consQueryNotMatch, $rootScope, $scope, $routeParams, $location, $filter, Seeker) {

        // init vars
        $rootScope.isLoadingScreen = true;
        var id = $routeParams.id;

//        // load specie list
//        if ( angular.isUndefined($rootScope.species) ) {
//            $rootScope.species = null;
//            Species.get().$promise.then( function(data) {
//                $rootScope.species = data;
//            });
//        }

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
apprisControllers.controller('DownloadsController', ['$rootScope', '$scope', 'Downloads',
    function($rootScope, $scope, Downloads) {
        // init pagination
        $scope.currentPage = 1;
        $scope.pageSize = 4;

        $scope.isArray = angular.isArray;

        // load all data files
        $scope.urlREADME = $rootScope.serverConf.hostWS + '/pub/README.txt';
        $scope.urlDOWNLOAD = $rootScope.serverConf.hostWS + '/pub';
        if ( $rootScope.serverConf.type == "gold" ) {
            $scope.baseUrlDownload = $rootScope.serverConf.hostWS + '/pub/current_release';
        }
        else {
            $scope.baseUrlDownload = $rootScope.serverConf.hostWS + '/pub/releases/' + $rootScope.downloads.version;
        }
        $scope.headers = [];
        $scope.datafiles = [];
        $scope.species = [];
        $scope.assemblies = [];
        $scope.datasets = [];
        if ( angular.isDefined($rootScope.downloads) ) {
            $scope.version = $rootScope.downloads.version;
            if (angular.isDefined($rootScope.downloads.headers) ) {
                $scope.headers = $rootScope.downloads.headers;
            }
            if (angular.isDefined($rootScope.downloads.datafiles) ) {
                $scope.datafiles = $rootScope.downloads.datafiles;
                angular.forEach($rootScope.downloads.datafiles, function(val) {
                    var specie = val['species'];
                    var assembly = val['assembly'];
                    var dataset = val['dataset'];
                    var path = val['path'];
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
        }
//        Downloads.get().$promise.then( function(data) {
//            $scope.version = data.version;
//            if (angular.isDefined(data.headers) ) {
//                $scope.headers = data.headers;
//            }
//            if (angular.isDefined(data.datafiles) ) {
//                $scope.datafiles = data.datafiles;
//                angular.forEach(data.datafiles, function(val) {
//                    var specie = val['species'];
//                    var assembly = val['assembly'];
//                    var dataset = val['dataset'];
//                    var path = val['path'];
//                    if ( $scope.species.indexOf(specie) == -1 ) {
//                        $scope.species.push(specie);
//                    }
//                    if ( $scope.assemblies.indexOf(assembly) == -1 ) {
//                        $scope.assemblies.push(assembly);
//                    }
//                    if ( $scope.datasets.indexOf(dataset) == -1 ) {
//                        $scope.datasets.push(dataset);
//                    }
//                });
//            }
//        });
    }
]);

apprisControllers.controller('NavTopController', ['$rootScope', function($rootScope) { } ]);

apprisControllers.controller('AboutController', ['$scope',
    function($scope) {
        $scope.logoCNIO = 'img/CNIO-logo_small.png';
        $scope.logoINB = 'img/INB-logo_small.png';
        $scope.linkCNIO = 'http://www.cnio.es';
    }
]);

apprisControllers.controller('ChangelogsController', ['$scope', 'Changelogs',
    function($scope, Changelogs) {
        // init pagination
        $scope.currentPage = 1;
        $scope.pageSize = 4;

        $scope.changelogs = Changelogs.query();
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

apprisControllers.controller('HelpController', ['$rootScope', '$scope', '$location', '$routeParams',
    function($rootScope, $scope, $location, $routeParams) {
        $scope.help = $routeParams.help;
        $scope.isActive = function (viewLocation) {
            return viewLocation === $location.path();
        };
    }
]);
