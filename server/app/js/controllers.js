'use strict';

/* Controllers */

var apprisControllers = angular.module('apprisControllers', []);

apprisControllers.run(function($rootScope, serverName, serverType, serverHost, serverHostWS, Server, Methods) {
    $rootScope.serverConf = {
        "name":    serverName,
        "type":    serverType,
        "host":    serverHost,
        "hostWS":  serverHostWS,
        "version": ""
    };
    $rootScope.server = Server.get();
    $rootScope.server.$promise.then( function(data) {
        $rootScope.serverConf.version = data.version;
        $rootScope.species = data.species;
        $rootScope.downloads = data.datafiles;
    });
    $rootScope.methods = Methods.get();
    $rootScope.isLoadingScreen = false;
})

/* The Gene Info Controllers */
apprisControllers.controller('GeneResultController', ['consUrlEnsembl', 'consUrlRefSeq', 'consUrlUniProt', '$rootScope', '$scope', '$routeParams', '$filter', 'Seeker',
    function(consUrlEnsembl, consUrlRefSeq, consUrlUniProt, $rootScope, $scope, $routeParams, $filter, Seeker) {

        // init vars
        var rawResults = null;
        $scope.limitString = 20;

        // get route parameters from type of mode
        // EXPORTER mode
        if ( angular.isDefined($routeParams.tid) && angular.isDefined($routeParams.species) && angular.isDefined($routeParams.id) ) {
            var queryData = {
                sp: $routeParams.species,
                id: $routeParams.id
            };
            if ( angular.isDefined($routeParams.as) ) { queryData.as = $routeParams.as }
            if ( angular.isDefined($routeParams.ds) ) { queryData.ds = $routeParams.ds }
            rawResults = Seeker.get(queryData);
        }
        // RUNNER mode
        else if ( angular.isDefined($routeParams.jobid) ) {
            rawResults = Seeker.get({
                id: $routeParams.jobid
            });
        }

        // create info
        if ( angular.isDefined(rawResults) ) {
            rawResults.then(function (data) {
                if ( angular.isDefined(data) && angular.isArray(data.match) && data.match.length > 0 ) {
                    var item = data.match[0]; // get the first
                    var specie_id = item.species;
                    var rst = [];
                    if ( angular.isDefined($rootScope.species[specie_id]) ) {
                        var speciesInfo = $rootScope.species[specie_id];
                        if ( angular.isDefined(item.namespace) && item.namespace == 'Gene_Id' && angular.isDefined($routeParams.sc) ) {
                            var link_source = '';
                            if ( $routeParams.sc == 'ensembl' ) {
                                link_source = consUrlEnsembl + '/' + specie_id + '/Gene/Summary?db=core;g=' + item.id
                            }
                            else if ( $routeParams.sc == 'refseq' ) {
                                link_source = consUrlRefSeq + '/' + 'gene/' + item.id
                            }
                            else if ( $routeParams.sc == 'uniprot' ) {
                                link_source = consUrlUniProt + '/' + 'uniprot/?query=' + item.id + '&fil=organism:' + specie_id;
                            }
                            if ( link_source != '' ) {
                                rst.push({
                                    "label": "Id",
                                    "value": item.id,
                                    "source": $routeParams.sc,
                                    "link_source": link_source
                                });
                            }
                        }
                        else {
                            rst.push({
                                "label": "Id",
                                "value": item.id
                            });
                        }
                        if ( angular.isDefined(item.dblink) && item.dblink !== null ) {
                            if ( angular.isDefined($filter('filter')(item.dblink, { "namespace": 'Ensembl_Gene_Id' })) && angular.isArray($filter('filter')(item.dblink, { "namespace": 'Ensembl_Gene_Id' })) && ($filter('filter')(item.dblink, { "namespace": 'Ensembl_Gene_Id' }).length > 0) ) {
                                var id = $filter('filter')(item.dblink, { "namespace": 'Ensembl_Gene_Id' }).map(function(e){ return e.id }).join(",");
                                rst.push({
                                    "label": "Ensembl Gene Id",
                                    "value": id,
                                    "source": "ensembl",
                                    "link_source":  consUrlEnsembl + '/' + specie_id + '/Gene/Summary?db=core;g=' + id
                                });
                            }
                            if ( angular.isDefined($filter('filter')(item.dblink, { "namespace": 'Refseq_Gene_Id' })) && angular.isArray($filter('filter')(item.dblink, { "namespace": 'Refseq_Gene_Id' })) && ($filter('filter')(item.dblink, { "namespace": 'Refseq_Gene_Id' }).length > 0) ) {
                                var id = $filter('filter')(item.dblink, { "namespace": 'Refseq_Gene_Id' }).map(function(e){ return e.id }).join(",");
                                rst.push({
                                    "label": "RefSeq Gene Id",
                                    "value": id,
                                    "source": "refseq",
                                    "link_source": consUrlRefSeq + '/' + 'gene/' + id
                                });
                            }
                            if ( angular.isDefined($filter('filter')(item.dblink, { "namespace": 'Uniprot_Gene_Id' })) && angular.isArray($filter('filter')(item.dblink, { "namespace": 'Uniprot_Gene_Id' })) && ($filter('filter')(item.dblink, { "namespace": 'Uniprot_Gene_Id' }).length > 0) ) {
                                var id = $filter('filter')(item.dblink, { "namespace": 'Uniprot_Gene_Id' }).map(function(e){ return e.id }).join(",");
                                rst.push({
                                    "label": "UniProt Gene Id",
                                    "value": id,
                                    "source": "uniprot",
                                    "link_source": consUrlUniProt + '/' + 'uniprot/?query=' + id + '&fil=organism:' + specie_id
                                });
                            }
                            if ( angular.isDefined($filter('filter')(item.dblink, { "namespace": 'External_Id' })) && angular.isArray($filter('filter')(item.dblink, { "namespace": 'External_Id' })) && ($filter('filter')(item.dblink, { "namespace": 'External_Id' }).length > 0) ) {
                                rst.push({
                                    "label": "Name",
                                    "value": $filter('filter')(item.dblink, { "namespace": 'External_Id' }).map(function(e){ return e.id }).join(",")
                                });
                            }
                        }
                        if ( angular.isDefined(item.biotype) && item.biotype !== null ) {
                            rst.push({
                                "label": "Biotype",
                                "value": item.biotype
                            });
                        }
                        if ( angular.isDefined(speciesInfo.scientific) && speciesInfo.scientific !== null ) {
                            rst.push({
                                "label": "Species",
                                "value": $rootScope.species[specie_id].scientific
                            });
                        }
                        if ( angular.isDefined(speciesInfo.assemblies) && angular.isArray(speciesInfo.assemblies) && (speciesInfo.assemblies.length > 0)) {
                            var assembly = $filter('filter')(speciesInfo.assemblies, { "id":  $routeParams.as });
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

/* The REST of Controllers */
apprisControllers.controller('DownloadsController', ['$rootScope', '$scope', '$filter',
    function($rootScope, $scope, $filter) {
        // init pagination
        $scope.currentPage = 1;
        $scope.pageSize = 4;

        // load all data files
        $scope.urlREADME = $rootScope.serverConf.hostWS + '/pub/README.md';
        $scope.urlDOWNLOAD = $rootScope.serverConf.hostWS + '/pub';
        if ( $rootScope.serverConf.type == "gold" ) {
            $scope.baseUrlDownload = $rootScope.serverConf.hostWS + '/pub/current_release/datafiles';
        }
        else {
            $scope.baseUrlDownload = $rootScope.serverConf.hostWS + '/pub/releases/' + $rootScope.serverConf.version + '/datafiles';
        }
        $scope.headers = [];
        $scope.downloads = [];
        $scope.species = [];
        $scope.assemblies = [];
        $scope.datasets = [];
        if ( angular.isDefined($rootScope.downloads) ) {
            $scope.headers = $rootScope.downloads;
            angular.forEach($rootScope.species, function(species, species_id) {
                angular.forEach(species.assemblies, function(assembly) {
                    angular.forEach(assembly.datasets, function(dataset, index) {
                        if ( angular.isDefined(dataset.datafiles) ) {
                            var ds = $filter('extractSourceName')(dataset.source);
                            if ( $scope.species.indexOf(species.common) == -1 ) {
                                $scope.species.push(species.common);
                            }
                            if ( $scope.assemblies.indexOf(assembly.name) == -1 ) {
                                $scope.assemblies.push(assembly.name);
                            }
                            if ( $scope.datasets.indexOf(ds) == -1 ) {
                                $scope.datasets.push(ds);
                            }
                            var path = species_id;
                            if ( index === 0 ) {
                                path += '/' + assembly.name;
                            }
                            else {
                                path += '/' + dataset.id;
                            }
                            var data = {
                                "species": species.common,
                                "assembly": assembly.name,
                                "dataset": ds,
                                "path": path,
                                "datafiles": []
                            };
                            angular.forEach(dataset.datafiles, function(datafile, key) {
                                data.datafiles[key] = datafile;
                            });
                            $scope.downloads.push(data);
                        }
                    });
                });
            });
        }
    }
]);

apprisControllers.controller('NavTopController', ['$scope',
    function($scope) {
        $scope.status = {
            isopen: false
        };

        $scope.toggleDropdown = function($event) {
            $event.preventDefault();
            $event.stopPropagation();
            $scope.status.isopen = !$scope.status.isopen;
        };
    }
]);

apprisControllers.controller('AboutController', ['$scope',
    function($scope) {
        $scope.logoCNIO = 'img/CNIO-logo_small.png';
        $scope.linkCNIO = 'http://www.cnio.es';
        $scope.logoINB = 'img/INB-logo_small.png';
        $scope.linkINB = 'http://www.inab.org';
        $scope.linkGENCODE = 'http://www.gencodegenes.org/';
        $scope.linkENCODE = 'https://www.encodeproject.org/';
        $scope.linkEnsembl = 'https://www.ensembl.org/';
    }
]);

apprisControllers.controller('ContactController', ['$scope',
    function($scope) {
        $scope.logoBlogger = 'img/blogger-ico.png';
        $scope.linkBlogger = 'http://appris-cnio.blogspot.com.es';
        $scope.logoTwitter = 'img/twitter-ico.png';
        $scope.linkTwitter = 'https://twitter.com/appris_cnio';
    }
]);

apprisControllers.controller('ChangelogsController', ['$scope',
    function($scope) {

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

apprisControllers.controller('ImptController', ['$rootScope', '$scope', '$location', '$routeParams',
    function($rootScope, $scope, $location, $routeParams) {
        $scope.help = $routeParams.help;
        $scope.isActive = function (viewLocation) {
            return viewLocation === $location.path();
        };
    }
]);
