/**
 *
 * Seeker - AngularJS module for seeker in APPRIS
 *
 */
var module = angular.module('appris.seeker', []);

/* CONTROLLERS */

module.controller('SeekerSimpleController', ['consPathSeeker', '$scope', '$location',
    function(consPathSeeker, $scope, $location) {

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
            $location.url(consPathSeeker+form.inQuery);
        };
    }
]);

module.controller('SeekerAdvancedController', ['consPathSeeker', '$scope', '$location', '$filter', 'Examples',
    function(consPathSeeker, $scope, $location, $filter, Examples) {

        // load examples
        if ( angular.isUndefined($scope.examples) ) {
            $scope.examples = Examples.query();
        }
        // submit search
        $scope.runJob = function (form) {
            var optQuery = '';
            // Validate
            // empty specie input
            if ( angular.isUndefined(form.gene) && (form.gene != '')) {
                $scope.genInvalid = true;
            }
            else {
                $scope.genInvalid = false;
            }
            // messages
            if ( $scope.genInvalid ) {
                return;
            }
            // create query
            var query = form.gene;
            if ( angular.isDefined(form.species) && angular.isObject(form.species) && form.species ) { optQuery += 'sp='+$filter('getSpeciesId')(form.species)+'&' }
            if ( angular.isDefined(form.dataset) && angular.isObject(form.dataset) && form.dataset ) { optQuery += 'ds='+form.dataset.id+'&' }
            if ( optQuery != '' ) {
                query +=  '?' + optQuery.replace(/&$/g,'');
            }
            // Redirect URL
            $location.url(consPathSeeker+query);
        };

        // fill with an example
        $scope.exSeekerController = function (inEx) {
            var ex = null;
            $scope.seekerAdvancedForm.species = null;
            $scope.seekerAdvancedForm.dataset = null;
            $scope.seekerAdvancedForm.gene = null;
            if ( inEx == 'gene_name_1' ) {
                ex = $scope.examples[2];
            }
            else if ( inEx == 'gene_name_2' ) {
                ex = $scope.examples[3];
            }
            else if ( inEx == 'gene_id_1' ) {
                ex = $scope.examples[4];
            }
            else if ( inEx == 'gene_id_2' ) {
                ex = $scope.examples[5];
            }
            else if ( inEx == 'gene_id_3' ) {
                ex = $scope.examples[6];
            }
            if ( angular.isObject(ex) ) {
                var sp_id = ex.species_id;
                var sc = ex.source;
                var as = ex.assembly;
                $scope.seekerAdvancedForm.gene = ex.gene;
                $scope.seekerAdvancedForm.species = $scope.species[sp_id];
                angular.forEach($scope.species[sp_id].assemblies, function(assembly) {
                    if ( assembly.id == as ) {
                        angular.forEach(assembly.datasets, function(dataset) {
                            if ( dataset.source.name == sc ) {
                                $scope.seekerAdvancedForm.dataset = dataset;
                            }
                        });
                    }
                });
            }
        };

    }
]);

module.controller('SeekerResultController', ['consQueryNotMatch', '$rootScope', '$scope', '$routeParams', '$location', '$filter', 'Seeker',
    function(consQueryNotMatch, $rootScope, $scope, $routeParams, $location, $filter, Seeker) {

        // init vars
        $rootScope.isLoadingScreen = true;
        var queryData = {
            id: $routeParams.id,
            sp: $routeParams.sp,
            ds: $routeParams.ds
        };
        $scope.limitString = 20;

        // create search summary
        var specieResult = {};
        var rawResults = Seeker.get(queryData);
        rawResults.then(function (data) {
            angular.forEach(data.match, function(item) {
                if (  angular.isUndefined(specieResult[item.species]) ) {
                    specieResult[item.species] = [];
                }
                var position = '-';
                if ( item.chr !== null && item.start !== null && item.end !== null ) {
                    position = item.chr+':'+item.start+'-'+item.end;
                }
                var dblink = {
                    External_Id:      '-',
                    Ensembl_Gene_Id:  '-',
                    Refseq_Gene_Id:   '-',
                    Uniprot_Gene_Id:  '-'
                };
                angular.forEach(item.dblink, function(item) {
                    if ( angular.isDefined(item.namespace) && item.namespace == 'External_Id' ) { dblink.External_Id += item.id+',' }
                    if ( angular.isDefined(item.namespace) && item.namespace == 'Ensembl_Gene_Id' ) { dblink.Ensembl_Gene_Id += item.id+',' }
                    if ( angular.isDefined(item.namespace) && item.namespace == 'Refseq_Gene_Id' ) { dblink.Refseq_Gene_Id += item.id+',' }
                    if ( angular.isDefined(item.namespace) && item.namespace == 'Uniprot_Gene_Id' ) { dblink.Uniprot_Gene_Id += item.id+',' }
                });
                if ( dblink.External_Id != '-' ) {
                    dblink.External_Id = dblink.External_Id.replace(/,$/g,'');
                    dblink.External_Id = dblink.External_Id.replace(/^-/g,'');
                }
                if ( dblink.Ensembl_Gene_Id != '-' ) {
                    dblink.Ensembl_Gene_Id = dblink.Ensembl_Gene_Id.replace(/,$/g,'');
                    dblink.Ensembl_Gene_Id = dblink.Ensembl_Gene_Id.replace(/^-/g,'');
                }
                if ( dblink.Refseq_Gene_Id != '-' ) {
                    dblink.Refseq_Gene_Id = dblink.Refseq_Gene_Id.replace(/,$/g,'');
                    dblink.Refseq_Gene_Id = dblink.Refseq_Gene_Id.replace(/^-/g,'');
                }
                if ( dblink.Uniprot_Gene_Id != '-' ) {
                    dblink.Uniprot_Gene_Id = dblink.Uniprot_Gene_Id.replace(/,$/g,'');
                    dblink.Uniprot_Gene_Id = dblink.Uniprot_Gene_Id.replace(/^-/g,'');
                }
                var speRst = {
                    "species":      item.species,
                    "assembly": {
                        "id": $filter('filter')($rootScope.species[item.species].assemblies, { "id": item.assembly })[0].id,
                        "name": $filter('filter')($rootScope.species[item.species].assemblies, { "id": item.assembly })[0].name
                    },
                    "source":       item.source,
                    "dataset":      item.dataset.replace(/\.([^$]*)$/g,''),
                    "label":        item.label,
                    "namespace":    item.namespace,
                    "id":           item.id,
                    "biotype":      item.biotype,
                    "dblink":       dblink,
                    "position":     position
                }
                specieResult[item.species].push(speRst);
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
                $location.url(consQueryNotMatch);
            }
            else if ( error.status >= 500 ) {
                $location.url(consQueryNotMatch);
            }
        });
    }
]);
