/**
 *
 * Runner - AngularJS module for status report in APPRIS
 *
 */
var module = angular.module('appris.runner', []);

/* CONTROLLERS */

module.controller('RunnerController', ['consPathServerStatus', '$rootScope', '$scope', '$location', '$filter', '$cacheFactory', 'Runner', 'Species', 'Methods', 'Examples',
    function(consPathServerStatus, $rootScope, $scope, $location, $filter, $cacheFactory, Runner, Species, Methods, Examples) {

        // init form var
        $rootScope.isLoadingScreen = false;
        $scope.runnerForm = {};
        $scope.runnerForm.species = {};
        $scope.runnerForm.methods = [];
        $scope.runnerForm.sequences = null;
        $scope.runnerForm.ensembl = null;
        $scope.runnerForm.gene = null;

//        // load list of species
//        if ( angular.isUndefined($rootScope.species) ) {
//            $rootScope.species = null;
//            Species.get().$promise.then( function(data) {
//                $rootScope.species = data;
//            });
//        }
//
//        // load list of methods
//        if ( angular.isUndefined($rootScope.methods) || (angular.isArray($rootScope.methods) && $rootScope.methods.length == 0 ) ) {
//            $rootScope.methods = null;
//            Methods.get().$promise.then( function(data) {
//                $rootScope.methods = data;
//            });
//        }

        // load examples
        if ( angular.isUndefined($scope.examples) ) {
            $scope.examples = Examples.query();
        }

        // submit runnerForm: run a Job
        $scope.runJob = function (form) {

            // Validate form

            // init flag values
            $scope.speInvalid = false;
            $scope.genInvalid = false;
            $scope.seqInvalid = false;
            $scope.dupInvalid = false;
            $scope.metInvalid = false;
            $scope.emailInvalid = false;
            $scope.formInvalid = false;


            // empty specie input
            if ( angular.isUndefined(form.species) ) {
                $scope.speInvalid = true;
            }
            // empty data input
            if ( angular.isUndefined(form.ensembl) && angular.isUndefined(form.gene) && angular.isUndefined(form.sequences) ) {
                $scope.inpInvalid = true;
            }
            else if ( ( angular.isUndefined(form.ensembl) || angular.isUndefined(form.gene)) && angular.isUndefined(form.sequences) ) {
                $scope.genInvalid = true;
            }
            else if ( angular.isUndefined(form.sequences) && angular.isUndefined(form.ensembl) && angular.isUndefined(form.gene) ) {
                $scope.seqInvalid = true;
            }
            // duplicate data input
            if ( angular.isObject(form.sequences) && angular.isObject(form.gene) ) {
                $scope.dupInvalid = true;
            }
            // empty data met
            if ( angular.isObject(form.methods) && !( $filter('isSelectedMethod')(form.methods) ) ) {
                $scope.metInvalid = true;
            }
            // invalid email if exits
            if ( angular.isDefined(form.email) && (form.email === '') ) {
                $scope.emailInvalid = true;
            }

            // messages
            if ( $scope.speInvalid || $scope.seqInvalid || $scope.genInvalid || $scope.dupInvalid || $scope.metInvalid ) {
                $scope.formInvalid = true;
                return;
            }

            // create postData
            var postData = {};
            if (    angular.isDefined(form.species) && angular.isObject(form.species) && form.species &&
                angular.isDefined(form.sequences) && form.sequences &&
                $filter('isSelectedMethod')(form.methods) )
            {
                postData.species = form.species.scientific;
                postData.sequences = form.sequences;
                postData.methods = $filter('selectedMethods')(form.methods).join(',');
                if ( angular.isDefined(form.email) ) { postData.email = form.email }
            }
            else if (   angular.isDefined(form.species) && angular.isObject(form.species) && form.species &&
                angular.isDefined(form.ensembl) && form.ensembl &&
                angular.isDefined(form.gene) && form.gene &&
                $filter('isSelectedMethod')(form.methods) )
            {
                postData.species = form.species.scientific;
                postData.id = form.gene;
                postData.e_version = form.ensembl.dataset;
                postData.methods = $filter('selectedMethods')(form.methods).join(',');
                if ( angular.isDefined(form.email) ) { postData.email = form.email }
            }
            else {
                form.$invalid = true;
                return;
            }

            // Submit job
            var newRunner = Runner.create(postData);

            // Trigger validation flag.
            $rootScope.isLoadingScreen = true;
            $scope.submitted = true;

            newRunner.then(function(data) {
                if ( angular.isDefined(data.jobid) ) {
                    var jobId = data.jobid;
                    $rootScope.isLoadingScreen = false;
                    $location.url(consPathServerStatus+jobId);
                }
                else {
                    $location.url(consPathServerStatus);
                }
            }, function (err) {
                //console.log('Request failed: '+JSON.stringify(err));
                //console.log('Request failed: '+err);
                $rootScope.isLoadingScreen = false;
                $location.url(consPathServerStatus);
            });
        };

        // fill with an example
        $scope.exRunnerController = function (inEx) {
            $scope.runnerForm.species = null;
            $scope.runnerForm.sequences = null;
            $scope.runnerForm.ensembl = null;
            $scope.runnerForm.gene = null;
            if ( inEx == 'seq' ) {
                var ex = $scope.examples[0];
                var spe_id = ex.species_id;
                $scope.runnerForm.species = $scope.species[spe_id];
                $scope.runnerForm.sequences = ex.sequences;
            }
            else if ( inEx == 'id' ) {
                var ex = $scope.examples[1];
                var spe_id = ex.species_id;
                $scope.runnerForm.species = $scope.species[spe_id];
                $scope.runnerForm.ensembl = $scope.species[spe_id].assemblies[0];
                $scope.runnerForm.gene = ex.gene;
            }
        };

    }
]);