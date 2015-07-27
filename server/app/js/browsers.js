/**
 *
 * Browsers - AngularJS module that controls the browsers
 *
 * Created by jmrodriguez on 4/10/15.
 *
 */

var module = angular.module('appris.browsers', []);

module.controller('BrowsersController', ['$scope', '$routeParams', '$filter', function($scope, $routeParams, $filter) {

    // Extract route parameters (global input parameters)
    $scope.species  = angular.isDefined($routeParams.species) ? $routeParams.species : null;
    $scope.tyid     = angular.isDefined($routeParams.tid) ? $routeParams.tid : null;
    $scope.gid      = angular.isDefined($routeParams.id) ? $routeParams.id : null;
    $scope.jobid    = angular.isDefined($routeParams.jobid) ? $routeParams.jobid : null;
    $scope.db       = angular.isDefined($routeParams.db) ? $routeParams.db : null;

    $scope.submitBrowsers = function(inids, inmethods, inalign) {
        // Create query for the browser.
        var ids = '';
        angular.forEach($filter('orderBy')(inids, 'transcript_id'), function(item) {
            if ( item.selected ) {
                ids += item.transcript_id+',';
            }
        });
        ids = ids.replace(/,$/g,'');
        var methods = 'appris,';
        angular.forEach(inmethods, function(item) {
            methods += item.name+',';
        });
        methods = methods.replace(/,$/g,'');

        // REST queries
        $scope.query = null;
        if ( $scope.species && $scope.tyid && $scope.gid ) { // EXPORTER mode
            $scope.query = {
                species: $scope.species,
                tid:     $scope.tyid,
                id:      $scope.gid,
                db:      $scope.db,
                ids:     ids,
                methods: methods,
                align:   inalign,
                format:  "json"
            };
        }
        else if ( $scope.jobid ) { // RUNNER mode
            $scope.query = {
                jobid:   $scope.jobid,
                ids:     ids,
                methods: methods,
                align:   inalign,
                format:  "json"
            };
        }

        $scope.browsers = [{
            id:         'annotation',
            title:      'Annotation Browser',
            refresh:    false,
            info:       '#infoBrowserAnnot',
            template:   'templates/apprisBrowserAnnot.tpl.html'
        },{
            id:         'sequence',
            title:      'Sequence Browser',
            refresh:    false,
            info:       '#infoBrowserSeq',
            template:   'templates/apprisBrowserSeq.tpl.html'
        }];
        if ( !$scope.isSeqRunner ) {
            $scope.browsers.push({
                id:         'genome',
                title:      'Genome Browser',
                refresh:    false,
                info:       '#infoBrowserGen',
                template:   'templates/apprisBrowserGen.tpl.html'
            });
        }

        $scope.browsers.activeTab = function() {
            return $scope.browsers.filter(function(tab){
                return tab.active;
            })[0];
        };
    };

    //$scope.$watchGroup(['resultAnnots', 'currentMethods'], function(newValues, oldValues, scope) { // For AngularJS 1.3
    $scope.$watchCollection('[resultAnnots, currentMethods]', function(newValues, oldValues, scope) { // For AngularJS 1.1.4
        if ( angular.isDefined(newValues[0]) && angular.isDefined(newValues[1]) ) {
            if ( angular.isArray(newValues[0]) && newValues[0].length > 0 && angular.isArray(newValues[1]) && newValues[1].length > 0 ) {
                $scope.submitBrowsers(newValues[0], newValues[1], $scope.browserDataAlign);
            }
        }
    });

}]);


/* FILTERS */

// Checks if at least one sequence is selected
apprisFilters.filter('isSeqSelected', function() {
    return function(input) {
        var selected = null;
        if ( angular.isDefined(input) && input.length > 0 ) {
            selected = false;
            angular.forEach(input, function(item) {
                if ( item.selected ) {
                    return true;
                }
            });
        }
        return selected;
    };
});