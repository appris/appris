/**
 *
 * genBrowser - AngularJS module for visualization of image of UCSC Genome Browser
 *
 */
var module = angular.module('appris.browser.gen', []);

/* CONTROLLERS */

module.controller('BrowserGenController', ['$scope', 'Sequencer', function($scope, Sequencer) {
    // Obtain query from parent controller
    $scope.query      = angular.isDefined($scope.$parent.query) ? $scope.$parent.query : null;

    // Get residues data
    $scope.query.operation = "genome";
    $scope.query.format = "html";
    var fGenome = Sequencer.get($scope.query);
    //var fGenome = Viewer.get($scope.query);
    fGenome.then(function (data) {
        $scope.genome =  data;
    }, function(error) {
        $scope.error = error;
    });

}]);


/* DIRECTIVES */

module.directive('browserGenTpl', ['$compile', function($compile) {

    return {
        restrict: 'E',
        scope: true,
        link: function(scope, element) {
            scope.$watch('genome', function(newValues) {
                if ( angular.isDefined(newValues) ) {
                    var h = '<div class="browser-gen-tpl">' + newValues + '</div>';
                    element.html($compile(h)(scope));
                    scope.browser.refresh = false;
                }
            });
            scope.$watch('error', function(newValue, oldValue, scope) {
                if ( angular.isDefined(newValue) ) {
                    var annots = newValue;
                    element.html($compile('<div class="browser error">' + annots + '</div>')(scope));
                    scope.browser.refresh = false;
                }
            });
            // by default
            scope.browser.refresh = true;
            element.html($compile('<div class="browser loading">Loading... If genomes fail to appear, please click on \'refresh browsers\'.</div>')(scope));
            // destroy element
            element.on('$destroy', function() {
            });
        }
    };
}]);
