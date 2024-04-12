'use strict';

/* Directives */

var apprisDirectives = angular.module('apprisDirectives', []);


apprisDirectives.directive('markdown', function() {
    var converter = new showdown.Converter();
    return {
        restrict: 'E',
        link: function(scope, element, attrs) {
            var htmlText = converter.makeHtml(element.text());
            element.html(htmlText);
        }
    }
});

/* TEMPLATES */

apprisDirectives.directive('navbarTopTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/navbarTop.tpl.html'
    };
});

apprisDirectives.directive('seekerTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/seeker.tpl.html'
    };
});

apprisDirectives.directive('toolsTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/tools.tpl.html'
    };
});

apprisDirectives.directive('toolsExtraTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/toolsExtra.tpl.html'
    };
});

apprisDirectives.directive('menuHelpTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/menuHelp.tpl.html'
    };
});

apprisDirectives.directive('alertTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/alert.tpl.html'
    };
});

apprisDirectives.directive('reportBrowsers', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/apprisBrowsers.tpl.html'
    };
});

/* PARTIAL PAGES */

apprisDirectives.directive('helpDatabase', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/database.html'
    };
});

apprisDirectives.directive('helpMethods', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/methods.html'
    };
});

apprisDirectives.directive('helpTrifid', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/trifid.html'
    };
});

apprisDirectives.directive('helpScores', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/scores.html'
    };
});

apprisDirectives.directive('helpImportantIsoforms', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/importantIsoforms.html'
    };
});

apprisDirectives.directive('helpApprisMane', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/apprisMane.html'
    };
});

apprisDirectives.directive('helpDbQuery', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/dbquery.html'
    };
});

apprisDirectives.directive('helpQueries', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/queries.html'
    };
});

apprisDirectives.directive('helpSeeker', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/seeker.html'
    };
});

apprisDirectives.directive('helpReport', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/report.html'
    };
});

apprisDirectives.directive('helpServer', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/server.html'
    };
});

apprisDirectives.directive('helpStatus', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/status.html'
    };
});

apprisDirectives.directive('helpServices', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/services.html'
    };
});

apprisDirectives.directive('helpSpecies', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/species.html'
    };
});

apprisDirectives.directive('helpAnnotations', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/annotations.html'
    };
});

apprisDirectives.directive('helpJobs', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/jobs.html'
    };
});

apprisDirectives.directive('helpFlags', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/flags.html'
    };
});

apprisDirectives.directive('helpIsoforms', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/isoforms.html'
    };
});

apprisDirectives.directive('helpBrowsers', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/browsers.html'
    };
});

apprisDirectives.directive('helpBrowserAnnot', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/browserAnnot.html'
    };
});

apprisDirectives.directive('helpBrowserSeq', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/browserSeq.html'
    };
});

apprisDirectives.directive('helpBrowserGen', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/browserGen.html'
    };
});

apprisDirectives.directive('helpBrowserForm', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/browserForm.html'
    };
});

apprisDirectives.directive('speciesList', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/species/speciesList.html'
    };
});

apprisDirectives.directive('seekerInfoGeneTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/seeker/infoGene.html'
    };
});

apprisDirectives.directive('seekerInfoInputsTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/seeker/infoInputs.html'
    };
});

apprisDirectives.directive('runnerInfoSpeciesTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/runner/infoSpecies.html'
    };
});

apprisDirectives.directive('runnerInfoInputsTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/runner/infoInputs.html'
    };
});

apprisDirectives.directive('runnerInfoMethodsTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/runner/infoMethods.html'
    };
});

apprisDirectives.directive('runnerInfoOptionalTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/runner/infoOptional.html'
    };
});

apprisDirectives.directive('resultInfoIsofTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/result/infoIsof.html'
    };
});
apprisDirectives.directive('resultInfoAnnotTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/result/infoAnnot.html'
    };
});

apprisDirectives.directive('resultInfoBrowserTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/result/infoBrowser.html'
    };
});

apprisDirectives.directive('resultInfoBrowserAnnotTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/result/infoBrowserAnnot.html'
    };
});

apprisDirectives.directive('resultInfoBrowserSeqTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/result/infoBrowserSeq.html'
    };
});

apprisDirectives.directive('resultInfoBrowserGenTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/result/infoBrowserGen.html'
    };
});

apprisDirectives.directive('helpDocker', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/docker.html'
    };
});

apprisDirectives.directive('helpPrevious', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/previous.html'
    };
});

apprisDirectives.directive('changeLogs', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/help/help/changelogs.html'
    };
});

apprisDirectives.directive('license', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/license.html'
    };
});

apprisDirectives.directive('menuImptTpl', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'templates/menuImpt.tpl.html'
    };
});

apprisDirectives.directive('tissueSpecific', function() {
    return {
        restrict: 'E',
        replace: true,
        templateUrl: 'partials/impt/tissues.html'
    };
});


/* CLEVER DIRECTIVES */

apprisDirectives.directive('validateRunnerForm',function(){
    return {
        require: "ngModel",
        link: function(scope, elm, attrs, Controller){
            var regex = /^\d$/;
            var validator = function(value) {
                //var floatValue = parseFloat(value);
                Controller.$setValidity('validNumber', regex.test(value));
                return value;
            };
            Controller.$parsers.unshift(validator);
            Controller.$formatters.unshift(validator);
        }
    };
});

apprisDirectives.directive('loadingScreen', function () {
    return {
        restrict: 'E',
        replace: true,
        template: '<div id="loading-screen" data-ng-show="isLoadingScreen"><div class="loading-text">Loading...</div></div>'
//        link: function (scope, element) {
//            scope.$watch('isLoadingScreen', function (val) {
//                if (val) {
//                    //$(element).show();
//                    scope.loading = true;
//                }
//                else {
//                    //$(element).hide();
//                    scope.loading = false;
//                }
//            });
//        }
    }
});

//apprisDirectives.directive('butterbar', ['$rootScope', function($rootScope) {
//    return {
//        link: function(scope, element, attrs) {
//            element.addClass('hide');
//            $rootScope.$on('$routeChangeStart', function() {
//                element.removeClass('hide');
//            });
//            $rootScope.$on('$routeChangeSuccess', function() {
//                element.addClass('hide');
//            });
//        }
//    };
//}]);

apprisDirectives.directive('spinJs', function() {
    return {
        restrict: 'AEC',
        replace: true,
        link: function(scope, element, attrs) {
            var spinner = new Spinner({
                className: 'spinner',
                lines: parseInt(attrs.lines, 10),
                length: parseInt(attrs.length,10),
                width: parseInt(attrs.width, 10),
                radius: parseInt(attrs.radius,10),
                position: attrs.position,
                top: attrs.top,
                left: attrs.left
            }).spin();
            element.append(spinner.el);
        }
    };
});

apprisDirectives.directive('focus', function() {
    return {
        link: function(scope, element, attrs) {
                element[0].focus();
        }
    };
});

apprisDirectives.directive('browserFrame', function ( $compile ) {
    return {
        restrict: 'AEC',
        scope: true,
        link: function ( scope, element, attrs) {
            scope.$watch('browserTab.frame', function (value) {
                if ( angular.isDefined(value) && value !== null ) {
                    if ( attrs.browserId === 'annotation' ) {
                        if ( angular.isArray(value) ) {
                            var accorHtml = '<accordion close-others="false">';
                            angular.forEach(value, function(item) {
                                var id = item.id;
                                var methods = item.methods;
                                var body = '<tabset>';
                                angular.forEach(methods, function(method) {
                                    var mid = method.id;
                                    var residues = method.residues;
                                    var isPos = true;
                                    var rows = '';
                                    angular.forEach(residues, function(res_value) {
                                        var cols = '';
                                        if ( angular.isDefined(res_value.start) && !angular.isDefined(res_value.end) ) {
                                            cols += '<td colspan="2">' + res_value.start + '</td>';
                                        }
                                        else {
                                            cols += '<td>' + res_value.start + '</td>';
                                        }
                                        if ( angular.isDefined(res_value.end) ) {
                                            isPos = false;
                                            cols += '<td>' + res_value.end + '</td>';
                                        }
                                        if ( angular.isDefined(res_value.annot) && !angular.isDefined(res_value.damaged)) {
                                            cols += '<td>' + res_value.annot + '</td>';
                                        }
                                        else {
                                            cols += '<td>' + res_value.annot + ' - damaged' +'</td>';
                                        }
                                        if ( cols !== '' ) { rows += '<tr>'+cols+'</tr>'; }

                                    });
                                    if ( rows != '' ) {
                                        var title = '';
                                        if ( isPos ) {
                                            title += '<tr><th colspan="2">residue</th><th>' + method.title + '</th></tr>';
                                        }
                                        else {
                                            title += '<tr><th>start</th><th>end</th><th>' + method.title + '</th></tr>';
                                        }
                                        var table = '<table class="annot-browser table table-bordered table-condensed">'+'<tbody>'+title+rows+'</tbody></table>';
                                        body += '<tab heading="'+scope.methods[mid].label+'">'+table+'</tab>';
                                    }
                                });
                                body += '</tabset>';
                                accorHtml += '<accordion-group class="accord">';
                                accorHtml += '<accordion-heading>'+
                                    id +
                                    '<i class="pull-right glyphicon" ng-class="{\'glyphicon-chevron-down\': browserTabAnnots, \'glyphicon-chevron-right\': !browserTabAnnots}"></i>'+
                                    '</accordion-heading>';
                                accorHtml += body;
                                accorHtml += '</accordion-group>';
                            });
                            accorHtml += '</accordion>';
                            element.html($compile('<div class="browser tpl annot">'+accorHtml+'</div>')(scope));
                        }
                        else {
                            element.html($compile('<div class="browser tpl annot">'+value+'</div>')(scope));
                        }
                    }
                    else if ( attrs.browserId === 'sequence' ) {
                        //element.html($compile('<div class="browser tpl seq">'+value+'</div>')(scope));
                        element.html($compile('<div class="browser tpl seq"><pre>'+value+'</pre></div>')(scope));
                    }
                    else if ( attrs.browserId === 'genome' ) {
                        element.html($compile('<div class="browser tpl gen">'+value+'</div>')(scope));
                    }
                }
                else { // loading step
                    element.html($compile('<div class="browser tpl">Loading...</div>')(scope));
                }
            });
        }
    };
});

apprisDirectives.directive('dropdownchoice', function() {
    return {
        link: function(scope, element) {
            element.on("click", function(e){
                if (e.target.id !== 'close' ) {
                    e.stopPropagation();
                }
            });
        }
    };
});
