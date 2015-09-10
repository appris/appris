'use strict';

/* App Module */

var apprisApp = angular.module('apprisApp', [

    'ngRoute',
    'ngAnimate',
    'ngSanitize',

    'apprisControllers',
    'apprisFilters',
    'apprisDirectives',
    'apprisServices',
    'apprisFilters',

    'appris.runner',
    'appris.status',
    'appris.report',
    'appris.browsers',
    'appris.browser.seq',
    'appris.browser.annot',
    'appris.browser.gen',

    'mgcrea.ngStrap.popover',
    'ui.bootstrap',

    'pagination',
    'checklist',
    'download'
]);

apprisApp.config(['$provide', function ($provide) {
    // SERVER INFO
/*
    // GOLD server - production
    $provide.value("serverName", '{APPRIS}');
    $provide.value("serverVersion", '2015_09.v10');
    $provide.value("serverType", 'gold');
    $provide.value("serverHost", 'http://appris.bioinfo.cnio.es');
    $provide.value("serverHostWS", 'http://apprisws.bioinfo.cnio.es');
    // BETA server - development
    $provide.value("serverName", '{APPRISdev}');
    $provide.value("serverVersion", '2015_09.v10');
    $provide.value("serverType", 'beta');
    $provide.value("serverHost", 'http://appris-dev.bioinfo.cnio.es');
    $provide.value("serverHostWS", 'http://apprisws-dev.bioinfo.cnio.es');
*/
    // ALPHA server - development in local
    $provide.value("serverName", '{APPRISloc}');
    $provide.value("serverVersion", '2015_09.v10');
    $provide.value("serverType", 'alpha');
    $provide.value("serverHost", 'http://local.es:3000');
    $provide.value("serverHostWS", 'http://local.es/ws');

    // CONSTANTS
    $provide.value("consUrlEnsembl", 'http://www.ensembl.org');
    $provide.value("consUrlFirestarligand", 'http://firedb.bioinfo.cnio.es/Php/ligand/index.html?id=');
    $provide.value("consUrlPDBligand", 'http://www.rcsb.org/pdb/ligand/ligandsummary.do?hetId=');
    $provide.value("consUrlPDBstructure", 'http://www.rcsb.org/pdb/explore/explore.do?structureId=');
    $provide.value("consUrlPfamfamily", 'http://pfam.xfam.org/family/');
    // constant paths
    $provide.value("consPageDatabase", '/database');
    $provide.value("consPageServer", '/server');
    $provide.value("consPathServerStatus", '/server/status/');
    $provide.value("consPathServerResult", '/server/result/');
    $provide.value("consPathSeeker", '/seeker/');
    $provide.value("consQueryNotMatch", '/query_not_match');
    $provide.value("consPageError", '/error');
    $provide.value("consPageNotFound", '/page_not_found');
}]);

apprisApp.config(['$compileProvider', function ($compileProvider) {
    $compileProvider.aHrefSanitizationWhitelist(/^\s*(https?|ftp|mailto|blob):/);
}]);


// Safely Prevent Template Caching in AngularJS
// AngularJS's $templateCache can be a pain in the ass. Sometimes we don't want templates to be cached. A quick Internet search to disable caching gives the following workaround:
// A solution is to tweak the above workaround so that new cache entries are removed on route change instead of indiscriminately removing all entries:
//app.run(function($rootScope, $templateCache) {
//    $rootScope.$on('$viewContentLoaded', function() {
//        $templateCache.removeAll();
//    });
//});
apprisApp.run(function($rootScope, $templateCache) {
    $rootScope.$on('$routeChangeStart', function(event, next, current) {
        if (typeof(current) !== 'undefined'){
            $templateCache.remove(current.templateUrl);
        }
    });
});

apprisApp.config(['$routeProvider', function ($routeProvider) {
//        console.log($httpProvider.defaults.headers);
//        $httpProvider.defaults.headers.post['Content-Type'] = 'application/json; charset=utf-8';
//        $httpProvider.defaults.headers.post['Accept'] = 'application/json, text/javascript';
//        $httpProvider.defaults.headers.post['Access-Control-Max-Age'] = '1728000';
//        $httpProvider.defaults.headers.common['Access-Control-Max-Age'] = '1728000';
//        $httpProvider.defaults.headers.common['Accept'] = 'application/json, text/javascript';
//        $httpProvider.defaults.headers.common['Content-Type'] = 'application/json; charset=utf-8';
//        $httpProvider.defaults.useXDomain = true;
//        delete $httpProvider.defaults.headers.common['X-Requested-With'];
//        console.log($httpProvider.defaults.headers);


//        $locationProvider.html5Mode(true);
//        $locationProvider.html5Mode('requireBase:true');

//        //initialize get if not there
//        if (!$httpProvider.defaults.headers.get) {
//            $httpProvider.defaults.headers.get = {};
//        }
//        //disable IE ajax request caching
//        $httpProvider.defaults.headers.get['If-Modified-Since'] = '0';

        $routeProvider.
            when('/species', {
//                controller: 'SpeciesController',
                templateUrl: 'partials/species.html'
            }).
            when('/server', {
                controller: 'RunnerController',
                templateUrl: 'partials/runner.html'
            }).
            when('/server/status/:jobid?', {
                controller: 'StatusController',
                templateUrl: 'partials/status.html'
            }).
            when('/server/result/:jobid', {
                //controller: 'ResultController',
                controller: 'ReportController',
                templateUrl: 'partials/report.html'
            }).
//            when('/seeker', {
//                redirectTo: '/exporter'
//            }).
            when('/database', {
//                controller: 'ExporterController',
                templateUrl: 'partials/database.html'
            }).
            when('/seeker/:id?', {
                controller: 'SeekerResultController',
                templateUrl: 'partials/seeker.html'
            }).
            when('/database/:tid/:species/:id/:methods?', {
                //controller: 'ResultController',
                controller: 'ReportController',
                templateUrl: 'partials/report.html'
            }).
            when('/tools', {
                templateUrl: 'partials/tools.html'
            }).
            when('/downloads', {
                templateUrl: 'partials/downloads.html'
            }).
            when('/help/:help*', {
                controller: 'HelpController',
                templateUrl: 'partials/help.html'
            }).
            when('/publications', {
                templateUrl: 'partials/publications.html'
            }).
            when('/about', {
                controller: 'AboutController',
                templateUrl: 'partials/about.html'
            }).
            when('/changelogs', {
                templateUrl: 'partials/changelogs.html'
            }).
            when('/', {
//                controller: 'ApprisController',
                templateUrl: 'partials/frontpage.html'
            }).
            when('/query_not_match', {
                templateUrl: 'templates/query_not_match.html'
            }).
            when('/page_not_found', {
                templateUrl: 'templates/page_not_found.html'
            }).
            when('/error', {
                templateUrl: 'templates/error.html'
            }).
            otherwise({
                templateUrl: 'templates/page_not_found.html'
            });
  }]);

