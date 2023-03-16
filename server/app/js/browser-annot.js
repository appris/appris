/**
 *
 * annotBrowser - AngularJS module for visualization of detailed annotations
 *
 */
var module = angular.module('appris.browser.annot', []);

/* CONTROLLERS */

module.controller('BrowserAnnotController', ['$scope', '$filter', 'Sequencer', function($scope, $filter, Sequencer) {
    // Obtain query from parent controller
    $scope.query      = angular.isDefined($scope.$parent.query) ? $scope.$parent.query : null;

    // Get residues data
    $scope.query.operation = "residues";
    var fResidues = Sequencer.get($scope.query);
    fResidues.then(function (data) {
        $scope.residues = data;
        $scope.methods =  $scope.$parent.methods;
    }, function(error) {
        $scope.error = error;
    });

}]);


/* DIRECTIVES */

module.directive('browserAnnotTpl', ['$compile', '$filter', 'consUrlFirestarligand', 'consUrlPDBstructure', 'consUrlPfamfamily', function($compile, $filter, consUrlFirestarligand, consUrlPDBstructure, consUrlPfamfamily) {

    function createAnnotReferences(residues, methods) {
        var elem = '';
        var charsTitle = 120;
        if ( angular.isArray(residues) ) {
            var accorHtml = '<accordion close-others="false">';
            angular.forEach(residues, function(item) {
                var id = item.id;
                var id_lbl = $filter('deleteSrcNames')(id);
                var body = '<tabset>';
                angular.forEach(item.methods, function(method) {
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
                        var annot = '-';
                        var pattern = '([^ ]*) \\(([^\\)]*)\\)';
                        if ( (mid === "functional_residue") && angular.isDefined(res_value.annot.match(pattern)) ) {
                            var match = res_value.annot.match(pattern);
                            if ( angular.isArray(match) && angular.isDefined(match[1]) && angular.isDefined(match[2]) ) {
                                annot = '';
                                angular.forEach(match[1].split(','), function(id,i) {
                                    var sc = match[2].split(',')[i];
                                    if ( id == "Cat_Site_Atl" ) {
                                        annot += id + " (" + sc + ") ";
                                    } else {
                                        annot += "<a href='"+consUrlFirestarligand+id+"' target='_blank'>" + id + "</a>(" + sc + ") ";
                                    }
                                });
                            }
                        }
                        else if ( (mid === "homologous_structure") && angular.isDefined(res_value.annot.match(pattern)) ) {
                            var match = res_value.annot.match(pattern);
                            if ( angular.isArray(match) && angular.isDefined(match[1]) && angular.isDefined(match[2]) ) {
                                var acc = match[1];
                                var eval = match[2];
                                var id = acc.replace(/_[A-Z0-9]{1,4}(?:_pdb_seq)?$/i,'');

                                const userKeyRegExp = /^AFM:/;
                                if ( userKeyRegExp.test(id) == true ) {
                                    id_clean = id.slice(4)
                                    annot = "<a href='"+ "https://alphafold.ebi.ac.uk/search/text/" + id_clean + "?suggested=true" +"' target='_blank'>"+acc+"</a>" + " (" + eval + ")";
                                }
                                else {
                                    annot = "<a href='"+consUrlPDBstructure+id+"' target='_blank'>"+acc+"</a>" + " (" + eval + ")";
                                }
                            }
                        }
                        else if ( (mid === "functional_domain") && angular.isDefined(res_value.annot.match(pattern)) ) {
                            var match = res_value.annot.match(pattern);
                            if ( angular.isArray(match) && angular.isDefined(match[1]) && angular.isDefined(match[2]) ) {
                                var acc = match[1];
                                var eval = match[2];
                                var id = acc.replace(/\.\d*$/g,'');
                                annot = "<a href='"+consUrlPfamfamily+id+"' target='_blank'>"+acc+"</a>" + " (" + eval + ")";
                            }
                        }
                        else {
                            annot = res_value.annot;
                        }

                        if ( angular.isDefined(res_value.annot) && !angular.isDefined(res_value.damaged)) {
                            cols += '<td>' + annot + '</td>';
                        }
                        else {
                            cols += '<td>' + annot + ' - damaged' +'</td>';
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
                        body += '<tab heading="'+methods[mid].label+'">'+table+'</tab>';
                    }
                });
                body += '</tabset>';
                var seqid = id_lbl;
                if ( id_lbl.length >= charsTitle) { seqid = id_lbl.substring(0,charsTitle-1) + '...' }
                accorHtml += '<accordion-group class="accord">';
                accorHtml += '<accordion-heading>'+
                    seqid +
                    '<i class="pull-right glyphicon" ng-class="{\'glyphicon-chevron-down\': browserTabAnnots, \'glyphicon-chevron-right\': !browserTabAnnots}"></i>'+
                    '</accordion-heading>';
                accorHtml += body;
                accorHtml += '</accordion-group>';
            });
            accorHtml += '</accordion>';
            elem = accorHtml;
        }
        else {
            elem = residues;
        }
        return elem;
    }

    return {
        restrict: 'E',
        scope: true,
        link: function(scope, element) {
            //scope.$watchGroup([''residues'], function(newValues, oldValues, scope) { // For AngularJS 1.3
            scope.$watchCollection('[residues, methods]', function(newValues, oldValues, scope) { // For AngularJS 1.1.4
                if ( angular.isDefined(newValues[0]) && angular.isDefined(newValues[1]) ) {
                    var annots = createAnnotReferences(newValues[0],newValues[1]);
                    var h = '<div class="browser-annot-tpl">' + annots + '</div>';
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
            element.html($compile('<div class="browser loading">Loading... If annotations fail to appear, please click on \'refresh browsers\'.</div>')(scope));
            // destroy element
            element.on('$destroy', function() {
            });
        }
    };
}]);
