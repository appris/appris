/**
 *
 * seqBrowser - AngularJS module for visualization of detailed annotations mapped onto the amino acid sequences
 *
 */
var module = angular.module('appris.browser.seq', []);

/* CONTROLLERS */

module.controller('BrowserSeqController', ['serverHostWS', 'consUrlFirestarligand', 'consUrlPDBstructure', 'consUrlPfamfamily', '$scope', '$filter', '$popover', 'Sequencer', function(serverHostWS, consUrlFirestarligand, consUrlPDBstructure, consUrlPfamfamily, $scope, $filter, $popover, Sequencer) {
    // Obtain query from parent controller
    $scope.query      = angular.isDefined($scope.$parent.query)   ? $scope.$parent.query : null;

    // Get reference sequences: sequences or align
    if ( angular.isDefined($scope.query.align) && $scope.query.align ) { $scope.query.operation = "align" }
    else { $scope.query.operation = "sequences" }
    var fSequences = Sequencer.get($scope.query);
    fSequences.then(function (data) {
        $scope.sequences = data;
    }, function(error) {
        $scope.error = error;
    });

    // Get cds data
    $scope.query.operation = "cds";
    var fCDS = Sequencer.get($scope.query);
    fCDS.then(function (data) {
        var report = $filter('filterCDSPerResidues')(data);
        $scope.resmetseqs = report.residues;
        $scope.metsequences = report.methods;
        $scope.checkmetseqs = [];
    }, function(error) {
        $scope.error = error;
    });

    // Get residues data
    $scope.query.operation = "residues";
    var fResidues = Sequencer.get($scope.query);
    fResidues.then(function (data) {
        var report = $filter('filterAnnotPerResidues')(data, $scope.$parent.methods);
        $scope.residues = report.residues;
        $scope.methods =  report.methods;
        $scope.checkmethods = [];
    }, function(error) {
        $scope.error = error;
    });

    // Get methods from given residue
    $scope.getMethodsFromRes = function(tid,ridx,mname) {
        $scope.query.operation = "residues";
        $scope.query.ids = tid;
        $scope.query.r = ridx;
        //$scope.query.methods = filterMethodsFromRes($scope.methods, mname);
        $scope.query.methods = mname;
        if ( angular.isDefined($scope.query.methods) && ($scope.query.methods != '') ) {
            var fResidues = Sequencer.get($scope.query);
            fResidues.then(function (data) {
                var id = tid + '_' + ridx;
                id = id.replace('.',"\\.");
                var e = angular.element( document.getElementById(id) );
                var asAServiceOptions = {
                    title:     'APPRIS annotations',
                    template:  'templates/popover.tpl.html',
                    placement: 'bottom',
                    autoClose: true,
                    trigger:   'manual'
                };
                var urlExporterWS = serverHostWS+'/rest/exporter';
                var urlRunnerResultWS = serverHostWS+'/rest/runner/result';
                    asAServiceOptions.content = createAnnotReferences($filter, $scope.query, data, urlExporterWS, urlRunnerResultWS, consUrlFirestarligand, consUrlPDBstructure, consUrlPfamfamily);
                var myPopover = $popover(e, asAServiceOptions);
                myPopover.$promise.then(myPopover.toggle);
            }, function(error) {
                $scope.error = error;
            });
        }
    };

}]);


/* DIRECTIVES */

module.directive('browserSeqTpl', ['$compile', '$filter', function($compile, $filter) {

    var resPerLine = 70;
    var charsTitle = 40;
    // NOTE: BE CAREFUL HARD-CORE for CORSAIR!!!
    function createSeqLegend(query) {
        var elem1 = '<h4>Annotations</h4>'+
                    '<span class="{{method.class}}" data-ng-repeat="(morder, method) in methods" data-ng-if=" (method.name != \'corsair\' && method.name != \'appris\') " >'+
                    '<input type="checkbox" id="checkbox-{{method.id}}" class="pseudo-checkbox sr-only" data-ng-model="checkmethods[method.id]" data-ng-change="printSeq(sequences,residues,checkmethods,resmetseqs,checkmetseqs)" >'+
                    '<label for="checkbox-{{method.id}}" class="pointer fancy-checkbox-label">{{method.desc}}</label>'+
                    '</span>';
        var elem2 = '<h4>Sequences</h4>'+
            '<span class="{{metseq.class}}" data-ng-repeat="metseq in metsequences">'+
            '<input type="checkbox" id="checkbox-{{metseq.id}}" class="pseudo-checkbox sr-only" data-ng-model="checkmetseqs[metseq.id]" data-ng-change="printSeq(sequences,residues,checkmethods,resmetseqs,checkmetseqs)" >'+
            '<label for="checkbox-{{metseq.id}}" class="pointer fancy-checkbox-label">{{metseq.desc}}</label>'+
            '</span>';
        var elem = '<div class="browser-seq-legend">';
        elem += '<h4>Highlight</h4>'+ elem1;
        if ( angular.isDefined(query.sc) && query.sc !== 'uniprot' && query.sc !== 'appris' ) {
            elem += elem2;
        }
        elem += '</div>';
        return elem;
    }

    function createSeqReferences(sequences, residues, methods, resmetseqs, metseqs) {
        var elem = '';
        if ( angular.isDefined(sequences.seqs) ) {
            elem = createSeqRef(sequences, residues, methods, resmetseqs, metseqs);
        }
        else if ( angular.isDefined(sequences.aln) ) {
            elem = createAlnRef(sequences, residues, methods, resmetseqs, metseqs);
        }
        return elem;
    }

    function splitSeqRef(sequences) {
        var split = [];
        angular.forEach(sequences, function(sequence) {
            var seq = [];
            var len = sequence.seq.length;
            var iseq = 0;
            angular.forEach(sequence.seq.match( new RegExp('.{0,'+resPerLine+'}', 'g')  ), function(seqline) {
                seq.push({
                    iseq: iseq,
                    seq: seqline
                });
                angular.forEach(seqline.split(''), function(s) {
                    if ( s !== '-' ) { iseq += 1 }
                });
            });
            split.push({
                id:   sequence.id,
                seqs: seq,
                len:  len
            });
        });
        return split;
    }

    function splitMatchAlnRef(sequences) {
        var split = [];
        var iseq = 0;
        angular.forEach(sequences.match( new RegExp('.{0,'+resPerLine+'}', 'g')  ), function(seqline) {
            split.push({
                iseq: iseq,
                seq: seqline
            });
            angular.forEach(seqline.split(''), function(s) {
                if ( s !== '-' ) { iseq += 1 }
            });
        });
        return split;
    }

    function createSeqRef(sequences, residues, methods, resmetseqs, metseqs) {
        var elem = '<browser-seq>';
        var splitseqs = splitSeqRef(sequences.seqs);
        if ( angular.isArray(splitseqs) && angular.isArray(splitseqs[0].seqs) ) {
            angular.forEach(splitseqs, function(spseq) {
                var tidx = spseq.id;
                var tidx_lbl = $filter('deleteSrcNames')(tidx);
                var splen = spseq.len;
                var spseqs = spseq.seqs;
                var seqline = '';
                angular.forEach(spseqs, function(spseq) {
                    if ( angular.isDefined(spseq.seq) && (spseq.seq !== '') ) {
                        var seq = spseq.seq;
                        var iseq = spseq.iseq;
                        angular.forEach(seq.split(''), function(s) {
                            if ( s !== '-' ) {
                                iseq += 1;
                                var mcls = '';
                                var mname = '';
                                if ( angular.isDefined(residues[tidx]) && angular.isDefined(residues[tidx][iseq]) && angular.isDefined(methods) ) {
                                    angular.forEach(residues[tidx][iseq], function(mres) {
                                        if ( methods[mres.id] ) {
                                            mcls += ' ' + mres.class + ' ';
                                            mname += mres.name + ',';
                                        }
                                    });
                                }
                                if ( angular.isDefined(resmetseqs[tidx]) && angular.isDefined(resmetseqs[tidx][iseq]) && angular.isDefined(metseqs) ) {
                                    var msres = resmetseqs[tidx][iseq]
                                    if ( metseqs[msres.id] ) {
                                        mcls += ' ' + msres.class + ' ';
                                        mname += msres.name + ',';
                                    }
                                }
                                if ( mcls != '' ) {
                                    var ridx = iseq;
                                    seqline += '<browser-seq-res id="' + tidx.replace('.',"\\.") + '_' + ridx + '" ridx="' + ridx + '" class="' + mcls + '" data-ng-click="getMethodsFromRes(\'' + tidx.replace('.',"\\.") + '\', \'' + ridx + '\', \'' + mname + '\' )">' + s + '</browser-seq-res>';
                                }
                                else {
                                    seqline += s;
                                }
                            }
                            else {
                                seqline += s;
                            }
                        });
                        seqline += '<br>';
                    }
                });
                var seqid = tidx_lbl;
                if ( tidx_lbl.length >= charsTitle) { seqid = tidx_lbl.substring(0,charsTitle-1) + '...' }
                elem += '<browser-seq-line '+tidx+'>';
                elem += '>' + seqid + '|' + splen +'<br>';
                elem += seqline;
                elem += '<br>';
                elem += '</browser-seq-line>';
            });
            elem += '<br>';
        }
        elem += '</browser-seq>';
        return elem;
    }

    function createAlnRef(sequences, residues, methods, resmetseqs, metseqs) {
        var elem = '<browser-seq>';
        var splitseqs = splitSeqRef(sequences.aln);
        var splitmatch = splitMatchAlnRef(sequences.match);
        if ( angular.isArray(splitseqs) && angular.isArray(splitseqs[0].seqs) ) {
            elem += '<browser-seq-line>';
            elem += Array(charsTitle+1).join(" ");
            for ( var i = 1; i <= resPerLine; i++ ) {
                if ( i === 1 ) { elem += i }
                else if ( i === resPerLine ) { elem += i }
                else if ( (i  % 10) === 0 ) { elem += '|' } else { elem += '.' }
            }
            elem += '<br>';
            elem += '</browser-seq-line>';
            var nlines = splitseqs[0].seqs.length;
            for ( var i = 0; i < nlines; i++ ) {
                angular.forEach(splitseqs, function(spseq) {
                    var tidx = spseq.id;
                    var tidx_lbl = $filter('deleteSrcNames')(tidx);
                    var spseqs = spseq.seqs;
                    if ( angular.isDefined(spseqs[i]) ) {
                        var spseq = spseqs[i];
                        if ( angular.isDefined(spseq.seq) && (spseq.seq !== '') ) {
                            var seq = spseq.seq;
                            var iseq = spseq.iseq;
                            var seqline = '';
                            angular.forEach(seq.split(''), function(s) {
                                if ( s !== '-' ) {
                                    iseq += 1;
                                    var mcls = '';
                                    var mname = '';
                                    if ( angular.isDefined(residues[tidx]) && angular.isDefined(residues[tidx][iseq]) && angular.isDefined(methods) ) {
                                        angular.forEach(residues[tidx][iseq], function(mres) {
                                            if ( methods[mres.id] ) {
                                                mcls += ' ' + mres.class + ' ';
                                                mname += mres.name + ',';
                                            }
                                        });
                                    }
                                    if ( angular.isDefined(resmetseqs[tidx]) && angular.isDefined(resmetseqs[tidx][iseq]) && angular.isDefined(metseqs) ) {
                                        var msres = resmetseqs[tidx][iseq]
                                        if ( metseqs[msres.id] ) {
                                            mcls += ' ' + msres.class + ' ';
                                            mname += msres.name + ',';
                                        }
                                    }
                                    if ( mcls != '' ) {
                                        var ridx = iseq;
                                        seqline += '<browser-seq-res id="' + tidx.replace('.',"\\.") + '_' + ridx + '" ridx="' + ridx + '" class="' + mcls + '" data-ng-click="getMethodsFromRes(\'' + tidx.replace('.',"\\.") + '\', \'' + ridx + '\', \'' + mname + '\' )">' + s + '</browser-seq-res>';
                                    }
                                    else {
                                        seqline += s;
                                    }
                                }
                                else {
                                    seqline += s;
                                }
                            });
                            var seqid = '';
                            if ( tidx_lbl.length >= charsTitle) { seqid = tidx_lbl.substring(0,charsTitle-4) + '...'}
                            else { seqid = tidx_lbl + Array(charsTitle - tidx_lbl.length ).join(" ") }
                            elem += '<browser-seq-line tidx="'+tidx+'">' + seqid + ' ' + seqline + '<br>' + '</browser-seq-line>';
                        }
                    }
                });
                var spmatch = splitmatch[i];
                if ( angular.isDefined(spmatch.seq) && (spmatch.seq !== '') ) {
                    var seqline = spmatch.seq;
                    var seqid = Array(charsTitle).join(" ");
                    elem += '<browser-seq-line match>' + seqid + ' ' + seqline + '<br>' + '</browser-seq-line>';
                }
                elem += '<br>';
            }
        }
        elem += '</browser-seq>';
        return elem;
    }

    return {
        restrict: 'E',

        link: function(scope, element) {
            //scope.$watchGroup(['sequences', 'residues'], function(newValues, oldValues, scope) { // For AngularJS 1.3
            scope.$watchCollection('[sequences, residues, checkmethods, resmetseqs, checkmetseqs]', function(newValues, oldValues, scope) { // For AngularJS 1.1.4
                scope.printSeq(newValues[0], newValues[1], newValues[2], newValues[3], newValues[4]);
            });
            scope.printSeq = function(sequences, residues, methods, resmetseqs, metseqs) {
                if ( angular.isDefined(sequences) && angular.isDefined(residues) && angular.isDefined(methods) && angular.isDefined(resmetseqs) && angular.isDefined(metseqs) ) {
                    var seqs   = createSeqReferences(sequences, residues, methods, resmetseqs, metseqs);
                    var legend = createSeqLegend(scope.query);
                    var h = '<div class="browser-seq-tpl">'+
                        '<div class="browser-seq-legends pull-right">'+legend+'</div>'+
                        '<div class="browser-seq-seqs pull-left">'+seqs+'</div>'+
                        '</div>';
                    scope.browser.refresh = true;
                    element.html($compile(h)(scope));
                    scope.browser.refresh = false;
                }
                else {
                    scope.browser.refresh = true;
                    element.html($compile('<div class="browser loading">Loading...</div>')(scope));
                }
            }
            scope.$watch('error', function(newValue, oldValue, scope) {
                if ( angular.isDefined(newValue) ) {
                    var annots = newValue;
                    element.html($compile('<div class="browser error">' + annots + '</div>')(scope));
                    scope.browser.refresh = false;
                }
            });
            // destroy element
            element.on('$destroy', function() {
            });

        }
    };
}]);

module.directive("browserSeq", function() {
    return {
        restrict: "E",
        link: function(scope, element){
            // destroy element
            element.on('$destroy', function() {
            });
        }
    }
});

/* FILTERS */

// Convert residues JSON to Hash Object where each key is the residue annotations
// Convert methods JSON to Hash Object where each key is the method annotations
module.filter('filterAnnotPerResidues', function(){
    return function(input, methods) {
        var outresidues = {};
        var outmethods = {};
        angular.forEach(input, function(item) {
            var tid     = angular.isDefined(item.id) ? item.id : null;
            var resreport  = {};
            angular.forEach(item.methods, function(tmethod) {
                var mid = tmethod.id;
                var mname = tmethod.name;
                var morder = methods[mid].order;
                var minclude = true;
                angular.forEach(tmethod.residues, function(tresidue) {
                    var rstart = parseInt(tresidue.start);
                    var rend   = angular.isDefined(tresidue.end) ? parseInt(tresidue.end) : parseInt(tresidue.start);
                    var rannot   = tresidue.annot;
                    var mclass = '';
                    if ( angular.isDefined(tresidue.type) ) {
                        var rtype   = tresidue.type;
                        mclass = methods[mid].type[rtype].class;
                    }
                    else { mclass = methods[mid].class }
                    if ( (rannot.indexOf('ALTERNATIVE') > -1) || (rannot.indexOf('MINOR') > -1) ) { minclude = false }
                    for (var r = rstart; r <= rend; r++) {
                        if ( !angular.isDefined(resreport[r]) ) {
                            resreport[r] = [{
                                title: tmethod.title,
                                label: tmethod.label,
                                id:    mid,
                                name:  mname,
                                class: mclass,
                                start: rstart,
                                end:   rend,
                                annot: rannot
                            }];
                        } else {
                            resreport[r].push({
                                title: tmethod.title,
                                label: tmethod.label,
                                id:    mid,
                                name:  mname,
                                class: mclass,
                                start: rstart,
                                end:   rend,
                                annot: rannot
                            });
                        }
                    }
                });
                if ( !angular.isDefined(outmethods[morder]) ) {
                    if ( minclude ) {
                        outmethods[morder] = methods[mid];
                    }
                }
            });
            if ( resreport ) {
                outresidues[tid] = resreport;
            }
        });
        return {
            residues:   outresidues,
            methods:    outmethods
        };
    }
});

module.filter('filterCDSPerResidues', function(){
    return function(input) {
        var outresidues = {};
        var outmethods = {};
        angular.forEach(input, function(item) {
            var tid     = angular.isDefined(item.id) ? item.id : null;
            var resreport  = {};
            var mid  = 'cds';
            var mlabel = 'CDS protein region';
            var alter = 0;
            angular.forEach(item.cds, function(tcds) {
                var rstart = parseInt(tcds.start);
                var rend   = parseInt(tcds.end);
                for (var r = rstart; r <= rend; r++) {
                    var mclass = 'cds_' + alter;
                    resreport[r] = {
                        id:    mid,
                        label: mlabel,
                        class: mclass
                    };
                }
                alter = ( alter == 0 )? 1 : 0;
            });
            if ( resreport ) {
                outresidues[tid] = resreport;
            }
            outmethods[0] = {
                id:    mid,
                label: mlabel,
                desc:  mlabel,
                class: mid
            };
        });
        return {
            residues:   outresidues,
            methods:    outmethods
        };
    }
});

/* COMMON functions */
function filterMethodsFromRes(methods, mclass) {
    var filmethods = '';
    angular.forEach(methods, function(method) {
        var mcls = '';
        if ( angular.isDefined(method.type) ) {
            var rtype   = method.type;
            mcls = ' ' + method.type[rtype].class + ' ';
        }
        else { mcls = ' ' + method.class + ' ' }
        var re = new RegExp(mcls, 'g');
        if ( mclass.match(re) ) {
            filmethods += method.name + ',';
        }
    });
    return filmethods;
}
function createAnnotReferences($filter, query, residues, urlExporter, urlRunnerRst, urlFirestar, urlPDB, urlPfam) {
    var elem = '';
    var charsTitle = 40;
    if ( angular.isArray(residues) ) {
        var accorHtml = '<accordion close-others="false">';
        angular.forEach(residues, function(item) {
            var id = item.id;
            var lbl = $filter('deleteSrcNames')(id);
            var id_lbl = lbl;
            if ( lbl.length >= charsTitle) { id_lbl = lbl.substring(0,charsTitle-1) }
            var body = '<tabset>';
            angular.forEach(item.methods, function(method) {
                var mid = method.id;
                var mname = method.name;
                var mlabel = method.label;
                var mresidues = method.residues;
                var isPos = true;
                var rows = '';
                var alink = '';
                angular.forEach(mresidues, function(mres) {
                    var cols = '';
                    if ( angular.isDefined(mres.start) && !angular.isDefined(mres.end) ) {
                        cols += '<td colspan="2">' + mres.start + '</td>';
                    }
                    else {
                        cols += '<td>' + mres.start + '</td>';
                    }
                    if ( angular.isDefined(mres.end) ) {
                        isPos = false;
                        cols += '<td>' + mres.end + '</td>';
                    }
                    var annot = '-';
                    var pattern = '([^ ]*) \\(([^\\)]*)\\)';
                    if ( (mid === "functional_residue") && angular.isDefined(mres.annot.match(pattern)) ) {
                        var match = mres.annot.match(pattern);
                        if ( angular.isArray(match) && angular.isDefined(match[1]) && angular.isDefined(match[2]) ) {
                            annot = '';
                            angular.forEach(match[1].split(','), function(id,i) {
                                var sc = match[2].split(',')[i];
                                annot += "<a href='"+urlFirestar+id+"' target='_blank'>" + id + "</a>(" + sc + ") ";
                            });
                        }
                    }
                    else if ( (mid === "homologous_structure") && angular.isDefined(mres.annot.match(pattern)) ) {
                        var match = mres.annot.match(pattern);
                        if ( angular.isArray(match) && angular.isDefined(match[1]) && angular.isDefined(match[2]) ) {
                            var acc = match[1];
                            var eval = match[2];
                            var id = acc.replace(/_[A-Z0-9]{1,4}$/i,'');
                            annot = "<a href='"+urlPDB+id+"' target='_blank'>"+acc+"</a>" + " (" + eval + ")";
                        }
                    }
                    else if ( (mid === "functional_domain") && angular.isDefined(mres.annot.match(pattern)) ) {
                        var match = mres.annot.match(pattern);
                        if ( angular.isArray(match) && angular.isDefined(match[1]) && angular.isDefined(match[2]) ) {
                            var acc = match[1];
                            var eval = match[2];
                            var id = acc;
                            annot = "<a href='"+urlPfam+id+"' target='_blank'>"+acc+"</a>" + " (" + eval + ")";

                        }
                    }
                    else {
                        annot = mres.annot;
                    }

                    if ( angular.isDefined(query.species) && angular.isDefined(query.id) && angular.isDefined(query.ids) ) {
                        var al = urlExporter + '/' + 'id/' + query.species + '/' + query.id + '?' + 'ids=' + encodeURIComponent(query.ids) + '&methods=' + mname + '&format=raw';
                        if ( angular.isDefined(query.as) && query.as !== null ) {
                            al += '&as=' + query.as;
                        }
                        if ( angular.isDefined(query.sc) && query.sc !== null ) {
                            al += '&sc=' + query.sc;
                        }
                        if ( angular.isDefined(query.ds) && query.ds !== null ) {
                            al += '&ds=' + query.ds;
                        }
                        // NOTE: HARD-CORE PROTEO does not print results in detail!!!
                        if ( mname != 'proteo' ) {
                            alink += '<a class="pull-right" href="' +  al + '" target="blank">In detail<i class="glyphicon glyphicon-new-window" style="text-decoration: none; margin-left: 3px"></a></i>';
                        }
                    }
                    else if ( angular.isDefined(query.jobid) ) {
                        var al = urlRunnerRst + '/' + query.jobid + '?' + 'ids=' + encodeURIComponent(query.ids) + '&methods=' + mname + '&format=raw';
                        alink += '<a class="pull-right" href="' +  al + '" target="blank">In detail<i class="glyphicon glyphicon-new-window" style="text-decoration: none; margin-left: 3px"></a></i>';
                    }

                    if ( angular.isDefined(mres.annot) && !angular.isDefined(mres.damaged)) {
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
                        title += '<tr><th colspan="2">residue</th><th>' + method.title + alink + '</th></tr>';
                    }
                    else {
                        title += '<tr><th>start</th><th>end</th><th>' + method.title + alink + '</th></tr>';
                    }
                    var table = '<table class="annot-browser table table-bordered table-condensed">'+'<tbody>' + title + rows + '</tbody></table>';
                    body += '<tab heading="'+mlabel+'">' + table + '</tab>';
                }
            });
            body += '</tabset>';
            accorHtml += '<accordion-group class="accord">';
            accorHtml += '<accordion-heading>'+
                id_lbl +
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