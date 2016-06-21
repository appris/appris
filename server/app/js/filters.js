'use strict';

/* Filters */

var apprisFilters = angular.module('apprisFilters', ['ngResource']).config(["$provide", function ($provide) {
    $provide.value("principal1", 'PRINCIPAL:1');
    $provide.value("principal2", 'PRINCIPAL:2');
    $provide.value("principal3", 'PRINCIPAL:3');
    $provide.value("principal4", 'PRINCIPAL:4');
    $provide.value("principal5", 'PRINCIPAL:5');
    $provide.value("alternative1", 'ALTERNATIVE:1');
    $provide.value("alternative2", 'ALTERNATIVE:2');
    $provide.value("minor", 'MINOR');
}]);;


/**
 * GENERAL FILTERS
 **/

// Replaces string
apprisFilters.filter('replace', function(){
    return function(input, org, repl){
        if ( !angular.isDefined(input) || input === null || input == '' ) return input;
        var filtered = input.replace(org,repl);
        return filtered;
    };
});

// Split string
apprisFilters.filter('split', function(){
    return function(input, splitChar, splitIndex) {
        // do some bounds checking here to ensure it has that index
        return input.split(splitChar)[splitIndex];
    }
});

// Create Sequence Object
apprisFilters.filter('seqIO', function(){
    return function(input) {
        var output = [];
        var inputs = input.split('>');
        angular.forEach(inputs, function(item) {
            var items = item.split('\n');
            var out = {
                id:     items[0],
                seq:    items[1]
            }
            output.push(out);
        });
        return output;
    }
});

// Checks if given 'key' exists within object
apprisFilters.filter('hasParameter', function() {
    return function(items, key) {
        if ( angular.isDefined(items) && angular.isDefined(items[key]) ) {
            return true;
        }
        else {
            return false;
        }
    };
});

// Sorts the object by given parameter
apprisFilters.filter('orderObjectBy', function(hasParameterFilter) {
    return function(items, field, reverse) {
        if (!angular.isObject(items)) return items;
        var filtered = [];
        angular.forEach(items, function(item) {
            if ( angular.isObject(item) && hasParameterFilter(item,field) ) {
                filtered.push(item) ;
            }
        });
        if ( angular.isDefined(items) ) {
            filtered.sort(function (a, b) {
                return (a[field] > b[field] ? 1 : -1);
            });
            if(reverse) filtered.reverse();
        }
        return filtered;
    };
});

// AngularJS default filter with the following expression:
// "person in people | filter: {name: $select.search, age: $select.search}"
// performs a AND between 'name: $select.search' and 'age: $select.search'. We want to perform a OR.
apprisFilters.filter('seqFilter', function() {
    return function(items, props) {
        var out = [];

        if (angular.isArray(items)) {
            items.forEach(function(item) {
                var itemMatches = false;

                var keys = Object.keys(props);
                for (var i = 0; i < keys.length; i++) {
                    var prop = keys[i];
                    var text = props[prop].toLowerCase();
                    if (item[prop].toString().toLowerCase().indexOf(text) !== -1) {
                        itemMatches = true;
                        break;
                    }
                }

                if (itemMatches) {
                    out.push(item);
                }
            });
        } else {
            // Let the output be the input untouched
            out = items;
        }

        return out;
    };
});

apprisFilters.filter('capitalize', function() {
    return function(input) {
        return (!!input) ? input.replace(/([^\W_]+[^\s-]*) */g, function(txt){return txt.charAt(0).toUpperCase() + txt.substr(1).toLowerCase();}) : '';
    }
});

apprisFilters.filter('delete_version', function() {
    return function(input) {
        return (!!input) ? input.replace(/\.([^$]*)$/g, '') : '';
    }
});


/**
 * SPECIFIC FILTERS
 **/

// For SPECIES PAGE: retrieves the 'popular genomes'
apprisFilters.filter('showHousesGenomes', function() {
    return function(items) {
        if (!angular.isObject(items)) return items;
        var filtered = [];
        angular.forEach(items, function(item) {
            if ( angular.isObject(item) && item.category === 'houses' ) {
                filtered.push(item) ;
            }
        });
        return filtered;
    };
});


/* For SERVER FORM */

// Retrieves the methods that are prepare to be exectured from methods.json (control: runner)
apprisFilters.filter('runnerMethods', function() {
    return function(items) {
        if (!angular.isObject(items)) return items;
        var filtered = [];
        angular.forEach(items, function(item, key) {
            if ( angular.isObject(item) ) {
                if ( key !== 'principal_isoform' ) {
                    if ( angular.isDefined(item.control) && item.control.indexOf('runner') >= 0 ) {
                        var dat = item;
                        dat.runner = true;
                        filtered.push(dat) ;
                    }
                }
            }
        });
        return filtered;
    };
});

// Checks if the select component of method has at least one selected element
apprisFilters.filter('isSelectedMethod', function() {
    return function(items, key) {
        if (!angular.isObject(items)) return false;
        var filtered = false;
        angular.forEach(items, function(item) {
            if ( angular.isDefined(item) && item === true ) {
                filtered = true;
            }
        });
        return filtered;
    };
});

// Create list of methods that has been selected
apprisFilters.filter('selectedMethods', function() {
    return function(items) {
        if (!angular.isObject(items)) return false;
        var filtered = [];
        angular.forEach(items, function(item, key) {
            if ( angular.isDefined(item) && item === true ) {
                filtered.push(key);
            }
        });
        filtered.push('appris');
        return filtered;
    };
});

// Extract the list of sources
apprisFilters.filter('extractSourceName', function(capitalizeFilter) {
    return function(item) {
        if (!angular.isObject(item)) return false;
        var filtered = "";
        if ( angular.isDefined(item.label) ) {
            filtered = item.label;
        }
        else {
            filtered = capitalizeFilter(item.name) + item.version;
        }
        return filtered;
    };
});

// Print suffix for datafiles
apprisFilters.filter('suffixDatafiles', function() {
    return function(item) {
        var item = item.toLocaleLowerCase();
        if ( item === "txt" ) { return item }
        else if ( item === "bed" || item == "gtf" ) { return item + '.gz' }
    };
});

/* For RESULT REPORT */

//// Creates report from results with genome input
//apprisFilters.filter('convertTransScoreObj', function() {
//    return function(input) {
//        var filtered = {};
//        var idList = [];
//        var selected = [];
//        angular.forEach(input, function(item) {
//            if ( angular.isDefined(item.annotation) ) {
//                var iTrans = item.transcript_id;
//                var sLabel = item.type;
//                var sScore = item.score;
//                if ( !(filtered[iTrans]) ) {
//                    filtered[iTrans] = {};
//                    filtered[iTrans]['transcript_id'] = iTrans;
//                    if ( angular.isDefined(item.transcript_name) ) {
//                        filtered[iTrans]['transcript_name'] = item.transcript_name;
//                    }
//                    if ( angular.isDefined(item.biotype) ) {
//                        filtered[iTrans]['biotype'] = item.biotype;
//                    }
//                    if ( angular.isDefined(item.length_aa) ) {
//                        filtered[iTrans]['length_aa'] = item.length_aa;
//                    }
//                    if ( angular.isDefined(item.ccds_id) ) {
//                        filtered[iTrans]['ccds_id'] = item.ccds_id;
//                    }
//                    if ( angular.isDefined(item.no_codons) ) {
//                        filtered[iTrans]['no_codons'] = item.no_codons;
//                    }
//                    var id = {}
//                    id['id'] = iTrans;
//                    if ( angular.isDefined(item.length_aa) ) {
//                        id['length_aa'] = item.length_aa;
//                    }
//                    idList.push(id);
//                }
//                if ( sLabel == "principal_isoform") {
//                    var sAnnot = '-';
////                    if ( angular.isDefined(item.annotation) ) {
////                        if ( item.annotation == "Principal Isoform" ) {
////                            sAnnot = principal_label;
////                        }
////                        else if ( item.annotation == "Possible Principal Isoform" ) {
////                            sAnnot = candidate_label;
////                        }
////                        else if ( item.annotation == "No Principal Isoform" ) {
////                            sAnnot = alternative_label;
////                        }
////                        filtered[iTrans][sLabel] = sAnnot;
////                    }
//                    if ( angular.isDefined(item.reliability) ) {
//                        sAnnot = item.reliability;
//                    }
//                    filtered[iTrans][sLabel] = sAnnot;
//                }
//                else if ( sLabel == "peptide_signal" || sLabel == "mitochondrial_signal" ) {
//                    filtered[iTrans][sLabel] = sAnnot;
//                }
//                else {
//                    filtered[iTrans][sLabel] = sScore;
//                }
//            }
//        });
//
//        Object.keys(filtered).forEach(function (key) {
//            var val = filtered[key];
//            selected.push(val);
//        });
//        return [idList, selected];
//    };
//});

//// Checks if the result is "SeqRUNNER mode"
//apprisFilters.filter('isSeqRunner', function() {
//    return function(input) {
//        var isSeqResult = null;
//        if ( angular.isDefined(input) && input.length > 0 ) {
//            isSeqResult = false;
//            var data = input[0];
//            if ( data.mode == "sequence" ) {
//                isSeqResult = true;
//            }
//        }
//        return isSeqResult;
//    };
//});

//// Creates report from results with seq input (SeqRUNNER mode)
//apprisFilters.filter('convertSeqRunnerObj', function() {
//    return function(input) {
//        var idList = [];
//        var selected = [];
//        if ( angular.isDefined(input) && input.length > 0 ) {
//            var data = input[0];
//            if ( data.appris && angular.isObject(data.appris) ) {
//                var dat = data.appris;
//                angular.forEach(dat, function(item, key) {
//                    if ( angular.isObject(item) ) {
//                        var filtered = {};
//                        var id = {}
//                        id['id'] = key;
//                        idList.push(id);
//                        filtered['transcript_id'] = key;
//                        if ( angular.isDefined(item.principal_isoform_signal) ) {
//                            var sAnnot = '-';
////                            var item2 = item.principal_isoform_signal;
////                            if ( item2 == "YES" ) {
////                                sAnnot = principal_label;
////                            }
////                            else if ( item2 == "UNKNOWN" ) {
////                                sAnnot = candidate_label;
////                            }
////                            else if ( item2 == "NO" ) {
////                                sAnnot = alternative_label;
////                            }
//                            if ( angular.isDefined(item.reliability) ) {
//                                sAnnot = item.reliability;
//                            }
//                            filtered['principal_isoform'] = sAnnot;
//                        }
//                        //if ( angular.isDefined(item.reliability) ) {
//                        //    filtered['reliability'] = item.reliability;
//                        //}
//                        if ( angular.isDefined(item.functional_residues_score) ) {
//                            filtered['functional_residue'] = item.functional_residues_score;
//                        }
//                        if ( angular.isDefined(item.homologous_structure_score) ) {
//                            filtered['homologous_structure'] = item.homologous_structure_score;
//                        }
//                        if ( angular.isDefined(item.vertebrate_conservation_score) ) {
//                            filtered['vertebrate_conservation'] = item.vertebrate_conservation_score;
//                        }
//                        if ( angular.isDefined(item.domain_score) ) {
//                            filtered['functional_domain'] = item.domain_score;
//                        }
//                        if ( angular.isDefined(item.transmembrane_helices_score) ) {
//                            filtered['transmembrane_signal'] = item.transmembrane_helices_score;
//                        }
//                        if ( angular.isDefined(item.peptide_signal) ) {
//                            filtered['signal_peptide_mitochondrial'] = item.peptide_signal;
//                        }
//                        if ( angular.isDefined(item.mitochondrial_signal) ) {
//                            filtered['signal_peptide_mitochondrial'] += '/'+item.mitochondrial_signal;
//                        }
//
////                        angular.forEach(item, function(item2, key2) {
////                            if ( key2 == "principal_isoform_signal") {
////                                var sAnnot = '-';
////                                if ( item2 == "YES" ) {
////                                    sAnnot = principal_label;
////                                }
////                                else if ( item2 == "UNKNOWN" ) {
////                                    sAnnot = candidate_label;
////                                }
////                                else if ( item2 == "NO" ) {
////                                    sAnnot = alternative_label;
////                                }
////                                filtered['principal_isoform'] = sAnnot;
////                            }
////                            else if ( key2 == "reliability") {
////                                filtered['reliability'] = item2;
////                            }
////                            else if ( key2 == "functional_residues_score") {
////                                filtered['functional_residue'] = item2;
////                            }
////                            else if ( key2 == "homologous_structure_score") {
////                                filtered['homologous_structure'] = item2;
////                            }
////                            else if ( key2 == "vertebrate_conservation_score") {
////                                filtered['vertebrate_conservation'] = item2;
////                            }
////                            else if ( key2 == "domain_score") {
////                                filtered['functional_domain'] = item2;
////                            }
////                            else if ( key2 == "transmembrane_helices_score") {
////                                filtered['transmembrane_signal'] = item2;
////                            }
////                            else if ( key2 == "peptide_signal" || key2 == "mitochondrial_signal" ) {
////                                if ( angular.isUndefined(filtered.signal_peptide_mitochondrial) ) {
////                                    filtered['signal_peptide_mitochondrial'] = item2;
////                                }
////                                else {
////                                    filtered['signal_peptide_mitochondrial'] += '/'+item2;
////                                }
////                            }
////                        });
//                        selected.push(filtered);
//                    }
//                });
//            }
//        }
//        return [idList, selected];
//    };
//});
//
//// Replaces the appris labels for "css class"
//apprisFilters.filter('activeAnnotClass', function(principal1, principal2, principal3, principal4, principal5, alternative1, alternative2) {
//    return function(input, type){
//        var filtered = '';
//        if ( angular.isDefined(input) ) {
//            if ( (input == principal1) ||  (input == principal2) || (input == principal3) || (input == principal4) || (input == principal5) ) {
//                if ( angular.isDefined(type) && type == 'label' ) { filtered = "success" }
//                else { filtered = "principal" }
//            }
//            else if ( (input == alternative1) || (input == alternative2) ) {
//                if ( angular.isDefined(type) && type == 'label' ) { filtered = "warning" }
//                else { filtered = "candidate" }
//            }
//        }
//        return filtered;
//    };
//});

