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

// Delete source names from string
apprisFilters.filter('deleteSrcNames', function(){
    return function(input) {
        return (!!input) ? input.replace(/[ensembl:|refseq:|uniprot:]/g, '').replace(/^\++/g, '').replace(/\++$/g,'').replace(/\++/, '+') : '';
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
                filtered.push(item);
            }
        });
        return filtered;
    };
});

// From species Object (coming from server.json): Create species id
apprisFilters.filter('getSpeciesId', function() {
    return function(input) {
        var speciesId = null;
        if ( angular.isDefined(input) && angular.isObject(input) && input ) {
            speciesId = input.scientific.toLowerCase().replace(/ /g,'_');
        }
        return speciesId;
    };
});

// From species Object (coming from server.json): Create list of dataset object
apprisFilters.filter('getDatasets', function() {
    return function(input) {
        if (!angular.isObject(input)) return input;
        var filtered = [];
        angular.forEach(input, function(items) {
            angular.forEach(items.datasets, function(item) {
                if ( angular.isObject(item) ) {
                    filtered.push(item);
                }
            });
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
        else { return item + '.gz' }
    };
});
