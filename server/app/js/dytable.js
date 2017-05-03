//var app = angular.module('plunker', []);
//
//app.controller('MainCtrl', function($scope) {
//    $scope.items = [{
//        id: 1,
//        name: "BMW",
//        test: "Test ge",
//        country: "Germany"
//    }, {
//        id: 2,
//        name: "Honda",
//        test: "Test jp",
//        country: "Japan"
//    }, {
//        id: 3,
//        name: "Samsung",
//        test: "Test ko",
//        country: "Korea"
//    },{
//        id: 4,
//        name: "KK",
//        test: "Test sp",
//        country: "Spain"
//    }];
//
//    $scope.columns = [{
//        id: "column2",
//        title: "Manufacturer",
//        directive: "secondcolumn",
//        visible: true
//    }, {
//        id: "column1",
//        title: "ID",
//        directive: "firstcolumn",
//        visible: true
//    }, {
//        id: "column3",
//        title: "Country",
//        directive: "thirdcolumn",
//        visible: false
//    }, {
//        id: "column4",
//        title: "test",
//        directive: "fourcolumn",
//        visible: true
//    }];
//
//    $scope.shuffleColumnOrder = function() {
//        $scope.columns = $scope.columns.sort(function() {
//            return .5 - Math.random();
//        });
//    }
//});

/**
 * table - AngularJS module for dynamic columns ina angularjs datatable.
 *
 */
var module = angular.module('dytable', []);


app.directive('item', function($compile) {
    function createTDElement(directive) {
        var table = angular.element('<table><tr><td ' + directive + '></td></tr></table>');
        return table.find('td');
    }

    function render(element, scope) {
        var column, html, i;
        for (i = 0; i < scope.columns.length; i++) {
            column = scope.columns[i];
            if (column.visible) {
                html = $compile(createTDElement(column.directive))(scope);
                element.append(html);
            }
        }

    }

    return {
        restrict: 'A',
        scope: {
            item: "=",
            columns: "="
        },
        controller: function($scope, $element) {
            $scope.$watch(function() {
                return $scope.columns;
            }, function(newvalue, oldvalue) {
                if (newvalue !== oldvalue) {
                    $element.children().remove();
                    render($element, $scope);
                    $compile($element.contents())($scope);
                }
            }, true);
        },
        compile: function() {
            return function(scope, element) {
                render(element, scope);
            }

        }
    };

});

app.directive("firstcolumn", function() {
    return {
        restrict: 'A',
        template: '{{item.id}}'
    }
});

app.directive("secondcolumn", function() {
    return {
        restrict: 'A',
        template: '{{item.name}}'
    }
});

app.directive("thirdcolumn", function() {
    return {
        restrict: 'A',
        template: '{{item.country}}'
    }
});

app.directive("fourcolumn", function() {
    return {
        restrict: 'A',
        template: '{{item.test}}'
    }
});