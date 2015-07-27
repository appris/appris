/**
 *
 * Status - AngularJS module for status report in APPRIS
 *
 */
var module = angular.module('appris.status', []);

/* CONTROLLERS */

module.controller('StatusController', ['consPathServerResult', '$scope', '$location', '$routeParams', '$timeout', 'Status',
    function(consPathServerResult, $scope, $location, $routeParams, $timeout, Status) {
        $scope.num_queries = 0;
        // recursive method that retrieves status
        $scope.control = function(jobid, data){
            if ( !angular.isDefined(data) ) {
                $scope.status = {
                    type:   'danger',
                    status: 'ERROR',
                    header: "An error occurred attempting the job: "+jobid,
                    log: "The status is not defined",
                    progress: 0
                };
            }
            else {
                if ( (data.status == 'NOT_FOUND' && $scope.num_queries > 2) ){
                    $scope.status = {
                        type:   'danger',
                        status: 'ERROR',
                        header: "An error occurred attempting the job: "+jobid,
                        log: data.log,
                        progress: 0
                    };
                }
                else if ( data.status == 'FINISHED' || data.status == 'FAILURE' || data.status == 'ERROR' ){
                    $scope.status = {
                        status: data.status,
                        log:    data.log,
                        progress: 100
                    };
                    if ( data.status == 'FINISHED' ) {
                        $scope.status.type = 'success';
                        $scope.status.header = "Your job has finished succesfully";
                        $location.url(consPathServerResult+jobid);
                    }
                    else if ( data.status == 'FAILURE' ) {
                        $scope.status.type = 'warning';
                        $scope.status.header = "Your job has finished with some errors";
                        $location.url(consPathServerResult+jobid);
                    }
                    else if ( data.status == 'ERROR' ) {
                        $scope.status.type = 'danger';
                        $scope.status.header = "Your job has failed";
                    }
                }
                else if ( data.status == 'PENDING' || data.status == 'RUNNING' || (data.status == 'NOT_FOUND' && $scope.num_queries <= 2) ){
                    var status = data.status;
                    var log = data.log;
                    if ( data.status == 'NOT_FOUND' && $scope.num_queries <= 2 ) { log = "The job is waiting to be processed"; }
                    else if ( log === "" ) { log = "Your job is running"; }
                    var progress = data.progress*100;
                    $scope.status = {
                        type:   'info',
                        header: "Your job is running to be processed... Please be patient",
                        status: status,
                        log:    log,
                        progress: progress
                    };
                    $scope.num_queries += 1;
                    Status.get({ jobid: jobid },
                        function(dat) {
                            $timeout(function() { $scope.control(jobid, dat); }, 10000);
                        },
                        function(err) {
                            $scope.status = {
                                type: 'danger',
                                status: 'ERROR',
                                header:    "Failed to received: "+jobid,
                                progress: 0
                            };
                        }
                    );
                }
            }
        };

        // init status control
        if ( angular.isUndefined($routeParams.jobid) ) {
            $scope.status = {
                type: 'danger',
                status: "ERROR",
                header: "An error occurred attempting the job",
                log: "The job is not defined",
                progress: 0
            };
        }
        else {
            var jobid = $routeParams.jobid;
            $scope.status = {
                type: 'info',
                status: "PENDING",
                header: "Your job is waiting to be processed... Please be patient",
                log: "The job is waiting to be processed",
                progress: 0
            };
            $scope.control(jobid, $scope.status);
        }
    }
]);