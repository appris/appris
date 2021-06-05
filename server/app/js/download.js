// directive allows to provide a function to be executed to get data to be downloaded
// attributes:
// download-response - Required. Function to get data. It must return a promise. It must be declared on the $scope.
// download-success - Optional. Function to be executed if download-response function was successfully resolved. It must be declared on the $scope.
// download-error - Optional. Function to be executed if download-response function return a promise rejection. It must be declared on the $scope.
// download-mime - Optional. provide a mime type of data being downloaded. Defaulted to "text/plain"
// download-name - Optional. name of the file to download. Defaulted to "download.txt"
// download-backup-url - in case browser does not support dynamic download, this url will be called to get the file

var module = angular.module('download', []);

module.directive('downloadResponse', [ '$parse', function ($parse) {

    function saveMethod1(type, data, filename, contentType) {
        var blob = new Blob([data], { type: contentType });

        // Support for saveBlob method (Currently only implemented in Internet Explorer as msSaveBlob, other extension in case of future adoption)
        // TODO: Maybe try this method using FileSaver.js
        if( navigator.msSaveBlob ) {
            navigator.msSaveBlob(blob, filename);
            console.log("SaveBlob Success");
        }
        else {
            // Try using other saveBlob implementations, if available
            var saveBlob = navigator.webkitSaveBlob || navigator.mozSaveBlob || navigator.saveBlob;
            if(saveBlob === undefined) throw "saveBlob is not supported. Falling back to the next method";
            saveBlob(blob, filename);
            console.log("SaveBlob Success");
        }
    }

    function saveMethod2(type, data, filename, contentType) {
        // Get the blob url creator
        var urlCreator = window.URL || window.webkitURL || window.mozURL || window.msURL;
        if (urlCreator) {
            // Try to use a download link
            var link = document.createElement("a");
            var url;
            if (type === 'download') {
                // Prepare a blob URL
                var blob = new Blob([data], { type: contentType });
                url = urlCreator.createObjectURL(blob);
                link.setAttribute("href", url);

                // Set the download attribute (Supported in Chrome 14+ / Firefox 20+)
                link.setAttribute("download", filename);

                // Simulate clicking the download link
                var event = document.createEvent('MouseEvents');
                event.initMouseEvent('click', true, true, window, 1, 0, 0, 0, 0, false, false, false, false, 0, null);
                link.dispatchEvent(event);

                console.log("Download link Success");
            }
            else if (type === 'link') {
                // Prepare a blob URL
                // Use application/octet-stream when using window.location to force download
                var blob = new Blob([data], { type: contentType });
                url = urlCreator.createObjectURL(blob);
//                window.location = url;
                window.open(url, '_blank');

                console.log("window.location Success");
            }
        } else {
            throw 'UrlCreator not supported. Falling back to the next method';
        }
    }

    function saveMethod3(attrs) {
        if (attrs.downloadBackupUrl && attrs.downloadBackupUrl != '') {
            console.log('opening ' + attrs.downloadBackupUrl);
            window.open('https://' + document.domain + attrs.downloadBackupUrl, '_blank');
        } else {
            throw 'Could not download a file using any of the available methods. Also you did not provide a backup download link. No more bullets left...';
        }
    }

    function convertJSONToString(data) {
        return JSON.stringify(data, undefined, 2);
    }

    return {
        restrict: 'A',
        scope: false,
        link:function (scope, elm, attrs) {
            var getDataHandler = $parse(attrs.downloadResponse);

            elm.on('click', function() {
                var promise = getDataHandler(scope);
                promise.then(
                    function (data) {
                        if (attrs.downloadSuccess && attrs.downloadSuccess != '') {
                            var successHandler = $parse(attrs.downloadSuccess);
                            successHandler(scope);
                        }

                        var filename = attrs.downloadName || 'download.txt';
                        var contentType = attrs.downloadMime || 'text/plain';
                        var type = attrs.downloadType || "link";

                        if ( contentType === 'application/json' ) {
                            data = convertJSONToString(data);
                        }

                        try {
                            saveMethod1(type, data, filename, contentType);
                            return;
                        } catch (e) {
                            console.log(e);
                            try {
                                saveMethod2(type, data, filename, contentType);
                                return;
                            } catch (e) {
                                console.log(e);
                                try {
                                    saveMethod3(attrs);
                                    return;
                                } catch (e) {
                                    throw e;
                                }
                                throw e;
                            }
                            throw e;
                        }
                    },
                    function(data) {
                        if (attrs.downloadError && attrs.downloadError != '') {
                            var errorHandler = $parse(attrs.downloadError);
                            errorHandler(scope);
                        }
                    }
                );
            });

        }
    };

}]);
