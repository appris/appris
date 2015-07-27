/**
 * Created by jmrodriguez on 4/4/14.
 */

exports.config = {
    seleniumAddress: 'http://localhost:4444/wd/hub',

    // To run plain JS files, uncomment the following line:
    specs: [
        '../test/e2e/*.js'
    ],

    // Set capabilities.
    capabilities: {
        'browserName': 'firefox'
    },

//    multiCapabilities: [{
//        'browserName': 'chrome'
//    }, {
//        'browserName': 'firefox'
//    }],

    jasmineNodeOpts: {
        showColors: true,
        defaultTimeoutInterval: 30000
    }
}
