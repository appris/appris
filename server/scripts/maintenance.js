/**
 * Created by jmrodriguez on 8/6/15.
 */
var express = require('express');
var maintenance = require('maintenance');
var app = express();

//app.set('port', process.env.PORT || 3000);
app.get('/', function (req, res) {
    console.log(req.url);
    res.send(200);
});

var options = {
    current: true,                      // current state, default **false**
    httpEndpoint: true,                 // expose http endpoint for hot-switch, default **false**,
    url: '/app/mt',                     // if `httpEndpoint` is on, customize endpoint url, default **'/maintenance'**
    accessKey: 'xx4zUU8Cyy7',           // token that client send to authorize, if not defined `access_key` is not used
    view: 'myview.html',                // view to render on maintenance, default **'maintenance.html'**
    api: '/api',                        // for rest API, species root URL to apply, default **undefined**
    status: 503,                        // status code for response, default **503**
    message: 'will be back'             // response message, default **'sorry, we are on maintenance'**
};

maintenance(app, options);