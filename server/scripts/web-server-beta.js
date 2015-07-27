#!/usr/bin/node

/**
 * Module dependencies.
 */

var express = require('express');
var http = require('http');
var path = require('path');
var app = express();
//var routes = require('./routes');
//var user = require('./routes/user');


/**
 * All Environments
 */
app.set('port', process.env.PORT || 3001);
app.use(express.favicon());
app.use(express.logger('dev'));
app.use(express.json());
app.use(express.urlencoded());
app.use(express.methodOverride());
app.use(express.cookieParser('your secret here'));
app.use(express.session());
app.use(app.router);
app.use(express.static(path.join(__dirname, path.normalize('../app_beta'))));
//app.use(express.static(path.join(__dirname, path.normalize('..'))));

/**
 * CORS Support in my Node.js web app written with Express
 */
app.all('/*', function(req, res, next) {
    var oneDay = 86400000;
    res.header('Access-Control-Allow-Origin', '*');
    res.header('Access-Control-Allow-Headers', 'Content-Type, X-Requested-With');
    res.header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS');
    res.header('Access-Control-Allow-Headers', 'X-Requested-With, Content-Type, X-XSRF-TOKEN');
    res.header('Cache-Control', 'max-age=' + oneDay);
    next();
});


// development only
if ('development' == app.get('env')) {
  app.use(express.errorHandler());
}

//app.get('/', routes.index);
//app.get('/users', user.list);

http.createServer(app).listen(app.get('port'), function(){
  console.log('Express server listening on port: '+ app.get('port'));
});
