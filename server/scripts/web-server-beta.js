#!/usr/bin/node

/**
 * Module dependencies.
 */

var express = require('express');
var http = require('http');
var path = require('path');
var util = require('util');
var app = express();
//var request = require('request');
//var routes = require('./routes');
//var user = require('./routes/user');

/**
 * All Environments
 */
app.set('port', process.env.PORT || 3001);
app.use( express.favicon(path.join(__dirname, path.normalize('../app/img/favicon-dev.ico'))) );
app.use( express.logger('dev') );
app.use( express.json() );
app.use( express.urlencoded() );
app.use( express.methodOverride() );
app.use( express.cookieParser('your secret here') );
app.use( express.session() );
app.use( app.router);
app.use( express.static(path.join(__dirname, path.normalize('../app'))) );

/**
 * CORS Support in my Node.js web app written with Express
 */
app.all('/*', function(req, res, next) {
    res.header('Access-Control-Allow-Origin', '*');
    res.header('Access-Control-Allow-Headers', 'Content-Type, X-Requested-With');
    res.header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS');
    res.header('Access-Control-Allow-Headers', 'X-Requested-With, Content-Type, X-XSRF-TOKEN');
    //var oneDay = 86400000;
    //res.header('Cache-Control', 'max-age=' + oneDay);
    res.header("Cache-Control", 'max-age=0, no-cache, no-store, must-revalidate');
    res.header("Pragma", 'no-cache');
    res.header("Expires", 0);
    next();
});

// development only
if ('development' == app.get('env')) {
    app.use(express.errorHandler());
}

/**
 * Create SERVER
 */
http.createServer(app).listen(app.get('port'), function(){
    console.log('Express server listening on port: '+ app.get('port'));
});

/**
 * Redirect to old web site
 */

// redirect from old download files to new ones
app.get('/download/data/:species/:file', function(req, res){
//    console.log("REQUEST: " + util.inspect(req, false, null));
//    console.log("SPECIE: "+req.params.species);
//    console.log("FILE: "+req.params.file);
    res.redirect('http://apprisws-dev.bioinfo.cnio.es/download/data/' + req.params.species + '/'+req.params.file);
});
app.get('/download/README.txt', function(req, res){
    res.redirect('http://apprisws-dev.bioinfo.cnio.es/download/README.txt');
});

// if we don't the specie
//app.get(/^\/download\/data\//, function(req, res){
app.get('/download/*', function(req, res){
    console.log("REQUEST_PATH: " + req.path);
    res.redirect('http://apprisws-dev.bioinfo.cnio.es' + req.path);
});

// redirect from old report site to new one.
app.get('/report.html', function(req, res){
//    res.status(500).json({ error: 'message' })
//    console.log("PARAMS: " + util.inspect(req.params, false, null));
    res.redirect('http://appris-dev.bioinfo.cnio.es/#/database' + '/id' + '/'+req.query.specie+ '/'+req.query.id );
});
