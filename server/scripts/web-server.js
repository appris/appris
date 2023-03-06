#!/usr/bin/node

/**
 * Module dependencies.
 */

var express = require('express');
var http = require('http');
var https = require('https');
var path = require('path');
var util = require('util');
var fs = require('fs');
var app = express();
//var request = require('request');
//var routes = require('./routes');
//var user = require('./routes/user');

/**
* All Environments
*/
app.set('port', process.env.PORT || 3000);
app.use( express.favicon(path.join(__dirname, path.normalize('../app/img/favicon.ico'))) );
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


// Certificate
//const privateKey = fs.readFileSync('/etc/letsencrypt/live/appris-tools.org/privkey.pem', 'utf8');
//const certificate = fs.readFileSync('/etc/letsencrypt/live/appris-tools.org/cert.pem', 'utf8');
//const ca = fs.readFileSync('/etc/letsencrypt/live/appris-tools.org/chain.pem', 'utf8');
//const credentials = {
//	key: privateKey,
//	cert: certificate,
//	ca: ca
//};

// Starting both http & https servers
//const httpServer = http.createServer(app);
//const httpsServer = https.createServer(credentials, app);

//httpServer.listen(80, () => {
//	console.log('HTTP Server running on port 80');
//});

//httpsServer.listen(443, () => {
//	console.log('HTTPS Server running on port 443');
//});


/**
 * Redirect http to https
 */
//app.get('*', function(req,res) {
//  console.log("REDIRECT to HTTPS");
//  console.log('https://appris-tools.org' + req.url);
//  if (req.headers['x-forwarded-proto'] != 'https') {
//    res.redirect('https://appris-tools.org' + req.url);
//console.log('https://' + req.headers.host + req.url);
//  }
//console.log("****");
//});

// set up a route to redirect http to https
//http.get('*', function(req, res) {  
//  console.log("REDIRECT to HTTPS");
//  console.log(req.headers.host);
//  console.log(req.url);
//
//  res.redirect('https://appris-tools.org' + req.url);
//});

// have it listen on 8080
// http.listen(8080);

/**
 * Redirect to old web site
 */

// redirect from old download files to new ones
//app.get('/download/data/:species/:file', function(req, res){
//    res.redirect('https://apprisws.bioinfo.cnio.es/download/data/' + req.params.species + '/'+req.params.file);
//});
//app.get('/download/README.txt', function(req, res){
//    res.redirect('https://apprisws.bioinfo.cnio.es/download/README.txt');
//});

// if we don't the specie
//app.get('/download/*', function(req, res){
//    console.log("REQUEST_PATH: " + req.path);
//    res.redirect('https://apprisws.bioinfo.cnio.es' + req.path);
//});

// redirect from old report site to new one.
//app.get('*', function(req, res){
//  console.log("REQUEST_PATH: " + req.path);
//  console.log(req);
//});
//app.get('http://appris-tools.org/.well-known/acme-challenge/*', function(req, res){
//  console.log("REQUEST_PATH: " + req.path);
//  res.redirect('http://apprisws.bioinfo.cnio.es' + req.path);
//});

// redirect from old report site to new one.
//app.get('/report.html', function(req, res){
//    res.redirect('https://appris.bioinfo.cnio.es/#/database' + '/id' + '/'+req.query.specie+ '/'+req.query.id +"?sc=ensembl" );
//});

