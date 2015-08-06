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
app.set('port', process.env.PORT || 3000);
app.use(express.static(path.join(__dirname, path.normalize('../app_mnt'))));

/**
 * Create SERVER
 */

http.createServer(app).listen(app.get('port'), function(){
    console.log('Express server listening on port: '+ app.get('port'));
});

