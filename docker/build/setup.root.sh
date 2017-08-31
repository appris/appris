#!/bin/bash

echo "-- init mysql" 
/etc/init.d/mysql start

echo "-- set the permissions for the 'Working directory'" 
setfacl -R -d -m g::rwx,o::rwx /home/appris/ws

