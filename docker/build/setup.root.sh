#!/bin/bash

echo "-- init mysql" 
/etc/init.d/mysql start

if [ ! -z "$(ls -A /home/appris/ws)" ]; then
  echo "-- set the permissions for the 'Working directory'"
  setfacl -R -d -m g::rwx,o::rwx /home/appris/ws
fi

