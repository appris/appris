#!/bin/bash

# Init MySQL
/etc/init.d/mysql start

# APPRIS env
export APPRIS_HOME="/opt/appris"
source ${APPRIS_HOME}/conf/apprisrc
source ${APPRIS_HOME}/conf/apprisrc.WS

