#/usr/bin/env bash

# Env vars of APPRIS web server
source ${PWD}/../../.bashrc
source ${PWD}/../../conf/apprisrc
source ${PWD}/../../conf/apprisrc.WS

# Then, go!
/opt/perlbrew/bin/perlbrew exec -q --with perl-5.18.2 \
     perl viewer.pl "$@"
