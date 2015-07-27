#!/bin/bash

# First we need to load environmental variables in order to use perl
source /etc/profile

# Method name
METHOD="getGTF"

# Input parameters
ID=${1}
SPECIE=${2}
CONF=${3}

# Get directory
PRGDIR="$(dirname "$0")"
case "$PRGDIR" in
	/*)
		;;
	*)
		PRGDIR="$PWD/$PRGDIR"
		;;
esac

# Dir declarations
REST_DIR="${PRGDIR}/../features/${ID}"
DATA_DIR="${REST_DIR}/data"
LOG_DIR="${REST_DIR}/logs"

# Log if applied
TRUE=1
FALSE=0
LOG_ENABLED=$TRUE
if [ "$LOG_ENABLED" -eq $TRUE ]; then
	LOG_FILENAME="${LOG_DIR}/${METHOD}.log"
	LOG_PARAMETERS="--loglevel=DEBUG --logfile=${LOG_FILENAME}"
fi


# print trace
echo "================="
echo "perl ${METHOD}.pl"
echo "	--conf=${CONF}"
echo "	--id=${ID}"
echo "	--species='${SPECIE}'"
echo "	--out-data=${DATA_DIR}/${ID}.gtf"
echo "	--out-transcripts=${DATA_DIR}/${ID}.transc.fa"
echo "	--out-translations=${DATA_DIR}/${ID}.transl.fa" 
echo "	${LOG_PARAMETERS}"

perl ${METHOD}.pl \
	--conf=${CONF} \
	--id=${ID} \
	--species="${SPECIE}" \
	--out-data=${DATA_DIR}/${ID}.gtf \
	--out-transcripts=${DATA_DIR}/${ID}.transc.fa \
	--out-translations=${DATA_DIR}/${ID}.transl.fa \
	${LOG_PARAMETERS}
 
echo "================="

 
 
