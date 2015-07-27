#!/bin/bash

# First we need to load environmental variables in order to use perl
source /etc/profile

# Method name
METHOD="getEComparaAlign"

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
echo "	--species='${SPECIE}'"
echo "	--data=${DATA_DIR}/${ID}.gtf"
echo "	--transcripts=${DATA_DIR}/${ID}.transc.fa"
echo "	--translations=${DATA_DIR}/${ID}.transl.fa"
echo "	--output_dir=${DATA_DIR}/" 
echo "	${LOG_PARAMETERS}"

perl ${METHOD}.pl \
	--conf=${CONF} \
	--species="${SPECIE}" \
	--data=${DATA_DIR}/${ID}.gtf \
	--transcripts=${DATA_DIR}/${ID}.transc.fa \
	--translations=${DATA_DIR}/${ID}.transl.fa \
	--output_dir=${DATA_DIR}/ \
	${LOG_PARAMETERS}

echo "================="

 
 
