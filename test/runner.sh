#!/bin/bash

# -----------------------------------------------------------------------------
# Load environmental variables
# If you do not trust the path, configure below:
SYSTEM=`uname -s`
if [ "$SYSTEM" = Darwin ]
then
	source /etc/bashrc
	source /etc/profile
	source ${HOME}/.bash_profile
elif [ "$SYSTEM" = Linux ]
then
	source /etc/profile
	source /etc/bash.bashrc
	source ${HOME}/.bashrc
fi

VER=0.1b
VERDATE="1-Apr-2013"

FILENAME=`basename $0`
FILENAME="${FILENAME%.*}"
DIRNAME=`dirname $0`

GENE_LIST=""
METHODS=""
OUTPUT_DIR=""
LOG_LEVEL=""


# -----------------------------------------------------------------------------
# Prepare environment for tests

echo "-- test environment"
export APPRIS_TEST_DIR="${APPRIS_HOME}/test"

export APPRIS_SAMPLES_DIR="${APPRIS_TEST_DIR}/samples"


# -----------------------------------------------------------------------------
# Clear old logs
echo "-- clear old logs"
find ${APPRIS_HOME}/test/logs -type f -name "*.log" -exec rm {} \;


# -----------------------------------------------------------------------------
# Prepare environment for Human-Gencode tests

NUM_TEST=1
TEST_SOURCE="TST.easy"
TEST_METHODS="firestar,matador3d,spade,corsair,thump,crash,appris"
TEST_SAMPLES="${APPRIS_SAMPLES_DIR}/easy.txt"
echo "-- run test ${NUM_TEST}"
echo "apprisall -c ${TEST_SOURCE} -g ${TEST_SAMPLES} -m ${TEST_METHODS} -l debug"
apprisall -c ${TEST_SOURCE} -g ${TEST_SAMPLES} -m ${TEST_METHODS} -l debug

NUM_TEST=$((${NUM_TEST}+1))
TEST_SOURCE="Test.1"
TEST_METHODS="spade,appris"
TEST_SAMPLES="${APPRIS_SAMPLES_DIR}/spade.txt"
#echo "-- run test ${NUM_TEST}"
#echo "apprisall -c ${TEST_SOURCE} -g ${TEST_SAMPLES} -m ${TEST_METHODS} -l debug"
#apprisall -c ${TEST_SOURCE} -g ${TEST_SAMPLES} -m ${TEST_METHODS} -l debug

NUM_TEST=$((${NUM_TEST}+1))
TEST_SOURCE="Test.1"
TEST_METHODS="appris"
TEST_SAMPLES="${APPRIS_SAMPLES_DIR}/appris.txt"
#echo "-- run test ${NUM_TEST}"
#echo "apprisall -c ${TEST_SOURCE} -g ${TEST_SAMPLES} -m ${TEST_METHODS} -l debug"
#apprisall -c ${TEST_SOURCE} -g ${TEST_SAMPLES} -m ${TEST_METHODS} -l debug











# DEPRECATED but MAYBE IT WOULD BE USEFUL
# -----------------------------------------------------------------------------
# Load external functions
source ${DIRNAME}/../scripts/bin/appris_env

# -----------------------------------------------------------------------------
# Prepare environment for battery test

#export APPRIS_TEST_DIR="${APPRIS_HOME}/test"

#export APPRIS_SAMPLES_DIR="${APPRIS_TEST_DIR}/samples"

#export APPRIS_ANNOT_DIR="${APPRIS_TEST_DIR}/annotations"

#export APPRIS_SCRIPTS_DIR="${APPRIS_HOME}/scripts"

#METHOD_SCRIPT="apprisall"
#SCRIPT="${APPRIS_SCRIPTS_DIR}/${METHOD_SCRIPT}.pl"


#export APPRIS_SPECIE="Homo sapiens"
#load_appris_specie_env "Hsap"

# print trace
#echo "================="
#echo "perl ${SCRIPT}"
#echo " --species="${APPRIS_SPECIE}""
#echo " --outpath=${APPRIS_ANNOT_DIR}"
#echo " ${GENCODE_PARAMETERS}"
#echo " ${ENSEMBL_PARAMETERS}"
#echo " ${GENE_LIST_PARAMETERS}"
#echo " ${METHOD_PARAMETERS}"
#echo " ${T_ALIGN_PARAMETERS}"
#echo " ${CACHED_PARAMETERS}"
#echo " ${CLUSTER_PARAMETERS}"
#echo " ${LOG_PARAMETERS}"

# run method
#perl ${SCRIPT} \
#	--species="${APPRIS_SPECIE}" \
#	--outpath=${APPRIS_ANNOT_DIR} \
#	${GENCODE_PARAMETERS} \
#	${ENSEMBL_PARAMETERS} \
#	${GENE_LIST_PARAMETERS} \
#	${METHOD_PARAMETERS} \
#	${T_ALIGN_PARAMETERS} \
#	${CACHED_PARAMETERS} \
#	${CLUSTER_PARAMETERS} \
#	${LOG_PARAMETERS} \
#echo "================="
	