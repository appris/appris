##########################################################
# SETTING UP ENVIRONMENT VARIABLES AND PERL LIBRARY PATH #
##########################################################

# -----------------------------------------------#
# Overwrite vars to run APPRIS in the web server #
# -----------------------------------------------#

export APPRIS_WSERVER_HOME="${APPRIS_HOME}/ws"

export APPRIS_WORKSPACE="${APPRIS_HOME}/annotations/wserver"

export APPRIS_WS_NAME="wserver"

export APPRIS_WSERVER_DOWNLOAD_DIR="${APPRIS_WSERVER_HOME}/download"

export APPRIS_WSERVER_DOWNLOAD_DATA_DIR="${APPRIS_WSERVER_DOWNLOAD_DIR}/data"

export APPRIS_NUM_PARALLEL_PROC=1

#export APPRIS_SCRIPTS_CLUSTER_INI="${APPRIS_SCRIPTS_CONF_DIR}/cluster.ini"

export APPRIS_WSERVER_SCRIPTS_DB_INI="${APPRIS_WSERVER_CONF_DIR}/apprisdb.ini"

export APPRIS_WSERVER_PIPELINE_STRUCTURE="firestar,matador3d,spade,corsair,thump,crash,appris"

export APPRIS_WSERVER_PIPELINE_STRUCTURE_SEQ="sequences,residues"

export APPRIS_WSERVER_PIPELINE_STRUCTURE_VIEW="align,genome,svg"

export APPRIS_WSERVER_OUTPUT_TYPE_SEQ="aa"

export APPRIS_WSERVER_OUTPUT_FORMAT="json"

export APPRIS_WSERVER_OUTPUT_FORMAT_SEQ="fasta"

export APPRIS_WSERVER_UCSC_URL="http://genome-euro.ucsc.edu/cgi-bin/hgTracks"

export APPRIS_WSERVER_UCSC_RENDER_URL="http://genome.ucsc.edu/cgi-bin/hgRenderTracks"

export APPRIS_WSERVER_OUTPUT_ENCODING="text/plain"

export APPRIS_WSERVER_REST_URL="http://apprisws.bioinfo.cnio.es/rest"

export APPRIS_WSERVER_REPORT_URL="http://appris.bioinfo.cnio.es/#/server/result"
