#!/bin/sh
#===============================================================================
#title        :cromwell-server-start.sh
#description  :Boots up a cromwell server
#author       :Leland Taylor
#date         :05/20/2019
#version      :0.1
#usage        :sh romwell-server-start.sh [input]
#input        :CROMWELL_WDL (required)
#               - WDL file of the workflow
#              CROMWELL_INPUT (required)
#               - Input paramters
#              CROMWELL_ROOT (optional)
#               - Output of cromwel jobs
#              CROMWELL_LOG (optional)
#               - Output of cromwel logs
#              CROMWELL_MEMORY (optional)
#               - Memory for Cromwell
#              CROMWELL_OPT (optional)
#               - Additional workflow options
#              CROMWELL_BACKEND (optional)
#               - Backend configuration
#              CROMWELL_CONC_WF (optional)
#               - Number of simultaneous workflows
#              CROMWELL_CONC_JL (optional)
#               - Maximum number of jobs running at once on slurm backend
#===============================================================================

#---------------------------- Init ---------------------------------------------
# get the dir of this script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
#echo ${SCRIPT_DIR}
#-------------------------------------------------------------------------------

#---------------------------- Get and check user input -------------------------
CROMWELL_WDL=$1
CROMWELL_INPUT=$2
CROMWELL_ROOT=${3-$(pwd)}
CROMWELL_LOG=${4-$(pwd)}
CROMWELL_MEMORY=${5-"30"}
CROMWELL_OPT=${6-"${SCRIPT_DIR}/../workflow_opts/slurm.json"}
CROMWELL_BACKEND=${7-"${SCRIPT_DIR}/../backend_with_db.conf"} # with MySQL
#CROMWELL_BACKEND=${7-"${SCRIPT_DIR}/../backend.conf"} # without MySQL
CROMWELL_CONC_WF=${8-"50"}
CROMWELL_CONC_JL=${9-"25"}

# check required files are present
if [ -z ${CROMWELL_WDL} ]; then
        echo "ERROR: CROMWELL_WDL required"
        exit
fi
if [ -z ${CROMWELL_INPUT} ]; then
        echo "ERROR: CROMWELL_INPUT required"
        exit
fi
if [ -z ${CROMWELL_ROOT} ]; then
        echo "ERROR: CROMWELL_ROOT required"
        exit
fi
if [ -z ${CROMWELL_LOG} ]; then
        echo "ERROR: CROMWELL_LOG required"
        exit
fi
if [ -z ${CROMWELL_MEMORY} ]; then
        echo "ERROR: CROMWELL_MEMORY required"
        exit
fi
if [ -z ${CROMWELL_OPT} ]; then
        echo "ERROR: CROMWELL_OPT required"
        exit
fi
if [ -z ${CROMWELL_BACKEND} ]; then
        echo "ERROR: CROMWELL_BACKEND required"
        exit
fi
if [ -z ${CROMWELL_CONC_WF} ]; then
        echo "ERROR: CROMWELL_CONC_WF required"
        exit
fi
if [ -z ${CROMWELL_CONC_JL} ]; then
        echo "ERROR: CROMWELL_CONC_JL required"
        exit
fi

echo CROMWELL_WDL       = "${CROMWELL_WDL}"
echo CROMWELL_INPUT     = "${CROMWELL_INPUT}"
echo CROMWELL_ROOT      = "${CROMWELL_ROOT}"
echo CROMWELL_LOG       = "${CROMWELL_LOG}"
echo CROMWELL_MEMORY    = "${CROMWELL_MEMORY} (Gb)"
echo CROMWELL_OPT       = "${CROMWELL_OPT}"
echo CROMWELL_BACKEND   = "${CROMWELL_BACKEND}"
echo CROMWELL_CONC_WF   = "${CROMWELL_CONC_WF}"
echo CROMWELL_CONC_JL   = "${CROMWELL_CONC_JL}"
echo "\n"
#-------------------------------------------------------------------------------

#--------------------------- Main ----------------------------------------------
cromwell \
    -Xmx${CROMWELL_MEMORY}g \
    -Dconfig.file=${CROMWELL_BACKEND} \
    -Dbackend.default=slurm \
    -Dbackend.providers.slurm.config.root=${CROMWELL_ROOT} \
    -Dworkflow-options.workflow-log-dir=${CROMWELL_LOG} \
    -Dsystem.input-read-limits.lines=400000000 \
    -Dsystem.input-read-limits.tsv=400000000 \
    -Dsystem.input-read-limits.json=400000000 \
    -Dsystem.max-concurrent-workflows=${CROMWELL_CONC_WF} \
    -Dbackend.providers.slurm.config.concurrent-job-limit=${CROMWELL_CONC_JL} \
    -Dcall-caching.enabled=false \
    run ${CROMWELL_WDL} \
    --inputs ${CROMWELL_INPUT} \
    --options ${CROMWELL_OPT} \
    --metadata-output metadata.json
#-------------------------------------------------------------------------------
