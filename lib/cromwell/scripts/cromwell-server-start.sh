#!/bin/sh
#===============================================================================
#title        :cromwell-server-start.sh
#description  :Boots up a cromwell server
#author       :Leland Taylor
#date         :05/20/2019
#version      :0.1
#usage        :sh romwell-server-start.sh [input]
#input        :CROMWELL_BACKEND (optional)
#               - Backend configuration
#              CROMWELL_ROOT (optional)
#               - Output of cromwel jobs
#              CROMWELL_LOG (optional)
#               - Output of cromwel logs
#              CROMWELL_MEMORY (optional)
#               - Memory for Cromwell
#===============================================================================

#---------------------------- Init ---------------------------------------------
# get the dir of this script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
#echo ${SCRIPT_DIR}
#-------------------------------------------------------------------------------

#---------------------------- Get and check user input -------------------------
CROMWELL_ROOT=${1-$(pwd)}
CROMWELL_LOG=${2-$(pwd)}
CROMWELL_MEMORY=${3-"20"}
CROMWELL_BACKEND=${4-"${SCRIPT_DIR}/../backend_with_db.conf"}
#CROMWELL_BACKEND=${4-"${SCRIPT_DIR}/../backend.conf"}

# check required files are present
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
if [ -z ${CROMWELL_BACKEND} ]; then
        echo "ERROR: CROMWELL_BACKEND required"
        exit
fi

echo CROMWELL_ROOT      = "${CROMWELL_ROOT}"
echo CROMWELL_LOG       = "${CROMWELL_LOG}"
echo CROMWELL_MEMORY    = "${CROMWELL_MEMORY} (Gb)"
echo CROMWELL_BACKEND   = "${CROMWELL_BACKEND}"
echo "\n"
#-------------------------------------------------------------------------------

#--------------------------- Main ----------------------------------------------
echo $(hostname -f)
cromwell \
    -Xmx${CROMWELL_MEMORY}g \
    -Dconfig.file=${CROMWELL_BACKEND} \
    -Dbackend.default=slurm \
    -Dbackend.providers.slurm.config.root=${CROMWELL_ROOT} \
    -Dworkflow-options.workflow-log-dir=${CROMWELL_LOG} \
    -Dsystem.input-read-limits.lines=400000000 \
    -Dsystem.input-read-limits.tsv=400000000 \
    -Dsystem.input-read-limits.json=400000000 \
    -Dsystem.workflow-restart=false \
    -Dsystem.max-concurrent-workflows=250 \
    -Dbackend.providers.slurm.config.concurrent-job-limit=400 \
    -Dcall-caching.enabled=false \
    server
#-------------------------------------------------------------------------------
