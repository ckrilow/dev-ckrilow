#!/bin/sh
#===============================================================================
#title        :cromwell-server-slurm-submit.sh
#description  :Submits a job to a cromwell instance
#author       :Leland Taylor
#date         :05/20/2019
#version      :0.1
#usage        :sh cromwell-server-slurm-submit.sh [input]
#input        :CROMWELL_IP (required)
#               - IP adress of cromwell server instance
#              CROMWELL_WDL (required)
#               - WDL file of the workflow
#              CROMWELL_INPUT (required)
#               - Input paramters
#              CROMWELL_OPT (required)
#               - Additional workflow options
#===============================================================================

#---------------------------- Get and check user input -------------------------
CROMWELL_IP=$1
CROMWELL_WDL=$2
CROMWELL_INPUT=$3
CROMWELL_OPT=$4
CROMWELL_WF_DEPS_IN=$5

# check required files are present
if [ -z ${CROMWELL_IP} ]; then
    echo "ERROR: CROMWELL_IP required"
    exit
fi
if [ -z ${CROMWELL_WDL} ]; then
    echo "ERROR: CROMWELL_WDL required"
    exit
fi
if [ -z ${CROMWELL_INPUT} ]; then
    echo "ERROR: CROMWELL_INPUT required"
    exit
fi
if [ -z ${CROMWELL_OPT} ]; then
    echo "ERROR: CROMWELL_OPT required"
    exit
fi
if [ -z ${CROMWELL_WF_DEPS_IN} ]; then
    CROMWELL_WF_DEPS=""
else
    CROMWELL_WF_DEPS="-F workflowDependencies=@${CROMWELL_WF_DEPS_IN}"
fi

echo CROMWELL_IP        = "${CROMWELL_IP}"
echo CROMWELL_WDL       = "${CROMWELL_WDL}"
echo CROMWELL_INPUT     = "${CROMWELL_INPUT}"
echo CROMWELL_OPT       = "${CROMWELL_OPT}"
echo CROMWELL_WF_DEPS   = "${CROMWELL_WF_DEPS}"
echo "\n"
#-------------------------------------------------------------------------------

#--------------------------- Main ----------------------------------------------
WF_ID=$(
    curl -X POST \
        --header "Accept: application/json" \
        -v "${CROMWELL_IP}:8000/api/workflows/v1" \
        -F workflowSource=@${CROMWELL_WDL} \
        -F workflowInputs=@${CROMWELL_INPUT} \
        -F workflowOptions=@${CROMWELL_OPT} \
        ${CROMWELL_WF_DEPS} \
    | jq '.id' \
    | sed s/\"//g
)

echo "\n\nLinks to workflow information:"
echo "\thttp://${CROMWELL_IP}:8000/api/workflows/v1/${WF_ID}/timing"
echo "\thttp://${CROMWELL_IP}:8000/api/workflows/v1/${WF_ID}/status"
echo "\thttp://${CROMWELL_IP}:8000/api/workflows/v1/${WF_ID}/outputs"
echo "\thttp://${CROMWELL_IP}:8000/api/workflows/v1/${WF_ID}/metadata"
#-------------------------------------------------------------------------------
