#!/bin/sh

# --------------------------- Get and check user input -------------------------
IMG=${1-"docker://letaylor/sc_qc_cluster:latest"}

# this will require file1 to be present
if [ -z $IMG ]; then
    echo "ERROR: IMG required"
    exit
fi
# Add singularity to path if not already there
if [ -z singularity ]; then
   echo "adding singularity to path"
   export PATH="/software/singularity-v3.5.3/bin:${PATH}"
   export SINGULARITY_CACHEDIR="/software/team152/nextflow/cache_singularity"
fi
# ------------------------------------------------------------------------------


# --------------------------- Main ---------------------------------------------
ROOT_DIR=$(pwd)
BUILD_DIR="/lustre/scratch119/humgen/projects/sc-eqtl-ibd/env/cache_singularity"
cd ${BUILD_DIR}

# NOTE: one may need to login to an interactive shell for the build to finish.
singularity pull ${IMG}
IMG_FILE=$(basename ${IMG} | sed s/':'/'_'/g)

# NOTE: /software/team152/nextflow/cache_singularity/ is only accessible
# on a head node.
mv ${BUILD_DIR}/${IMG_FILE}.sif /software/team152/nextflow/cache_singularity/${IMG_FILE}.img
cd ${ROOT_DIR}
# ------------------------------------------------------------------------------
