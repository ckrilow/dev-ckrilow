#!/bin/sh

# --------------------------- Get and check user input -------------------------
# shub://letaylor/docker-letaylor:sc_qc_cluster
IMG=${1-"docker://letaylor/sc_qc_cluster:latest"}

# this will require file1 to be present
if [ -z $IMG ]; then
    echo "ERROR: IMG required"
    exit
fi
# Add singularity to path if not already there
if [ -z "$(command -v singularity)" ]; then
   echo "adding singularity to path"
   export PATH="/software/singularity-v3.5.3/bin:${PATH}"
fi
# ------------------------------------------------------------------------------


# --------------------------- Main ---------------------------------------------
# Not in /software becuase no write privileges in interactive shell
#export SINGULARITY_CACHEDIR="${HOME}/.singularity"
export SINGULARITY_CACHEDIR="/software/team152/nextflow/cache_singularity"
#export SINGULARITY_CACHEDIR="/lustre/scratch119/humgen/projects/sc-eqtl-ibd/env/cache_singularity"
#singularity cache clean

ROOT_DIR=$(pwd)
#BUILD_DIR="/lustre/scratch119/humgen/projects/sc-eqtl-ibd/env/cache_singularity"
BUILD_DIR="/tmp"
cd ${BUILD_DIR}

IMG_FILE=$(basename ${IMG} | sed s/':'/'_'/g)
rm -f ${BUILD_DIR}/${IMG_FILE}.sif
#rm -f /software/team152/nextflow/cache_singularity/${IMG_FILE}.img

# NOTE: one may need to login to an interactive shell for the build to finish.
singularity pull ${IMG}

# NOTE: /software/team152/nextflow/cache_singularity/ is only accessible
# on a head node.
cp ${BUILD_DIR}/${IMG_FILE}.sif ${HOME}/${IMG_FILE}.img
cp ${BUILD_DIR}/${IMG_FILE}.sif /lustre/scratch119/humgen/projects/sc-eqtl-ibd/env/cache_singularity/${IMG_FILE}.img
#cp ${BUILD_DIR}/${IMG_FILE}.sif /software/team152/nextflow/cache_singularity/${IMG_FILE}.img
# cd ${ROOT_DIR}

# singularity exec /software/team152/nextflow/cache_singularity/${IMG_FILE}.img echo "Image is built!"
# ------------------------------------------------------------------------------
