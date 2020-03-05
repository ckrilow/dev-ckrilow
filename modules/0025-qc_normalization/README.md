
# Description
The methods used in this module are described in `docs/methods.pdf`.

Insert brief description of module. This is a Cromwell workflow example.



# Quickstart

Quickstart for deploying this pipeline on a local cluster.

## 1. Set up the environment

Install most of the packages via conda:
```bash
# the repo directory
REPO_MODULE="${HOME}/repo/path/to/this/module"

# install environment using Conda
conda env create --name example_module_1 --file ${REPO_MODULE}/environment.yml

# activate the new Conda environment
source activate example_module_1

# to update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```


## 2. Prepare the input files

The Cromwell pipeline takes as input:
1. File1.
2. File2.

Prepare the File1.
```bash
python "${REPO_MODULE}/scripts/prepare_file1_script.py" \
    --input_file input-file1.tsv > input-file1.tsv
```

## 3. Set up and run Cromwell

* **Note**: When scattering for many jobs without a MySQL backend, you will likely an out of memory error for Java (java.lang.OutOfMemoryError: Java heap space). To set up a MySQL ask your systems administrator to make a GCP MySQL instance.
* **Note:** If Cromwell does not shutdown gracefully and you are using a MySQL database for logging, you may need to delete and rebuild MySQL database.
* **Note**: If Cromwell ever stalls at the MySQL access stage, it is likely because the underlying Cromwell database (e.g., cromwell_db) within MySQL is corrupt. Delete the Cromwell database within MySQL. Create a new Cromwell database within MySQL. Re-run Cromwell.


Change to a temporary runtime dir:
```bash
mkdir -p /tmp/${USER}/example_module_1
cd /tmp/${USER}/example_module_1
```

Set up Cromwell:
```bash
# get the dir with the cromwell src wdl file
CROMWELL_SRC=${REPO_MODULE}
# update PATH to include the scripts
export PATH="${PATH}:${CROMWELL_SRC}/scripts"

# the cromwell libs
CROMWELL_LIB="${HOME}/my_workflow/lib/cromwell"
```

Set the input parameters for your Cromwell workflow. In the `cromwell-slurm-run.sh` script below we assume the input parameters are in `params-input.json`. A description of each parameter can be found in `example_module_1.wdl`.

The below commands will generate an *empty* input file that describes all of the various input parameters and their default values. You do not need to keep all of these parameters as many of them have defaults defined in the wdl script. See the provided example of an edited file: [example_runtime_setup/params-input.json](example_runtime_setup/params-input.json).
```bash
less -S example_module_1.wdl

womtool inputs example_module_1.wdl > params-input.json
vim params-input.json
```


In order to use the backend with a MySQL database, you will need to set the `url` (which points to the MySQL database), `user` (MySQL database user id), and `password` (MySQL database user password) parameters in `${CROMWELL_LIB}/backend_with_db.conf`.
```bash
# configure the backend
CROMWELL_BACKEND="$(pwd)/my_backend_with_db.conf"

# copy the general backend because backend_with_db.conf assumes backed.conf is
# in the same directory as backend_with_db.conf
cp ${CROMWELL_LIB}/backend.conf backend.conf # without MySQL

# edit the url, user, and password parameters in this file
cp ${CROMWELL_LIB}/backend_with_db.conf my_backend_with_db.conf # with MySQL
vim my_backend_with_db.conf
```

Check that you can connect to the MySQL database using MySQL directly:
```bash
mysql -h <ip_address> -u <user> -p
```



Run Cromwell using slurm on a local cluster with a MySQL backend. The provided Cromwell script sets some useful default parameters.
```bash
# make directories for Cromwell output on a local cluster
mkdir -p /tmp/${USER}/cromwell/cromwell-executions
mkdir -p /tmp/${USER}/cromwell/cromwell-workflow-logs

# get the dir with the cromwell src wdl file
CROMWELL_SRC=${REPO_MODULE}
# update PATH to include the scripts
export PATH="${PATH}:${CROMWELL_SRC}/scripts"

# set key parameters for cromwell-slurm-run.sh
CROMWELL_WDL="${CROMWELL_SRC}/example_module_1.wdl"
CROMWELL_INPUT="params-input.json"
CROMWELL_OPT="${CROMWELL_LIB}/workflow_opts/slurm.json"
CROMWELL_BACKEND="${CROMWELL_BACKEND}" # NOTE: must be an absolute path
CROMWELL_CONC_WF="100" # Number of simultaneous workflows
CROMWELL_CONC_JL="100" # Maximum number of jobs running at once on slurm backend

# run Cromwell giving 15 Gb of memory to the JVM
CMD="${CROMWELL_LIB}/scripts/cromwell-slurm-run.sh \
   ${CROMWELL_WDL} \
   ${CROMWELL_INPUT} \
   \"/scratch/${USER}/cromwell/cromwell-executions\" \
   \"/scratch/${USER}/cromwell/cromwell-workflow-logs\" \
   15 \
   ${CROMWELL_OPT} \
   ${CROMWELL_BACKEND} \
   ${CROMWELL_CONC_WF} \
   ${CROMWELL_CONC_JL}"
echo ${CMD}
```

See example runtime setup in the [example_runtime_setup](example_runtime_setup) dir.



# Pipeline

## Core pipeline scripts

* `example_module_1.wdl`: Performs the core analysis.

# Additional notes

## Analysis 1 of example_module_1

* Additional notes.


## Cromwell

### Run Cromwell on GCP:

If you are using GCP for pipeline execution you need to perform the following steps (but using credentials specific to your setup). For more details on these steps see GCP and the Cromwell website.

```bash
# by default the authentication is saved to:
# ${HOME}/.config/gcloud/application_default_credentials.json
gcloud auth application-default login --no-launch-browser

# the below method is preferable to simply exporting (in some cases simply
# exporting did not work properly)
GCP_SERVICE_ID="" # something@your-project.iam.gserviceaccount.com
GCP_SERVICE_JSON=${USER}/.config/gcloud/something.json
gcloud auth activate-service-account ${GCP_SERVICE_ID} \
    --key-file ${GCP_SERVICE_JSON}

# Note that in addition setting up the gsutil as described in the tutorials,
# you must also export the below flag
export GOOGLE_APPLICATION_CREDENTIALS=${HOME}/.config/gcloud/${USER}-ukbb-service-credentials.json

# set an account
gcloud auth list
gcloud config set account ${GCP_SERVICE_ID}
```

Edit the Cromwell options file `${CROMWELL_LIB}/workflow_opts/google.json` and add the following variables: `default_runtime_attributes` and `user_service_account_json`, and optionally `monitoring_image` (to keep track of resources consumed).
```bash
# pass default_runtime_attributes and user_service_account_json
# optionally add monitoring_image to keep track of resources consumed
cp "${CROMWELL_LIB}/workflow_opts/google.json" cromwell-opts-gcp.json
CROMWELL_OPT="cromwell-opts-gcp.json"
vim cromwell-opts-gcp.json
```

Edit the `${CROMWELL_LIB}/backend.conf` backend file to properly setup the cloud backend.
```bash
# configure the backend
CROMWELL_BACKEND="$(pwd)/my_backend_with_db.conf"

# copy the general backend because backend_with_db.conf assumes backed.conf is
# in the same directory as backend_with_db.conf
cp ${CROMWELL_LIB}/backend.conf backend.conf # without MySQL

# edit the url, user, and password parameters in this file
cp ${CROMWELL_LIB}/backend_with_db.conf my_backend_with_db.conf # with MySQL
vim my_backend_with_db.conf
```

Now run:
```bash
# specify the root where Cromwell will write results
CROMWELL_ROOT="gs://my_cloud_bucket/tmp"

cromwell \
    -Dconfig.file=${CROMWELL_BACKEND} \
    -Dbackend.default=google \
    -Dbackend.providers.google.config.project=${CROMWELL_ROOT} \
    -Dbackend.providers.google.config.root=${CROMWELL_ROOT} \
    run example_module_1.wdl \
    --inputs inputs-gcp.json  \
    --options ${CROMWELL_OPT} \
    --metadata-output metadata.json
```

Run Cromwell with call caching. With caching enabled, Cromwell writes the location of the output file from a call. If the inputs are the same, doesn't re-run the call, but rather just copies or links to the output (depending on the filesystems::gcs::caching::duplication settings). If copy is enabled, all of the files in the shard dir of the previous run will be copied to the current shard. If reference is enabled, a `call_caching_placeholder.txt` file is created that points to the previous run's shard dir. See more [here](https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching/).
```bash
cromwell \
    -Dconfig.file=${CROMWELL_BACKEND} \
    -Dbackend.default=google \
    -Dbackend.providers.google.config.project=${CROMWELL_ROOT} \
    -Dbackend.providers.google.config.root=${CROMWELL_ROOT} \
    -Dcall-caching.enabled=true \
    run example_module_1.wdl \
    --inputs inputs-gcp.json \
    --options ${CROMWELL_OPT} \
    --metadata-output metadata.json
```

Get a cost estimate of the Cromwell jobs on GCP:
```bash
python "${CROMWELL_LIBS}/cromwell-accountant/cost.py" < metadata.json
```

### Run Cromwell on slurm:
```bash
CROMWELL_BACKEND="${CROMWELL_LIB}/backend_with_db.conf"
CROMWELL_OPT="${CROMWELL_LIB}/workflow_opts/slurm.json"

# set a directory for a bunch of temp files
CROMWELL_ROOT="/tmp/${USER}/cromwell/cromwell-executions"
CROMWELL_LOGS="/tmp/${USER}/cromwell/cromwell-workflow-logs"

# increase the byte count limit for the proper input function
cromwell \
    -Dconfig.file=${CROMWELL_BACKEND} \
    -Dbackend.default=slurm \
    -Dbackend.providers.slurm.config.root=${CROMWELL_ROOT} \
    -Dworkflow-options.workflow-log-dir=${CROMWELL_LOGS} \
    -Dsystem.max-concurrent-workflows=10 \
    run ${CROMWELL_SRC}/example_module_1.wdl \
    --inputs inputs-local.json \
    --options ${CROMWELL_OPT} \
    --metadata-output metadata.json
```

Note that if running cromwell over a long list of files, you may need to increase the input file size limit. Below is an example. In such cases using a MySQL backend with caching is strongly recommended (otherwise the memory requirements will be huge).
```bash
INPUT_FILE="files.txt"
# get the byte count of the input file
wc -c < ${INPUT_FILE}

# increase the byte count limit for the proper input function (in this case tsv
# for read_tsv function in cromwell)
cromwell \
    -Dconfig.file=${CROMWELL_BACKEND} \
    -Dbackend.default=slurm \
    -Dbackend.providers.slurm.config.root=${CROMWELL_ROOT} \
    -Dworkflow-options.workflow-log-dir=${CROMWELL_LOGS} \
    -Dsystem.input-read-limits.tsv=1974965 \
    -Dsystem.max-concurrent-workflows=10 \
    run ${CROMWELL_SRC}/example_module_1.wdl \
    --inputs inputs-local.json \
    --options ${CROMWELL_OPT} \
    --metadata-output metadata.json
```
