#!/usr/bin/env nextflow

VERSION = "0.0.1"

def help_message() {
    log.info"""
    ================================================================
     single cell qc and normalization ${VERSION}
    ================================================================

    Runs basic single cell qc and normalization

    Usage:
    nextflow run main.nf

    Options:
        --file_paths_10x  Tab-delimited file containing experiment_id and
                          path_data_10xformat columns.

        --file_metadata   Tab-delimited file containing sample metadata.

        --output_dir      Directory name to save results to. (Defaults to
                          'nf-qc_normalization')

    Profiles:
        standard          local execution
        docker            local execution with docker [TODO]
        singularity       local execution with singularity [TODO]
        lsf               lsf cluster execution [TODO:with singularity]
    """.stripIndent()
}

// initialize script parameters
params.help = false
if (params.help){
    help_message()
    exit 0
}

params.output_dir = "nf-qc_normalization"

// example Channel init
// Channel
//     .fromPath( params.file_paths_10x )
//     .println()

// label 'big_mem' - NOTE to self.

process merge_10x_samples {
    tag { name }
    cache false
    scratch true
    cpus 1
    memory '50 GB'
    echo true

    publishDir path: "${params.output_dir}/sc_df-not_normalized",
               mode: 'copy',
               overwrite: 'true'

    input:
    file(params.file_paths_10x)
    file(params.file_metadata)

    output:
    file("sc_df.rds.gz") into merge_10x_samples_file

    script:
    """
    make_seurat_obj_from_10x_files.R \
        --files_merge ${params.file_paths_10x} \
        --metadata_file ${params.file_metadata} \
        --out_file sc_df \
        --verbose
    """
}


workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
