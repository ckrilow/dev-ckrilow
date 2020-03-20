#!/usr/bin/env nextflow

VERSION = "0.0.1" // do not edit, controlled by bumpversion


// set default parameters
params.output_dir    = "nf-qc_normalization"
params.help          = false


// startup messge - either with help message or the parameters supplied
def help_message() {
    log.info """
    ============================================================================
     single cell qc and normalization ~ v${VERSION}
    ============================================================================

    Runs basic single cell qc and normalization

    Usage:
    nextflow run main.nf -profile <local|lsf> [options]

    Options:
        --file_paths_10x   Tab-delimited file containing experiment_id and
                           path_data_10xformat columns.

        --file_metadata    Tab-delimited file containing sample metadata.

        --output_dir       Directory name to save results to. (Defaults to
                           'nf-qc_normalization')

    Profiles:
        local              local execution
        lsf                lsf cluster execution
        local_singularity  local execution with singularity [TODO]
        lsf_singularity    lsf cluster execution with singularity [TODO]
    """.stripIndent()
}

if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell qc and normalization ~ v${VERSION}
    ============================================================================
    file_paths_10x                : ${params.file_paths_10x}
    file_metadata                 : ${params.file_metadata}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
    // a dictionary way to accomplish the text above
    // def summary = [:]
    // summary['file_paths_10x'] = params.file_paths_10x
    // log.info summary.collect { k,v -> "${k.padRight(20)} : $v" }.join("\n")
}


// initalize Channels
// example Channel init
// Channel
//     .fromPath( params.file_paths_10x )
//     .println()

// label 'big_mem' - NOTE to self.


process merge_10x_samples {
    // Takes a list of raw 10x files and merges them into one Seurat object.
    // NOTE: once normalization is set, it would be faster to normalize per
    //     sample and then merge
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 1             // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir path: "${params.output_dir}/sc_df",
               mode: "copy",
               overwrite: "true"

    input:
    file(params.file_paths_10x)
    file(params.file_metadata)

    output:
    file("adata.h5") into results_merge_10x_samples

    script:
    """
    scanpy_merge-dev.py \
        --samplesheetdata_file ${params.file_paths_10x} \
        --metadata_file ${params.file_metadata}
    """
    // """
    // make_seurat_obj_from_10x_files.R \
    //     --files_merge ${params.file_paths_10x} \
    //     --metadata_file ${params.file_metadata} \
    //     --out_file sc_df \
    //     --verbose
    // """
}


process normalize_and_pca {
    // Takes AnnData object and nomalizes across samples.
    // NOTE: once normalization is set, it would be faster to normalize per
    //     sample and then merge
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir path: "${params.output_dir}/sc_df",
               mode: "copy",
               overwrite: "true"

    input:
    file(in_file) from results_merge_10x_samples

    output:
    file("adata-normalized_pca.h5") into normalize_and_pca_andata_file
    file("adata-metadata.tsv.gz") into normalize_and_pca_metadata_file
    file("adata-pcs.tsv.gz") into normalize_and_pca_pc_file

    script:
    """
    scanpy_normalize_pca.py \
        --h5_anndata ${in_file}
    """
    // """
    // normalize_seurat_obj.R \
    //     --file ${in_file} \
    //     --metadata_split_column sanger_sample_id \
    //     --normalization_method LogNormalize \
    //     --out_file sc_df \
    //     --integrate \
    //     --verbose
    // """
}


process harmony {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir path: "${params.output_dir}/sc_df",
               mode: "copy",
               overwrite: "true"

    input:
    file(in_file_pcs) from normalize_and_pca_pc_file
    file(in_file_meta) from normalize_and_pca_metadata_file

    output:
    file("adata-pcs-harmony.tsv.gz") into harmony_file

    script:
    """
    harmony_process_pcs.R \
        --pca_file ${in_file_pcs} \
        --metadata_file ${in_file_meta} \
        --metadata_columns sanger_sample_id \
        --n_pcs 15
    """
}


workflow.onComplete {
    // executed after workflow finishes
    // ------------------------------------------------------------------------
    log.info """\n
    ----------------------------------------------------------------------------
     pipeline execution summary
    ----------------------------------------------------------------------------
    Completed         : ${workflow.complete}
    Duration          : ${workflow.duration}
    Success           : ${workflow.success}
    Work directory    : ${workflow.workDir}
    Exit status       : ${workflow.exitStatus}
    """.stripIndent()
}
