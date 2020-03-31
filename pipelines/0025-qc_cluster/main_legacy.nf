#!/usr/bin/env nextflow

VERSION = "0.0.1" // do not edit, controlled by bumpversion


// set default parameters
params.output_dir    = "nf-qc_cluster"
params.help          = false
params.scale_vars_to_regress = ['', 'total_counts,age']
params.umap_colors_quantitative = ['age']
params.colors_categorical = ['sanger_sample_id,sex,leiden']


// startup messge - either with help message or the parameters supplied
def help_message() {
    log.info """
    ============================================================================
     single cell qc and clustering ~ v${VERSION}
    ============================================================================

    Runs basic single cell qc and clustering

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Options:
        --file_paths_10x   Tab-delimited file containing experiment_id and
                           path_data_10xformat columns.

        --file_metadata    Tab-delimited file containing sample metadata.

        --output_dir       Directory name to save results to. (Defaults to
                           'nf-qc_cluster')

    Profiles:
        local              local execution
        lsf                lsf cluster execution
        local_singularity  local execution with singularity [TODO]
        lsf_singularity    lsf cluster execution with singularity [TODO]

    Params file:
        Additional parameters for runtime maybe passed via a yaml file. See
        example in example_runtime_setup/params.yml
    """.stripIndent()
}


if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell qc and clustering ~ v${VERSION}
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

// Channel: variables to regress out prior to scaling.
Channel
    .fromList( params.scale_vars_to_regress )
    .set { scale__vars2regress }


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

    publishDir  path: "${outdir}",
                mode: "copy",
                overwrite: "true"

    input:
    file(params.file_paths_10x)
    file(params.file_metadata)

    output:
    val(outdir) into merge__outDir
    file("adata.h5") into merge__annData

    script:
    outdir = "${params.output_dir}"
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
    // Takes annData object and nomalizes across samples.
    // NOTE: once normalization is set, it would be faster to normalize per
    //     sample and then merge
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    // Publish the adata to dir the specifies the parameters used
    // Publish the PCs to a sub-directory
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("pcs.tsv.gz")) {
                       "reduced_dims-pca/${filename}"
                    } else {
                        "${filename}"
                    }
                },
                mode: "copy",
                overwrite: "true"

    input:
    val(outdir_prev) from merge__outDir
    file(infile) from merge__annData
    each vars2regress from scale__vars2regress

    output:
    tuple(
        val(outdir),
        file("adata-pcs.tsv.gz"),
        file("adata-normalized_pca.h5"),
        file("adata-metadata.tsv.gz")
    ) into normalize__results
    file("adata-normalized_pca.h5") into normalize__annData
    file("adata-metadata.tsv.gz") into normalize__metaData
    // file("adata-metadata.tsv.gz") into normalize__metaData

    script:
    outdir = "${outdir_prev}/normalize.total_count"
    if (vars2regress == "") {
        outdir = "${outdir}-scale.vars_to_regress=none"
        """
        scanpy_normalize_pca.py \
            --h5_anndata ${infile}
        """
    } else {
        outdir = "${outdir}-scale.vars_to_regress=${vars2regress}"
        """
        scanpy_normalize_pca.py \
            --h5_anndata ${infile} \
            --vars_to_regress ${vars2regress}
        """
    }
    // """
    // normalize_seurat_obj.R \
    //     --file ${infile} \
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

    publishDir  path: "${outdir}",
                mode: "copy",
                overwrite: "true"

    input:
    tuple(
        val(outdir_prev),
        file(infile_reducedDims),
        file(ignore_0),
        file(infile__metaData)
    ) from normalize__results
    file(test) from normalize__annData

    output:
    tuple(
        val(outdir),
        file("adata-pcs-harmony.tsv.gz")
    ) into harmony__results
    // val(outdir) into harmony__outDir
    // file("adata-pcs-harmony.tsv.gz") into harmony__results

    script:
    outdir = "${outdir_prev}/reduced_dims-harmony"
    """
    harmony_process_pcs.R \
        --pca_file ${infile_reducedDims} \
        --metadata_file ${infile__metaData} \
        --metadata_columns sanger_sample_id \
        --n_pcs 15
    """
}


process cluster {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir  path: "${outdir}",
                mode: "copy",
                overwrite: "true"

    input:
    tuple(
        val(outdir_prev),
        file(infile_reducedDims),
    ) from harmony__results
    file(infile_annData) from normalize__annData
    // tuple(
    //     val(ignore_0),
    //     file(ignore_1),
    //     file(infile_annData),
    //     file(ignore_2)
    // ) from normalize__results

    output:
    tuple val(outdir), file("test-clustered.tsv.gz") into cluster__results
    tuple val(outdir), file("test-clustered.h5") into cluster__annData
    file("*.pdf") optional true
    file("*.png") optional true

    script:
    outdir = "${outdir_prev}/cluster"
    """
    scanpy_cluster.py \
        --h5_anndata ${infile_annData} \
        --tsv_pcs ${infile_reducedDims} \
        --cluster_method leiden \
        --number_pcs 15 \
        --resolution 1 \
        --output_file test-clustered
    """
    // --output_file
}


// Channel: variables to regress out prior to scaling.
// Channel
//     .fromList( normalize__pcs, harmony__results )
//     .set { umap__reducedDims }
//
// process umap {
//     // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
//     // runs Harmony
//     // ------------------------------------------------------------------------
//     //tag { output_dir }
//     //cache false        // cache results from run
//     scratch false      // use tmp directory
//     echo true          // echo output from script
//     cpus 16            // cpu requirements
//     //memory '50 GB'   // memory requirements
//
//     publishDir  path: "${outdir}",
//                 mode: "copy",
//                 overwrite: "true"
//
//     input:
//     val(outdir_prev) from cluster__outDir
//     file(infile_annData) from cluster__annData
//     each file(infile_reducedDims) from umap__reducedDims
//
//     output:
//     file("*.png") into umap_plots
//
//     script:
//     outdir = "${outdir_prev}"
//     """
//     scanpy_umap.py \
//         --h5_anndata ${infile_annData} \
//         --tsv_pcs ${infile_reducedDims} \
//         --colors_quantitative age \
//         --colors_categorical sanger_sample_id,sex,leiden \
//         --number_pcs 15
//     """
//     // --output_file
// }


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
