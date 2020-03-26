#!/usr/bin/env nextflow

nextflow.preview.dsl=2

VERSION = "0.0.1" // do not edit, controlled by bumpversion


// include modules
include {
    merge_samples;
    normalize_and_pca;
    subset_pcs;
    harmony;
} from "./modules/core.nf"
include {
    umap;
    umap as umap__harmony;
    wf__cluster;
    wf__cluster as wf__cluster_harmony;
} from "./modules/cluster.nf"


// set default parameters
params.output_dir    = "nf-qc_cluster"
params.help          = false
params.reduced_dims__vars_to_regress = ['', 'total_counts,age']
params.reduced_dims__n_pcs = [3, 15]
// params.harmony__metadata_columns = ['sanger_sample_id', 'sanger_sample_id,bead_version']
// params.harmony__thetas = ['0.5', '1,0.5']
params.harmony__metadata_columns = ['sanger_sample_id,bead_version']
params.harmony__thetas = ['1,0.5']
params.cluster__resolutions = [1.0]
params.umap__colors_quantitative = 'age'
params.umap__colors_categorical = 'sanger_sample_id,sex'


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


// Initalize Channels.
// Channel: example init
// Channel
//     .fromPath( params.file_paths_10x )
//     .println()
// Channel: variables to regress out prior to scaling.
reduced_dims__vars_to_regress = Channel
    .fromList(params.reduced_dims__vars_to_regress)


// Run the workflow
workflow {
    main:
        // Merge the samples, perform cell + gene filtering, add metadata
        merge_samples(
            params.output_dir,
            params.file_paths_10x,
            params.file_metadata
        )
        // Normalize, regress (optional), scale, and calculate PCs
        normalize_and_pca(
            params.output_dir,
            merge_samples.out.anndata,
            reduced_dims__vars_to_regress
        )
        // Subset PCs to those for anlaysis
        subset_pcs(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            params.reduced_dims__n_pcs
        )
        // "Correct" PCs using Harmony
        harmony(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            params.reduced_dims__n_pcs,
            params.harmony__metadata_columns,
            params.harmony__thetas
        )
        // Make UMAPs of the reduced dimensions
        umap(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            subset_pcs.out.reduced_dims,
            params.umap__colors_quantitative,
            params.umap__colors_categorical
        )
        umap__harmony(
            harmony.out.outdir,
            harmony.out.anndata,
            harmony.out.reduced_dims,
            params.umap__colors_quantitative,
            params.umap__colors_categorical
        )
        // Cluster the results, varying the resolution.
        // Also, generate UMAPs of the results.
        wf__cluster(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            subset_pcs.out.metadata,
            subset_pcs.out.pcs,
            subset_pcs.out.reduced_dims,
            params.cluster__resolutions
        )
        wf__cluster_harmony(
            harmony.out.outdir,
            harmony.out.anndata,
            harmony.out.metadata,
            harmony.out.pcs,
            harmony.out.reduced_dims,
            params.cluster__resolutions
        )
    // NOTE: One could do publishing in the workflow like so, however
    //       that will not allow one to build the directory structure
    //       depending on the input data call. Therefore, we use publishDir
    //       within a process.
    // publish:
    //     merge_samples.out.anndata to: "${params.output_dir}",
    //         mode: "copy",
    //         overwrite: "true"
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
