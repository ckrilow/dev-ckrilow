#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    merge_samples;
    normalize_and_pca;
    subset_pcs;
    harmony;
} from "./modules/core.nf"
include {
    convert_seurat;
    umap;
    umap as umap__harmony;
    wf__cluster;
    wf__cluster as wf__cluster_harmony;
} from "./modules/cluster.nf"


// Set default parameters.
params.output_dir    = "nf-qc_cluster"
params.help          = false
// NOTE: The default parameters below were chosen to show the flexiblity of
//       this pipeline. They were not chosen because these are the values one
//       should use for final analysis.
// Default parameters for reduced dimension calculations.
params.reduced_dims = [
    vars_to_regress: [value: ["", "total_counts,age"]],
    n_dims: [value: [15, 30]]
]
// Default parameters for harmony.
params.harmony = [
    variables_and_thetas: [value: [
        [variable: "sanger_sample_id", theta: "1.0"],
        [variable: "sanger_sample_id,bead_version", theta: "1.0,0.2"]
    ]]
]
// Default parameters for cluster calculations.
params.cluster = [
    methods: [value: ["leiden", "louvain"]],
    resolutions: [value: [1.0, 3.0]],
]
// Default parameters for cluster marker gene calculations.
params.cluster_marker = [
    methods: [value: ["wilcoxon", "logreg"]]
]
// Default parameters for umap calculations.
params.umap = [
    colors_quantitative: [value: "age"],
    colors_categorical: [value: "sanger_sample_id,sex"],
]


// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell qc and clustering ~ v${VERSION}
    ============================================================================

    Runs basic single cell qc and clustering

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --file_paths_10x    Tab-delimited file containing experiment_id and
                            path_data_10xformat columns.

        --file_metadata     Tab-delimited file containing sample metadata.

        --file_sample_qc    YAML file containing sample quality control
                            filters.

        --variable_genes_exclude
                            Tab-delimited file with genes to exclude from
                            highly variable gene list. Must contain
                            ensembl_gene_id column. If no filter, then pass an
                            empty file.

    Other arguments:
        --output_dir        Directory name to save results to. (Defaults to
                            'nf-qc_cluster')

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml

    Profiles:
        local               local execution
        lsf                 lsf cluster execution
        local_singularity   local execution with singularity [TODO]
        lsf_singularity     lsf cluster execution with singularity [TODO]
    """.stripIndent()
}


// Boot message - either help message or the parameters supplied.
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
    file_sample_qc                : ${params.file_sample_qc}
    variable_genes_exclude        : ${params.variable_genes_exclude}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
    // A dictionary way to accomplish the text above.
    // def summary = [:]
    // summary['file_paths_10x'] = params.file_paths_10x
    // log.info summary.collect { k,v -> "${k.padRight(20)} : $v" }.join("\n")
}


// Initalize Channels.
// Channel: example init
// Channel
//     .fromPath( params.file_paths_10x )
//     .println()
// Channel: required files
// file_paths_10x = Channel
//     .fromPath(params.file_paths_10x)
// file_metadata = Channel
//     .fromPath(params.file_metadata)
// file_sample_qc = Channel
//     .fromPath(params.file_sample_qc)
// Channel: variables to regress out prior to scaling.
// reduced_dims__vars_to_regress = Channel
//     .fromList(params.reduced_dims__vars_to_regress)
// harmony__variables_and_thetas = Channel
//     .from(params.harmony.variables_and_thetas.value)


// Run the workflow.
workflow {
    main:
        // Merge the samples, perform cell + gene filtering, add metadata.
        merge_samples(
            params.output_dir,
            params.file_paths_10x,
            params.file_metadata,
            params.file_sample_qc
        )
        // Normalize, regress (optional), scale, and calculate PCs
        normalize_and_pca(
            params.output_dir,
            merge_samples.out.anndata,
            params.variable_genes_exclude,
            params.reduced_dims.vars_to_regress.value
        )
        // Make Seurat dataframes of the normalized anndata
        convert_seurat(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata
        )
        // Subset PCs to those for anlaysis
        subset_pcs(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            params.reduced_dims.n_dims.value
        )
        // "Correct" PCs using Harmony
        harmony(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            params.reduced_dims.n_dims.value,
            params.harmony.variables_and_thetas.value
        )
        // Make UMAPs of the reduced dimensions
        umap(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            subset_pcs.out.reduced_dims,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        umap__harmony(
            harmony.out.outdir,
            harmony.out.anndata,
            harmony.out.reduced_dims,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        // Cluster the results, varying the resolution.
        // Also, generate UMAPs of the results.
        wf__cluster(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            subset_pcs.out.metadata,
            subset_pcs.out.pcs,
            subset_pcs.out.reduced_dims,
            params.cluster.methods.value,
            params.cluster.resolutions.value,
            params.cluster_marker.methods.value
        )
        wf__cluster_harmony(
            harmony.out.outdir,
            harmony.out.anndata,
            harmony.out.metadata,
            harmony.out.pcs,
            harmony.out.reduced_dims,
            params.cluster.methods.value,
            params.cluster.resolutions.value,
            params.cluster_marker.methods.value
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
