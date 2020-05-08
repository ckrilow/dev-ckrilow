#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    wf__multiplet;
} from "./modules/multiplet.nf"
include {
    merge_samples;
    plot_qc;
    normalize_and_pca;
    subset_pcs;
    harmony;
    bbknn;
    lisi;
} from "./modules/core.nf"
include {
    convert_seurat;
    wf__cluster;
    wf__cluster as wf__cluster_harmony;
    wf__cluster as wf__cluster_bbknn;
} from "./modules/cluster.nf"
include {
    wf__umap;
    wf__umap as wf__umap_harmony;
    wf__umap as wf__umap_bbknn;
    // umap_calculate_and_plot;
    // umap_calculate_and_plot as umap_calculate_and_plot__harmony;
} from "./modules/umap.nf"


// Set default parameters.
params.output_dir           = "nf-qc_cluster"
params.help                 = false
params.run_multiplet        = false
params.file_sample_qc       = "no_file__file_sample_qc"
params.file_cellmetadata    = "no_file__file_cellmetadata"
params.genes_exclude_hvg    = "no_file__genes_exclude_hvg"
params.genes_score          = "no_file__genes_score"
// NOTE: The default parameters below were chosen to show the flexiblity of
//       this pipeline. They were not chosen because these are the values one
//       should use for final analysis.
// Default parameters for qc plots.
params.plots_qc = [
    facet_columns: [value: ["experiment_id"]],
    variable_columns_distribution_plots: [value: [
        "total_counts,pct_counts_mito_gene"
    ]],
]
// Default parameters for reduced dimension calculations.
params.reduced_dims = [
    vars_to_regress: [value: ["", "total_counts,age"]],
    n_dims: [value: [15, 30]]
]
// Default parameters for harmony.
params.harmony = [
    variables_and_thetas: [value: [
        [variable: "experiment_id", theta: "1.0"],
        [variable: "experiment_id,bead_version", theta: "1.0,0.2"]
    ]]
]
// Default parameters for lisi
params.lisi = [
    variables: [value: ["experiment_id,bead_version"]]
]
// Default parameters for cluster calculations.
params.cluster = [
    number_neighbors: [value: [15, 20]],
    methods: [value: ["leiden", "louvain"]],
    resolutions: [value: [1.0, 3.0]]
]
// Default parameters for cluster resolution validation.
params.cluster_validate_resolution = [
    sparsity: [value: [0.25, 0.1]],
    train_size_cells: [value: [10000, -1]],
    //number_cells: [value: [10000, -1]],
    //train_size_fraction: [value: [0.33]]
]
// Default parameters for cluster marker gene calculations.
params.cluster_marker = [
    methods: [value: ["wilcoxon", "logreg"]]
]
// Default parameters for umap calculations.
params.umap = [
    colors_quantitative: [value: "age"],
    colors_categorical: [value: "experiment_id,sex"],
    n_neighbors: [value: [15, 30]],
    umap_init: [value: ["X_pca", "spectral"]],
    umap_min_dist: [value: [0.5]],
    umap_spread: [value: [1.0, 5.0]]
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

    Optional arguments:
        --file_sample_qc    YAML file containing sample quality control
                            filters.

        --file_cellmetadata Tab-delimited file containing experiment_id and
                            data_path_cellmetadata columns.

        --genes_exclude_hvg
                            Tab-delimited file with genes to exclude from
                            highly variable gene list. Must contain
                            ensembl_gene_id column. If no filter, then pass an
                            empty file.

        --genes_score
                            Tab-delimited file with genes to use to score
                            cells. Must contain ensembl_gene_id and score_id
                            columns.  If one score_id == "cell_cycle", then
                            requires a grouping_id column with "G2/M" and "S".
                            If no filter, then pass an empty file.

        --output_dir        Directory name to save results to. (Defaults to
                            'nf-qc_cluster')

        --run_multiplet     Flag to run multiplet analysis. Output from
                            multiplet analysis will be added to the final
                            AnnData object. The flag only works if
                            file_cellmetadata is not set.

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml.

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
    file_cellmetadata             : ${params.file_cellmetadata}
    genes_exclude_hvg             : ${params.genes_exclude_hvg}
    genes_score                   : ${params.genes_score}
    output_dir (output folder)    : ${params.output_dir}
    run_multiplet                 : ${params.run_multiplet}
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
// Channel
//     .fromPath(params.file_paths_10x)
//     .splitCsv(header: true, sep: "\t", by: 1)
//     .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format))}
//     .view()
channel__file_paths_10x = Channel
    .fromPath(params.file_paths_10x)
    .splitCsv(header: true, sep: "\t", by: 1)
    .map{row -> tuple(
        row.experiment_id,
        file("${row.data_path_10x_format}/barcodes.tsv.gz"),
        file("${row.data_path_10x_format}/features.tsv.gz"),
        file("${row.data_path_10x_format}/matrix.mtx.gz")
    )}
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
        // Optionally run multiplet filters.
        if (
            params.run_multiplet &
            params.file_cellmetadata == "no_file__file_cellmetadata"
        ) {
            wf__multiplet(
                params.output_dir,
                channel__file_paths_10x
            )
            // NOTE: file__cellmetadata is already defined as path, so no need
            // to call file or path.
            file_cellmetadata = wf__multiplet.out.file__cellmetadata
        } else {
            // For some reason cannot use path here, must use file.
            file_cellmetadata = file(params.file_cellmetadata)
        }
        // Merge the samples, perform cell + gene filtering, add metadata.
        file_sample_qc = file(params.file_sample_qc)
        merge_samples(
            params.output_dir,
            params.file_paths_10x,
            params.file_metadata,
            file_sample_qc,
            file_cellmetadata,
            "sanger_sample_id"
        )
        // Make QC plots of the merged data.
        plot_qc(
            params.output_dir,
            merge_samples.out.anndata,
            params.plots_qc.facet_columns.value,
            params.plots_qc.variable_columns_distribution_plots.value
        )
        // Normalize, regress (optional), scale, and calculate PCs.
        genes_exclude_hvg = file(params.genes_exclude_hvg)
        genes_score = file(params.genes_score)
        normalize_and_pca(
            params.output_dir,
            merge_samples.out.anndata,
            genes_exclude_hvg,
            genes_score,
            params.reduced_dims.vars_to_regress.value
        )
        // Make Seurat dataframes of the normalized anndata
        // convert_seurat(
        //     normalize_and_pca.out.outdir,
        //     normalize_and_pca.out.anndata
        // )
        // Subset PCs to those for anlaysis
        subset_pcs(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            normalize_and_pca.out.param_details,
            params.reduced_dims.n_dims.value
        )
        // "Correct" PCs using Harmony
        harmony(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            normalize_and_pca.out.param_details,
            params.reduced_dims.n_dims.value,
            params.harmony.variables_and_thetas.value
        )
        // Run BBKNN
        bbknn(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            normalize_and_pca.out.param_details,
            params.reduced_dims.n_dims.value,
            'experiment_id'
        )
        // TODO: There is a bug below where lisi will be called for each
        // normalize_and_pca call. It just means there will be some duplicate
        // output files in each normalize_and_pca dir and a bit of wasted CPU.
        lisi_input = subset_pcs.out.reduced_dims_params.collect().mix(
            harmony.out.reduced_dims_params.collect()
        ).mix(
            bbknn.out.reduced_dims_params.collect()
        )
        lisi(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.metadata,
            params.lisi.variables.value,
            lisi_input.collect()
        )
        // Scatter-gather UMAP plots
        wf__umap(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            // subset_pcs.out.metadata,
            // subset_pcs.out.pcs,
            subset_pcs.out.reduced_dims,
            params.umap.n_neighbors.value,
            params.umap.umap_init.value,
            params.umap.umap_min_dist.value,
            params.umap.umap_spread.value,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        wf__umap_harmony(
            harmony.out.outdir,
            harmony.out.anndata,
            // harmony.out.metadata,
            // harmony.out.pcs,
            harmony.out.reduced_dims,
            params.umap.n_neighbors.value,
            params.umap.umap_init.value,
            params.umap.umap_min_dist.value,
            params.umap.umap_spread.value,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        // NOTE: for BBKNN, we specifically pass the PCs to the reduced dims
        ///      slot not the UMAPS.
        // NOTE: for BBKNN n_neighbors is not needed since already calculated
        wf__umap_bbknn(
            bbknn.out.outdir,
            bbknn.out.anndata,
            // bbknn.out.metadata,
            bbknn.out.pcs,
            // bbknn.out.reduced_dims,
            ['-1'],  // params.cluster.number_neighbors.value,
            params.umap.umap_init.value,
            params.umap.umap_min_dist.value,
            params.umap.umap_spread.value,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        // Make UMAPs of the reduced dimensions - no scatter gather
        // umap_calculate_and_plot(
        //     subset_pcs.out.outdir,
        //     subset_pcs.out.anndata,
        //     subset_pcs.out.reduced_dims,
        //     params.umap.colors_quantitative.value,
        //     params.umap.colors_categorical.value,
        //     params.umap.n_neighbors.value,
        //     params.umap.umap_init.value,
        //     params.umap.umap_min_dist.value,
        //     params.umap.umap_spread.value
        // )
        // umap_calculate_and_plot__harmony(
        //     harmony.out.outdir,
        //     harmony.out.anndata,
        //     harmony.out.reduced_dims,
        //     params.umap.colors_quantitative.value,
        //     params.umap.colors_categorical.value,
        //     params.umap.n_neighbors.value,
        //     params.umap.umap_init.value,
        //     params.umap.umap_min_dist.value,
        //     params.umap.umap_spread.value
        // )
        // Cluster the results, varying the resolution.
        // Also, generate UMAPs of the results.
        wf__cluster(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            subset_pcs.out.metadata,
            subset_pcs.out.pcs,
            subset_pcs.out.reduced_dims,
            params.cluster.number_neighbors.value,
            params.cluster.methods.value,
            params.cluster.resolutions.value,
            params.cluster_validate_resolution.sparsity.value,
            params.cluster_validate_resolution.train_size_cells.value,
            // params.cluster_validate_resolution.number_cells.value,
            // params.cluster_validate_resolution.train_size_fraction.value,
            params.cluster_marker.methods.value,
            params.umap.n_neighbors.value,
            params.umap.umap_init.value,
            params.umap.umap_min_dist.value,
            params.umap.umap_spread.value
        )
        wf__cluster_harmony(
            harmony.out.outdir,
            harmony.out.anndata,
            harmony.out.metadata,
            harmony.out.pcs,
            harmony.out.reduced_dims,
            params.cluster.number_neighbors.value,
            params.cluster.methods.value,
            params.cluster.resolutions.value,
            params.cluster_validate_resolution.sparsity.value,
            params.cluster_validate_resolution.train_size_cells.value,
            // params.cluster_validate_resolution.number_cells.value,
            // params.cluster_validate_resolution.train_size_fraction.value,
            params.cluster_marker.methods.value,
            params.umap.n_neighbors.value,
            params.umap.umap_init.value,
            params.umap.umap_min_dist.value,
            params.umap.umap_spread.value
        )
        wf__cluster_bbknn(
            bbknn.out.outdir,
            bbknn.out.anndata,
            bbknn.out.metadata,
            bbknn.out.pcs,
            bbknn.out.reduced_dims,
            ['-1'],  // params.cluster.number_neighbors.value,
            params.cluster.methods.value,
            params.cluster.resolutions.value,
            params.cluster_validate_resolution.sparsity.value,
            params.cluster_validate_resolution.train_size_cells.value,
            // params.cluster_validate_resolution.number_cells.value,
            // params.cluster_validate_resolution.train_size_fraction.value,
            params.cluster_marker.methods.value,
            ['-1'],  // params.umap.n_neighbors.value,
            params.umap.umap_init.value,
            params.umap.umap_min_dist.value,
            params.umap.umap_spread.value
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
