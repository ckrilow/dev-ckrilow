#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    wf__multiplet;
} from "./modules/multiplet.nf"
include {
    prep_merge_samples;
    merge_samples;
    plot_predicted_sex;
    plot_qc;
    estimate_pca_elbow;
    normalize_and_pca;
    subset_pcs;
    plot_pcs;
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
// Default key to add in metadata
params.metadata_key_column = [
    value: "experiment_id"
]
// Default parameters for qc plots.
params.plots_qc = [
    facet_columns: [value: ["experiment_id"]],
    variable_columns_distribution_plots: [value: [
        "total_counts,pct_counts_mito_gene"
    ]]
]
// Default parameters for reduced dimension calculations.
// run_downstream_analysis: If false don't run clustering or umaps
params.reduced_dims = [
    run_downstream_analysis: false,
    vars_to_regress: [value: ["", "total_counts,pct_counts_mito_gene"]],
    n_dims: [value: [15, 30]]
]
// Default parameters for harmony.
params.harmony = [
    run_process: false,
    variables_and_thetas: [value: [
        [variable: "experiment_id", theta: "1.0"],
        [variable: "experiment_id,phase", theta: "1.0,0.2"]
    ]]
]
// Default parameters for bbknn
params.bbknn = [
    run_process: false,
    batch_variable: [value: ["experiment_id"]]
]
// Default parameters for lisi
params.lisi = [
    run_process: false,
    variables: [value: ["experiment_id,phase"]]
]
// Default parameters for cluster calculations.
params.cluster = [
    number_neighbors: [value: [15, 20]],
    methods: [value: ["leiden"]],
    resolutions: [value: [1.0, 3.0]],
    variables_boxplot: [value: ["total_reads,pct_counts_mito_gene"]],
    known_markers: [run_process: false]
]
// Default parameters for cluster resolution validation.
params.cluster_validate_resolution = [
    sparsity: [value: [0.0001]],
    train_size_cells: [value: [-1]]
]
// Default parameters for cluster marker gene calculations.
params.cluster_marker = [
    methods: [value: ["wilcoxon"]]
]
// Default parameters for umap calculations.
params.umap = [
    run_process: true,
    colors_quantitative: [value: "total_counts"],
    colors_categorical: [value: "experiment_id,phase"],
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
//n_experiments = file(params.file_paths_10x).countLines()

// Initialize known markers channel
// cluster__known_markers is a list of tsv files, first serialize
// the array then run plot_known_markers
if (params.cluster.known_markers.run_process) {
    channel__cluster__known_markers = Channel
        .fromList(params.cluster.known_markers.value)
        .map{row -> tuple(row.file_id, file(row.file))}
} else {
    channel__cluster__known_markers = tuple('', '')
}


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
        prep_merge_samples(channel__file_paths_10x)
        merge_samples(
            params.output_dir,
            params.file_paths_10x,
            params.file_metadata,
            file_sample_qc,
            file_cellmetadata,
            params.metadata_key_column.value,
            prep_merge_samples.out.barcodes.collect(),
            prep_merge_samples.out.features.collect(),
            prep_merge_samples.out.matrix.collect()
        )
        // Predict sex from gene expression and check against phenotypes.
        plot_predicted_sex(
            params.output_dir,
            merge_samples.out.anndata
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
        // Estimate number of PCs to use using eblow from variance explained
        estimate_pca_elbow(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata
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
        if (params.harmony.run_process) {
            harmony(
                normalize_and_pca.out.outdir,
                normalize_and_pca.out.anndata,
                normalize_and_pca.out.metadata,
                normalize_and_pca.out.pcs,
                normalize_and_pca.out.param_details,
                params.reduced_dims.n_dims.value,
                params.harmony.variables_and_thetas.value
            )
        }
        // Run BBKNN
        if (params.bbknn.run_process) {
            bbknn(
                normalize_and_pca.out.outdir,
                normalize_and_pca.out.anndata,
                normalize_and_pca.out.metadata,
                normalize_and_pca.out.pcs,
                normalize_and_pca.out.param_details,
                params.reduced_dims.n_dims.value,
                params.bbknn.batch_variable.value
            )
        }
        // TODO: There is a bug below where lisi will be called for each
        // normalize_and_pca call. It just means there will be some duplicate
        // output files in each normalize_and_pca dir and a bit of wasted CPU.
        if (params.lisi.run_process) {
            lisi_input = subset_pcs.out.reduced_dims_params.collect()
            if (params.harmony.run_process) {
                lisi_input = lisi_input.mix(
                    harmony.out.reduced_dims_params.collect()
                )
            }
            if (params.bbknn.run_process) {
                lisi_input = lisi_input.mix(
                    bbknn.out.reduced_dims_params.collect()
                )
            }
            lisi(
                normalize_and_pca.out.outdir,
                normalize_and_pca.out.metadata,
                params.lisi.variables.value,
                lisi_input.collect()
            )
        }
        // Scatter-gather UMAP plots
        if (
            params.reduced_dims.run_downstream_analysis &
            params.umap.run_process
        ) {
            wf__umap(
                subset_pcs.out.outdir,
                subset_pcs.out.anndata,
                subset_pcs.out.metadata,
                subset_pcs.out.pcs,
                subset_pcs.out.reduced_dims,
                "False",
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
            // Use the data with all of the umaps calculated for downstream
            // clustering, so that we have all of the umap dims in adata.
            cluster_subset_pcs__outdir = wf__umap.out.outdir
            cluster_subset_pcs__anndata = wf__umap.out.anndata
            cluster_subset_pcs__metadata = wf__umap.out.metadata
            cluster_subset_pcs__pcs = wf__umap.out.pcs
            cluster_subset_pcs__reduced_dims = wf__umap.out.reduced_dims
        } else if (params.reduced_dims.run_downstream_analysis) {
            // If running downstream analysis and no umaps, set input for
            // downstream analysis
            cluster_subset_pcs__outdir = subset_pcs.out.outdir
            cluster_subset_pcs__anndata = subset_pcs.out.anndata
            cluster_subset_pcs__metadata = subset_pcs.out.metadata
            cluster_subset_pcs__pcs = subset_pcs.out.pcs
            cluster_subset_pcs__reduced_dims = subset_pcs.out.reduced_dims
        }
        if (params.harmony.run_process & params.umap.run_process) {
            wf__umap_harmony(
                harmony.out.outdir,
                harmony.out.anndata,
                harmony.out.metadata,
                harmony.out.pcs,
                harmony.out.reduced_dims,
                "False",
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
            // Use the data with all of the umaps calculated for downstream
            // clustering, so that we have all of the umap dims in adata.
            cluster_harmony__outdir = wf__umap_harmony.out.outdir
            cluster_harmony__anndata = wf__umap_harmony.out.anndata
            cluster_harmony__metadata = wf__umap_harmony.out.metadata
            cluster_harmony__pcs = wf__umap_harmony.out.pcs
            cluster_harmony__reduced_dims = wf__umap_harmony.out.reduced_dims
        } else if (params.harmony.run_process) {
            cluster_harmony__outdir = harmony.out.outdir
            cluster_harmony__anndata = harmony.out.anndata
            cluster_harmony__metadata = harmony.out.metadata
            cluster_harmony__pcs = harmony.out.pcs
            cluster_harmony__reduced_dims = harmony.out.reduced_dims
        }
        // NOTE: for BBKNN, we specifically pass the PCs to the reduced dims
        ///      slot not the UMAPS.
        // NOTE: for BBKNN n_neighbors is not needed since already calculated
        if (params.bbknn.run_process & params.umap.run_process) {
            wf__umap_bbknn(
                bbknn.out.outdir,
                bbknn.out.anndata,
                bbknn.out.metadata,
                bbknn.out.pcs,
                bbknn.out.reduced_dims,
                "True",  // Don't look at the reduced_dims parameter
                ["-1"],  // params.cluster.number_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
            // Use the data with all of the umaps calculated for downstream
            // clustering, so that we have all of the umap dims in adata.
            cluster_bbknn__outdir = wf__umap_bbknn.out.outdir
            cluster_bbknn__anndata = wf__umap_bbknn.out.anndata
            cluster_bbknn__metadata = wf__umap_bbknn.out.metadata
            cluster_bbknn__pcs = wf__umap_bbknn.out.pcs
            cluster_bbknn__reduced_dims = wf__umap_bbknn.out.reduced_dims
        } else if (params.bbknn.run_process) {
            cluster_bbknn__outdir = bbknn.out.outdir
            cluster_bbknn__anndata = bbknn.out.anndata
            cluster_bbknn__metadata = bbknn.out.metadata
            cluster_bbknn__pcs = bbknn.out.pcs
            cluster_bbknn__reduced_dims = bbknn.out.reduced_dims
        }
        // START LEGACY CODE ----------------------------------------------
        // NOTE: Legacy code due to it being hard to compare
        // Make UMAPs of the reduced dimensions - no scatter gather.
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
        // END LEGACY CODE ------------------------------------------------
        // Cluster the results, varying the resolution.
        // Also, generate UMAPs of the results.i
        if (params.reduced_dims.run_downstream_analysis) {
            wf__cluster(
                cluster_subset_pcs__outdir,
                cluster_subset_pcs__anndata,
                cluster_subset_pcs__metadata,
                cluster_subset_pcs__pcs,
                cluster_subset_pcs__reduced_dims,
                "False",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value
            )
            plot_pcs(
                cluster_subset_pcs__outdir,
                cluster_subset_pcs__anndata,
                params.reduced_dims.n_dims.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
        }
        if (params.harmony.run_process) {
            wf__cluster_harmony(
                cluster_harmony__outdir,
                cluster_harmony__anndata,
                cluster_harmony__metadata,
                cluster_harmony__pcs,
                cluster_harmony__reduced_dims,
                "False",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value
            )
        }
        if (params.bbknn.run_process) {
            wf__cluster_bbknn(
                cluster_bbknn__outdir,
                cluster_bbknn__anndata,
                cluster_bbknn__metadata,
                cluster_bbknn__pcs,
                cluster_bbknn__reduced_dims,
                "True",  // use_pcs_as_reduced_dims
                ["-1"],  // params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                ["-1"],  // params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value
            )
        }
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
