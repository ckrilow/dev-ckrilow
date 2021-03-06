#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    cellbender__rb__get_input_cells;
    cellbender__remove_background;
    cellbender__gather_qc_input;
} from "./modules/core.nf"


// Set default parameters.
params.output_dir           = "nf-preprocessing"
params.help                 = false
// Default parameters for cellbender remove_background
params.cellbender_rb = [
    total_ndroplets_subtract_factor: [value: [65000]],
    epochs: [value: [200]],
    learning_rate: [value: [0.001, 0.0001]]
]

// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell preprocessing ~ v${VERSION}
    ============================================================================

    Runs basic single cell preprocessing

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --file_paths_10x    Tab-delimited file containing experiment_id and
                            path_data_10xformat columns.

    Optional arguments:
        --output_dir        Directory name to save results to. (Defaults to
                            '${params.output_dir}')

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml.

    Profiles:
        local               local execution
        lsf                 lsf cluster execution
    """.stripIndent()
}


// Boot message - either help message or the parameters supplied.
if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell preprocessing ~ v${VERSION}
    ============================================================================
    file_paths_10x                : ${params.file_paths_10x}
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
        file("${row.data_path_10x_format}/matrix.mtx.gz"),
        row.ncells_expected,
        row.ndroplets_include_cellbender
    )}
//n_experiments = file(params.file_paths_10x).countLines()


// Run the workflow.
workflow {
    main:
        // Prep the data for cellbender
        cellbender__rb__get_input_cells(
            params.output_dir,
            channel__file_paths_10x,
            100, // lower_bound_cell_estimate
            10, // lower_bound_total_droplets_included
            params.cellbender_rb.total_ndroplets_subtract_factor.value
        )
        // Correct counts matrix to remove ambient RNA
        cellbender__remove_background(
            cellbender__rb__get_input_cells.out.outdir,
            cellbender__rb__get_input_cells.out.experiment_data,
            cellbender__rb__get_input_cells.out.expected_cells,
            cellbender__rb__get_input_cells.out.total_droplets_include,
            params.cellbender_rb.epochs.value,
            params.cellbender_rb.learning_rate.value
        )
        gather_qc_input = cellbender__remove_background.out.filtered_10x.groupTuple(by: 2)
            .reduce([:]) { map, tuple ->
                def map_key = "epoch_" + tuple[3][0] + "_lr_" + tuple[4][0]
                def key_list = map[map_key]
                if (!key_list) {
                    key_list = [[tuple[0][0], tuple[2], tuple[1][0]]]
                } else {
                    key_list.add([tuple[0][0], tuple[2], tuple[1][0]])
                }
                map[map_key] = key_list
                return(map)
            }
            .flatMap()
            .map {  entry ->
                combined_data = [entry.key, [], [], []]
                entry.value.each {
                    combined_data[1].add(it[0])
                    combined_data[2].add(it[1])
                    combined_data[3].add(it[2])
                }
                return(combined_data)
            }

        cellbender__gather_qc_input(
            params.output_dir,
            gather_qc_input
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
