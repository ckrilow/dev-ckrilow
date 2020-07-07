#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    cellphonedb__make_database;
    cellphonedb__run_analysis;
    cellphonedb__plot_dot;
    cellphonedb__plot_heatmap;
} from "./modules/core.nf"


// Set default parameters.
params.output_dir           = "nf-cell_cell_interaction"
params.help                 = false
// Default parameters for cellbender remove_background
params.cellphonedb = [
    database: [value: ''],
    database_version: [value: ''],
    iterations: [value: [1000, 10000]],
    threshold: [value: [1]]
]

// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     cell cell interactions ~ v${VERSION}
    ============================================================================

    Runs cell cell interaction analysis from single cell data.

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --h5ad              Scanpy anndata file.

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
     cell cell interactions ~ v${VERSION}
    ============================================================================
    adata_file                    : ${params.adata_file}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
    // A dictionary way to accomplish the text above.
    // def summary = [:]
    // summary['file_paths_10x'] = params.file_paths_10x
    // log.info summary.collect { k,v -> "${k.padRight(20)} : $v" }.join("\n")
}


// Run the workflow.
workflow {
    main:
        // If we do not have a database file, then attempt to download one
        if (params.cellphonedb.database.value == "") {
            make_cellphonedb_database(
                params.output_dir,
                params.cellphonedb.database_version.value
            )
            file_db = make_cellphonedb_database.out.database
        } else {
            file_db = file(params.cellphonedb.database.value)
        }
        // Now run cellphonedb
        cellphonedb__run_analysis(
            params.output_dir,
            file(params.h5ad),
            file_db,
            params.cellphonedb.iterations.value,
            params.cellphonedb.threshold.value
        )
        // Make plots of the results
        cellphonedb__plot_dot(
            cellphonedb__run_analysis.out.outdir,
            cellphonedb__run_analysis.out.pvalues,
            cellphonedb__run_analysis.out.means
        )
        // cellphonedb__plot_heatmap(
        //     cellphonedb__run_analysis.out.outdir,
        //     cellphonedb__run_analysis.out.pvalues,
        //     0.05
        // )
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
