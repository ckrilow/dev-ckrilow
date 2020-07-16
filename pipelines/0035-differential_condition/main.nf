#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    get_cell_label_list;
    run_diffxpy;
} from "./modules/core.nf"


// Set default parameters.
params.output_dir           = "nf-differential_condition"
params.help                 = false
// Default parameters for cellbender remove_background
params.anndata_cell_label = [value: 'cluster']
params.differential_condition = [value: ['cd_status']]


// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell differential condition ~ v${VERSION}
    ============================================================================

    Runs basic single cell preprocessing

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --file_anndata      Anndata file with cell type labels.

    Optional arguments:
        --output_dir        Directory name to save results to. (Defaults to
                            '${params.output_dir}')

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml.

    Profiles:
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
     single cell differential condition ~ v${VERSION}
    ============================================================================
    file_anndata                  : ${params.file_anndata}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
}


// Initalize Channels.
// anndata = Channel
//     .fromPath(params.file_anndata)


// Run the workflow.
workflow {
    main:
        // Get a list of all of the cell types
        get_cell_label_list(
            params.file_anndata,
            params.anndata_cell_label.value
        )
        // For each cell type run differential expression
        cell_labels = get_cell_label_list.out.cell_labels
            .splitCsv(header: false, sep: ',')
        run_diffxpy(
            params.output_dir,
            params.file_anndata,
            params.differential_expression.condition.value,
            params.differential_expression.covariates.value,
            params.anndata_cell_label.value,
            '15',
            params.differential_expression.diffxpy_method.value
        )
        // TODO: for each condition... for each covariate ... for each method
        // merge the output across all clusters 
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
