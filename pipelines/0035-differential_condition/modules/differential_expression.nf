#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process get_cell_label_list {
    // Get all of the cell labels in an anndata file
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    input:
        path(anndata)
        val(anndata_cell_label)

    output:
        path("cell_labels.csv", emit: cell_labels)

    script:
        runid = random_hex(16)
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "get_cell_label_list: ${process_info}"
        013-get_cell_label_list.py \
            --h5_anndata ${anndata} \
            --cell_label ${anndata_cell_label}
        """
}


process run_diffxpy {
    // Run diffxpy
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    label 'gpu'          // use GPU
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(cell_label_column)
        each cell_label_analyse
        each condition_column
        each covariate_columns
        each method

    output:
        val(outdir, emit: outdir)
        path("${outfile}-de_results.tsv.gz", emit: results)
        path("${outfile}-de_results_obj.joblib.gz", emit: results_obj)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        outfile = "test"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_diffxpy: ${process_info}"
        rm -fr plots
        015-run_diffxpy.py \
            --h5_anndata ${anndata} \
            --condition_column ${condition_column} \
            --covariate_columns ${covariate_columns} \
            --cell_label_column ${cell_label_column} \
            --cell_label_analyse ${cell_label_analyse} \
            --method ${method} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


workflow wf__differential_expression {
    take:
        outdir
        anndata
        anndata_cell_label
        condition
        covariates
        diffxpy_method
    main:
        // Get a list of all of the cell types
        get_cell_label_list(
            anndata,
            anndata_cell_label
        )
        // For each cell type compute differential expression for that cell
        // type
        cell_labels = get_cell_label_list.out.cell_labels
            .splitCsv(header: false, sep: ',')
            // .map{row -> tuple(
            //     row[0]
            // )}
        cell_labels.view()
        // Run diffxpy with all combinations of conditions (e.g., sex,
        // disease status), covariates (e.g., size_factors, age),
        // methods (e.g., wald)
        run_diffxpy(
            outdir,
            anndata,
            anndata_cell_label,
            // '1',  // just run on first cluster for development
            cell_labels,  // run for all clusters for run time
            condition,
            covariates,
            diffxpy_method
        )
        // TODO
        // For each condition... merge results across covariates, methods, and
        // celltypes
    emit:
        cell_labels = get_cell_label_list.out.cell_labels
}
