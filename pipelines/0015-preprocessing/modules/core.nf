#!/usr/bin/env nextflow

// NOTE: label 'big_mem' may be useful at some point


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process cellbender__rb__get_input_cells {
    // Calculates thresholds for input cells of cellbender__remove_background
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            val(ncells_expected),
            val(ndroplets_include_cellbender)
        )
        val(lower_bound_cell_estimate)
        val(lower_bound_total_droplets_included)
        each total_ndroplets_subtract_factor

    // TODO: see if possible to return value in expected_cells.txt and
    // total_droplets_included.txt rather than the file itself
    output:
        val(outdir, emit: outdir)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            val(ncells_expected),
            val(ndroplets_include_cellbender),
            emit: experiment_data
        )
        path(
            "${outfile}-expected_cells.txt",
            emit: expected_cells
        )
        path(
            "${outfile}-total_droplets_included.txt",
            emit: total_droplets_include
        )
        path("${outfile}-cell_estimate_cutoff.tsv.gz")
        path("${outfile}-total_droplets_cutoff.tsv.gz")
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/${experiment_id}"
        outdir = "${outdir}/umi_count_estimate"
        outdir = "${outdir}-ncells_expected__${ncells_expected}"
        outdir = "${outdir}-lb_cell_estimate__${lower_bound_cell_estimate}"
        outdir = "${outdir}-lb_total_droplets_include__${lower_bound_total_droplets_included}"
        outdir = "${outdir}-ndroplets_subtract_factor__${total_ndroplets_subtract_factor}"
        outfile = "umi_count_estimates"
        cmd__expected_ncells = ""
        if ("${ncells_expected}" != "NA") {
            cmd__expected_ncells = "--expected_ncells ${ncells_expected}"
        }
        cmd__droplets_include = ""
        if ("${ndroplets_include_cellbender}" != "NA") {
            cmd__droplets_include = "--total_ndroplets ${ndroplets_include_cellbender}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__remove_background__get_input_cells: ${process_info}"
        rm -fr plots
        mkdir txd_input
        ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
        ln --physical ${file_10x_features} txd_input/features.tsv.gz
        ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz
        015-get_estimates_from_umi_counts.py \
            --tenxdata_path txd_input \
            --output_file ${outfile} \
            --total_ndroplets_subtract_factor \
                ${total_ndroplets_subtract_factor} \
            --lower_bound_expected_ncells ${lower_bound_cell_estimate} \
            --lower_bound_total_droplets_included \
                ${lower_bound_total_droplets_included} \
            ${cmd__expected_ncells} \
            ${cmd__droplets_include}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process cellbender__remove_background {
    // Remove ambient RNA
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
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            val(ncells_expected),
            val(ndroplets_include_cellbender)
        )
        path(expected_cells)
        path(total_droplets_include)
        each epochs
        each learning_rate

    output:
        val(outdir, emit: outdir)
        path("${outfile}.h5", emit: h5)
        path("${outfile}_filtered.h5", emit: h5_filtered)
        path("${outfile}_cell_barcodes.csv", emit: barcodes)
        path("${outfile}.log", emit: log)
        path("${outfile}_filtered_10x_mtx/barcodes.tsv.gz", emit: tenx_barcodes)
        path("${outfile}_filtered_10x_mtx/features.tsv.gz", emit: tenx_features)
        path("${outfile}_filtered_10x_mtx/matrix.mtx.gz", emit: tenx_matrix)
        tuple(
            val(experiment_id),
            val(ncells_expected),
            val("${outdir}/${outfile}_filtered_10x_mtx"),
            val(epochs),
            val(learning_rate),
            emit: filtered_10x
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    //       USER may need to play with the traing fraction and learning rate
    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}" // /${experiment_id}"
        lr_string = "${learning_rate}".replaceAll("\\.", "pt")
        outfile = "cellbender"
        outfile = "${outfile}-epochs_${epochs}"
        outfile = "${outfile}-learningrate_${lr_string}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__remove_background: ${process_info}"
        rm -fr plots
        mkdir txd_input
        ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
        ln --physical ${file_10x_features} txd_input/features.tsv.gz
        ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz
        cellbender remove-background \
            --input txd_input \
            --output ${outfile} \
            --cuda \
            --expected-cells \$(cat ${expected_cells}) \
            --total-droplets-included \$(cat ${total_droplets_include}) \
            --model full \
            --z-dim 200 \
            --z-layers 1000 \
            --low-count-threshold 10 \
            --epochs ${epochs} \
            --empty-drop-training-fraction 0.5 \
            --learning-rate ${learning_rate}
        [ -f ${outfile} ] && mv ${outfile} ${outfile}.h5
        convert-cellbender_10xmtx.py \
            -h5 ${outfile}_filtered.h5 \
            -g background_removed \
            -od ${outfile}_filtered_10x_mtx
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process cellbender__gather_qc_input {
    // Prepare cell bender output for qc_cluster pipeline. For each epoch and
    // learning rate, gather 10x matrix files into format for qc_cluster
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(epoch_lr_key),
            val(experiment_ids),
            val(tenx_mtx_paths),
            val(ncells_expected)
        )

    output:
        val(outdir, emit: outdir)
        path(outfile, emit: formatted_input)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/qc_cluster_input_files"
        outfile = "${epoch_lr_key}.tsv"
        experiment_ids = experiment_ids.join(",")
        tenx_mtx_paths = tenx_mtx_paths.join(",")
        ncells_expected = ncells_expected.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__gather_qc_input: ${process_info}"
        echo "Epoch and Learning rate key: ${epoch_lr_key}"
        echo "Experiment ids: ${experiment_ids}"
        echo "Matrix paths: ${tenx_mtx_paths}"
        echo "Number of cells expected: ${ncells_expected}"
        035-prepare_qc_cluster_input.py \
            -id ${experiment_ids} \
            -dir ${tenx_mtx_paths} \
            -n ${ncells_expected} \
            -of ${outfile}
        """
}
