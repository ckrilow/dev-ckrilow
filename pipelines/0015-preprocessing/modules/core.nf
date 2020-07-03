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

    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            val(ncells_expected)
        )
        val(lower_bound_cell_estimate)
        val(lower_bound_total_droplets_included)

    // TODO: see if possible to return value in expected_cells.txt and
    // total_droplets_included.txt rather than the file itself
    output:
        val(outdir, emit: outdir)
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
        outfile = "cellbender-umi_count_estimates"
        outfile = "${outfile}-lower_bound_cell_estimate__${lower_bound_cell_estimate}"
        outfile = "${outfile}-lower_bound_total_droplets_included__${lower_bound_total_droplets_included}"
        cmd__expected_ncells = ""
        if ("${ncells_expected}" != "NA") {
            cmd__expected_ncells = "--expected_ncells ${ncells_expected}"
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
            --lower_bound_cell_estimate ${lower_bound_cell_estimate} \
            --lower_bound_total_droplets_included \
                ${lower_bound_total_droplets_included} \
            ${cmd__expected_ncells}
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

    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix),
            val(ncells_expected)
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
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    // TODO: convert out h5 file to barcodes.tsv.gz, features.tsv.gz and
    //           matrix.mtx.gz
    //       USER may need to play with the traing fraction and learning rate
    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}" // /${experiment_id}"
        lr_string = "${learning_rate}".replaceAll("\\.", "pt")
        outfile = "cellbender-epochs_${epochs}"
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
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
