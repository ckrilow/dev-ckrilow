#!/usr/bin/env nextflow

// NOTE: label 'big_mem' may be useful at some point


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
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
        val(epochs)
        val(learning_rate)

    output:
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
        outdir = "${outdir_prev}/${experiment_id}"
        lr_string = "${learning_rate}".replaceAll("\\.", "pt")
        outfile = "cellbender-epochs_${epochs}"
        outfile = "${outfile}-learningrate_${lr_string}"
        cmd__expected_ncells = ""
        if ("${ncells_expected}" != "NA") {
            cmd__expected_ncells = "--cmd__expected_ncells ${ncells_expected}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellbender__remove_background: ${process_info}"
        rm -fr plots
        mkdir input
        ln --physical ${file_10x_barcodes} input/barcodes.tsv.gz
        ln --physical ${file_10x_features} input/features.tsv.gz
        ln --physical ${file_10x_matrix} input/matrix.mtx.gz
        015-get_estimates_from_umi_counts.py \
            --tenxdata_path input \
            --output_file ${experiment_id} \
            --lower_bound_cell_estimate 100 \
            --lower_bound_total_droplets_included 10 \
            ${cmd__expected_ncells}
        cellbender remove-background \
            --input input
            --output ${outfile} \
            --cuda \
            --expected-cells $(cat ${experiment_id}-expected_cells.txt) \
            --total-droplets-included $(cat ${experiment_id}-total_droplets_included.txt) \
            --z-dim 200 \
            --z-layers 1000 \
            --low-count-threshold 10 \
            --epochs ${epochs} \
            --empty-drop-training-fraction 0.5 \
            --learning-rate ${learning_rate}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
