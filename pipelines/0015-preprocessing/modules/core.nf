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
            path(file_10x_matrix)
        )

    output:
        path("cellbender.h5", emit: h5)
        path("cellbender_filtered.h5", emit: h5_filtered)
        path("cellbender_cell_barcodes.csv", emit: barcodes)
        path("cellbender.log", emit: log)
        path("plots/*.png") optional true

    // TODO: (a) automatically set --expected-cells and
    //       --total-droplets-included based on UMI curve
    //       (b) convert out h5 file to barcodes.tsv.gz, features.tsv.gz and
    //           matrix.mtx.gz
    //       USER may need to play with the traing fraction and learning rate
    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
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
        cellbender remove-background \
            --input input
            --output cellbender \
            --cuda \
            --expected-cells 6036 \
            --total-droplets-included 250000 \
            --z-dim 200 \
            --z-layers 1000 \
            --low-count-threshold 5 \
            --epochs 200 \
            --empty-drop-training-fraction 0.25 \
            --learning-rate 0.0001
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
