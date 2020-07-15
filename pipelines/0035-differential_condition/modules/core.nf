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
        015-run_diffxpy.py
            --input txd_input \
            --condition_column disease_status \
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
