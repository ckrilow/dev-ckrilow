#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process run_scrublet {
    // Runs scrublet for each sample.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "copy",
                overwrite: "true"

    //each smplid_and_datapath
    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )

    output:
        val(outdir, emit: outdir)
        val(experiment_id, emit: experiment_id)
        path("${runid}-${outfile}-scrublet.tsv.gz", emit: multiplet_calls)
        val("${outdir}/${outfile}-scrublet.tsv.gz",
            emit: multiplet_calls_published
        )
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/multiplet"
        outdir = "${outdir}.method=scrublet"
        outfile = "${experiment_id}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_scrublet: ${process_info}"
        TMP_DIR=\$(mktemp -d -p \$(pwd))
        ln --physical ${file_10x_barcodes} \$TMP_DIR
        ln --physical ${file_10x_features} \$TMP_DIR
        ln --physical ${file_10x_matrix} \$TMP_DIR
        run_scrublet.py \
            --tenxdata_dir \$TMP_DIR \
            --n_simulated_multiplet 10000 \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process make_cellmetadata_pipeline_input {
    // Makes a input tsv file of all of the
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(scrublet__experiment_id)
        path(scrublet__multiplet_calls)

    output:
        val(outdir, emit: outdir)
        path('file_cellmetadata.tsv', emit: file__cellmetadata)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "make_pipeline_input_file: ${process_info}"
        # Note: the default paste delim is tab
        paste ${scrublet__experiment_id} ${scrublet__multiplet_calls} \
            | awk 'BEGIN{print "experiment_id\tper_cell_metadata"}1' \
            > file_cellmetadata.tsv
        """
}


// workflow wf__multiplet {
//     take:
//         outdir
//         sample_id
//         paths_10x
//     main:
//         // Identify multiplets using scrublet.
//         run_scrublet(
//             outdir
//             sample_id
//             paths_10x
//         )
//         // Generate input file for merge based in multiplets
//         // make_pipeline_input_file(
//         //     run_scrublet.out.outdir,
//         //     run_scrublet.out.multiplet_annotations
//         // )
//     // Return merged input data file.
// }
