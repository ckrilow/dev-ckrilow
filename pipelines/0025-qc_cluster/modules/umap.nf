#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process umap_calculate {
    // UMAP from reduced_dims.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    // Don't publish these results as they are just temporary
    // publishDir  path: "${outdir}",
    //             saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    //             mode: "copy",
    //             overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        //path(file__metadata) // not really needed
        //path(file__pcs) // not really needed
        path(file__reduced_dims)
        each n_neighbors
        each umap_init
        each umap_min_dist
        each umap_spread

    output:
        val(outdir, emit: outdir)
        tuple(
            val("${in_file_id}"),
            file("${runid}-${outfile}.h5ad"),
            emit: anndata
        )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Get a unique identifier for the input file that will be later
        // used as a key to merge umap jobs on the same input.
        //
        // Here the unique identifier is file__reduced_dims, since the reduced
        // dims function does not change anndata but just writes a tsv file.
        //
        // If file__reduced_dims == file object, then use
        // file__reduced_dims.name...
        in_file_id = file__reduced_dims.toString().tokenize('-').get(0)
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "${outfile}-umap"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "umap_calculate: ${process_info}"
        umap_calculate.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --n_neighbors ${n_neighbors} \
            --umap_init ${umap_init} \
            --umap_min_dist ${umap_min_dist} \
            --umap_spread ${umap_spread} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-${outfile}
        """
        //--calculate_densities \
}



process umap_gather {
    // Merge UMAP from reduced_dims (reduce or gather).
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
        tuple(val(key), path(files__anndata))

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}.h5ad", emit: anndata)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "adata-umap"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        files__anndata = files__anndata.join(',')
        """
        echo "umap_gather: ${process_info}"
        umap_gather.py \
            --h5_anndata_list ${files__anndata} \
            --output_file ${runid}-${outfile}
        """
}


process umap_plot_swarm {
    // Plot UMAPs.
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
        path(file__anndata)
        val(colors_quantitative)
        val(colors_categorical)
        val(drop_legend_n)

    output:
        path("plots/*.png")
        path("plots/*.pdf") optional true
        // path("plots/*.svg") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "umap"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative ${colors_quantitative}"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical ${colors_categorical}"
        }
        // drop_legend_n = "-1"
        // if (cmd__colors_cat.contains("experiment_id")) {
        //     drop_legend_n = "8"
        // }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "umap: ${process_info}"
        umap_plot.py \
            --h5_anndata ${file__anndata} \
            --number_cpu ${task.cpus} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --drop_legend_n ${drop_legend_n} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process umap_calculate_and_plot {
    // UMAP from reduced_dims.
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
        path(file__anndata)
        path(file__reduced_dims)
        val(colors_quantitative)
        val(colors_categorical)
        each n_neighbors
        each umap_init
        each umap_min_dist
        each umap_spread

    output:
        path("plots/*.png")
        path("plots/*.pdf") optional true
        // path("plots/*.svg") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "umap"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative ${colors_quantitative}"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical ${colors_categorical}"
        }
        drop_legend_n = "-1"
        if (cmd__colors_cat.contains("experiment_id")) {
            drop_legend_n = "8"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "umap_calculate_and_plot: ${process_info}"
        umap_calculate_and_plot.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --n_neighbors ${n_neighbors} \
            --umap_init ${umap_init} \
            --umap_min_dist ${umap_min_dist} \
            --umap_spread ${umap_spread} \
            --number_cpu ${task.cpus} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --drop_legend_n ${drop_legend_n} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
        // --calculate_densities
}


workflow wf__umap {
    take:
        outdir
        anndata
        // metadata
        // pcs
        reduced_dims
        n_neighbors
        umap_init
        umap_min_dist
        umap_spread
        colors_quantitative
        colors_categorical
    main:
        // Calculate UMAPs with various settings
        umap_calculate(
            outdir,
            anndata,
            // metadata,
            // pcs,
            reduced_dims,
            n_neighbors,
            umap_init,
            umap_min_dist,
            umap_spread
        )
        // Gather step.
        // Gather by tuple ... if we just to a collect, then will get all
        // umap_calculate calls, not split by reduced_dims. See link below:
        // http://nextflow-io.github.io/patterns/index.html#_process_outputs_into_groups
        umap_gather(
            outdir,
            umap_calculate.out.anndata.groupTuple()
        )
        // Make plots
        umap_plot_swarm(
            umap_gather.out.outdir,
            umap_gather.out.anndata,
            // cluster.out.reduced_dims,
            colors_quantitative,
            colors_categorical,
            '8'
        )
    emit:
        // Return merged input data file.
        outdir = umap_gather.out.outdir
        file__cellmetadata = umap_gather.out.anndata
}
