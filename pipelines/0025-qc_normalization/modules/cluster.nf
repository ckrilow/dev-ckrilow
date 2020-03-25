
process cluster {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir  path: "${outdir}",
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__reduced_dims)

    output:
        val(outdir, emit: outdir)
        path("test-clustered.h5", emit: anndata)
        path("test-clustered.tsv.gz", emit: clusters)
        path("*.pdf") optional true
        path("*.png") optional true

    script:
        outdir = "${outdir_prev}/cluster"
        """
        scanpy_cluster.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --cluster_method leiden \
            --number_pcs 15 \
            --resolution 1 \
            --output_file test-clustered
        """
}


process umap {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir  path: "${outdir}",
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__reduced_dims)

    output:
        path("*.png")

    script:
        outdir = "${outdir_prev}"
        """
        scanpy_umap.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --colors_quantitative age \
            --colors_categorical sanger_sample_id,sex,leiden \
            --number_pcs 15
        """
}


workflow wf_cluster {
    take:
        outdir
        anndata
        reduced_dims
    main:
        // Cluster the results, varying the resolution.
        cluster(
            outdir,
            anndata,
            reduced_dims
        )
        // Generate UMAPs of the results.
        umap(
            cluster.out.outdir,
            cluster.out.anndata,
            reduced_dims
        )
}
