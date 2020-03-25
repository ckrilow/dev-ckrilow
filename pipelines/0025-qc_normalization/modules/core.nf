
process merge_samples {
    // Takes a list of raw 10x files and merges them into one Seurat object.
    // NOTE: once normalization is set, it would be faster to normalize per
    //     sample and then merge
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 1             // cpu requirements
    //memory '50 GB'   // memory requirements

    publishDir  path: "${outdir}",
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file_paths_10x)
        path(file_metadata)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("adata.h5", emit: anndata)

    script:
        outdir = "${outdir_prev}"
        """
        scanpy_merge-dev.py \
            --samplesheetdata_file ${file_paths_10x} \
            --metadata_file ${file_metadata}
        """
}


process normalize_and_pca {
    // Takes annData object and nomalizes across samples.
    // NOTE: once normalization is set, it would be faster to normalize per
    //     sample and then merge
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script
    cpus 16            // cpu requirements
    //memory '50 GB'   // memory requirements

    // Publish the adata to dir the specifies the parameters used
    // Publish the PCs to a sub-directory
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("pcs.tsv.gz")) {
                       "reduced_dims-pca/${filename}"
                    } else {
                        "${filename}"
                    }
                },
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each vars_to_regress

    output:
        val(outdir, emit: outdir)
        path("adata-normalized_pca.h5", emit: anndata)
        path("adata-metadata.tsv.gz", emit: metadata)
        val(outdir_pca, emit: outdir_pca)
        path("adata-pcs.tsv.gz", emit: pcs)

    script:
        outdir = "${outdir_prev}/normalize.total_count"
        if (vars_to_regress == "") {
            outdir = "${outdir}-scale.vars_to_regress=none"
            outdir_pca = "${outdir}/reduced_dims-pca"
            """
            scanpy_normalize_pca.py \
                --h5_anndata ${file__anndata}
            """
        } else {
            outdir = "${outdir}-scale.vars_to_regress=${vars_to_regress}"
            outdir_pca = "${outdir}/reduced_dims-pca"
            """
            scanpy_normalize_pca.py \
                --h5_anndata ${file__anndata} \
                --vars_to_regress ${vars_to_regress}
            """
        }
}


process harmony {
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
        path(file__meta_data)
        path(file__reduced_dims)

    output:
        val(outdir, emit: outdir)
        path("adata-pcs-harmony.tsv.gz", emit: reduced_dims)

    script:
        outdir = "${outdir_prev}/reduced_dims-harmony"
        """
        harmony_process_pcs.R \
            --pca_file ${file__reduced_dims} \
            --metadata_file ${file__meta_data} \
            --metadata_columns sanger_sample_id \
            --n_pcs 15
        """
}
