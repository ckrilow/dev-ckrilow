#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


// NOTE: This raises erroneous warning.
// WARN: There's no process matching config selector: umap
include {
    umap_calculate_and_plot as umap_calculate_and_plot__cluster;
} from "./umap.nf"


process cluster {
    // Clusters results.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        path(file__reduced_dims)
        each number_neighbors
        each method
        each resolution
        // tuple(val(outdir_prev), path(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}-clustered.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path("${runid}-${outfile}-clustered.tsv.gz", emit: clusters)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        resolution_str = "${resolution}" //.replaceAll("\\.", "pt")
        outdir = "${outdir_prev}/cluster"
        outdir = "${outdir}.number_neighbors=${number_neighbors}"
        outdir = "${outdir}.method=${method}"
        outdir = "${outdir}.resolution=${resolution_str}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cluster: ${process_info}"
        0055-scanpy_cluster.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --number_neighbors ${number_neighbors} \
            --cluster_method ${method} \
            --resolution ${resolution} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-${outfile}-clustered
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process cluster_validate_resolution {
    // Validate the resolution for clusters.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        path(file__reduced_dims)
        path(file__clusters)
        each sparsity
        each train_size_cells
        // each number_cells_downsample
        // each train_size_fraction

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path(file__clusters, emit: clusters)
        path("${runid}-${outfile}-lr_model.joblib.gz", emit: model)
        path(
            "${runid}-${outfile}-test_result.tsv.gz",
            emit: model_test_result
        )
        path(
            "${runid}-${outfile}-lr_coef.tsv.gz",
            emit: model_coefficient
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/validate_resolution"
        // outdir = "${outdir}.method=${method}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        // Add downsampling information.
        n_cells_downsample_txt = "none"
        // cmd__dask = "--dask_scale 500"
        // if (number_cells_downsample > 0) { // If downsample cells no dask
        //     n_cells_downsample_txt = "${number_cells_downsample}"
        //     cmd__dask = ""
        // }
        outfile = "${outfile}-n_cells_downsample=${n_cells_downsample_txt}"
        // Add sparsity information.
        sparsity_txt = "${sparsity}".replaceAll("\\.", "pt")
        outfile = "${outfile}-sparsity=${sparsity_txt}"
        // Add training cell count size.
        train_size_cells_txt = "none"
        cmd__dask = "--dask_scale 500"
        cmd__train_cells = ""
        if (train_size_cells > 0) {
            train_size_cells_txt = "${train_size_cells}"
            cmd__train_cells = "--train_size_cells ${train_size_cells}"
            cmd__dask = ""
        }
        outfile = "${outfile}-train_size_cells=${train_size_cells}"
        // train_size_txt = "${train_size_fraction}".replaceAll("\\.", "pt")
        // outfile = "${outfile}-train_size_fraction=${train_size_txt}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cluster_validate_resolution: ${process_info}"
        0057-scanpy_cluster_validate_resolution.py \
            --h5_anndata ${file__anndata} \
            --sparsity ${sparsity} \
            ${cmd__dask} \
            ${cmd__train_cells} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
        // --number_cells ${number_cells_downsample} \
        // --train_size_fraction ${train_size_fraction} \
}


process cluster_markers {
    // Find markers for clusters.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        path(file__reduced_dims)
        path(file__clusters)
        each method

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path(file__clusters, emit: clusters)
        path(
            "${runid}-${outfile}-cluster_markers.tsv.gz",
            emit: cluster_markers
        )
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/cluster_markers"
        outdir = "${outdir}.method=${method}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cluster: ${process_info}"
        0057-scanpy_cluster_markers.py \
            --h5_anndata ${file__anndata} \
            --rank_genes_method ${method} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process convert_seurat {
    // Converts anndata h5 file to a Seurat data object.
    // TODO: automatically add reduced_dims to Seurat data object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}.rds.gz", emit: seurat_data)
        path("${outdir_relative}/*", emit: matrix_data)
        // tuple(path("${outdir_relative}/*"), emit: matrix_data)

    script:
        runid = random_hex(16)
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outdir_relative = "${runid}-matrices-${outfile}"
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "convert_seurat: ${process_info}"
        outdir_relative_full_path=\$(pwd)"/${outdir_relative}"
        convert-anndata_10x.py \
            --h5_anndata ${file__anndata} \
            --output_dir ${outdir_relative}
        ln -s \${outdir_relative_full_path}/matrix-X.mtx.gz ${outdir_relative}/matrix.mtx.gz
        convert-10x_seurat.R \
            --in_dir \${outdir_relative_full_path} \
            --metadata_file \${outdir_relative_full_path}/metadata-barcodes.tsv.gz \
            --count_matrix_file \${outdir_relative_full_path}/matrix-counts.mtx.gz \
            --out_file ${runid}-${outfile}
        unlink ${outdir_relative}/matrix.mtx.gz
        """
}


workflow wf__cluster {
    take:
        outdir
        anndata
        metadata
        pcs
        reduced_dims
        cluster__number_neighbors
        cluster__methods
        cluster__resolutions
        cluster_validate_resolution__sparsity
        cluster_validate_resolution__train_size_cells
        // cluster_validate_resolution__number_cells
        // cluster_validate_resolution__test_size
        cluster_marker__methods
        n_neighbors
        umap_init
        umap_min_dist
        umap_spread
    main:
        // Cluster the results, varying the resolution.
        cluster(
            outdir,
            anndata,
            metadata,
            pcs,
            reduced_dims,
            cluster__number_neighbors,
            cluster__methods,
            cluster__resolutions
        )
        // Validate the resolution
        cluster_validate_resolution(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.metadata,
            cluster.out.pcs,
            cluster.out.reduced_dims,
            cluster.out.clusters,
            cluster_validate_resolution__sparsity,
            cluster_validate_resolution__train_size_cells
            // cluster_validate_resolution__number_cells,
            // cluster_validate_resolution__test_size
        )
        // Make Seurat dataframes of the clustered anndata
        // convert_seurat(
        //     cluster.out.outdir,
        //     cluster.out.anndata
        // )
        // Generate UMAPs of the results.
        umap_calculate_and_plot__cluster(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.reduced_dims,
            '',
            'cluster',
            n_neighbors,
            umap_init,
            umap_min_dist,
            umap_spread
        )
        // Find marker genes for clusters
        cluster_markers(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.metadata,
            cluster.out.pcs,
            cluster.out.reduced_dims,
            cluster.out.clusters,
            cluster_marker__methods
        )
}
