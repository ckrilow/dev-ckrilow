#!/usr/bin/env nextflow

// NOTE: label 'big_mem' may be useful at some point


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process merge_samples {
    // Takes a list of raw 10x files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "copy",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file_paths_10x)
        path(file_metadata)
        path(file_params)
        path(file_cellmetadata)
        val(metadata_key)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("${runid}-adata.h5ad", emit: anndata)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        cmd__params = ""
        if (file_params.name != "no_file__file_sample_qc") {
            cmd__params = "--params_yaml ${file_params}"
        }
        cmd__cellmetadata = ""
        if (file_cellmetadata.name != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_samples: ${process_info}"
        0025-scanpy_merge.py \
            --tenxdata_file ${file_paths_10x} \
            --sample_metadata_file ${file_metadata} \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-adata \
            ${cmd__params} \
            ${cmd__cellmetadata}
        """
}


process plot_qc {
    // Takes annData object, generates basic qc plots
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
        each facet_columns
        each variable_columns_distribution_plots

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        // Append run_id to output file.
        outfile = "${runid}-${outfile}"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_qc: ${process_info}"
        plot_qc_umi_nfeature_mt.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            ${cmd__facet_columns}
        plot_anndata_distribution_across_cells.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            --variable_columns ${variable_columns_distribution_plots} \
            ${cmd__facet_columns}
        plot_anndata_distribution_across_cells.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            --variable_columns ${variable_columns_distribution_plots} \
            --ecdf \
            ${cmd__facet_columns}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process normalize_and_pca {
    // Takes annData object, nomalizes across samples, calculates PCs.
    // NOTE: Once normalization is set, it would be faster to normalize per
    //       sample and then merge.
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
        path(file__genes_exclude_hvg)
        path(file__genes_score)
        each vars_to_regress

    output:
        val(outdir, emit: outdir)
        path("${runid}-adata-normalized_pca.h5ad", emit: anndata)
        path("${runid}-adata-metadata.tsv.gz", emit: metadata)
        path("${runid}-adata-pcs.tsv.gz", emit: pcs)
        path(
            "${runid}-adata-normalized_pca-counts.h5ad",
            emit: anndata_filtered_counts
        )
        path("plots/*.pdf")
        path("plots/*.png") optional true
        // tuple(
        //     val(outdir),
        //     path("${runid}-adata-normalized_pca.h5ad"),
        //     path("${runid}-adata-metadata.tsv.gz"),
        //     path("${runid}-adata-pcs.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/normalize=total_count"
        // Add any variables we are regressing to the output dir.
        if (vars_to_regress == "") {
            outdir = "${outdir}.vars_to_regress=none"
            cmd__vars_to_regress = ""
        } else {
            outdir = "${outdir}.vars_to_regress=${vars_to_regress}"
            cmd__vars_to_regress = "--vars_to_regress ${vars_to_regress}"
        }
        // Add details on the genes we are exlcuding from hgv list.
        file_vge = "${file__genes_exclude_hvg.getSimpleName()}"
        outdir = "${outdir}.hvg_exclude=${file_vge}"
        // Add details on the scores we are using.
        file_score = "${file__genes_score.getSimpleName()}"
        outdir = "${outdir}.scores=${file_score}"
        // Customize command for optional files.
        cmd__genes_exclude_hvg = ""
        if (file__genes_exclude_hvg.name != "no_file__genes_exclude_hvg") {
            cmd__genes_exclude_hvg = "--variable_genes_exclude ${file__genes_exclude_hvg}"
        }
        cmd__genes_score = ""
        if (file__genes_score.name != "no_file__genes_score") {
            cmd__genes_score = "--score_genes ${file__genes_score}"
        }
        // Basic details on the run.
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "normalize_pca: ${process_info}"
        0035-scanpy_normalize_pca.py \
            --h5_anndata ${file__anndata} \
            --output_file ${runid}-adata \
            --number_cpu ${task.cpus} \
            ${cmd__vars_to_regress} \
            ${cmd__genes_exclude_hvg} \
            ${cmd__genes_score}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
        // Old version with bash evaluation of optional commands
        //
        // echo "normalize_pca: ${process_info}"
        // # If there are entries in the variable_genes_exclude file, add it to
        // # the call.
        // cmd__vg_exclude="--variable_genes_exclude ${file__genes_exclude_hvg}"
        // val=\$(cat ${file__genes_exclude_hvg} | wc -l)
        // if [ \$val -eq 0 ]; then cmd__vg_exclude=""; fi
        // # If there are entries in the score_genes file, add it to the call.
        // cmd__score_genes="--score_genes ${file__genes_score}"
        // val=\$(cat ${file__genes_score} | wc -l)
        // if [ \$val -eq 0 ]; then cmd__score_genes=""; fi
        // 0035-scanpy_normalize_pca.py \
        //     --h5_anndata ${file__anndata} \
        //     --output_file ${runid}-adata \
        //     --number_cpu ${task.cpus} \
        //     ${cmd__vars_to_regress} \
        //     \${cmd__vg_exclude} \
        //     \${cmd__score_genes}
        // mkdir plots
        // mv *pdf plots/ 2>/dev/null || true
        // mv *png plots/ 2>/dev/null || true
}


process subset_pcs {
    // Takes PCs (rows = cell barcodes) and subsets down to a specified number.
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
        each n_pcs

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("${runid}-reduced_dims.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/reduced_dims-pca"
        outdir = "${outdir}.n_pcs=${n_pcs}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "subset_pcs: ${process_info}"
        0045-subset_pca_file.py \
            --tsv_pcs ${file__pcs} \
            --number_pcs ${n_pcs} \
            --output_file ${runid}-reduced_dims
        """
}


process harmony {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo true          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
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
        each n_pcs
        each variables_and_thetas

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("${runid}-reduced_dims.tsv.gz"),
        //     val(n_pcs),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/reduced_dims-harmony"
        outdir = "${outdir}.n_pcs=${n_pcs}"
        outdir = "${outdir}.variables=${variables_and_thetas.variable}"
        theta_str = "${variables_and_thetas.theta}" // .replaceAll("\\.", "pt")
        outdir = "${outdir}.thetas=${theta_str}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "harmony: ${process_info}"
        0045-harmony_process_pcs.py \
            --pca_file ${file__pcs} \
            --metadata_file ${file__metadata} \
            --metadata_columns ${variables_and_thetas.variable} \
            --theta ${variables_and_thetas.theta} \
            --n_pcs ${n_pcs} \
            --out_file ${runid}-reduced_dims
        """
        // 0045-harmony_process_pcs.R \
        //     --pca_file ${file__pcs} \
        //     --metadata_file ${file__metadata} \
        //     --metadata_columns ${variables_and_thetas.variable} \
        //     --theta ${variables_and_thetas.theta} \
        //     --n_pcs ${n_pcs} \
        //     --out_file ${runid}-reduced_dims
}
