#!/usr/bin/env nextflow

// NOTE: label 'big_mem' may be useful at some point


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process merge_samples {
    // Takes a list of raw 10x files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode         // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
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
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

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
        rm -fr plots
        0025-scanpy_merge.py \
            --tenxdata_file ${file_paths_10x} \
            --sample_metadata_file ${file_metadata} \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-adata \
            ${cmd__params} \
            ${cmd__cellmetadata}
        0026-plot_filtered_cells.py \
            --tsv_file ${runid}-adata-cell_filtered_per_experiment.tsv.gz \
            --output_file ${runid}-adata-cell_filtered_per_experiment
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_predicted_sex {
    // Takes annData object, plots the predicted sex fron gene expression
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Append run_id to output file.
        outfile = "${runid}-scatterplot-sex_sample_swap_check"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_predicted_sex: ${process_info}"
        rm -fr plots
        0028-plot_predicted_sex.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_qc {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
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
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Append run_id to output file.
        outfile = "${runid}-${outfile}"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        // Run distribution across cells if a value is specified
        cmd__anndataobs = ""
        cmd__anndataobs_ecdf = ""
        if (variable_columns_distribution_plots != "") {
            cmd__anndataobs = "plot_anndataobs_distribution_across_cells.py"
            cmd__anndataobs = "${cmd__anndataobs} --h5_anndata ${file__anndata}"
            cmd__anndataobs = "${cmd__anndataobs} --output_file ${outfile}"
            cmd__anndataobs = "${cmd__anndataobs} --variable_columns ${variable_columns_distribution_plots}"
            cmd__anndataobs = "${cmd__anndataobs} ${cmd__facet_columns}"
            cmd__anndataobs_ecdf = "${cmd__anndataobs} --ecdf"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_qc: ${process_info}"
        rm -fr plots
        plot_qc_umi_nfeature_mt.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            ${cmd__facet_columns}
        ${cmd__anndataobs}
        ${cmd__anndataobs_ecdf}
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
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
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
        val("${param_details}", emit: param_details)
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
        // Add any variables we are regressing to the output dir.
        param_details="vars_to_regress=none"
        if (vars_to_regress == "") {
            cmd__vars_to_regress = ""
        } else {
            param_details = "vars_to_regress=${vars_to_regress}"
            cmd__vars_to_regress = "--vars_to_regress ${vars_to_regress}"
        }
        outdir = "${outdir_prev}/normalize=total_count.${param_details}"
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
        rm -fr plots
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


process estimate_pca_elbow {
    // Takes annData object, estiamtes the elbow in PC var explained.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}.tsv.gz", emit: pca_elbow_estimate)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        outfile = "${outfile}-knee"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "estimate_pca_elbow: ${process_info}"
        rm -fr plots
        0030-estimate_pca_elbow.py \
            --h5_anndata ${file__anndata} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process subset_pcs {
    // Takes PCs (rows = cell barcodes) and subsets down to a specified number.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("${runid}-reduced_dims.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        param_details = "${params__pcs}-pca"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "subset_pcs: ${process_info}"
        0045-subset_pca_file.py \
            --tsv_pcs ${file__pcs} \
            --number_pcs ${n_pcs} \
            --output_file ${runid}-reduced_dims
        cp ${runid}-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        """
}


process harmony {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs
        each variables_and_thetas

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("${runid}-reduced_dims.tsv.gz"),
        //     val(n_pcs),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        param_details = "${params__pcs}-harmony"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        param_details = "${param_details}.variables=${variables_and_thetas.variable}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        outdir = "${outdir}.thetas=${variables_and_thetas.theta}"
        theta_str = "${variables_and_thetas.theta}".replaceAll("\\.", "pt")
        param_details = "${param_details}.thetas=${theta_str}"
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
        cp ${runid}-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        """
        // NOTE: below code for harmony in R
        // 0045-harmony_process_pcs.R \
        //     --pca_file ${file__pcs} \
        //     --metadata_file ${file__metadata} \
        //     --metadata_columns ${variables_and_thetas.variable} \
        //     --theta ${variables_and_thetas.theta} \
        //     --n_pcs ${n_pcs} \
        //     --out_file ${runid}-reduced_dims
}


process bbknn {
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs
        each batch_var

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}-bbknn.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )

    script:
        runid = random_hex(16)
        param_details = "${params__pcs}-bbknn"
        param_details = "${param_details}.batch=${batch_var}"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "bbknn: ${process_info}"
        0045-bbknn.py \
            --h5_anndata ${file__anndata} \
            --batch_key ${batch_var} \
            --n_pcs ${n_pcs} \
            --output_file ${runid}-${outfile}-bbknn
        cp ${runid}-${outfile}-bbknn-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        mv ${runid}-${outfile}-bbknn-reduced_dims.tsv.gz \
            ${runid}-reduced_dims.tsv.gz
        """
}


process lisi {
    // Takes a list of reduced_dims and calculates lisi
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

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
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__metadata)
        val(variables)
        file(file__reduced_dims)
        //tuple(val(label__reduced_dims), file(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path(file__metadata, emit: metadata)
        path("${runid}-${outfile}-lisi.tsv.gz", emit: clusters)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__metadata job.
        outfile = "${file__metadata}".minus(".tsv.gz")
            .split("-").drop(1).join("-")
        file__reduced_dims = file__reduced_dims.join("::")
        label__reduced_dims = file__reduced_dims
            .replaceAll("reduced_dims-", "")
            .replaceAll(".tsv.gz", "")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "lisi: ${process_info}"
        rm -fr plots
        0047-lisi.py \
            --reduced_dims_tsv ${file__reduced_dims} \
            --reduced_dims_tsv_labels ${label__reduced_dims} \
            --metadata_tsv ${file__metadata} \
            --metadata_columns ${variables} \
            --perplexity 30 \
            --output_file ${runid}-${outfile}-lisi
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
