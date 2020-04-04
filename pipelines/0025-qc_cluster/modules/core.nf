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

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("${runid}-adata.h5", emit: anndata)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        """
        echo "merge_samples: ${process_info}"
        0025-scanpy_merge.py \
            --tenxdata_file ${file_paths_10x} \
            --metadata_file ${file_metadata} \
            --metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key "sanger_sample_id" \
            --params_yaml ${file_params} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-adata
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
        path(file__variable_genes_exclude)
        each vars_to_regress

    output:
        val(outdir, emit: outdir)
        path("${runid}-adata-normalized_pca.h5", emit: anndata)
        path("${runid}-adata-metadata.tsv.gz", emit: metadata)
        path("${runid}-adata-pcs.tsv.gz", emit: pcs)
        path("*.pdf") optional true
        path("*.png") optional true
        // tuple(
        //     val(outdir),
        //     path("${runid}-adata-normalized_pca.h5"),
        //     path("${runid}-adata-metadata.tsv.gz"),
        //     path("${runid}-adata-pcs.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/normalize.total_count"
        if (vars_to_regress == "") {
            outdir = "${outdir}-scale.vars_to_regress=none"
            cmd__vars_to_regress = ""
        } else {
            outdir = "${outdir}-scale.vars_to_regress=${vars_to_regress}"
            cmd__vars_to_regress = "--vars_to_regress ${vars_to_regress}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "normalize_pca: ${process_info}"
        cmd__vg_exclude="--variable_genes_exclude ${file__variable_genes_exclude}"
        val=\$(cat ${file__variable_genes_exclude} | wc -l)
        if [ \$val -eq 0 ]; then cmd__vg_exclude=""; fi
        0035-scanpy_normalize_pca.py \
            --h5_anndata ${file__anndata} \
            --output_file ${runid}-adata \
            --number_cpu ${task.cpus} \
            ${cmd__vars_to_regress} \
            \${cmd__vg_exclude}
        """
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
                    if (filename.endsWith("normalized_pca.h5")) {
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
        outdir = "${outdir}-n_pcs=${n_pcs}"
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
                    if (filename.endsWith("normalized_pca.h5")) {
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
        outdir = "${outdir}-n_pcs=${n_pcs}"
        outdir = "${outdir}-variables=${variables_and_thetas.variable}"
        outdir = "${outdir}-thetas=${variables_and_thetas.theta}"
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
