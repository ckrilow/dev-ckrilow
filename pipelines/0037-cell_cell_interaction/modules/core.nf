#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process cellphonedb__make_database {
    // Generates cellphonedb database
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        val(version)

    output:
        path(".cpdb/releases/*/cellphone.db", emit: database)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Set default version to latest databse
        if ("${version}" == "") {
            version = "latest"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellphonedb__make_database: ${process_info}"
        cellphonedb database download --version ${version}
        """
}


process cellphonedb__run_analysis {
    // Runs cellphonedb analysis.
    // 1. Takes adata file and converts it into cellphone db input
    // 2. Runs cellphonedb
    // 3. Plots cellphonedb results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(h5_anndata)
        path(database)
        each iterations
        each threshold

    output:
        val(outdir, emit: outdir)
        path(
            "pvalues.tsv",
            emit: tsv_pvalues
        )
        path(
            "means.tsv",
            emit: tsv_means
        )
        path(
            "significant_means.tsv",
            emit: tsv_significant_means
        )

    // NOTE: DLT - I don't know if this should be statistical_analysis or just
    // analysis
    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        outfile = "cellphonedb"
        outfile = "${outfile}-iterations_${iterations}"
        thresh_string = "${threshold}".replaceAll("\\.", "pt")
        outfile = "${outfile}-threshold_${thresh_string}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellphonedb__run_analysis: ${process_info}"
        0015-h5ad_to_cellphonedb.py \
            --h5_anndata ${h5_anndata} \
            --output_file cellphonedb
        cellphonedb method statistical_analysis \
            test_meta.txt \
            test_counts.txt \
            --counts-data ensembl \
            --project-name ${outfile} \
            --database ${database} \
            --iterations ${iterations} \
            --threshold ${threshold} \
            --result-precision 3 \
            --output-path \$(pwd) \
            --output-format tsv \
            --pvalue 0.05 \
            --debug-seed 0 \
            --threads ${task.cpus} \
            --verbose
        """
}


process cellphonedb__plot_dot {
    // Dotplot of cellphonedb results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(pvalues)
        path(means)

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellphonedb__plot_dot: ${process_info}"
        rm -fr plots
        cellphonedb plot dot_plot \
            --pvalues-path ${pvalues} \
            --means-path ${means} \
            --output-path \$(pwd) \
            --output-name dotplot.pdf \
            --verbose
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process cellphonedb__plot_heatmap {
    // Heatmap of cellphonedb results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(pvalues)
        val(pvalue_cutoff)

    output:
        val(outdir, emit: outdir)
        path(
            "count_network-pvcutoff__${pvalue_cut_string}.tsv",
            emit: count_network_tsv
        )
        path(
            "interactions_count-pvcutoff__${pvalue_cut_string}.tsv",
            emit: interactions_count_tsv
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        pvalue_cut_string = "${pvalue_cutoff}".replaceAll("\\.", "pt")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cellphonedb__plot_dot: ${process_info}"
        rm -fr plots
        cellphonedb plot heatmap_plot \
            --pvalues-path ${pvalues} \
            --output-path \$(pwd) \
            --count-name heatmap_count-pvcutoff__${pvalue_cut_string}.pdf \
            --log-name heatmap_log_count-pvcutoff__${pvalue_cut_string}.pdf \
            --count-network-name count_network-pvcutoff__${pvalue_cut_string}.tsv \
            --interaction-count-name interactions_count-pvcutoff__${pvalue_cut_string}.tsv \
            --pvalue ${pvalue_cutoff}
            --verbose
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
