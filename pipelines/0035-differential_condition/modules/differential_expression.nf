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
    //label 'gpu'          // use GPU
    label 'long_job'
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(cell_label_column)
        each cell_label_analyse
        each model
        // each condition_column
        // each covariate_columns
        each method

    output:
        val(outdir, emit: outdir)
        tuple(
            val(runid), // Need random hex to control grouping
            val(condition_column),
            val(cell_label_analyse),
            val(covariate_columns),
            val(method),
            path("${outfile}-de_results.tsv.gz"),
            val(outdir),
            emit: results
        )
        path("${outfile}-de_results_obj.joblib.gz")
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        condition_column = model.variable
        // Sort out covariates
        covariate_columns_discrete = model.covariate_discrete
        covariate_columns_continuous = model.covariate_continuous
        covariate_columns = ""  // list of all covariates used in model
        cmd__covar = ""
        if (covariate_columns_discrete != "") {  // add disc cov call
            cmd__covar = "${cmd__covar} --covariate_columns_discrete ${covariate_columns_discrete}"
            covariate_columns = "${covariate_columns_discrete},"
        }
        if (covariate_columns_continuous != "") {  // add contin cov call
            cmd__covar = "${cmd__covar} --covariate_columns_continuous ${covariate_columns_continuous}"
            covariate_columns = "${covariate_columns}${covariate_columns_continuous}"
        }
        if (covariate_columns.endsWith(",")) {
            covariate_columns = covariate_columns.substring(
                0,
                covariate_columns.length() - 1
            )
        }
        if (covariate_columns == "") {
            covariate_columns = "none"
        }
        // cell_label_analyse comes in array-format.
        cell_label_analyse = cell_label_analyse[0] // Get first element.
        outdir = "${outdir_prev}/${condition_column}/diffxpy/"
        outdir = "${outdir}cell_label=${cell_label_analyse}"
        outdir = "${outdir}_covariates=${covariate_columns}"
        outdir = "${outdir}_method=${method}"
        outfile = "cell_label__${cell_label_analyse}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_diffxpy: ${process_info}"
        rm -fr plots
        015-run_diffxpy.py \
            --h5_anndata ${anndata} \
            --condition_column ${condition_column} \
            --cell_label_column ${cell_label_column} \
            --cell_label_analyse ${cell_label_analyse} \
            --method ${method} \
            --output_file ${outfile} \
            ${cmd__covar}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process run_mast {
    // Run MAST
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(cell_label_column)
        each cell_label_analyse
        each model
        // each condition_column
        // each covariate_columns
        each method

    output:
        val(outdir, emit: outdir)
        tuple(
            val(runid), //need random hex to control grouping
            val(condition_column),
            val(cell_label_analyse),
            val(covariate_columns),
            val(method),
            path("${outfile}-de_results.tsv.gz"),
            val(outdir),
            emit: results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        condition_column = model.variable
        // Sort out covariates
        covariate_columns_discrete = model.covariate_discrete
        covariate_columns_continuous = model.covariate_continuous
        covariate_columns = ""  // list of all covariates used in model
        cmd__covar = ""
        if (covariate_columns_discrete != "") {  // add disc cov call
            cmd__covar = " --covariate_columns_discrete ${covariate_columns_discrete}"
            covariate_columns = "${covariate_columns_discrete},"
        }
        if (covariate_columns_continuous != "") {  // add contin cov call
            cmd__covar = "${cmd__covar} --covariate_columns_continuous ${covariate_columns_continuous}"
            covariate_columns = "${covariate_columns}${covariate_columns_continuous}"
        }
        if (covariate_columns.endsWith(",")) {
            covariate_columns = covariate_columns.substring(
                0,
                covariate_columns.length() - 1
            )
        }
        if (covariate_columns == "") {
            covariate_columns = "none"
        }
        // cell_label_analyse comes in array-format.
        cell_label_analyse = cell_label_analyse[0] // Get first element.
        outdir = "${outdir_prev}/${condition_column}/mast/"
        outdir = "${outdir}cell_label=${cell_label_analyse}"
        outdir = "${outdir}_covariates=${covariate_columns}"
        outdir = "${outdir}_method=${method}"
        outfile = "cell_label__${cell_label_analyse}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_mast: ${process_info}"
        rm -fr plots
        016-prepare_mast_input.py \
            --h5ad_file ${anndata} \
            --condition_column ${condition_column} \
            --cell_label_column ${cell_label_column} \
            --cell_label_analyse ${cell_label_analyse} \
            --output_dir mast_input \
            ${cmd__covar}
        016-run_mast.R \
            --mtx_dir mast_input \
            --cell_label_column ${cell_label_column} \
            --cell_label_analysed ${cell_label_analyse} \
            --condition_column ${condition_column} \
            --method ${method} \
            --out_file ${outfile} \
            --cores_available ${task.cpus} \
            ${cmd__covar}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process merge_dataframes {
    // Merge resulting dataframes from diffxpy
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            val(result_keys),
            val(result_paths)
        )

    output:
        val(outdir, emit: outdir)
        tuple(
            val(condition),
            path("${outfile}-de_results.tsv.gz"),
            emit: merged_results
        )

    script:
        runid = random_hex(16)
        outdir = outdir_prev
        outfile = "${condition}_merged_results"
        result_keys = result_keys.join(",")
        result_paths = result_paths.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_dataframes: ${process_info}"
        019-merge_de_dataframes.py \
            --dataframe_keys ${result_keys} \
            --dataframe_paths ${result_paths} \
            --output_file ${outfile}
        """
}


process plot_de_results {
    // Generate plots from the merged data frames to evaluate
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )
        each mean_expression_filter

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = outdir_prev
        outfile = "${condition}-merged"
        // script automatically adds expression filter
        // outfile = "${condition}-mean_expr_filt__${mean_expression_filter}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_de_results: ${process_info}"
        rm -fr plots
        019-compare_de_results.py \
            --dataframe ${merged_df} \
            --columns_to_compare de_method,covariates_passed \
            --mean_expression_filter 0.1 \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process run_rna_enrich {
    // Run RNA Enrich for each DE result
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label_analyse),
            val(covariate_columns),
            val(method),
            path(de_results),
            val(outdir_prev)
        )
        each rna_enrich_config

    output:
        val(outdir, emit: outdir)
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label_analyse),
            val(covariate_columns),
            val(method),
            path("${outfile}-rna_enrich_results.tsv"),
            val(outdir),
            emit: results
        )
        path("${outfile}-ensembl_entrez_mapping.tsv.gz", emit: mappings)
        path("plots/*.jpg") optional true

    script:
        runid = random_hex(16)
        // Inputs besides `runid` come in array format
        cell_label_analyse = cell_label_analyse[0]
        outdir = outdir_prev[0]
        outdir = "${outdir}/rna_enrich-"
        outdir = "${outdir}pval=${rna_enrich_config.pval_cutoff}"
        outdir = "${outdir}_database=${rna_enrich_config.database}"
        outdir = "${outdir}_min_counts=${rna_enrich_config.min_counts}/"
        outfile = "cell_label__${cell_label_analyse}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_de_results: ${process_info}"
        rm -fr plots
        031-run_rna_enrich.R \
            --de_file ${de_results} \
            --rna_enrich_script $baseDir/bin/legacy/rna_enrich.R \
            --min_gene_to_test 10 \
            --max_gene_to_test 99999 \
            --p_val_cutoff ${rna_enrich_config.pval_cutoff} \
            --database ${rna_enrich_config.database} \
            --min_counts ${rna_enrich_config.min_counts} \
            --output_file ${outfile} \
            --verbose
        mkdir plots
        mv *jpg plots/ 2>/dev/null || true
        """
}

workflow wf__differential_expression {
    take:
        outdir
        anndata
        anndata_cell_label
        model
        diffxpy_method_config
        mast_method_config
        plot_config
        rna_enrich_config
    main:
        // Get a list of all of the cell types
        get_cell_label_list(
            anndata,
            anndata_cell_label
        )
        // For each cell type compute differential expression for that cell
        // type
        cell_labels = get_cell_label_list.out.cell_labels
            .splitCsv(header: false, sep: ',')
            // .map{row -> tuple(
            //     row[0]
            // )}
        // Run diffxpy with all combinations of conditions (e.g., sex,
        // disease status), covariates (e.g., size_factors, age),
        // methods (e.g., wald)
        if (diffxpy_method_config.run_process) {
            run_diffxpy(
                outdir,
                anndata,
                anndata_cell_label,
                // '1',  // just run on first cluster for development
                cell_labels,  // run for all clusters for run time
                model,
                diffxpy_method_config.value
            )
        }
        // Run MAST with all combinations of conditions (e.g., sex,
        // disease status), covariates (e.g., size_factors, age),
        // methods (e.g., wald)
        if (mast_method_config.run_process) {
            run_mast(
                outdir,
                anndata,
                anndata_cell_label,
                // '1',  // just run on first cluster for development
                cell_labels,  // run for all clusters for run time
                model,
                mast_method_config.value
            )
        }
        // Combine results of all of the models (e.g., different methods like
        // diffxpy and MAST, different covariate sets)
        if (diffxpy_method_config.run_process & mast_method_config.run_process) {
            de_results = run_diffxpy.out.results.groupTuple(by: 0)
                .concat(run_mast.out.results.groupTuple(by: 0))
        } else if (diffxpy_method_config.run_process) {
            de_results = run_diffxpy.out.results.groupTuple(by: 0)
        } else if (mast_method_config.run_process) {
            de_results = run_mast.out.results.groupTuple(by: 0)
        }
        de_results_merged = de_results
            .reduce([:]) { map, tuple ->
                def dataframe_key = "cell_label=" + tuple[2][0]
                dataframe_key += "::covariates=" + tuple[3][0].replaceAll(
                    ",",
                    "-"
                )
                dataframe_key += "::method=" + tuple[4][0]

                def map_key = tuple[1][0] // structure map by condition
                def key_list = map[map_key]
                if (!key_list) {
                    key_list = [[dataframe_key, tuple[5][0]]]
                } else {
                    key_list.add([dataframe_key, tuple[5][0]])
                }
                map[map_key] = key_list
                return(map)
            }
            .flatMap()
            .map {  entry ->
                combined_data = [entry.key, [], []]
                entry.value.each {
                    combined_data[1].add(it[0])
                    combined_data[2].add(it[1])
                }
                return(combined_data)
            }
        merge_dataframes(
            outdir,
            de_results_merged
        )
        // Basic plots of the differential expression results across all models
        plot_de_results(
            outdir,
            merge_dataframes.out.merged_results,
            plot_config.mean_expression_filter.value
        )

        // Run RNA_enrich on DE results
        run_rna_enrich(
            de_results,
            rna_enrich_config
        )
        // TODO: combine and compare all of the enrichment analysis results

    emit:
        cell_labels = get_cell_label_list.out.cell_labels
}
