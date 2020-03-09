#!/usr/bin/env Rscript

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

SCRIPT_NAME <- "normalize_seurat_obj.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("future"))

# need to increase the RAM
# see https://satijalab.org/seurat/v3.0/future_vignette.html
future::plan("multiprocess", workers = 16)
options(future.globals.maxSize = 100000 * 1024^2)


#' Normalises each sample in merged Seurat data
#'
#' \code{normalize_seurat_obj_list} normalises each sample in Seurat dataframe
#' that has multiple samples.
#'
#' @param sc_df Seurat data object.
#'     Seurat data object with multiple experiments / samples concatenated.
#' @param method String.
#'     Normalisation method (LogNormalize | SCTransform).
#' @param filter_variable_genes_regex String.
#'     Regex that will remove variable genes. Example:
#'     "^HLA-|^IG[HJKL]|^RNA|^MT|^RP" removes HLA, immunoglobulin, RNA, MT,
#'     and RP genes based on HUGO gene names.
#' @param assay String.
#'     Assay slot to use in the Seurat data object.
#' @param sctransform_variable_regress List.
#'     Variables in meta.data to regress out of assay during SCTransform.
#'     NOTE: if merging across multiple samples and you want to explicitly
#'     model a batch variable that spans multiple samples, use harmony by
#'     itself or after seurat integration.
#'
#' @return List.
#'     List of multiple Seurat data objects that have been normalised.
#'
#' @importFrom Seurat SplitObject
#' @importFrom Seurat SCTransform
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ExpMean
#' @importFrom Seurat LogVMR
#' @importFrom stringr str_detect
#' @export
normalize_seurat_obj_list <- function(
        seurat_list,
        method = "LogNormalize",
        filter_variable_genes_regex = "",
        assay = "RNA",
        sctransform_variable_regress = c(),
        verbose = TRUE
    ) {
    #seurat_list <- Seurat::SplitObject(sc_df, split.by = grouping_var_id)
    for (i in 1:length(seurat_list)) {
        cat("[normalize_seurat_obj_list]\t", i, "/", length(seurat_list), "\n")

        if (method == "SCTransform") {
            assay_id <- "SCT"
            # performs normalization, variance stablization, and feature
            # selection
            #
            # See SCTransform tab from
            # https://satijalab.org/seurat/v3.1/integration.html
            #
            # https://github.com/satijalab/seurat/issues/1957
            # The SCTransform workflow results in Pearson residuals, which you
            # can think of as similar to what you get after ScaleData in the
            # standard workflow, which proceeds: NormalizeDate, ScaleData,
            # FindVariableFeatures
            #
            # OK to ignore: Warning message in theta.ml(y = y, mu = fit$fitted):
            # https://github.com/ChristophH/sctransform/issues/25 and
            # https://github.com/satijalab/seurat/issues/1378
            # These warnings are showing that there are some genes for which it
            # is hard to reliably estimate theta (presumably because of very few
            # non-zero observations). Usually we don't worry about these
            # warnings too much, since we regularize the parameters in a later
            # step, thus averaging out uncertainty of individual gene
            # parameters.
            #
            # NOTE on vars.to.regress
            # https://github.com/ChristophH/sctransform/issues/41
            withCallingHandlers({
                seurat_list[[i]] <- Seurat::SCTransform(
                    seurat_list[[i]],
                    do.correct.umi = TRUE,
                    vars.to.regress = sctransform_variable_regress,
                    variable.features.n = 3000,
                    do.scale = TRUE,
                    do.center = TRUE,
                    seed.use = 1448145, # 1448145 is the detault
                    verbose = verbose
                )
            }, warning = function(w) {
                if (startsWith(
                    conditionMessage(w),
                    "iteration limit reached"
                )) {
                    invokeRestart("muffleWarning")
                } else if (
                    stringr::str_detect(conditionMessage(w), "sqrt")
                ) { # TODO: fix catch
                    invokeRestart("muffleWarning")
                }
            }) # end CallingHandlers
        } else if (method == "LogNormalize") {
            assay_id <- "RNA"
            # the standard LogNormalize workflow is:
            #     NormalizeData, FindVariableFeatures, ScaleData
            # example:
            # https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html

            # https://satijalab.org/seurat/v3.1/immune_alignment.html
            # LogNormalize can be done on a per sample basis and merged
            # to integrate data with different exposures.
            seurat_list[[i]] <- Seurat::NormalizeData(
                seurat_list[[i]],
                normalization.method = "LogNormalize",
                scale.factor = 10000, # default is 10000
                verbose = verbose
            )
            seurat_list[[i]] <- Seurat::FindVariableFeatures(
                seurat_list[[i]],
                selection.method = "vst",
                mean.function = Seurat::ExpMean,
                dispersion.function = Seurat::LogVMR,
                nfeatures = 3000,
                verbose = verbose
            )
        } else {
            stop("[normalize_seurat_obj_list] ERROR")
        } # end method == "SCTransform"

        # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO
        # gene names
        #filter_variable_genes_regex <- "^HLA-|^IG[HJKL]|^RNA|^MT|^RP"
        if (filter_variable_genes_regex != "") {
            var_genes <- seurat_list[[i]]@assays[[assay_id]]@var.features
            cat(
                "[normalize_seurat_obj_list] ",
                "number var genes prior to filters:\t",
                length(var_genes),
                "\n"
            )
            var_genes <- grep(
                filter_variable_genes_regex,
                var_genes,
                invert = TRUE,
                value = TRUE
            )
            # NOTE: does not remove features, just drops them from the
            #       variable features used for clustering
            seurat_list[[i]]@assays[[assay_id]]@var.features <- intersect(
                var_genes,
                rownames(seurat_list[[i]]@assays[[assay_id]]@data)
            )
            cat(
                "[normalize_seurat_obj_list] ",
                "number var genes after filters:\t",
                length(var_genes),
                "\n"
            )
        } # end  if (filter_variable_genes_regex != "") {
    } # end for each sample loop

    return(seurat_list)
}



#' Integrates a list of Seurat data objects.
#'
#' \code{integrate_seurat_obj_list} integrates a list of Seurat data objects.
#'
#' @param seurat_list List.
#'     List of Seurat data objects.
#' @param method_norm String.
#'     Normalisation method. Valid options: c("LogNormalize", "SCTransform").
#' @param method_integrate String.
#'     Integration method. Valid options: c("seurat", "harmony",
#'     "seurat-harmony"). If "seurat" is in the integration method, the data is
#'     integrated using Seurat 3 which performs canonical correlation analysis
#'     and then mutual nearest neighbors to identify anchors between datasets.
#'     If "seurat-harmony" then integrate using Seurat 3 then run
#'     harmony on the resulting reduced dimensions. Seurat-harmony would be
#'     particularly useful to integrate samples (using Seurat) and then to
#'     further remove any batch effects from the data using harmony (e.g., if
#'     some samples were processed on one chip and other samples on another).
#'     If one wanted to merge the data without using any correction, then
#'     run in "harmony" mode and use the "pca" reductions in the output. All
#'     harmony adjusted reductions will be labeled as "harmony" followed by the
#'     theta values used.
#' @param n_integration_pcs Numeric.
#'     Number of PCs for data integration (for Seurat::FindIntegrationAnchors
#'     and Seurat::FindIntegrationAnchor calls).
#' @param n_pcs_dim_reduction Numeric.
#'     Number of PCs to calculate in the final merged data. If 0, then defaults
#'     to 10% of the total number of cells. If this value is < 200, then will
#'     override and set n_pcs_dim_reduction to 200.
#' @param scale_variable_regress List.
#'     Variables in meta.data to regress out of assay in Seurat::ScaleData
#'     function call. Passed to Seurat::ScaleData.
#' @param scale_split_by List
#'     Name of variable in metadata of Seurat data object or a vector or
#'     factor defining grouping of cells for Seurat::ScaleData function call.
#'     Passed to Seurat::ScaleData.
#' @param run_jackstraw Logical.
#'     Estimate number of PCs using jackstraw.
#' @param harmony_covariates List.
#'     Name of variables in metadata of Seurat data object to use as covariates
#'     in harmony analysis.
#' @param harmony_theta Data.frame
#'     Data.frame of theta values to use for each harmony covariate. Columns =
#'     covariate id and rows = different theta combinations. Example:
#'     data.frame("rep_id" = c(0.1, 0.1), "smpl_id" = c(0.1, 4.5)) would run
#'     (i) iteration 1: where theta for rep_id = 0.1 and smpl_id = 0.1 and
#'     (ii) iteration 2: where heta for rep_id = 0.1 and smpl_id = 4.5. If
#'     the data.frame is empty, will default to a theta of 2 for all covariates.
#' @param verbose Logical.
#'     If TRUE, print extra information to stdout.
#'
#' @return Seurat data object.
#'     An integrated Seurat data object with the below possible values in the
#'     "reductions" slot depending on the call parameters:
#'     <"seurat3__"|null><pca|harmony_[underscore seperated list of thetas that
#'     were used for covariates]>. For instance if the user specified
#'     "seurat-harmony" for method_integrate and theta of 0.3 and 1, the
#'     resulting reductions would be seurat__pca (normal pcs without harmony)
#'     and seurat__harmony_0pt3_1 (normal pcs adjusted by harmony with two
#'     covariates whose corresponding thetas are 0.3 and 1).
#'
#' @import Seurat
#' @importFrom harmony RunHarmony
#' @export
integrate_seurat_obj_list <- function(
        seurat_list,
        method_norm = "LogNormalize",
        method_integrate = "seurat",
        n_integration_pcs = 30,
        n_pcs_dim_reduction = 0,
        scale_variable_regress = c(),
        scale_split_by = NULL,
        run_jackstraw = FALSE,
        harmony_covariates = c(),
        harmony_theta = data.frame(),
        verbose = TRUE
    ) {

    if (!(method_norm %in% c("LogNormalize", "SCTransform"))) {
        stop("[integrate_seurat_obj_list]:\tERROR in method_norm")
    }
    if (!(method_integrate %in% c("seurat", "harmony", "seurat-harmony"))) {
        stop("[integrate_seurat_obj_list]:\tERROR in method_integrate.")
    }
    if ((method_integrate %in% c("harmony", "seurat-harmony") &
        length(harmony_covariates) == 0)
    ) {
        stop("[integrate_seurat_obj_list]:\tERROR specify harmony_covariates.")
    }

    # init the assay id that will be used
    # if we don't use suerat integration, then valid options are RNA or SCT
    # otherwise, should be seurat-harmony.
    assay_id <- "RNA"
    if (method_norm == "SCTransform") {
        assay_id <- "SCT"
    }
    if (method_integrate %in% c("seurat", "seurat-harmony")) {
        assay_id <- "integrated"
    }

    # set the first part of the reduction id
    reduction_id_1 <- ""
    if (method_integrate %in% c("seurat", "seurat-harmony")) {
        reduction_id_1 <- "seurat3__"
    }
    reduction_pca_id <- paste0(reduction_id_1, "pca")

    # Select features for downstream integration
    # Choose the features to use when integrating multiple datasets.
    # This function ranks features by the number of datasets they appear
    # in, breaking ties by the median rank across datasets. It returns
    # the highest features by this ranking.
    var_features <- Seurat::SelectIntegrationFeatures(
        object.list = seurat_list,
        nfeatures = 3000, # default = 2000
        verbose = verbose
    )

    if (method_integrate == "harmony") {
        # get the union of the variable features for each sample
        # this is old... better to use Seurat::SelectIntegrationFeatures
        # var_features <- unique(unlist(lapply(seurat_list,
        #     FUN = function(x) {x@assays[[assay_id]]@var.features}
        # )))

        # Merge a list of seurat objects
        # NOTE: information apart from the normalized data (e.g., variable
        # features, scaled data will be lost)
        seurat_df <- merge_seurat_list(seurat_list)

        # Set the variable features of the merged slot
        seurat_df@assays[[assay_id]]@var.features <- var_features

        # scale the data
        # NOTE: harmony seems to work better if the data is centered and
        # scaled which effectively means one weights every gene the same...
        # otherwise gene weightings would change as function of magnitude and
        # variance of gene expression
        seurat_df <- Seurat::ScaleData(
            seurat_df,
            vars.to.regress = scale_variable_regress,
            model.use = "linear", # linear, negbinom, poisson
            do.scale = TRUE,
            do.center = TRUE,
            split.by = scale_split_by,
            assay = assay_id, # should be integrated
            verbose = verbose
        )

    } else if (method_integrate %in% c("seurat", "seurat-harmony")) {
        if (verbose) {
            cat("[integrate_seurat_obj_list]:\t",
                "integrating data using seurat.\n"
            )
        }

        if (method_norm == "SCTransform") {
            integrate_method_norm <- "SCT"
            # PrepSCTIntegration, which ensures that all necessary Pearson
            # residuals have been calculated.
            seurat_list <- Seurat::PrepSCTIntegration(
                object.list = seurat_list,
                anchor.features = var_features,
                verbose = verbose
            )
        } else if (method_norm == "LogNormalize") {
            integrate_method_norm <- "LogNormalize"
        } else {
            stop("[integrate_seurat_obj_list]:\tERROR method_norm")
        } # end if (method_norm == "SCTransform") {

        # Identify anchors
        # Commands are identical to the standard workflow, but
        # make sure to set normalization.method = 'SCT':
        anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurat_list,
            anchor.features = var_features,
            normalization.method = integrate_method_norm,
            dims = 1:n_integration_pcs,
            verbose = verbose
        )
        if (verbose) {
            cat("[integrate_seurat_obj_list]:\t",
                "completed Seurat::FindIntegrationAnchor.\n"
            )
        }

        # Integrate samples
        # NOTE: seurat_df@assays[["integrated"]]@var.features == var_features
        seurat_df <- Seurat::IntegrateData(
            anchorset = anchors,
            normalization.method = integrate_method_norm,
            dims = 1:n_integration_pcs,
            verbose = verbose
        )
        if (verbose) {
            cat("[integrate_seurat_obj_list]:\t",
                "completed Seurat::IntegrateData.\n"
            )
        }

        # After running Seurat::IntegrateData, the assay defaults to integrated
        # however, we can also set this explicity to make that clear
        # NOTE: throws error if the assay does not exist
        Seurat::DefaultAssay(object = seurat_df) <- "integrated"

        if (method_norm == "LogNormalize") {
            # scale the data
            seurat_df <- Seurat::ScaleData(
                seurat_df,
                vars.to.regress = scale_variable_regress,
                model.use = "linear", # linear, negbinom, poisson
                do.scale = TRUE,
                do.center = TRUE,
                split.by = scale_split_by,
                assay = assay_id, # should be integrated
                verbose = verbose
            )
        }
    } else {
        stop("[integrate_seurat_obj_list]:\tERROR invalid method_integrate")
    }# end if (method_integrate == "harmony") {

    # for max n_pcs of integrated use 10% of total cells
    if (n_pcs_dim_reduction == 0) {
        n_pcs_dim_reduction <- round(nrow(seurat_df@meta.data) * 0.1)
    }
    if (n_pcs_dim_reduction < 200) {
        n_pcs_dim_reduction <- 200
    }
    if (verbose) {
        cat("[integrate_seurat_obj_list]:\t",
            "n_pcs_dim_reduction = ", n_pcs_dim_reduction, "\n"
        )
    }

    # PCA
    # warning message about n pc is not an issue
    #     see https://github.com/satijalab/seurat/issues/1249
    seurat_df <- Seurat::RunPCA(
        object = seurat_df,
        assay = assay_id,
        npcs = n_pcs_dim_reduction,
        verbose = verbose,
        features = seurat_df@assays[[assay_id]]@var.features
    )
    seurat_df@reductions[[reduction_pca_id]] <- seurat_df@reductions[["pca"]]
    seurat_df@reductions[["pca"]] <- NULL

    if (run_jackstraw) {
        if (verbose) {
            cat("[integrate_seurat_obj_list]:\trunning jackstraw\n")
        }
        seurat_df <- Seurat::JackStraw(
            object = seurat_df,
            reduction = reduction_pca_id,
            num.replicate = 100,
            dims = n_pcs_dim_reduction
        )
        seurat_df <- Seurat::ScoreJackStraw(
            object = seurat_df,
            reduction = reduction_pca_id,
            dims = 1:n_pcs_dim_reduction
        )
        n_sig_pcs <- sum(
            seurat_df@reductions$pca@jackstraw$overall.p.values[
                , "Score"
            ] < 0.05
        )
        if (verbose) {
            cat("[integrate_seurat_obj_list]:\tn_sig_pcs (jackstraw) = ",
                n_sig_pcs,
                "\n"
            )
        }
    }

    # if harmony, then we need to update the PC embeddings
    if (method_integrate %in% c("harmony", "seurat-harmony")) {
        if (verbose) {
            cat("[integrate_seurat_obj_list]:\t",
                "integrating data using harmony.\n"
            )
        }
        if ((length(harmony_covariates) > 0) &
            (ncol(harmony_theta) != length(harmony_covariates))) {
            cat("[integrate_seurat_obj_list]:\t",
                "setting harmony_theta to 2 for all covariates.\n"
            )
            # set harmony_theta to default value of 2
            harmony_theta <- data.frame(
                t(rep(2, length(harmony_covariates)))
            )
            colnames(harmony_theta) <- harmony_covariates
        }
        col_check <- unlist(lapply(harmony_covariates, FUN = function(x) {
            return(x %in% colnames(harmony_theta))
        }))
        if (!all(col_check)) {
            stop("[integrate_seurat_obj_list]:\t",
                "ERROR: harmony covariate names and columns of theta matrix",
                 " do not match."
            )
        }
        # One could also integrate using Harmony. The advantage of Harmony is
        # that one can also include covariates. Note that Harmony acts on PC
        # embeddings so one does not split the data into each sample prior to
        # integration
        for (i in 1:nrow(harmony_theta)) {
            reduction_id <- gsub(
                "\\.",
                "pt",
                paste(c("harmony", harmony_theta[i,]), collapse = "_")
            )
            reduction_id <- paste0(reduction_id_1, reduction_id)
            if (verbose) {
                cat("[integrate_seurat_obj_list]:",
                    "harmony:\t", reduction_id, "\n"
                )
            }
            seurat_df <- harmony::RunHarmony(
                object = seurat_df,
                reduction = reduction_pca_id,
                group.by.vars = harmony_covariates,
                theta = unlist(harmony_theta[1,harmony_covariates]),
                dims.use = NULL, # if null use max # of pcs,
                #reduction.save = reduction_id, # where to save the embeddings
                assay.use = assay_id, # uses the normalized RNA data
                epsilon.cluster = -Inf, # Set to -Inf to never stop early.
                epsilon.harmony = -Inf # Set to -Inf to never stop early.
            )
            # to fix an error in harmony reduction.save
            seurat_df@reductions[[reduction_id]] <- seurat_df@reductions[[
              "harmony"
            ]]
            seurat_df@reductions[["harmony"]] <- NULL
        } # end harmony covariates itration
    } # end if (method_integrate == "harmony") {

    return(seurat_df)
} # end integrate_list



#' Merges a list of Seurat data objects.
#'
#' \code{merge_seurat_list} normalises each sample in Seurat dataframe
#' that has multiple samples.
#'
#' @param seurat_list List.
#'     List of multiple Seurat data objects to be merged.
#' @param list_ids List.
#'     List of character ids for each element of seurat_list. This will be
#'     used to rename cells with Seurat::RenameCells add.cell.id.
#'
#' @return Seurat data object.
#'     A Seurat data object that has been merged from the data.
#'
#' @importFrom Seurat RenameCells
#' @export
merge_seurat_obj_list <- function(seurat_list, list_ids = c()) {
    sc_df <- seurat_list[[1]]
    if (length(list_ids) > 0) {
        sc_df <- Seurat::RenameCells(
            object = sc_df,
            add.cell.id = names(list_ids)[[1]],
            for.merge = TRUE
        )
    }

    for (i in seq(2, length(seurat_list))) {
        if (length(list_ids) > 0) {
            seurat_obj <- Seurat::RenameCells(
                object = seurat_list[[i]],
                add.cell.id = names(list_ids)[[i]],
                for.merge = TRUE
            )
        } else {
            seurat_obj <- seurat_list[[i]]
        }
        # The merge will not preserve reductions, graphs, logged commands,
        # or feature-level metadata that were present in the original objects.
        sc_df <- merge(
            sc_df,
            seurat_obj,
            merge.data = TRUE
        )
    } # end for loop

    return(sc_df)
}



#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("-f", "--file"),
            type = "character",
            help = paste0(
                "Rds of a Seurat object."
            )
        ),

        optparse::make_option(c("--metadata_split_column"),
            type = "character",
            default = "sanger_sample_id",
            help = paste0(
                "Metadata column to split by for merge (e.g., orig.ident).",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--normalization_method"),
            type = "character",
            default = "LogNormalize",
            help = paste0(
                "Normalization method. Valid options LogNormalize or SCT for",
                " SCTransform.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--integrate"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Perform Seurat 3 integration. If FALSE, then data is merged,",
                " as cell readouts are concatenated.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--out_file"),
            type = "character",
            default = "sc_dat",
            help = paste0(
                "Name (and possibly path) of output file. Will have rds.gz",
                " appended to it.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--verbose"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Verbose mode (write extra info to std.err).",
                " [default: %default]"
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Reads in a Seurat object, performs normalization and optionally",
            " performs integration."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
    # if you also have non-boolean optional args:
    getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for (item in parserObj@options) {
            optionStrings <- append(optionStrings,
                                    c(item@short_flag, item@long_flag))
        }
        optionStrings
    }
    optStrings <- getOptionStrings(parser)
    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    # read in the data
    sc_df_tmp <- readRDS(param[["file"]])
    # seurat_list <- Seurat::SplitObject(
    #     sc_df,
    #     split.by = param[["metadata_split_column"]]
    # )

    # Run the proper normalization
    seurat_list <- normalize_seurat_obj_list(
        seurat_list = Seurat::SplitObject(
            sc_df_tmp,
            split.by = param[["metadata_split_column"]]
        ),
        method = param[["normalization_method"]],
        filter_variable_genes_regex = "",
    )

    # Either integrate the final data using
    # * Seurat 3 integration
    # * perform no integration - just concatenating data across samples (merge)
    #
    # Note this also runs PC calculation
    if (param[["integrate"]]) {
        sc_df <- integrate_seurat_obj_list(
            seurat_list,
            method_norm = param[["normalization_method"]],
            method_integrate = "seurat",
            n_integration_pcs = 30,
            verbose = TRUE,
            harmony_covariates = param[["metadata_split_column"]]
        )
    } else {
        sc_df <- merge_seurat_obj_list(
            seurat_list
        )
    }

    # save final dataframe
    saveRDS(
        sc_df,
        file = paste0(param[["out_file"]], ".rds.gz"),
        compress = TRUE
    )

    return(0)
}

main <- function() {
    # run analysis
    run_time <- system.time(df_results <- command_line_interface())
    message(paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours."
    ))
    return(0)
}



#' Quick Uniform Manifold Approximation and Projection (UMAP) plots.
#'
#' \code{plot_umap_devonly}: UMAP plots of Seurat data object.
#'
#' @param seurat_df Seurat data object.
#'     Seurat data object.
#' @param n_pcs_analysis Integer.
#'     Number of PCs to use for projections. If NULL, calculates number of
#'     PCs using the inflection point of the knee plot.
#' @param plot_knee Logical.
#'     If TRUE, make the knee plot.
#' @param reduction_id String.
#'     Id of the reduction variable in the Seurat data object to plot.
#' @param factor_ids List.
#'     List of factor columns in meta.data to color cells by.
#' @param numeric_ids List.
#'     List of numeric columns in meta.data to color cells by.
#' @param title String.
#'     Plot title.
#' @param verbose Logical.
#'     If TRUE, print extra information to stdout.
#'
#' @return List.
#'     List of ggplot2 objects.
#'
#' @import ggplot2
#' @importFrom Seurat RunUMAP
#' @importFrom inflection uik
#' @export
plot_umap_devonly <- function(
        seurat_df,
        n_pcs_analysis = NULL,
        plot_knee = FALSE,
        reduction_id = "pca",
        factor_ids = c("sample_id"),
        numeric_ids = c("nCount_RNA"),
        title = "",
        verbose = TRUE
    ) {

    plt_list <- list()
    plt_i <- 1

    if (is.null(n_pcs_analysis)) {
        # Calculate knee (recommended): http:/dx.doi.org/10.2139/ssrn.3043076
        stdev_rank <- rank(seurat_df@reductions[['pca']]@stdev*-1)
        n_pcs_analysis <- inflection::uik(
            stdev_rank,
            log10(seurat_df@reductions[['pca']]@stdev)
        )
    }
    if (verbose) {
        cat("n_pcs_analysis:\t", n_pcs_analysis, "\n")
    }

    if (plot_knee) {
        plt_data <- data.frame(
            x = stdev_rank,
            y = seurat_df@reductions[['pca']]@stdev
        )
        plt <- ggplot2::ggplot(plt_data, ggplot2::aes(
            x = x,
            y = y
        ))
        plt <- plt + ggplot2::theme_bw(base_size = 12)
        plt <- plt + ggplot2::geom_point()
        plt <- plt + ggplot2::labs(
            x = "Principal component",
            y = "Standard deviation",
            title = title
        )
        plt <- plt + ggplot2::geom_vline(
            xintercept = n_pcs_analysis,
            linetype = "dashed"
        )
        plt_list[[plt_i]] <- plt
        plt_i <- plt_i + 1
    }

    # run UMAP
    seurat_df <- Seurat::RunUMAP(
        seurat_df,
        reduction = reduction_id,
        dims = 1:n_pcs_analysis
    )

    if (length(factor_ids) > 0) {
        for(i in 1:length(factor_ids)) {
            if (verbose) {
                cat("umap plot for:\t", factor_ids[[i]], "\n")
            }
            # plt_list[[plt_i]] <- Seurat::DimPlot(
            #     seurat_df,
            #     reduction = "umap",
            #     pt.size = 0.1,
            #     alpha = 0.5,
            #     label = FALSE,
            #     group.by = factor_ids[[i]]
            # )
            plt_data <- data.frame(
                x = seurat_df@reductions[["umap"]]@cell.embeddings[,1],
                y = seurat_df@reductions[["umap"]]@cell.embeddings[,2],
                color = seurat_df@meta.data[[factor_ids[[i]]]]
            )
            plt <- ggplot2::ggplot(plt_data, ggplot2::aes(
                x = x,
                y = y,
                color = color
            ))
            plt <- plt + ggplot2::theme_bw(base_size = 12)
            plt <- plt + ggplot2::geom_point(alpha = 0.5, size = 0.1)
            plt <- plt + ggplot2::labs(
                x = "UMAP 1",
                y = "UMAP 2",
                color = factor_ids[[i]],
                title = title
            )
            plt <- plt + ggplot2::guides(
                color = ggplot2::guide_legend(
                    override.aes = list(alpha = 1, size = 1)
                )
            )
            plt_list[[plt_i]] <- plt
            plt_i <- plt_i + 1
        }
    }

    if (length(numeric_ids) > 0) {
        for(i in 1:length(numeric_ids)) {
            if (verbose) {
                cat("umap plot for:\t", numeric_ids[[i]], "\n")
            }
            #plt_list[[plt_i]] <- Seurat::FeaturePlot(
            #    object = seurat_df,
            #    reduction = "umap",
            #    pt.size = 0.1,
            #    features = numeric_ids[[i]]
            #)
            plt_data <- data.frame(
                x = seurat_df@reductions[["umap"]]@cell.embeddings[,1],
                y = seurat_df@reductions[["umap"]]@cell.embeddings[,2],
                color = as.numeric(seurat_df@meta.data[[numeric_ids[[i]]]])
            )
            plt <- ggplot2::ggplot(plt_data, ggplot2::aes(
                x = x,
                y = y,
                color = color
            ))
            plt <- plt + ggplot2::theme_bw(base_size = 12)
            plt <- plt + ggplot2::geom_point(alpha = 0.5, size = 0.1)
            plt <- plt + ggplot2::labs(
                x = "UMAP 1",
                y = "UMAP 2",
                color = numeric_ids[[i]],
                title = title
            )
            plt <- plt + ggplot2::guides(
                color = ggplot2::guide_legend(
                    override.aes = list(alpha = 1, size = 1)
                )
            )
            plt <- plt + ggplot2::scale_colour_gradientn(
                colors = wesanderson::wes_palette(
                    "Zissou1",
                    100,
                    type = "continuous"
                ),
                guide = "colorbar"
            )
            plt_list[[plt_i]] <- plt
            plt_i <- plt_i + 1
        }
    }

    return(plt_list)
}



# code for development
dev <- function() {
    require(Seurat)
    require(tibble)
    require(reshape2)
    require(dplyr)
    if (!require(SeuratData)) {
        devtools::install_github('satijalab/seurat-data')
        require(SeuratData)
        SeuratData::InstallData("panc8.SeuratData")
    }
    utils::data("panc8")

    filt <- grepl("indrop", panc8@meta.data[["dataset"]])
    df_dev <- subset(
        panc8,
        cells = colnames(panc8@assays[['RNA']]@counts)[filt]
    )
    df_dev@meta.data[["sample_id"]] <- df_dev@meta.data[["orig.ident"]]


    method_norm <- "LogNormalize" # SCTransform or LogNormalize
    seurat_list <- Seurat::SplitObject(df_dev, split.by = "sample_id")
    seurat_list <- normalize_seurat_obj_list(
        seurat_list,
        filter_variable_genes_regex = "^HLA-|^IG[HJKL]|^RNA|^MT|^RP",
        #sctransform_variable_regress = c("sample_id"),
        method = method_norm # SCTransform or LogNormalize
    )
    harmony_covariates <- c("replicate", "sample_id")
    # harmony_theta ncol match harmony_covariates
    harmony_theta <- data.frame(
      "replicate" = c(0.1, 0.1),
      "sample_id" = c(0.1, 4.5)
    )
    df_dev3 <- integrate_seurat_obj_list(
        seurat_list,
        method_norm = method_norm,
        method_integrate = "seurat-harmony", # seurat, harmony, seurat-harmony
        n_integration_pcs = 30,
        n_pcs_dim_reduction = 200,
        scale_variable_regress = c(),
        scale_split_by = NULL,
        run_jackstraw = FALSE,
        harmony_covariates = harmony_covariates,
        harmony_theta = harmony_theta,
        verbose = TRUE
    )
    print(names(df_dev3@reductions))
    n_pcs_max <- ncol(df_dev3@reductions[["seurat3__pca"]]@cell.embeddings)

    plt_list <- list()
    plt_i <- 1
    pdf(file = paste0("dev.pdf"), height = 10, width = 10)

    # add pca plots
    plts <- plot_umap_devonly(
        df_dev3,
        reduction_id = "seurat3__pca",
        factor_ids = c("sample_id"),
        numeric_ids = c(),
        title = "seurat3__pca"
    )
    for (i in 1:length(plts)) {
        print(plts[[i]])
        plt_list[[plt_i]] <- plts[[i]]
        plt_i <- plt_i + 1
    }

    # add harmony embeddings if available
    harmony_embeddings <- unlist(
        grep("harmony", names(df_dev3@reductions), value = TRUE)
    )
    if (length(harmony_embeddings) > 0) {
        for (i in 1:length(harmony_embeddings)) {
            cat("harmony:\t", harmony_embeddings[[i]], "\n")
            plts <- plot_umap_devonly(
                df_dev3,
                reduction_id = harmony_embeddings[[i]],
                factor_ids = c("sample_id"),
                numeric_ids = c(),
                title = harmony_embeddings[[i]]
            )
            for (i in 1:length(plts)) {
                print(plts[[i]])
                plt_list[[plt_i]] <- plts[[i]]
                plt_i <- plt_i + 1
            }
        }
    }

    # Local inverse simpsons index (LISI) is a measurement of mixing
    # can apply to data labelled by sample (integration) or batch or cell type.
    # * If cell type that don't want mixed, value near 1 is good. In this
    # case LISI represents the number of cells from each cell type in a
    # neighborhood.
    # * If sample then a value close to the number of samples is good. In this
    # case LISI represents the number of cells from each sample in a
    # neighborhood.
    #
    # NOTE: need to run require(Rcpp) to make this work properly
    lisi_res_list <- list()
    for (i in names(df_dev3@reductions)) {
        lisi_res_list[[i]] <- lisi::compute_lisi(
            df_dev3@reductions[[i]]@cell.embeddings,
            df_dev3@meta.data,
            harmony_covariates # variables to compute LISI for
        ) %>%
        tibble::rownames_to_column(var = "cell_id") %>%
        reshape2::melt(
            id.vars = c("cell_id"),
            variable.name = "label",
            value.name = "lisi"
        )
        lisi_res_list[[i]]$reduction <- i
    }
    lisi_res <- data.table::rbindlist(lisi_res_list)

    plt <- ggplot2::ggplot(lisi_res, ggplot2::aes(
        x = lisi,
        fill = reduction
    ))
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::geom_density(alpha = 0.25)
    plt <- plt + ggplot2::facet_grid(label ~ .)
    plt <- plt + ggplot2::labs(
        x = "LISI",
        #y = "Cumulative density",
        color = "Reduction",
        title = ""
    )
    print(plt)

    plt <- ggplot2::ggplot(lisi_res, ggplot2::aes(
        x = lisi,
        color = reduction
    ))
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::stat_ecdf(alpha = 0.5)
    plt <- plt + ggplot2::facet_grid(label ~ .)
    plt <- plt + ggplot2::labs(
        x = "LISI",
        y = "Cumulative density",
        color = "Reduction",
        title = ""
    )
    print(plt)
    plt_list[[plt_i]] <- plt
    plt_i <- plt_i + 1

    # all plots in matrix
    if (length(plt_list) > 0) {
        n_col <- floor(sqrt(length(plt_list)))
        do.call(gridExtra::grid.arrange, c(plt_list, ncol = n_col))
    }

    dev.off()
    print(paste0(getwd(), "/dev.pdf"))
}



# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
