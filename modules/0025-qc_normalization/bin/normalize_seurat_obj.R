#!/usr/bin/env Rscript

SCRIPT_NAME <- "normalize_seurat_obj.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("future"))

# need to increase the RAM
# see https://satijalab.org/seurat/v3.0/future_vignette.html
future::plan("multiprocess", workers = 16)
options(future.globals.maxSize = 100000 * 1024^2)


#' Normalizes a list of Seurat objects
#'
#' \code{normalize_seurat_obj_list} takes a list of Seurat objects and 
#' normalizes.
#'
#' @param seurat_list List of Seurat objects.
#'     List of Seurat objects - see Seurat::SplitObject. 
#' @param method String.
#'     Normalization method. Valid options: LogNormalize or SCT. 
#'     Default = LogNormalize.
#'     
#' @return List of Seurat objects.
#'     List of Seurat objects with normalization performed.
#'
#' @importFrom Seurat SCTransform
#' @importFrom Seurat NormalizeData
#' @importFrom stringr str_detect
#' @export
normalize_seurat_obj_list <- function(
        seurat_list,
        method = "LogNormalize",
        verbose = TRUE
    ) {
    for (i in 1:length(seurat_list)) {
        if (verbose) {
            cat(
                "[normalize_seurat_obj_list]: ",
                i, "/", length(seurat_list), "\n"
            )   
        }
        
        if (method == "SCT") {
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
            withCallingHandlers({
                seurat_list[[i]] <- Seurat::SCTransform(
                    seurat_list[[i]],
                    do.correct.umi = TRUE,
                    #vars.to.regress = c("mitoProportionCounts"),
                    variable.features.n = 3000,
                    do.scale = TRUE,
                    do.center = TRUE,
                    seed.use = 1448145,
                    verbose = FALSE
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
            # https://satijalab.org/seurat/v3.1/immune_alignment.html
            # LogNormalize can be done on a per sample basis and merged
            # to integrate data with different exposures.
            seurat_list[[i]] <- Seurat::NormalizeData(
                seurat_list[[i]],
                normalization.method = "LogNormalize",
                scale.factor = 10000, # default is 10000
                verbose = FALSE
            )
        } else {
            stop("ERROR")
        } # end method == "SCTransform"
    } # end for
    
    return(seurat_list)
}


#' Normalizes a list of Seurat objects
#'
#' \code{integrate_seurat_obj_list} takes a list of Seurat objects and 
#' performs Seurat 3 integration
#'
#' @param seurat_list List of Seurat objects.
#'     List of Seurat objects that have been normalized.
#' @param method String.
#'     Normalization method - this changes the proper way to call Seurat 
#'     integration. Valid options: LogNormalize or SCT. 
#'     Default = LogNormalize.
#' @param n_integration_pcs Integer.
#'     Number of PCs to use for integration.
#'     Default = 30.
#'
#' @return List of Seurat objects.
#'     List of Seurat objects with normalization performed.
#'
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ExpMean
#' @importFrom Seurat LogVMR
#' @importFrom Seurat SelectIntegrationFeatures
#' @importFrom Seurat PrepSCTIntegration
#' @importFrom Seurat FindIntegrationAnchors
#' @importFrom Seurat IntegrateData
#' @export
integrate_seurat_obj_list <- function(
        seurat_list,
        method = "LogNormalize",
        filter_variable_genes_regex = "",
        n_integration_pcs = 30,
        verbose = TRUE
    ) {
        
    if (method == "SCT") {
        # Next, select features for downstream integration, and run
        # PrepSCTIntegration, which ensures that all necessary Pearson
        # residuals have been calculated.
        if (verbose) {
            cat("[integrate_seurat_obj_list]: ",
                "selecting features for downstream integration\n"
            )
        }
        seurat_list_features <- Seurat::SelectIntegrationFeatures(
            object.list = seurat_list,
            nfeatures = 3000, # default = 2000
            verbose = FALSE
        )
        seurat_list <- Seurat::PrepSCTIntegration(
            object.list = seurat_list,
            anchor.features = seurat_list_features,
            verbose = TRUE
        )
        
        # Identify anchors
        # Commands are identical to the standard workflow, but
        # make sure to set normalization.method = 'SCT':
        anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurat_list,
            anchor.features = seurat_list_features,
            normalization.method = "SCT",
            dims = 1:n_integration_pcs,
            verbose = FALSE
        )
        
        # Integrate samples
        sc_df_integrated <- Seurat::IntegrateData(
            anchorset = anchors,
            normalization.method = "SCT",
            dims = 1:n_integration_pcs,
            verbose = FALSE
        )
    } else if (method == "LogNormalize") {
        for (i in 1:length(seurat_list)) {
            if (verbose) {
                cat(
                    "[integrate_seurat_obj_list]: ",
                    i, "/", length(seurat_list), "\n"
                )  
            }
            seurat_list[[i]] <- Seurat::FindVariableFeatures(
                seurat_list[[i]],
                selection.method = "vst",
                mean.function = Seurat::ExpMean,
                dispersion.function = Seurat::LogVMR,
                nfeatures = 3000,
                verbose = FALSE
            )
            
            # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO
            # gene names
            #var_regex <- '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
            if (filter_variable_genes_regex != "") {
                var_genes <- seurat_list[[i]]@assays[[assay]]@var.features
                if (verbose) {
                    cat(
                        "Number var genes prior to filters:\t",
                        length(var_genes),
                        "\n"
                    )
                }
                var_genes <- grep(
                    filter_variable_genes_regex,
                    var_genes,
                    invert = TRUE,
                    value = TRUE
                )
                # NOTE: does not remove features, just drops them from the
                #       variable features used for clustering
                seurat_list[[i]]@assays[[assay]]@var.features <- intersect(
                    var_genes,
                    rownames(seurat_list[[i]]@assays[[assay]]@data)
                )
                cat(
                    "Number var genes after filters:\t",
                    length(var_genes),
                    "\n"
                )
            } # end filter_variable_genes_regex
        } # end seurat_list iteration
        
        if (verbose) {
            cat("Running integration anchors\n")
        }
        anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurat_list,
            normalization.method = "LogNormalize",
            dims = 1:n_integration_pcs,
            verbose = FALSE
        )
        if (verbose) {
            cat("found anchors\n")
        }
        
        # Integrate samples
        sc_df_integrated <- Seurat::IntegrateData(
            anchorset = anchors,
            normalization.method = "LogNormalize",
            # dims = 1:n_integration_pcs,
            # verbose = FALSE
        )
        if (verbose) {
            cat("done with integration\n")
        }    
        Seurat::DefaultAssay(object = sc_df_integrated) <- "integrated"
        
        # Now scale the data
        if (verbose) {
            cat("scaling the data\n")
        }
        sc_df_integrated <- Seurat::ScaleData(
            sc_df_integrated,
            vars.to.regress = c(
                "nCount_RNA"
                #"orig.ident",
                #"active.ident",
                #"S.Score",
                #"G2M.Score",
                #"mitoProportionCounts"
            ),
            model.use = "linear", # linear, negbinom, poisson
            do.scale = TRUE,
            do.center = TRUE
        )
    } else {
        stop("ERROR")
    }
    
    return(sc_df_integrated)
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
            "Reads in a Seurat object and performs normalization."
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
        method = param[["normalization_method"]]
    )
    
    # Either integrate the final data using 
    # * Seurat 3 integration
    # * perform no integration - just concatenating data across samples (merge)
    if (param[["integrate"]]) {
        sc_df <- integrate_seurat_obj_list(
            seurat_list,
            method = param[["normalization_method"]],
            filter_variable_genes_regex = "",
            n_integration_pcs = 30,
            verbose = TRUE
        )
    } else {
        stop()
    }
    
    # calculate PCs on data
    # if (verbose) {
    #     cat("calculating n_pcs_dim_reduction\n")
    # }
    # # for max n_pcs of integrated use 10% of total cells
    # n_pcs_dim_reduction <- round(nrow(sc_df_integrated@meta.data) * 0.1)
    # # if (n_pcs_dim_reduction > 200) {
    # #     n_pcs_dim_reduction <- 200
    # # }
    # n_pcs_dim_reduction <- 200
    # if (verbose) {
    #     cat("n_pcs_dim_reduction:\t", n_pcs_dim_reduction, "\n")
    # }
    # PCA
    # warning message about n pc is not an issue
    #     see https://github.com/satijalab/seurat/issues/1249
    # sc_df_integrated <- Seurat::RunPCA(
    #     object = sc_df_integrated,
    #     npcs = n_pcs_dim_reduction,
    #     verbose = FALSE
    # )
    # run_jackstraw <- FALSE
    # if (run_jackstraw) {
    #     if (verbose) {
    #         cat("running jackstraw\n")
    #     }
    #     sc_df_integrated <- Seurat::JackStraw(
    #         object = sc_df_integrated,
    #         num.replicate = 100,
    #         dims = n_pcs_dim_reduction
    #     )
    #     sc_df_integrated <- Seurat::ScoreJackStraw(
    #         object = sc_df_integrated,
    #         dims = 1:n_pcs_dim_reduction
    #     )
    #     n_sig_pcs <- sum(
    #         sc_df_integrated@reductions$pca@jackstraw$overall.p.values[
    #             , "Score"
    #             ] < 0.05
    #     )
    #     cat("n_sig_pcs (jackstraw):\t", n_sig_pcs, "\n")
    # }
    
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


# code for development
dev <- function() {
    sc_df <- readRDS("sc_df.rds.gz")
}

# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
