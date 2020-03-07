#!/usr/bin/env Rscript

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

SCRIPT_NAME <- "plot_umap.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("inflection"))
suppressPackageStartupMessages(library("Seurat"))



#' Quick Uniform Manifold Approximation and Projection (UMAP) plots.
#'
#' \code{plot_umap}: UMAP plots of Seurat data object.
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
plot_umap <- function(
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



#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("-f", "--files_merge"),
            type = "character",
            help = paste0(
                "Tab-delimited input file, col1 = 'experiment_id',",
                " col2 = 'path_data_10xformat'."
            )
        ),

        optparse::make_option(c("--metadata_file"),
            type = "character",
            help = paste0(
              "Tab-delimited metadata file, must have a column labelled",
              " sanger_sample_id that maps to path_data_10xformat."
            )
        ),

        optparse::make_option(c("--out_file"),
            type = "character",
            default = "sc_df",
            help = paste0(
                "Name (and possibly path) of output file. Will have rds.gz",
                " appended to it. ",
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
            "Makes Seurat data object with metadata on samples and cells."
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
    files_df <- read.csv(param[["files_merge"]], sep = "\t")

    metadata_df <- read.csv(param[["metadata_file"]], sep = "\t")
    rownames(metadata_df) <- as.character(metadata_df$sanger_sample_id)
    # add metadata info on time to chromium
    metadata_df$time_to_chromium_processing <- chron::times(
            metadata_df[["chromium_time"]]
        ) - chron::times(metadata_df[["collection_time"]])

    # format files for merge seurat dataframe script
    files <- structure(
        as.character(files_df$path_data_10xformat),
        names = as.character(files_df$experiment_id)
    )

    # merge the data into a seurat data object
    store_genes_as_ensembl_ids <- FALSE
    sc_df <- load_10x_data_seurat(
        files = files,
        metadata_df = metadata_df,
        store_genes_as_ensembl_ids = store_genes_as_ensembl_ids,
        min_avg_counts = 0,
        min_n_cells_exprsing_gene = 0,
        seurat_min_cells = 0,
        seurat_min_features = 0,
        verbose = param[["verbose"]]
    )

    # remove ERCCs
    remove_erccs <- TRUE
    if (remove_erccs) {
        filt <- grepl("^ERCC", rownames(sc_df))
        if (param[["verbose"]]) {
            cat("Dropping ERCCs from data. Summary of ERCC counts:\n")
            print(summary(Matrix::colSums(sc_df@assays$RNA@counts[filt,])))
        }
        sc_df <- sc_df[!filt,]
    }

    # update basic metrics
    # sc_df@meta.data[["DetFeatures"]] <- Matrix::colSums(
    #     sc_df@assays$RNA@counts != 0
    # )
    # sc_df@meta.data[["TotalCounts"]] <- Matrix::colSums(
    #     sc_df@assays$RNA@counts
    # )

    # add additional info: calculate fraction MT counts
    if (!store_genes_as_ensembl_ids) {
        gene_filter <- grep(
            pattern = "^MT-",
            x = rownames(sc_df@assays[['RNA']]@counts),
            value = TRUE
        )
        sc_df@meta.data[["mitoProportionCounts"]] <- Matrix::colSums(
            sc_df@assays[['RNA']]@counts[gene_filter, ]
        ) / Matrix::colSums(sc_df@assays[['RNA']]@counts)
    }

    # TODO convert to single cell experiment and annotate cell cycle cells with
    # cyclone

    # TODO run scrublet and annotate cells with likelihood of >1 cells being
    # present.

    # save final un-normalized dataframe\
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

    # NOTE: normalization should be run per sample, this is just for demo
    # purposes only
    df_dev <- Seurat::NormalizeData(
        df_dev,
        normalization.method = "LogNormalize",
        scale.factor = 10000 # default is 10000
    )
    df_dev <- Seurat::FindVariableFeatures(
        df_dev,
        selection.method = "vst",
        mean.function = Seurat::ExpMean,
        dispersion.function = Seurat::LogVMR,
        nfeatures = 3000
    )
    df_dev <- Seurat::ScaleData(
        df_dev,
        model.use = "linear", # linear, negbinom, poisson
        do.scale = TRUE,
        do.center = TRUE,
        assay = "RNA"
    )
    # for max n_pcs use 10% of total cells
    #n_pcs_dim_reduction <- round(nrow(df_dev@meta.data) * 0.1)
    n_pcs_dim_reduction <- 25
    # PCA
    # warning message about n pc is not an issue
    #     see https://github.com/satijalab/seurat/issues/1249
    df_dev <- Seurat::RunPCA(
        object = df_dev,
        assay = "RNA",
        npcs = n_pcs_dim_reduction,
        features = df_dev@assays[["RNA"]]@var.features
    )

    # add pca plots
    plts <- plot_umap(
        df_dev,
        reduction_id = "pca",
        factor_ids = c("sample_id"),
        numeric_ids = c(),
        title = "pca"
    )

    pdf(file = paste0("dev.pdf"), height = 10, width = 10)
    for (i in 1:length(plts)) {
        print(plts[[i]])
    }
    dev.off()
    print(paste0(getwd(), "/dev.pdf"))
}



# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    dev()
    #main()
}
