#!/usr/bin/env Rscript

SCRIPT_NAME <- "make_seurat_obj_from_10x_files.R"

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(chron))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))


#' Reads a list of files into a Seurat object
#'
#' \code{load_10x_data_seurat} takes a list of files and returns a merged
#' Seurat object. Note that one should still collapse samples to a joint
#' reference space.
#'
#' @param files List.
#'     List where name is the sample name and value is a matrix dir from 10x.
#'     Each name of this list should be unique.
#' @param metadata_df Data.frame.
#'     Row names of dataframe should correspond to sample names of used to
#'     identify each file (i.e., all(names(files) %in% rownames(metadata_df)))
#' @param store_genes_as_ensembl_ids Logical.
#'     If true, genes stored as ensemble ids.
#' @param min_avg_counts Numeric.
#'     For each sample, filter features for features that have an average number
#'     of counts >min_avg_counts across cells (in that sample).
#' @param min_n_cells_exprsing_gene Numeric.
#'     For each sample, filter features for features that are expressed in
#'     >min_n_cells_exprsing_gene number of cells (in that sample).
#' @param seurat_min_cells Numeric.
#'     For each sample, include features detected in at least this many cells.
#'     Will subset the counts matrix as well. To reintroduce excluded features,
#'     create a new object with a lower cutoff.
#'     Passed to Seurat::CreateSeuratObject for each sample.
#' @param seurat_min_features Numeric.
#'     For each sample, include cells where at least this many features are
#'     detected. Passed to Seurat::CreateSeuratObject for each sample.
#'
#' @return Seurat data object.
#'     Merged Seurat data object.
#'
#' @importFrom Seurat Read10X
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat RenameCells
#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#' @export
load_10x_data_seurat <- function(
        files,
        metadata_df = NULL,
        store_genes_as_ensembl_ids = FALSE,
        min_avg_counts = 0.01,
        min_n_cells_exprsing_gene = 10,
        seurat_min_cells = 3,
        seurat_min_features = 100,
        verbose = TRUE
    ) {

    # read in the data:
    # 1. barcodes.tsv:
    #    file with the cell IDs, representing all cells quantified
    #    Barcodes are listed in the order of data presented in the
    #    matrix file (i.e. these are the column names)
    # 2. genes.tsv
    #    file with the gene IDs, representing all genes quantified
    #    The order of these genes corresponds to the order of the
    #    rows in the matrix file (i.e. these are the row names).
    # 3. matrix.mtx
    #    matrix of counts per gene for every cell
    #    The rows are associated with the gene IDs above and
    #    columns correspond to the cellular barcodes
    gene_column <- 1
    if (store_genes_as_ensembl_ids) {
        gene_column <- 2
    }
    sc_df_10x <- Seurat::Read10X(
        data.dir = files[[1]],
        gene.column = gene_column
    )

    if (min_avg_counts > 0 & min_n_cells_exprsing_gene > 0) {
        filt <- (Matrix::rowMeans(sc_df_10x) > min_avg_counts) &
            (Matrix::rowSums(sc_df_10x > 0) > min_n_cells_exprsing_gene)
        if (verbose) {
           cat(
                "[", names(files)[[1]], "]",
                "remaining features after filter on the union of ",
                "min_avg_counts and min_n_cells_exprsing_gene filters:\t",
                sum(!filt),
                "\n"
            )
        }
        sc_df_10x <- sc_df_10x[filt,]
    } else {
        if (min_avg_counts > 0) {
            filt <- (Matrix::rowMeans(sc_df_10x) > min_avg_counts)
            if (verbose) {
                cat(
                    "[", names(files)[[1]], "]",
                    "remaining features after filter on",
                    " min_avg_counts:\t",
                    sum(!filt),
                    "\n"
                )
            }
            sc_df_10x <- sc_df_10x[filt,]
        }
        if (min_n_cells_exprsing_gene > 0) {
            filt <- (Matrix::rowSums(sc_df_10x > 0) > min_n_cells_exprsing_gene)
            if (verbose) {
                cat(
                    "[", names(files)[[1]], "]",
                    "remaining features after filter on",
                    " min_n_cells_exprsing_gene:\t",
                    sum(!filt),
                    "\n"
                )
            }
            sc_df_10x <- sc_df_10x[filt,]
        }
    }

    # Initialize the Seurat object with the raw (non-normalized data)
    # * Keep all genes expressed in >= 3 cells
    # * Keep all cells with at least 100 detected genes
    #seurat_min_cells <- 3
    #seurat_min_features <- 100
    # Default metatdata
    # * orig.ident: set this to sample identity
    # * nCount_RNA: number of UMIs per cell (row)
    # * nFeature_RNA: number of genes detected per cell (row)
    sc_df <- Seurat::CreateSeuratObject(
        counts = sc_df_10x,
        min.cells = seurat_min_cells,
        min.features = seurat_min_features,
        project = names(files)[[1]], # saved to @meta.data[['orig.ident']]
        assay = "RNA"
    )

    # Rename the cellids so that we can merge
    sc_df <- Seurat::RenameCells(
        object = sc_df,
        add.cell.id = names(files)[[1]],
        for.merge = TRUE
    )


    # Add the other samples
    for (i in seq(2, length(files))){
        seurat_data <- Seurat::Read10X(data.dir = files[[i]])

        if (min_avg_counts > 0 & min_n_cells_exprsing_gene > 0) {
            filt <- (Matrix::rowMeans(seurat_data) > min_avg_counts) &
                (Matrix::rowSums(seurat_data > 0) > min_n_cells_exprsing_gene)
            if (verbose) {
                cat(
                    "[", names(files)[[i]], "]",
                    "remaining features after filter on the union of ",
                    " min_avg_counts and min_n_cells_exprsing_gene:\t",
                    sum(!filt),
                    "\n"
                )
            }
            seurat_data <- seurat_data[filt,]
        } else {
            if (min_avg_counts > 0) {
                filt <- (Matrix::rowMeans(seurat_data) > min_avg_counts)
                if (verbose) {
                    cat(
                        "[", names(files)[[i]], "]",
                        "remaining features after filter on",
                        " min_avg_counts:\t",
                        sum(!filt),
                        "\n"
                    )
                }
                seurat_data <- seurat_data[filt,]
            }
            if (min_n_cells_exprsing_gene > 0) {
                filt <- (
                    Matrix::rowSums(seurat_data > 0) > min_n_cells_exprsing_gene
                )
                if (verbose) {
                    cat(
                        "[", names(files)[[i]], "]",
                        "remaining features after filter on",
                        " min_n_cells_exprsing_gene:\t",
                        sum(!filt),
                        "\n"
                    )
                }
                seurat_data <- seurat_data[filt,]
            }
        }

        seurat_obj <- Seurat::CreateSeuratObject(
            counts = seurat_data,
            min.cells = seurat_min_cells,
            min.features = seurat_min_features,
            project = names(files)[[i]], # saved to @meta.data[['orig.ident']]
            assay = "RNA"
        )
        seurat_obj <- Seurat::RenameCells(
            object = seurat_obj,
            add.cell.id = names(files)[[i]],
            for.merge = TRUE
        )

        # Grow the sc_df... add.cell.id = append sample id to barcode
        sc_df <- merge(sc_df, seurat_obj)

        if (verbose) {
            cat(
                "[", names(files)[[i]], "]",
                "dims after merge:\t",
                dim(sc_df),
                "\n"
            )
        }
    }

    # if we have metadata information add it.
    if (!is.null(metadata_df)) {
        # not sure why this has to be vectors, but this is the only way that
        # it would work
        for (col in colnames(metadata_df)) {
            tmp_array <- metadata_df[[col]]
            names(tmp_array) <- rownames(metadata_df)
            sc_df <- AddMetaData(
                object = sc_df,
                metadata = tmp_array[sc_df@meta.data[["orig.ident"]]],
                #metadata = tmp_array[sc_df@active.ident],
                col.name = col
            )
        }

        # check to make sure the metadata in the Seurat object is correct
        check_order <- identical(
            sc_df@meta.data[["orig.ident"]],
            as.character(Seurat::Idents(sc_df))
        )
        if (!check_order) {
            stop("ERROR in adding metadata to Seurat object")
        }
        for (col in colnames(metadata_df)) {
            tmp_array <- metadata_df[[col]]
            names(tmp_array) <- rownames(metadata_df)
            tmp_df <- unique(sc_df@meta.data[,c("orig.ident", col)])
            tmp_array2 <- tmp_df[[col]]
            names(tmp_array2) <- tmp_df[["orig.ident"]]
            if (!identical(tmp_array2, tmp_array[tmp_df[["orig.ident"]]])) {
                stop("ERROR in adding metadata to Seurat object")
            }
        }
    }

    return(sc_df)
} # end load_10x_data_seurat



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
            default = "sc_dat",
            help = paste0(
                "Name (and possibly path) of output file. Will have rds.gz",
                " appended to it. ",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--verbose"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Verbose mode (write extra info to std.err).",
                " [default %default]"
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
    # load demo data from Seurat
    # UPDATE: for some reason bug here in Seurat code
    #
    # require(Seurat)
    # require(tibble)
    # require(reshape2)
    # require(dplyr)
    # if (!require(SeuratData)) {
    #     devtools::install_github('satijalab/seurat-data@v0.2.1')
    #     require(SeuratData)
    # }
    #
    # SeuratData::InstallData("panc8.SeuratData")
    # utils::data("panc8")
    # filt <- grepl("indrop", panc8@meta.data[["dataset"]])
    # df_dev <- subset(
    #     panc8,
    #     cells = colnames(panc8@assays[['RNA']]@counts)[filt]
    # )
    # df_dev@meta.data[["sample_id"]] <- df_dev@meta.data[["orig.ident"]]
    # seurat_df <- df_dev

    # load demo data from Hemberg lab
    # require(SingleCellExperiment)
    # # load Baron data
    # # https://hemberg-lab.github.io/scRNA.seq.datasets/human/pancreas/
    # df_dev <- readRDS(gzcon(url(
    #     "https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/baron-human.rds"
    # )))

    # write Seurat object to 10x format
    # DropletUtils::write10xCounts(
    #     x = seurat.object@assays$RNA@counts,
    #     path = output.path
    # )


    # load actual data from Anderson lab to play around with
    # f <- "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/leland_dev/input_files.tsv"
    # files_df <- read.csv(f, sep = "\t")
    #
    # meta_data_file <- "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/minimal_samples_spreadsheet.tsv"
    # meta_df <- read.csv(meta_data_file, sep = "\t")
    # rownames(meta_df) <- as.character(meta_df$sanger_sample_id)
    #
    # # format files for merge seurat dataframe script
    # files <- structure(
    #     as.character(files_df$path_data_10xformat),
    #     names = as.character(files_df$experiment_id)
    # )
    #
    # sc_df <- load_10x_data_seurat(
    #     files = files,
    #     metadata_df = metadata_df,
    #     min_n_cells_exprsing_gene = 0,
    #     seurat_min_cells = 0,
    #     seurat_min_features = 0
    # )
}

# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
