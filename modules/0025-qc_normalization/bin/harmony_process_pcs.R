#!/usr/bin/env Rscript

SCRIPT_NAME <- "harmony_process_pcs.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("harmony"))


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
#' @param verbose Logical.
#'     Write extra output to std.out. Default = TRUE.
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
        optparse::make_option(c("-f", "--pca_file"),
            type = "character",
            help = paste0(
                "Tab-delimited input file of PCs. Columns: 'cell_barcode'",
                " followed by PC1,PC2,etc."
            )
        ),

        optparse::make_option(c("--metadata_file"),
            type = "character",
            help = paste0(
              "Tab-delimited metadata file, must have a column labelled",
              " cell_barcode that maps to pca_file."
            )
        ),
        
        optparse::make_option(c("--metadata_columns"),
            type = "character",
            help = paste0(
              "Comma separated string of columns to use in metadata_file."
            )
        ),
        
        optparse::make_option(c("--theta"),
            type = "character",
            default = "",
            help = paste0(
                "Comma separated string of theta values (corresponding to",
                " metadata_columns). If '' then sets theta to 2 for all", 
                " columns. Larger values of theta result in more diverse", 
                " clusters.",
                " [default: %default]"
            )
        ),
        
        optparse::make_option(c("--n_pcs"),
            type = "numeric",
            default = 0,
            help = paste0(
                "Number of PCs to use from pca_file. If 0 then use all.",
                " [default: %default]"
            )
        ),
        
        
        optparse::make_option(c("--out_file"),
            type = "character",
            default = "",
            help = paste0(
                "Name (and possibly path) of output file. Will have tsv.gz",
                " appended to it. If '' then add '-harmony.tsv.gz' to",
                " pca_file.",
                " [default: %default]"
            )
        )
        
        # optparse::make_option(c("--verbose"),
        #     type = "logical",
        #     action = "store_true",
        #     default = FALSE,
        #     help = paste0(
        #         "Verbose mode (write extra info to std.err).",
        #         " [default: %default]"
        #     )
        # )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Runs Harmony using a PC and metadata file from scanpy."
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
    
    # Read in the PCA file.
    f_pca <- param[["pca_file"]]
    mtx_pca <- data.table::fread(
        cmd = paste("gunzip -c", f_pca), 
        sep = "\t", 
        header = TRUE, 
        stringsAsFactors = FALSE
    )
    if (!("cell_barcode" %in% colnames(mtx_pca))) {
        stop("Missing cell_barcode column.")
    }
    mtx_pca <- as.matrix(mtx_pca, rownames = "cell_barcode")
    
    # Check the input pcs. 
    n_pcs <- param[["n_pcs"]]
    if (n_pcs > length(colnames(mtx_pca))) {
        stop("Invalid n_pcs value.")
    } else if (n_pcs < 0) {
        stop("Invalid n_pcs vaue.")
    } else if (n_pcs == 0) { 
        n_pcs <- length(colnames(mtx_pca))
    }
 
    
    # Get the metadata_file columns that we want to adjust with Harmony.
    metadata_columns <- unlist(strsplit(param[["metadata_columns"]], ","))
    
    
    # Read in the metadata file.
    f_meta <- param[["metadata_file"]]
    df_meta <- data.frame(data.table::fread(
        cmd = paste("gunzip -c", f_meta), 
        sep = "\t", 
        header = TRUE, 
        stringsAsFactors = FALSE
    ))
    if (!("cell_barcode" %in% colnames(df_meta))) {
        stop("Missing cell_barcode column.")
    }
    rownames(df_meta) <- df_meta[["cell_barcode"]]
    df_meta <- df_meta[rownames(mtx_pca), metadata_columns]
    
    
    # Get the theta values for each column (if none, set to 2 for all columns).
    theta <- rep(2, length(metadata_columns))
    if (param[["theta"]] != '') {
        theta <- unlist(strsplit(param[["theta"]]))
    }
    
    
    harmony_embeddings <- harmony::HarmonyMatrix(
        data_mat = mtx_pca[, seq(1, n_pcs)],
        meta_data = df_meta,
        do_pca = FALSE,
        #npcs = n_pcs, # only relevant if do_pca == TRUE
        verbose = TRUE,
        # epsilon.harmony = -Inf, # Set to -Inf to never stop early.
        vars_use = metadata_columns,
        theta = theta
    )
    colnames(harmony_embeddings) <- gsub(
        "PC", 
        "harmony",
        colnames(harmony_embeddings)
    )
    out_col_order <- colnames(harmony_embeddings)
    out_col_order <- c("cell_barcode", out_col_order)
    harmony_embeddings <- as.data.frame(harmony_embeddings)
    harmony_embeddings[["cell_barcode"]] <- rownames(harmony_embeddings)
    
    
    base <- param[["out_file"]]
    if (base == "") {
        base <- paste0(
            gsub(".tsv.gz", "", param[["pca_file"]]),
            "-harmony"
        )
    }
    gzfh <- gzfile(paste(base, ".tsv.gz", sep = ""), "w", compression = 9)
    write.table(
        harmony_embeddings[out_col_order], 
        gzfh,
        row.names = FALSE, 
        col.names = TRUE,
        quote = FALSE, 
        sep = "\t", 
        na = ""
    )
    close(gzfh)

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
    f_pca <- "adata-pcs.tsv.gz"
    mtx_pca <- data.table::fread(
        cmd = paste("gunzip -c", f_pca), 
        sep = "\t", 
        header = TRUE, 
        stringsAsFactors = FALSE
    )
    mtx_pca <- as.matrix(mtx_pca, rownames = "cell_barcode")
    
    f_meta <- "adata-metadata.tsv.gz"
    df_meta <- data.frame(data.table::fread(
        cmd = paste("gunzip -c", f_meta), 
        sep = "\t", 
        header = TRUE, 
        stringsAsFactors = FALSE
    ))
    rownames(df_meta) <- df_meta[["cell_barcode"]]
    df_meta <- df_meta[rownames(mtx_pca), c("sanger_sample_id", "bead_version")]
    
    n_pcs <- 10
    harmony_embeddings <- harmony::HarmonyMatrix(
        data_mat = mtx_pca[, seq(1, n_pcs)],
        meta_data = df_meta,
        do_pca = FALSE,
        #npcs = n_pcs,
        verbose = TRUE,
        # epsilon.harmony = -Inf, # Set to -Inf to never stop early.
        # vars_use = c("sanger_sample_id"),
        # theta = c(1)
        vars_use = c("sanger_sample_id", "bead_version"),
        theta = c(1, 0.2)
    )
    colnames(harmony_embeddings) <- gsub(
        "PC", 
        "harmony",
        colnames(harmony_embeddings)
    )
    out_col_order <- colnames(harmony_embeddings)
    out_col_order <- c("cell_barcode", out_col_order)
    harmony_embeddings <- as.data.frame(harmony_embeddings)
    harmony_embeddings[["cell_barcode"]] <- rownames(harmony_embeddings)
    
    base <- "harmony_dev"
    gzfh = gzfile(paste(base, ".tsv.gz", sep = ""), "w", compression = 9)
    write.table(
        harmony_embeddings[out_col_order], 
        gzfh,
        row.names = FALSE, 
        col.names = TRUE,
        quote = TRUE, 
        sep = "\t", 
        na = ""
    )
    close(gzfh)
}

# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
