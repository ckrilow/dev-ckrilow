#!/usr/bin/env Rscript

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Seurat"))


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
#' @importFrom data.table fread
#' @importFrom Seurat Read10X
#' @importFrom Seurat CreateSeuratObject
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("--file_paths_10x"),
            type = "character",
            help = paste0(
                "Input file_paths_10x.tsv (all samples in the data freeze)."
            )
        ),

        optparse::make_option(c("--cellranger_irods_locations"),
            type = "character",
            help = paste0(
              "Locations of each file in irods (csv)."
            )
        ),

        optparse::make_option(c("--metadata"),
            type = "character",
            help = paste0(
              "Metadata file."
            )
        ),

        optparse::make_option(c("--metadata_dict"),
            type = "character",
            help = paste0(
              "Metadata definitions."
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Makes data for OTAR partner sharing."
        )
    )

    # A hack to fix a bug in optparse that won't let you use positional args
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

    # Read in the sample sheet
    df_samples <- read.csv(param[["file_paths_10x"]], sep = "\t")

    # Read in cellranger_irods_locations
    df_cellranger <- read.csv(
        param[["cellranger_irods_locations"]],
        sep = ","
    )
    rownames(df_cellranger) <- df_cellranger$sanger_sample_id

    # Make output file with irods data location of the files to release.
    df_cellranger_out <- df_cellranger[df_samples$experiment_id, ]
    file_or_dir_endings_irods <- c(
        "cellranger_filtered" = "filtered_feature_bc_matrix/",
        "cellranger_unfiltered" = "raw_feature_bc_matrix/",
        "bam_file" = "possorted_genome_bam.bam",
        "bam_index" = "possorted_genome_bam.bai"

    )
    for (i in names(file_or_dir_endings_irods)) {
        df_cellranger_out[[i]] <- paste0(
            df_cellranger_out$cellranger_irods_path,
            "/",
            file_or_dir_endings_irods[[i]]
        )
    }
    df_cellranger_out$cellranger_irods_path <- NULL
    write.table(
        df_cellranger_out,
        "data_release-irods_locations.tsv",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t",
        na = ""
    )



    # Read in the sample metadata
    df_metadata <- read.csv(
        param[["metadata"]],
        sep = "\t"
    )
    rownames(df_metadata) <- df_metadata$sanger_sample_id

    # Read in the sample metatadictionary
    df_metadata_dict <- read.csv(
        param[["metadata_dict"]],
        sep = "\t"
    )

    # Subset the sample metadata
    # EDIT: instead of subsetting here, just keep the
    # cols_sample_metadata <- c(
    #     "date_of_sample",
    #     "date_of_plate_submission",
    #     "patient_id",
    #     "sanger_sample_id",
    #     "sex",
    #     "age",
    #     "biopsy_type",
    #     "disease_status",
    #     "inflammation_status",
    #     "cd_type",
    #     "time_to_chromium_processing",
    #     "total_reads"
    #     # smoking_status
    #     # experimentalist
    #     # epithelial_immune_ratio
    #     # protocol
    #     # collagenase
    #     # inhibitor
    #     # DNase
    #     # ack_lysis_buffer
    #     # frozen_processed
    #     # chip_well_position
    #     # enzyme_lot
    #     # bead_version
    #     # bead_lot
    #     # chip_version
    #     # chip_lot
    #     # id_run
    #     # lane
    #     # library_id
    #     # library_type
    #     # study
    #     # study_id
    # )
    cols_sample_metadata <- df_metadata_dict$column_name
    df_metadata <- df_metadata[df_samples$experiment_id, cols_sample_metadata]
    write.table(
        df_metadata,
        "data_release-metadata.tsv",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )

    write.table(
        df_metadata_dict,
        "data_release-metadata_column_description.tsv",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
}


main <- function() {
    command_line_interface()
}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    main()
}
