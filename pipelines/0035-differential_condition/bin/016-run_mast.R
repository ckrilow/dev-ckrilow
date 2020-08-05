#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-i", "--mtx_dir"),
                        type = "character",
                        default = "matrix_dir",
                        help = "Directory containing input files for MAST"
  ),

  optparse::make_option(c("-l", "--cell_label_column"),
                        type = "character",
                        default = "cluster",
                        help = "Column to use for cell label."
  ),

  optparse::make_option(c("-a", "--cell_label_analysed"),
                        type = "character",
                        default = "",
                        help = "Cell label."
  ),

  optparse::make_option(c("-c", "--condition_column"),
                        type = "character",
                        default = "condition",
                        help = "Column to use as condition."
  ),

  optparse::make_option(c("-d", "--covariate_columns_discrete"),
                        type = "character",
                        default = "",
                        help = "Discrete covariates to include in the model."
  ),

  optparse::make_option(c("-z", "--covariate_columns_continuous"),
                        type = "character",
                        default = "",
                        help = "Continuous covariates to include in the model."
  ),

  optparse::make_option(c("-m", "--method"),
                        type = "character",
                        default = "bayesglm",
                        help = "MAST method to use to model DE."
  ),

  optparse::make_option(c("-o", "--out_file"),
                        type = "character",
                        default = "",
                        help = "Base output name."
  ),

  optparse::make_option(c("-n", "--cores_available"),
                        type = "integer",
                        default = 1,
                        help = "Number of cores to use."
  ),

  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = TRUE,
                        help = ""
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Calculates differentially expressed genes using MAST."
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
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
################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
################################################################################

################################ Functions #####################################

cast_covariates <- function(df,
                            cols,
                            cast_func,
                            cast_func_description,
                            verbose) {
  if (verbose) {
    print(sprintf("Casting columns to be %s...", cast_func_description))
  }
  for (col in cols) {
    df[col] <- cast_func(df[[col]])
  }
  return(df)
}

## Code is lifted from
## https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
run_MAST <- function(matrix,
                     cell_data,
                     feature_data,
                     condition_col,
                     covariates,
                     de_method = "bayesglm",
                     verbose=TRUE) {
  # Format data into sca object
  sca <- MAST::FromMatrix(exprsArray = matrix,
                          cData = cell_data,
                          fData = feature_data)

  ## Get formula
  formula_str <- sprintf("~ %s", condition_col)
  if (length(covariates) != 0) {
    formula_str <- sprintf("%s + %s",
                           formula_str,
                           paste(covariates, collapse = " + "))
  }
  formula <- formula(formula_str)

  if (verbose) {
    print(sprintf("Calculating differential expression using the formula: %s",
                  formula_str))
    print(sprintf("Calculating differential expression using the method: %s",
                  de_method))
    print(sprintf(paste("Randomly selecting value in the `%s` column to",
                        "use as base condition in fitting the model."),
                  condition_col))
  }

  # Select a condition to use as the 'reference' condition. Can pass this
  # in through parameters later.
  condition_data <- factor(SummarizedExperiment::colData(sca)[[condition_col]])
  reference_condition <- levels(condition_data)[1]
  condition_data <- relevel(condition_data, reference_condition)
  SummarizedExperiment::colData(sca)[condition_col] <- condition_data

  if (verbose) {
    print(sprintf("Selected the value '%s' as the reference condition.",
                  reference_condition))
    print("Fitting the model...")
  }

  zlm_fit <- MAST::zlm(formula,
                       sca,
                       method = de_method,
                       silent = FALSE,
                       parallel = TRUE)
  if (verbose) {
    print("Done fitting the model.")
  }

  test_var <- paste(c(condition_col,
                      setdiff(levels(condition_data), reference_condition)),
                    collapse = "")
  if (verbose) {
    print(sprintf("Performing Log-Ratio Test for the variable %s...",
                  test_var))
  }
  results <- summary(zlm_fit, doLRT = test_var)
  results_dt <- results$datatable

  # Structure results into interpretable format
  fcHurdle <- merge(results_dt[contrast == test_var & component =='H',
                               .(primerid, `Pr(>Chisq)`)], #hurdle P values
                    results_dt[contrast == test_var & component =='logFC',
                               .(primerid, coef, ci.hi, ci.lo)],
                    by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle <- merge(fcHurdle, as.data.frame(mcols(sca)), by='primerid')

  # Add important metadata
  fcHurdle$de_method <- sprintf('mast-%s', de_method)
  fcHurdle$condition <- condition_col
  fcHurdle$reference_condition <- reference_condition
  fcHurdle$coef_condition <- setdiff(levels(condition_data),
                                     reference_condition)

  ## Order by FDR
  fcHurdle <- fcHurdle[order(fcHurdle$fdr, decreasing = F)]
  return(fcHurdle)
}

plot_volcano_plot <- function(df,
                              fc_col,
                              p_val_col,
                              fc_threshold = 1.8,
                              p_val_threshold = 0.05) {
  fc_threshold <- log2(fc_threshold)
  df <- df[!(is.na(df[[fc_col]]))  & !(is.na(df[[p_val_col]])),]
  df$significant <- apply(df, 1, function(x) {
    return(abs(as.numeric(x[[fc_col]])) >= fc_threshold &
             abs(as.numeric(x[[p_val_col]])) <= p_val_threshold)
  })
  df$significant <- factor(df$significant,
                           levels = c(TRUE, FALSE),
                           labels = c("True", "False"))

  df$neg_log10 <- -log10(df[[p_val_col]])
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = fc_col,
                                                  y = "neg_log10",
                                                  color = "significant")) +
    ggplot2::geom_point(size = .5) +
    ggplot2::labs(x = "log2(FC)",
                  y = "-log10(p-value)",
                  color = sprintf("Log2(FC)>=%s and\nFDR<=%s",
                                  round(fc_threshold, digits = 2),
                                  p_val_threshold)) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("#CF9400", "Black"))
  return(plot)
}


plot_ma_plot <- function(df,
                         mean_expr_col,
                         fc_col,
                         mean_expr_threshold = 0.0001,
                         fc_threshold = 1.8) {
  fc_threshold <- log2(fc_threshold)
  df <- df[!(is.na(df[[mean_expr_col]]))  & !(is.na(df[[fc_col]])),]
  df$significant <- apply(df, 1, function(x) {
    return(abs(as.numeric(x[[mean_expr_col]])) >= mean_expr_threshold &
             abs(as.numeric(x[[fc_col]])) >= fc_threshold)
  })
  df$significant <- factor(df$significant,
                           levels = c(TRUE, FALSE),
                           labels = c("True", "False"))

  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x=mean_expr_col,
                                                  y=fc_col,
                                                  color="significant")) +
    ggplot2::geom_point(size=.5) +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Mean Expression",
                  y="log2(FC)",
                  color=sprintf("Mean Expr>=%s and\nLog2(FC)>=%s",
                                mean_expr_threshold,
                                round(fc_threshold, digits=2))) +
    ggplot2::scale_color_manual(values=c("#CF9400", "Black"))
  return(plot)
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file_base <- arguments$options$out_file

# Re-set options to allow multicore
old <- options(stringsAsFactors = FALSE,
               mc.cores=arguments$options$cores_available)
on.exit(options(old), add = TRUE)

# Read data in, cast to correct data types
if (verbose) {
  print("Reading in the data...")
}
mtx_file_dir = arguments$options$mtx_dir
logtpm_matrix <- as(Matrix::readMM(sprintf("%s/log1p_cp10k/matrix.mtx.gz",
                                           mtx_file_dir)),
                    "matrix")
logtpm_features <- read.csv(sprintf("%s/log1p_cp10k/features.tsv.gz",
                                    mtx_file_dir),
                            sep ="\t",
                            header=F,
                            col.names = c("gene_id", "gene_symbol", "type"),
                            row.names = 1)
logtpm_barcodes <- read.csv(sprintf("%s/log1p_cp10k/barcodes.tsv.gz",
                                    mtx_file_dir),
                            sep ="\t",
                            header=F)$V1
logtpm_metadata <- read.csv(sprintf("%s/log1p_cp10k/cell_metadata.tsv.gz",
                                    mtx_file_dir),
                            sep ="\t",
                            header=T,
                            row.names=1)
logtpm_metadata <- logtpm_metadata[logtpm_barcodes, , drop=F]
rownames(logtpm_matrix) <- rownames(logtpm_features)
colnames(logtpm_matrix) <- rownames(logtpm_metadata)

# Variable casting
discrete_covs <- strsplit(x=arguments$options$covariate_columns_discrete,
                          split=",",
                          fixed=TRUE)[[1]]
logtpm_metadata <- cast_covariates(logtpm_metadata, discrete_covs,
                                   as.character,
                                   "characters",
                                   verbose)
continuous_covs <- strsplit(x=arguments$options$covariate_columns_continuous,
                            split=",",
                            fixed=TRUE)[[1]]
logtpm_metadata <- cast_covariates(logtpm_metadata,
                                   continuous_covs,
                                   as.numeric,
                                   "numeric",
                                   verbose)

## Filter all covariates with a single value
## MAST throws an error if not
covariates_passed <- c(discrete_covs, continuous_covs)
if (length(covariates_passed) > 0) {
    covariates <- covariates_passed[sapply(covariates_passed, function(x) {
      remove <- length(unique(logtpm_metadata[[x]])) <= 1
      if (verbose && remove) {
        print(sprintf(
            "Covariate `%s` only has one value, removing from MAST list.",
            x
        ))
      }
      return(!remove)
    })]
} else {
    covariates <- covariates_passed
}

# Run MAST
de_results <- run_MAST(matrix = logtpm_matrix,
                       cell_data = logtpm_metadata,
                       feature_data = logtpm_features,
                       condition_col = arguments$options$condition_column,
                       covariates = covariates,
                       de_method = arguments$options$method,
                       verbose = verbose)

# Fit dataframe to match diffxpy
de_results$covariates_passed <- paste(covariates_passed, collapse = ",")
de_results$covariates <- paste(covariates, collapse = ",")
de_results$cell_label_column <- arguments$options$cell_label_column
de_results$cell_label_analysed <- arguments$options$cell_label_analysed
names(de_results)[names(de_results) == "primerid"] <- "gene"
names(de_results)[names(de_results) == "Pr(>Chisq)"] <- "pval"
de_results$log2fc <- de_results$coef

# Add mean expression from counts data
counts_matrix <- as(Matrix::readMM(sprintf("%s/counts/matrix.mtx.gz",
                                           mtx_file_dir)),
                    "matrix")
counts_features <- read.csv(sprintf("%s/counts/features.tsv.gz", mtx_file_dir),
                            sep ="\t",
                            header=F,
                            col.names = c("gene_id", "gene_symbol", "type"),
                            row.names = 1)
counts_barcodes <- read.csv(sprintf("%s/counts/barcodes.tsv.gz", mtx_file_dir),
                            sep ="\t",
                            header=F)$V1
rownames(counts_matrix) <- rownames(counts_features)
colnames(counts_matrix) <- rownames(counts_barcodes)
de_results$mean <- Matrix::rowMeans(counts_matrix[de_results$gene,])

## Save result
if (verbose) {
  print("Wriiting DE results...")
}
gz_file <- gzfile(sprintf("%s-de_results.tsv.gz", output_file_base),
                  "w",
                  compression = 9)
write.table(x= de_results,
            file = gz_file,
            sep="\t",
            col.names=T,
            row.names=F,
            quote=F)
close(gz_file)

## Plot results
vol_plot <- plot_volcano_plot(de_results,
                              'log2fc',
                              'fdr',
                              1.8,
                              0.05)
png(file = sprintf("%s-plot_volcano.png", output_file_base))
print(vol_plot)
dev.off()

ma_plot <- plot_ma_plot(de_results,
                        'mean',
                        'log2fc',
                        0.0001,
                        1.8)
png(file = sprintf("%s-plot_ma.png", output_file_base))
print(ma_plot)
dev.off()

if (verbose) {
  print("Done.")
}

################################################################################
