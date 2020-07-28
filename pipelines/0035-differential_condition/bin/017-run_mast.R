#!/usr/bin/env Rscript

plot_settings <- function(plt,
                          base_size,
                          plot_marg_top,
                          text_size_axis_labels,
                          text_size_legend_labels) {
  plt <- plt +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.margin = margin(
        t = plot_marg_top, r = 45, b = 44, l = 5, unit = "pt"
      ),
      axis.title = ggplot2::element_text(size = text_size_axis_labels),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.spacing.x = grid::unit(7, 'pt'), # spacing btween key and text
      legend.title = ggplot2::element_blank(),
      legend.title.align = 0,
      legend.text = ggplot2::element_text(size = text_size_legend_labels),
      legend.justification = "center",
      legend.box.margin = ggplot2::margin( # no right margin, crank up left
        t = 0, r = 1, b = 0, l = 1, "pt"
      ),
      legend.box.spacing = grid::unit(1, "pt"), # spacing btwn plotting area & legend box
      legend.key.width = unit(2,"line"), # distance between point and label
      legend.key.height = grid::unit(1, "pt") # spacing btwn lines of legend
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 9)))
  return(plt)
}

plot_volcano_plot <- function(df,
                              fc_column,
                              p_val_column,
                              x_axis_threshold,
                              y_axis_threshold,
                              title) {
  # df$significant <- apply(df, 1, function(x) {
  #   return(abs(as.numeric(x[x_axis])) >= x_axis_threshold & abs(as.numeric(x[y_axis])) >= y_axis_threshold)
  # })
  # 
  # df$significant <- factor(df$significant, levels=c(TRUE, FALSE), labels=c("Significant", "Insignificant"))
  
  # plot <- ggplot2::ggplot(df, ggplot2::aes_string(x=x_axis, y=y_axis, color="significant")) +
  #   ggplot2::geom_point() +
  #   ggplot2::scale_x_continuous(trans='log10') +
  #   ggplot2::scale_y_continuous(trans='log10') +
  #   # ggplot2::geom_hline(yintercept = y_axis_threshold, linetype="dashed", color="red") +
  #   # ggplot2::geom_vline(xintercept = x_axis_threshold, linetype="dashed", color="red") +
  #   # ggplot2::geom_vline(xintercept = -x_axis_threshold, linetype="dashed", color="red") +
  #   ggplot2::labs(x="Log(FC)", y="Log(p-value)", color="Significance") +
  #   ggplot2::scale_color_manual(values=c("Red", "Black"))
  
  df$log10p <- -log10(df[[p_val_column]])
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x=fc_column, y="log10p")) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::labs(x="log2(FC)", y="-log10(p-value)") +
    ggplot2::theme_bw()
  return(plot)
}


plot_ma_plot <- function(df,
                         mean_expr_col,
                         fc_col,
                         x_axis_threshold,
                         y_axis_threshold,
                         title) {
  # df$significant <- apply(df, 1, function(x) {
  #   return(abs(as.numeric(x[x_axis])) >= x_axis_threshold & abs(as.numeric(x[y_axis])) >= y_axis_threshold)
  # })
  # 
  # df$significant <- factor(df$significant, levels=c(TRUE, FALSE), labels=c("Significant", "Insignificant"))
  
  # plot <- ggplot2::ggplot(df, ggplot2::aes_string(x=mean_expr_column, y=p_val_col, color="significant")) +
  #   ggplot2::geom_point() +
  #   ggplot2::geom_hline(yintercept = y_axis_threshold, linetype="dashed", color="red") +
  #   ggplot2::geom_vline(xintercept = x_axis_threshold, linetype="dashed", color="red") +
  #   ggplot2::geom_vline(xintercept = -x_axis_threshold, linetype="dashed", color="red") +
  #   ggplot2::labs(x=x_axis_title, y=y_axis_title, title=title, color="Significance") +
  #   ggplot2::scale_color_manual(values=c("Red", "Black"))
  
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x=mean_expr_column, y=fc_col)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::labs(x="log10(Mean Expression)", y="log2(FC)") 
  return(plot)
}

## Code is lifted from https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
run_MAST <- function(mtx_file_dir,
                     cell_label_column,
                     cell_label_analysed,
                     condition_col,
                     covariates_list_str,
                     de_method = "bayesglm",
                     output_file_base,
                     verbose=TRUE) {
  
  ## Load libraries
  library(MAST)
  library(Matrix)
  library(ggplot2)
  
  if (verbose) {
    print("Reading in the data...")
  }
  
  logtpm_matrix <- Matrix::readMM(sprintf("%s/log1p_cp10k/matrix.mtx.gz", mtx_file_dir))
  logtpm_matrix <- as(logtpm_matrix, "matrix")
  
  features <- read.csv(sprintf("%s/log1p_cp10k/features.tsv.gz", mtx_file_dir), sep ="\t", header=F, col.names = c("gene_id", "gene_symbol", "type"), row.names = 1)
  barcodes <- read.csv(sprintf("%s/log1p_cp10k/barcodes.tsv.gz", mtx_file_dir), sep ="\t", header=F)$V1
  metadata <- read.csv(sprintf("%s/log1p_cp10k/cell_metadata.tsv.gz", mtx_file_dir), sep ="\t", header=T, row.names=1)
  metadata <- metadata[barcodes, , drop=F]
  
  # Format metadata as correct type -- R should already do this, but double check
  metadata_type <- read.csv(sprintf("%s/log1p_cp10k/covariate_type.tsv.gz", mtx_file_dir), sep ="\t", header=T)
  for (i in seq(1, nrow(metadata_type))) {
    var <- metadata_type[i,'variable']
    type <- metadata_type[i,'type']
    
    if (type == 'categorical' && !is.character(metadata[[var]])) {
      if (verbose) {
        print(sprintf("The categorical variable '%s' is not categorical. Attempting to force categorical.", var))
      }
      metadata[var] <- as.character(metadata[[var]])
    } else if (type == 'continuous' && is.character(metadata[[var]])) { ## only check if it's character
      if (verbose) {
        print(sprintf("The continuous variable '%s' is not numeric. Attempting to force numeric.", var))
      }
      metadata[var] <- as.numeric(metadata[[var]])
    }
  }
  
  # Set dimnames
  rownames(logtpm_matrix) <- rownames(features)
  colnames(logtpm_matrix) <- rownames(metadata)
  
  sca <- MAST::FromMatrix(exprsArray=logtpm_matrix, cData = metadata, fData=features)
  
  ## Get formula
  covariates <- strsplit(x=covariates_list_str, split=",", fixed=TRUE)[[1]]
  formula_str <- sprintf("~ %s + %s", condition_col, paste(covariates, collapse = " + "))
  formula <- formula(formula_str)
  
  if (verbose) {
    print(sprintf("Calculating differential expression using the formula: %s", formula_str))
    print(sprintf("Calculating differential expression using the method: %s", de_method))
    print(sprintf("Randomly selecting value in the `%s` column to use as base condition in fitting the model.", condition_col))
  }
  
  condition_data <- factor(colData(sca)[[condition_col]])
  reference_condition <- levels(condition_data)[1]
  condition_data <- relevel(condition_data, reference_condition)
  colData(sca)[condition_col] <- condition_data
  
  if (verbose) {
    print(sprintf("Selected the value '%s' as the reference condition.", reference_condition))
    print("Fitting the model...")
  }
  
  zlm_fit <- MAST::zlm(formula, sca, method = de_method, silent=FALSE, parallel = TRUE)
  
  if (verbose) {
    print("Done fitting the model.")
  }
  
  test_var <- paste(c(condition_col, setdiff(levels(condition_data), reference_condition)), collapse = "")
  if (verbose) {
    print(sprintf("Performing Log-Ratio Test for the variable %s...", test_var))
  }
  results <- summary(zlm_fit, doLRT=test_var)
  results_dt <- results$datatable
  
  # Structure results into interpretable format
  fcHurdle <- merge(results_dt[contrast==test_var & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    results_dt[contrast==test_var & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle <- merge(fcHurdle, as.data.frame(mcols(sca)), by='primerid')
  fcHurdle$de_software <- "MAST"
  fcHurdle$de_method <- de_method
  fcHurdle$condition <- condition_col
  fcHurdle$reference_condition <- reference_condition
  fcHurdle$test_condition <- setdiff(levels(condition_data), reference_condition)
  fcHurdle$covariates <- covariates_list_str
  fcHurdle$cell_label_column <- cell_label_column
  fcHurdle$cell_label_analysed <- cell_label_analysed
  
  ## Add mean expression from counts data
  counts_matrix <- Matrix::readMM(sprintf("%s/counts/matrix.mtx.gz", mtx_file_dir))
  counts_matrix <- as(counts_matrix, "matrix")
  counts_features <- read.csv(sprintf("%s/counts/features.tsv.gz", mtx_file_dir), sep ="\t", header=F, col.names = c("gene_id", "gene_symbol", "type"), row.names = 1)
  counts_barcodes <- read.csv(sprintf("%s/counts/barcodes.tsv.gz", mtx_file_dir), sep ="\t", header=F)$V1
  rownames(counts_matrix) <- rownames(counts_features)
  colnames(counts_matrix) <- rownames(counts_barcodes)
  fcHurdle$meanExpr <- Matrix::rowMeans(counts_matrix[fcHurdle$primerid,])
  
  ## Order by FDR
  fcHurdle <- fcHurdle[order(fcHurdle$fdr, decreasing=F)]
  
  vol_plot <- plot_volcano_plot(fcHurdle,
                                'coef',
                                'Pr(>Chisq)',
                                0,
                                0,
                                "Dud")
  vol_plot <- plot_settings(vol_plot,  46, 20, 30, 30)
  png(file = sprintf("%s-plot_volcano.png", output_file_base),
      height = cm(10.5),
      width = cm(13.5)
  )
  print(vol_plot)
  dev.off()
  
  gz_file <- gzfile(sprintf("%s-de_results.tsv.gz", output_file_base), "w", compression = 9)
  write.table(x= fcHurdle,
              file = gz_file,
              sep="\t",
              col.names=T,
              row.names=F,
              quote=F)
  close(gz_file)
  
  if (verbose) {
    print("Done.")
  }
}

# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # disable strings as factors, but re-enable upon exit
  ## Need to set the global allowance higher for the amount of data
  ## Expect cores at end 
  old <- options(stringsAsFactors = FALSE, future.globals.maxSize= 5368709120, mc.cores=as.numeric(args[8]))
  on.exit(options(old), add = TRUE)
  
  run_MAST(args[1], args[2], args[3], args[4], args[5], args[6], args[7])
}
