#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-i", "--de_file"),
                        type = "character",
                        default = "de_file",
                        help = "TSV containing DE genes"
  ),
  
  # optparse::make_option(c("-g", "--grouping_column"),
  #                       type = "character",
  #                       default = "df_key",
  #                       help = "Column to use to group the DE results"
  # ),
  
  optparse::make_option(c("-s", "--rna_enrich_script"),
                        type = "character",
                        default = "rna_enrich.r",
                        help = "Path to RNA Enrich script."
  ),
  
  optparse::make_option(c("-m", "--min_gene_to_test"),
                        type = "integer",
                        default = 10,
                        help = "Minimum number of unique gene IDs in a GO term."
  ),
  
  optparse::make_option(c("-n", "--max_gene_to_test"),
                        type = "integer",
                        default = 99999,
                        help = "Maximum number of unique gene IDs in a GO term."
  ),
  
  optparse::make_option(c("-p", "--p_val_cutoff"),
                        type = "integer",
                        default = 1,
                        help = "Entrez gene IDs in each category with 
                        p-values<sig.cutoff will be returned. Default to use
                        everything."
  ),
  
  optparse::make_option(c("-d", "--database"),
                        type = "character",
                        default = "KEGG",
                        help = "database to be tested- choices are 'GO', 
                        'KEGG', 'Cytoband','custom'"
  ),
  
  optparse::make_option(c("-c", "--min_counts"),
                        type = "integer",
                        default = 0,
                        help = "Minimum avg expression for genes to be included
                        in the analysis"
  ),
  
  optparse::make_option(c("-o", "--output_file"),
                        type = "character",
                        default = "",
                        help = "Base output name."
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
suppressPackageStartupMessages(source(arguments$options$rna_enrich_script))
suppressPackageStartupMessages(library(org.Hs.eg.db))
################################################################################

################################ Functions #####################################

retrieve_entrez_ids <- function(de_genes) {
  ensembl_entrez_map <- as.list(org.Hs.egENSEMBL2EG)
  de_genes$entrez_ids <- sapply(de_genes$gene, function(x) {
    ids <- ensembl_entrez_map[[x]] 
    ## sort and pick first one
    ## Maybe we should include both?
    ids_numeric <- sort(as.numeric(ids), decreasing=F)
    return(ids_numeric[1])
  })
  return(de_genes)
}

################################################################################


######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file_base <- arguments$options$output_file

# Read in data
de_genes <- read.csv(arguments$options$de_file,
                     sep="\t",
                     header=T)

# Need to grab Entrez IDs and filter to those with Entrez IDs
if (verbose) {
  print("Mapping Ensembl IDs to Entrez IDs...")
}
de_genes <- retrieve_entrez_ids(de_genes)

# Save mappings
gz_file <- gzfile(sprintf("%s-ensembl_entrez_mapping.tsv.gz", output_file_base),
                  "w", 
                  compression = 9)
write.table(x= de_genes[,c('gene', 'entrez_ids')],
            file = gz_file,
            sep="\t",
            col.names=T,
            row.names=F,
            quote=F)
close(gz_file)


# Run RNA Enrich
enrich_results <- rna_enrich(sigvals = de_genes$pval,
                             geneids = de_genes$entrez_ids,
                             avg_readcount = de_genes$mean,
                             species = "hsa",
                             direction = de_genes$coef,
                             min.g = arguments$options$min_gene_to_test,
                             max.g = arguments$options$max_gene_to_test,
                             sig.cutoff = arguments$options$p_val_cutoff,
                             database = arguments$options$database,
                             read_lim = arguments$options$min_counts,
                             conceptList = NULL, ## TODO: add support
                             nullsetList = NULL, ## TODO: add support
                             plot_file = sprintf('%s-rna_enrich_plot.jpg',
                                                 output_file_base),
                             plot_height = 400,
                             plot_width = 400,
                             results_file = sprintf('%s-rna_enrich_results.tsv',
                                                    output_file_base))

################################################################################



