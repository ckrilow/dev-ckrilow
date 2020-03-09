library("tidyverse")

dat <- read.csv(
    "/Users/lt9/repo/scrna_cellranger/sync_status/samples_metainfo.tsv",
    sep = "\t"
)

dat <- subset(dat, sample_status == "Sequenced + Cell Ranger")
dat <- subset(dat, biopsy_type_original %in% c("ti", "neoti"))
#dat <- subset(dat, disease_status_original %in% c("cd"))
dat <- subset(dat, protocol %in% c("tissue_v2", "blood_final"))

table(dat[c("biopsy_type", "disease_status")])


dat$experiment_id <- dat$sanger_sample_id
dat$file_path <- paste0(
    "/home/ubuntu/data/scrna_cellranger/results/iget_cellranger/full_data/",
    dat$biopsy_type_original,
    "/",
    dat$disease_status_original,
    "/",
    dat$sanger_sample_id,
    "/filtered_feature_bc_matrix"
)


dat <- dat %>% dplyr::select(experiment_id, sanger_sample_id, file_path)

write.table(
    dat,
    "input_sample_files.tsv",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    sep = "\t",
    na = ""
)
