#Example usage

#Rscript save_input_files.R --data_path "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data" --meta_path "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/sync_status/samples_metainfo.tsv" --disease_status "All" --tissue "All"

#Rscript save_input_files.R --data_path "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data" --meta_path "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/sync_status/samples_metainfo.tsv" --disease_status "All" --tissue "TI"

#Rscript save_input_files.R --data_path "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data" --meta_path "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/sync_status/samples_metainfo.tsv" --disease_status "All" --tissue "Rectum"


library("optparse")


option_list <- list(
  make_option(c("--data_path"), type = "character"),
  make_option(c("--meta_path"), type = "character"),
  make_option(c("--output"),type = "character"),
  make_option(c("--disease_status"),type = "character", default="All"),
  make_option(c("--tissue"),type = "character", default="All" )
)

opt=parse_args(OptionParser(option_list=option_list))

  ########## ########## LOAD DATA

    info=read.csv(file = opt$meta_path, sep = "\t", stringsAsFactors = FALSE, colClasses=)
    col.times=read.csv(file = opt$meta_path, sep = "\t")[ ,c('chromium_time', 'collection_time')]
    info$time_to_chromium_processing <- chron::times(col.times[["chromium_time"]]) - chron::times(col.times[["collection_time"]])
    info=info[which(info$sample_status=="Sequenced + Cell Ranger"),]
    info=info[which(info$protocol=="tissue_v2"),]

    if(opt$tissue=="TI" | opt$tissue=="Rectum"){
      info=info[info$biopsy_type==opt$tissue,]
    } else {

    }

    if(opt$disease_status=="Crohn's disease" | opt$disease_status=="Healthy"){
      info=info[info$disease_status==opt$disease_status,]
    } else {

    }



    if(length(grep("Pilot", info$sanger_sample_id))>0){info=info[-grep("Pilot", info$sanger_sample_id),]}
    paths=list.files(opt$data_path, pattern="matrix.mtx.gz", recursive = TRUE, full.names = TRUE)
    paths=paths[grep("filtered_feature_bc_matrix", paths)]
    paths=sub("/matrix.mtx.gz", "", paths)

    df=data.frame("experiment_id"=sapply(strsplit(paths,"/"), "[[", length(strsplit(paths,"/")[[1]])-1), "data_path_10x_format"=paths)
    df=df[df$experiment_id %in% info$sanger_sample_id,]
    df=df[match(info$sanger_sample_id, df$experiment_id),]

    df$short_experiment_id <- gsub("_study_of_dissociation_methods_for_human_gut_tissues", "", df$experiment_id)
    df$short_experiment_id <- gsub("_Disease_Collection_Study", "", df$experiment_id)


    write.table(df, file='file_paths_10x.tsv', quote=FALSE, sep='\t', row.names=FALSE)
    #write.table(info, file='file_metadata.tsv', quote=FALSE, sep='\t', row.names=FALSE)
