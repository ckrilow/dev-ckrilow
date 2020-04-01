
  data.path="/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data"
  meta.path.file="/nfs/users/nfs_m/mk32/repo/scrna_cellranger/sync_status/samples_metainfo.tsv"
  output.path="/nfs/users/nfs_m/mk32/repo/sc_nextflow/studies/gut-freeze002"
  
  
  ########## ########## LOAD DATA
  status=c(0, "Healthy", "Crohn's disease")
  
  for(i in 1:length(status)){
    
    
    info=read.csv(file = meta.path.file, sep = "\t", stringsAsFactors = FALSE, colClasses=)
    col.times=read.csv(file = meta.path.file, sep = "\t")[ ,c('chromium_time', 'collection_time')]
    info$time_to_chromium_processing <- chron::times(col.times[["chromium_time"]]) - chron::times(col.times[["collection_time"]])
    info=info[which(info$sample_status=="Sequenced + Cell Ranger"),]
    info=info[which(info$protocol=="tissue_v2"),]
    info=info[info$biopsy_type=="TI",]
    if(status[i]>0){info=info[info$disease_status==status[i],]} 
    if(length(grep("Pilot", info$sanger_sample_id))>0){info=info[-grep("Pilot", info$sanger_sample_id),]}
    paths=list.files(data.path, pattern="matrix.mtx.gz", recursive = TRUE, full.names = TRUE)
    paths=paths[grep("filtered_feature_bc_matrix", paths)]
  
    
    df=data.frame("experiment_id"=sapply(strsplit(paths,"/"), "[[", length(strsplit(paths,"/")[[1]])-2), "data_path_10x_format"=paths)
    df=df[df$experiment_id %in% info$sanger_sample_id,]
    df=df[match(info$sanger_sample_id, df$experiment_id),] 

    df$short_experiment_id <- gsub("_study_of_dissociation_methods_for_human_gut_tissues", "", df$experiment_id)
    df$short_experiment_id <- gsub("_Disease_Collection_Study", "", df$experiment_id) 

    if (status[i]==0) {
      subdir="ti-cd_healthy"
    } else if (status[i]=="Healthy") {
      subdir="ti-healthy"
    } else {
      subdir="ti-cd"
    }
    
    
    out.path=file.path(output.path, subdir)
    system(paste0("mkdir -p ", out.path))
    
    write.table(df, file=file.path(out.path,'file_paths_10x.tsv'), quote=FALSE, sep='\t', row.names=FALSE)
    write.table(info, file=file.path(out.path,'file_metadata.tsv'), quote=FALSE, sep='\t', row.names=FALSE)

    
  }


  
  
  
  
  
  