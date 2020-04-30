#!/usr/bin/env Rscript

library(dmm)
library(Seurat)
library(ggplot2)
library(hdf5r)
library(optparse)

optionList <- list(
     optparse::make_option(c("-od", "--output_dir"),
                           type = "character",
                           help = paste0("Basename of output files,
                           assuming output in current working  directory."
                           )
     ),

     optparse::make_option(c("-of", "--output_file"),
                           type = "character",
                           help = paste0("Basename of output files,
                           assuming output in current working  directory."
                           )
     ),

    optparse::make_option(c("-m", "--maximum_de"),
                          type = "numeric",
                          help = paste0(
                            ""
                          )
    ),

    optparse::make_option(c("-a", "--auc_difference"),
                          type="float",
                          help = paste0(
                            ""
                          )
    )
    )

opt = parse_args(OptionParser(option_list=optionList))

output.path=opt$od
basename=opt$of
max_de=opt$m #Maximum differentially expressed genes
auc_diff=opt$a #Difference for AUC, truncate de genes
              #with: AUC < 0 + auc_diff and AUC > 1 - auc_diff


# Data
X=file.path(output.path, "_X.csv")
X <- H5File$new(file.path(output.path, "_X.h5"), mode="r")
M=X[['df']][1:X[['df']]$dims[1],1:X[['df']]$dims[2]]
cells=X[['cells']][1:X[['cells']]$dims]
genes=X[['genes']][1:X[['genes']]$dims]
df=as.data.frame(M)
colnames(df) <- cells
rownames(df) <- genes

# Metadata
obs.file=file.path(output.path, "obs.csv")
obs=read.csv(file = obs.file, sep = "\t",
                              stringsAsFactors = FALSE,
                              row.names = 1)

seur <- CreateSeuratObject(counts=df, assay = "RNA", meta.data = obs)
Idents(object=seur) <- factor(seur$leiden)

# Initialize

clusters_active=unfactor(seur@active.ident)
clusters_unique=sort(unique(unfactor(seur@active.ident)))

update_min=0
k=1

merging_progress=list()
merging_matrix=list()

merging_progress[[k]] <- clusters_active


while(update_min<max_de) {

  mmm=matrix(data=NA, nrow=length(clusters_unique),
                      ncol=length(clusters_unique))

  for(j in 1:length(clusters_unique)){
    for(i in 1:length(clusters_unique)){

      if(clusters_unique[j]==clusters_unique[i]){

        mmm[j,i]=NA

      } else {

        Y = FindMarkers(seur,
                        ident.1 = clusters_unique[j],
                        ident.2 = clusters_unique[i],
                        min.pct = 0.25,
                        random.seed = 1,
                        test.use = "roc")

        # matrix with DE genes between each pair of clusters
        mmm[j,i]=sum(Y$myAUC<auc_diff, Y$myAUC>(1-auc_diff))
      }

    }

  }

  df=data.frame(mmm)
  rownames(df) <- clusters_unique
  colnames(df) <- clusters_unique

  merging_matrix[[k]]=df

  update_min=min(mmm, na.rm = TRUE) #find minimum de genes
  cond=which(mmm == update_min, arr.ind = TRUE) #find clusters with min de gene

  if(length(cond) == 2){
    clusters_active=replace(clusters_active,
      clusters_active==clusters_unique[cond[1]],
      clusters_unique[cond[2]]) # merge clusters

  } else if (length(cond) > 2){
    clusters_active=replace(clusters_active,
      clusters_active==clusters_unique[cond[1,][1]],
      clusters_unique[cond[1,][2]])  # merge clusters

  }

  Idents(object=seur)=factor(clusters_active)
  clusters_unique.update=sort(unique(unfactor(seur@active.ident)))

  k=k+1

  clusters_unique=clusters_unique.update

  merging_progress[[k]]=clusters_active


}


# Save results
merging_progress=do.call(cbind,merging_progress)
colnames(merging_progress) <- paste0("merge_step_",0:(k-1))
rownames(merging_progress) <- colnames(seur)
write.table(merging_progress,
  file=file.path(output.path, basename),
  sep='\t')
