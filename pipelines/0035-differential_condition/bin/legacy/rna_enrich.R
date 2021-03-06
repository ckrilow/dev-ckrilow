#make changes from current LRpath/RNA-Enrich script and save to dropbox, add main figure files, supplement as pdf

rna_enrich = function(sigvals,geneids, avg_readcount, species, direction=NULL, min.g = 10,  max.g=NA, sig.cutoff=0.05, database="GOBP", read_lim=5,conceptList=NULL, nullsetList=NULL, plot_file="RNA-Enrich_plot.jpg",plot_height=480, plot_width=480,results_file="RNA-Enrich_results.txt") {
###########################################################################################
#  Function for RNA-Enrich: A cut-off free functional enrichment testing method 
#  for RNA-seq with improved detection power
#  Written by: Chee Lee, University of Michigan, 2015
###########################################################################################
##
##  This function uses modified random sets method with weights calculated from a spline to 
##  adjust for any relationship between read count and differential expression p-value to test for 
##  enriched (or depleted) biological categories in gene expression data. 
##
##  Please acknowledge your use of LRpath in publications by referencing:
##  Lee C, Patil S, Sartor MA. RNA-Enrich: A cut-off free functional enrichment testing 
##  method for RNA-seq with improved detection power. 
##
##  Inputs:
##  9 parameters: sigvals, geneids, avg_readcount,species, direction, min.g, max.g, sig.cutoff, database, read_lim
##  "sigvals" a vector of p-values, same length and order as "geneids"
##  "geneids" a vector of Entrez gene IDs, may contain duplicates and missing values, same length and order as sigvals.
##  "species"  Used only for KEGG, GO, and Cytoband.  human="hsa", mouse="mmu", rat="rno", "dme" for fly, 
##  	"dre" for zebrafish, "cel" for C.elegans, "sce" for yeast etc. (no default)
##  "direction"  Optional parameter for directional testing.  If desired, a vector, of same length and order as sigvals.
##		Values should be negative if down-regulated and positive if up-regulated. (e.g. -1 or 1) Only the sign is used.
##  "min.g"  The minimum number of unique gene IDs analyzed in category to be tested
##		default = 10
##  "max.g"  The maximum number of unique gene IDs analyzed in category to be tested
##		default = NA (99999)
##  "sig.cutoff" Entrez gene IDs in each category with p-values<sig.cutoff will be returned
##		default = 0.05
##  "database" database to be tested- choices are "GO", "KEGG", "Cytoband","custom", 
##		default = "GOBP"
##		"custom" for user provided gene sets
##		when "custom" is selected, conceptList is list of custom gene set IDs, 
##		nullsetList is list of genes connected to conceptList, 
##  	should be a 1 gene set ID to 1 gene ID relationship between conceptList and nullsetList.
##  "read_lim" Read limit of genes that will be included in analysis
##  	default = 5
##	"conceptList" List of concepts or gene sets, use only when database="custom"
##	"nullsetList" List of genes, use only when database="custom"
##  	conceptList and nullsetList should be a 1 gene set ID to 1 gene ID relationship between conceptList and nullsetList.
##  "plot_file" Name of file of output plot
##  	default = "RNA-Enrich_plot.jpg"
##  "plot_height" Height in pixels of plot
##  	default = 480 pixels
##  "plot_width" Width in pixels of plot
##  	default = 480 pixels
##  "results_file" Name of file of tab-delimited text output of results
##  	default = "RNA-Enrich_results.txt"
##
##  Outputs:
##  object is dataframe with the following columns:
##	1) GO, KEGG, or Concept ID	
##	2) GO, KEGG, or Concept term - name of category
##	3) database - Database of gene sets that was tested
##	4) n.genes - Number of unique Entrez Gene IDs in category
##	5) coeff - coefficient of slope (positive values indicate enrichment or up-regulation; negative values 
##		indicate depletion or down-regulation)
##  6) odds.ratio -  Odds ratio, as measure of strenth of enrichment (or depletion)
##  7) status - "Enriched" or "Depleted", or if direction is used "Up-regulated" or "Down-regulated"
## 	8) p.value	-  P-value that slope does not equal zero, i.e. that term is enriched (or depleted)
##	9) FDR	- False Discovery Rate (Benjamini & Hochberg, 1995)
##  10) sig.genes	- comma separated Entrez gene ids in category with p-value<"sig.cutoff"
##
##    Examples:
##    pvalues <- IBMT.results$IBMT.p[,1]
##    entrez.geneid <- unlist(mget(featureNames(affy.eset),env=hgu133aENTREZID))
##    GO.results <- LRpath(sigvals=pvalues, geneids=entrez.geneid, species="hsa", database="GO")
##    KEGG.results <- LRpath(sigvals=pvalues, geneids=entrez.geneid, database="KEGG",
##		species="hsa")
##
##########################################################################################



	rna_enrich_plot(sigvals=sigvals, geneids=geneids, avg_readcount=avg_readcount, filename=plot_file, plot_height=plot_height, plot_width=plot_width)

	de_df = data.frame(geneid = geneids, p.value = sigvals, readcount = avg_readcount)

	if(!is.null(direction)){
		de_df$direction = direction
		de_df = stats::aggregate(cbind(p.value,readcount,direction) ~ geneid, de_df, mean) 
		de_df$direction = ifelse(de_df$direction>0,1,-1)
	}else{ 	
		de_df = stats::aggregate(cbind(p.value,readcount) ~ geneid, de_df, mean)
	}
		
	
	
	de_df = na.omit(de_df)
	de_df = subset(de_df, geneid !='')
	
	if(!species %in% c('dme','sce')){
		de_df$geneid = suppressWarnings(as.numeric(as.character(de_df$geneid)))
		de_df = na.omit(de_df)
	}
library("annotate")
library("stats")
	
 # Get entire gene list using GO.db, for RNA-Enrich, must test each branch of GO separately for later weighting step to work correctly
		if (database %in% c('GOBP','GOCC','GOMF')) {
			 library("GO.db")
			 if (species=="hsa") {
				library("org.Hs.eg.db")
				xx<-as.list(org.Hs.egGO2ALLEGS)
				ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="mmu") {
			library("org.Mm.eg.db")
				xx<-as.list(org.Mm.egGO2ALLEGS)
				ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="rno") {
			library("org.Rn.eg.db")
			xx<-as.list(org.Rn.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="dme") {
			library("org.Dm.eg.db")
			xx<-as.list(org.Dm.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="dre") {
			library("org.Dr.eg.db")
			xx<-as.list(org.Dr.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="cel") {
			library("org.Ce.eg.db")
			xx<-as.list(org.Ce.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="sce") {
			library("org.Sc.sgd.db")
			xx<-as.list(org.Sc.sgdGO2ALLORFS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }

     ##Remove GO ids that aren't mapped to any Gene id
     xx<-xx[!is.na(xx)]

	## Determine which Entrez Gene IDs are annotated anywhere in GO
	## Gene ids not annotated in GO are left out of the environment
     if (species %in% c("sce")) {	
		nullset<-unique(unlist(xx)) 
	} else {
        nullset<-as.numeric(unique(unlist(xx)))
	}
}
#########  KEGG works for human, mouse, rat, and any other species that uses Entrez geneID
if (database=="KEGG") {	
	library(KEGG.db)
	xx<-as.list(KEGGPATHID2EXTID)
	xx<-xx[!is.na(xx)]
	xx<-xx[grep(species,names(xx))]

	## Determine which Entrez Gene IDs are annotated anywhere in KEGG
	if (species %in% c("dme","sce")) {						
	   nullset <-names(as.list(KEGGEXTID2PATHID))			
	} else								
	nullset<- as.numeric(names(as.list(KEGGEXTID2PATHID)))	
	nullset<- nullset[!is.na(nullset)]
}


#######  For this, we'll use data in species specific R packages, e.g. org.Hs.eg.db for human
	if (database=="Cytoband") {
		 if (species=="hsa") {
		library("org.Hs.eg.db")
			###  Get GO list of Entrez genes
		xx1<-as.list(org.Hs.egMAP)
		 }
		 else if (species!="hsa") {
			stop("Only human is supported in Cytoband")
		 }
	
	## Currently, list names are Entrez IDs, and values are cytobands.  Need to switch.
	xx1<-xx1[!is.na(xx1)]
	cyto.1<-as.vector(sapply(xx1,function(x) { x[1] }))  ## names(xx1) are the entrez IDs to match this
	xx<-list()
	for (i in 1:length(unique(cyto.1))) {
		xx[[i]]<- names(xx1)[cyto.1==unique(cyto.1)[i]]
	}
	names(xx)<-as.character(unique(cyto.1))

	## Get null set vector for Cytoband
	nullset<- as.numeric(names(xx1))
	nullset<- nullset[!is.na(nullset)]
}

if (database=="custom"){
if(!is.null(conceptList) & !is.null(nullsetList) ){
	xx = by(nullsetList, conceptList, list)
	xx = xx[!is.na(xx)]
	nullset = unique(as.numeric(unlist(xx)))
}else{
			stop("Need conceptList and nullsetList to proceed!")
}}
	
	
	
	catsizes<-sapply(xx,length);
	if(is.na(max.g)){max.g = max(catsizes)}
	xx<-xx[catsizes>=min.g & catsizes<=max.g];
	ncats<-length(xx);

	#make empty d dataframe to initiate gpw creation
	d_allgenes = data.frame(geneid=nullset[which(!nullset %in% de_df$geneid)],p.value=NA,readcount=0);
	if(!is.null(direction)){
		d_allgenes$direction = NA;
		d_allgenes$direction = as.numeric(d_allgenes$direction)
	};
	gpw = rbind(de_df,d_allgenes);


	# Add log10 p-values to gpw, add min p[!=0] to all pvals so can log transform later
	gpw$log10_p.value = (-1)*log10(gpw$p.value + min(gpw$p.value[gpw$p.value !=0],na.rm=T));

	# remove genes with read lim readcounts (default is 5) and log readcount
	if(read_lim==0){ 
		gpw = subset(gpw, readcount>read_lim);
	}else{
		gpw = subset(gpw, readcount>=read_lim);
	}
	gpw$log10_readcount = log10(gpw$readcount);

	# Sort by readcount 
	gpw = gpw[order(gpw$readcount,decreasing=T),];

	gpw = gpw[!is.na(gpw$log10_p.value),];
	# Create model - adjust for readcount
	model = "log10_p.value ~ s(log10_readcount,bs='cr')";

	# Compute binomial spline fit.
	library(mgcv); library(rms);
	fit = gam(as.formula(model),data=gpw,family="gaussian");

	# Compute weights for each gene, based on the predicted prob(DE) for each gene. 
	pP = fitted(fit);
	w0 = 1 / (pP/mean(gpw$log10_p.value,na.rm=T));
	w0 = w0 / mean(w0,na.rm=T);

	gpw$weight = w0;
	gpw$prob_P = pP;
	gpw$resid.dev = resid(fit,type="deviance");

	# change log pvals if directional test
	if (!is.null(direction)) {
		gpw$log10_p.value = gpw$log10_p.value * gpw$direction
	}

	cols = c("geneid","readcount","log10_readcount","p.value","log10_p.value","weight","prob_P","resid.dev");
	gpw = subset(gpw,select=cols);

	#now run test	

	lrm.fast = function(x,y) {
		fit = lrm.fit(x,y);
		vv = diag(fit$var);
		cof = fit$coef;
		z = cof/sqrt(vv);
		pval = pchisq(z^2,1,lower.tail=F);
		c(cof[2],pval[2]);
	};
	
	# Restrict our genes/weights/peaks to only those genes in the genesets. 
	gpw = subset(gpw,geneid %in% nullset);
	
	# Re-normalize weights so that mean is 1 (if any genes were removed from previous gpw step)
	gpw$weight = gpw$weight / mean(gpw$weight);
	
	r_pvals = c();
	r_effects = c();
	r_go_ids = c();
	r_go_genes = c();
	r_go_genes_num = c();
	r_go_genes_DE = c();
	r_go_genes_DE_num = c();
	r_go_genes_avg_length = c();
	for (i in 1:length(xx)) {
		
		go_id = names(xx)[i];
		
		if(species %in% c('dme','sce')){
			go_genes = as.character(xx[[i]]);
		}else{
			go_genes = as.numeric(as.character(xx[[i]]));
		}
		# Eliminate GO genes that aren't in the gpw. 
		go_genes = go_genes[go_genes %in% gpw$geneid];
		
		# A boolean vector of gene membership of all geneids in go_genes
		b_genes = gpw$geneid %in% go_genes;
		sg_go = gpw$log10_p.value[b_genes];
		wg_go = gpw$weight[b_genes];
		
		r_go_genes[i] = paste(go_genes,collapse=";");
		r_go_genes_num[i] = length(go_genes);
		
		r_go_genes_DE_num[i] = length(gpw$geneid[gpw$p.value<=sig.cutoff & b_genes==TRUE]);
		go_genes_DE = gpw$geneid[gpw$p.value<=sig.cutoff & b_genes==TRUE]; 
		if(length(go_genes_DE)==0){go_genes_DE=NA}
		r_go_genes_DE[i] = paste(go_genes_DE,collapse=", ");
		
		testm = cbind(y=as.numeric(b_genes), x=gpw$weight*gpw$log10_p.value);
		
		if(sum(testm[,'y']==1) == nrow(gpw) | sum(testm[,'y']) ==0){
		ep = c(NA,NA)
		}else{ep = lrm.fast(testm[,"x"], testm[,"y"])}
		
		r_effects[i] = ep[1];
		r_pvals[i] = ep[2];
		r_go_ids[i] = go_id;
		
		if ((i-ncats/10)>0&(i-ncats/10)<1.1) {
		print("10% categories finished.")
		};
		if ((i-ncats/5)>0&(i-ncats/5)<1.1) {
			print("20% categories finished.")
		};
		if ((i-ncats/2)>0&(i-ncats/2)<1.1) {
			print("50% categories finished.")
		};
		if ((i-9*ncats/10)>0&(i-9*ncats/10)<1.1) {
			print("90% categories finished.")
		};
	}
	
	fdr = p.adjust(r_pvals,method="BH");
  
	is_depleted = r_effects < 0; 
	is_enriched = r_effects > 0;
  
	enr = rep(NA,length(r_go_ids));
	if(!is.null(direction)){
		enr[is_depleted] = "down";
		enr[is_enriched] = "up";
	}else{
		enr[is_depleted] = "depleted";
		enr[is_enriched] = "enriched";
    }
 ##---annotate concept_id/r_go_ids with description, right now just have GO and KEGG, NAs for others
  concept_descriptions = NA
  if(database %in% c('GOBP','GOCC','GOMF')){
	concept_descriptions = Term(r_go_ids);
  };
  if(database=="KEGG") {
	KEGGterms<- as.list(KEGGPATHNAME2ID);
	kterms<- as.vector(names(KEGGterms));
	k.ids<- paste(species,unlist(KEGGterms),sep="");
	keggrows<- match(r_go_ids,k.ids);
	concept_descriptions<-kterms[keggrows];
  };
  if (database=="custom") {
	concept_descriptions = r_go_ids;
  };

	results = data.frame(
		"Concept ID"=r_go_ids,
		"Concept name" = concept_descriptions,
		"database" = database, 
		"n.genes"=r_go_genes_num, #changed 7/8, incorrect before
		"coeff"=r_effects,
		"odds.ratio"=exp(r_effects),
		"status"=enr,
		"p.value"=r_pvals,
		"FDR"=fdr,
		"sig.genes"=r_go_genes_DE,
		stringsAsFactors=F
	);
	results = results[ !is.na(results$p.value),];

	results = results[order(results$FDR,results$p.value),];
	write.table(results, file=results_file, sep="\t",row.names=F, quote=F)
	return(results);
}

rna_enrich_plot = function(sigvals,geneids,avg_readcount,filename="rna_enrich_plot.jpg", plot_height=480, plot_width=480){
	library(mgcv); library(rms);
	de_df = data.frame(geneid = geneids, p.value = sigvals, readcount = avg_readcount)
	de_df = na.omit(de_df)
	de_df = de_df[de_df$readcount>=5,]
	de_df = de_df[order(de_df$p.value),]
		
	rownames(de_df) = 1:nrow(de_df)
	de_df$group = ceiling(as.numeric(rownames(de_df))/25)
	gavgrc = as.numeric(by(de_df$readcount, de_df$group, mean))
	gavgp = as.numeric(by(de_df$p.value, de_df$group, mean))
	log10_gavgrc = log10(gavgrc)
	log10_gavgp = -log10(gavgp + min(gavgp[gavgp!=0]))

	de_df$log10_p.value = (-1)*log10(de_df$p.value + min(de_df$p.value[de_df$p.value !=0],na.rm=T))
	de_df$log10_readcount = log10(de_df$readcount)
	de_df = de_df[!is.na(de_df$log10_p.value),]
	model = "log10_p.value ~ s(log10_readcount,bs='cr')"
	  # Compute binomial spline fit.  
	fit = gam(as.formula(model),data=de_df,family="gaussian")
	  # Compute weights for each gene, based on the predicted prob(peak) for each gene. 
	pP = fitted(fit);
	 
	jpeg(filename, height=plot_height, width=plot_width)
	par(mar= c(5, 5, 4, 2) + 0.1)
	  
	plot(log10_gavgp ~ log10_gavgrc, ylab='-log10(p-value) - binned 25 genes', xlab='log10(avg read count) - binned 25 genes',pch=20, cex.axis=1.5, cex.lab=1.5) 
	lines(pP[order(de_df$log10_readcount)]~sort(de_df$log10_readcount),col='red',lwd=2)
	
	dev.off()
}
