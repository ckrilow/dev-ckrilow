
# DETAILS

* Project number: OTAR2-057
* Team Email: otar2057@opentargets.org
* Data freeze: v002


# SUMMARY

Single cell RNA-seq (scRNA-seq) of (A) terminal ileum biopsies (25 Crohn's disease and 8 healthy participants) and (B) rectum biopsies (9 healthy participants). This release includes: (1) 3' 10x scRNA-seq data processed with CellRanger v3.0.2 (genome reference = GRCh38) by Sanger DNA pipelines, (2) phenotype information on particpants, and (3) sample metadata.

Note that we include two CellRanger output datasets, labelled "filtered_feature_bc_matrix" and "raw_feature_bc_matrix". "Filtered" refers to data where empty droplets have been removed as part of standard CellRanger analysis. For the "raw_feature" data, empty droplets have not been removed.


# PROJECT OVERVIEW AND GOALS

The goal of this project to integrate data across all available and relevant sources in order to answer the following questions:

1) What genes have the highest probability of mediating GWAS signals and, therefore, to play a causal role in IBD?

2) What are the broader phenotypic consequences of perturbing these genes?

This effort will inform drug development efforts in IBD, as we will identify potential novel targets (list of causal genes), on-target safety concerns (PheWAS), additional indications (PheWAS), and early estimates of efficacy and therapeutic windows (dose-response modeling). Through detailed review of individual genes as potential drug targets, this work will also serve as an application-focused pilot for the Open Targets genetics portal. Examples pursued in this proposal will provide use cases to help guide the development of work flows, tools, and visualization strategies.

It is often difficult to draw insight into drug targets from GWAS studies because many of the signals lie outside of genes, in sequence that is presumed regulatory. Without knowledge of which genes are regulated by these non-coding variants, one cannot infer drug targets or uncover new disease insights. This is because the ability of GWAS to pinpoint IBD drug targets hinges on our capacity to correlate protein function with IBD risk, using genetic variants as a causal anchor. Thus, if we want to translate our significant IBD locus discovery efforts into biological understanding and drug targets, it is vital we create definitive genetic maps of gene regulation in disease relevant cells.

While a lot of effort has been put into identifying eQTLs in the peripheral immune system, there are currently no cell-type specific eQTLs maps for gut tissue. This could explain why only around 15% of IBD loci currently colocalise with a known eQTL. Thankfully, gut tissue is readily accessible through routine clinical care. This, combined with recent progress in cellular and genomic technologies, presents a unique opportunity for deriving novel genetically validated IBD drug targets. In this project we will perform single-cell RNA sequencing (scRNAseq) of healthy terminal ileum biopsies from 350 patients undergoing screening for colorectal cancer to create the definitive TI cellular map. We will genotype these individuals, identify genetic variants associated with gene-expression (eQTLs), and combine these with existing IBD GWAS data to identify genes, cell types and pathways causally underpinning the disease.
