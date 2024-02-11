install.packages(c("gsl","ggplot2","remotes","rmarkdown"))
# If you haven't install Bioconductor yet:
# install.packages("BiocManager")
# BiocManager::install()

BiocManager::install(c("GenomicRanges", "rtracklayer", "EnrichedHeatmap", "AnnotationHub", 
                       "ensembldb", "edgeR", "esATAC", "sechm", "motifmatchr","rGREAT",
                       "bsseq","DMRcate","data.table","InteractionSet","chromVAR","limma",
                       "universalmotif", "MotifDb", "TFBSTools", "Biostrings", "PWMEnrich",
                       "Rsubread","Rfastp"))
BiocManager::install("Bioconductor/BiocFileCache")
BiocManager::install("ETHZ-INS/epiwraps")
