install.packages(c("gsl","ggplot2","remotes"))
# If you haven't install Bioconductor yet:
# install.packages("BiocManager")
# BiocManager::install()

BiocManager::install(c("GenomicRanges", "rtracklayer", "EnrichedHeatmap", "AnnotationHub", 
                       "ensembldb", "edgeR", "esATAC", "sechm","genomation","Rsubread","Rfastp"))
BiocManager::install("ETHZ-INS/epiwraps")
