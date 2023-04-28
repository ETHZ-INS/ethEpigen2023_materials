# Assignment

* Download ATAC-seq peak counts in the hippocampus upon stress (subset of the original data, already in SummarizedExperiment format) :
  * https://ethz-ins.org/content/mouse_mm38_hippocampus.peakCounts.SE.rds
* Using this object, perform a chromVAR motif analysis, and run 2 differential motif accessibility analyses, respectively:
  * comparing stressed (denoted ‘FSS’ – forced swim stress) and control animals
  * comparing male and female animals
* For each analysis, report the top most significant motifs, plot a heatmap of the normalized accessibility scores across the samples for those motifs, and write a short paragraph interpreting the results.

Save your assignment in a R markdown named `assignment.Rmd`, render it, and push the html file to this folder in your github repository
