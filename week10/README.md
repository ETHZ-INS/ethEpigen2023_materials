# Assignment

* Download and decompress the following archive:
  * https://ethz-ins.org/content/w10.assignment.zip
  * This contains the bigwig files and peaks (bed) files for three TFs of the CREB family (all restricted to chr1; aligned against the hg38 genome)
* Use clustering and visualization to illustrate the relationship between the binding of the different proteins
* Use enrichment analysis (either GO or motif) on at least one of the clusters
* Write a paragraph describing your results

Save your assignment in a R markdown named `assignment.Rmd`, render it, and push the html file to this folder in your github repository

### Tip

To get a clearer picture, focus on high-confidence peaks from each factor to define the universe of regions, e.g.:

```
peaks <- list.files(pattern="bed$")
# we first import the peaks
peaks <- lapply(peaks, rtracklayer::import.bed)
# we'll focus on the high-quality peaks
peaks <- lapply(peaks, FUN=function(x) x[x$score>800])
# we get the union of non-redundant regions
regions <- reduce(unlist(GRangesList(peaks)))
```
