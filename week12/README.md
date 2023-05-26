# Assignment

* Choose a transcription factor (e.g. p300), and obtain peaks from ENCODE (ChIP-seq in a human context!)
* Isolate the peaks that are:
  - Between 2.5kb and 10kb from a TSS
  - More than 10kb from a TSS
* For each set of peaks:
  - Subset to those peaks that have a predicted distal target(s) using Salviato et al. (2021)
    - You can download a GRanges of those interactions at https://ethz-ins.org/content/hg38.SalviatoDistalEnhancerTargets.GR.rds 
  - Find the nearest TSS for each peak
  - In what proportion of the cases is the predicted target the closest gene?
* Hints:
  - you can use the annotateRegions function, as we did in week 4, to get the gene nearest to each peak
  - beware not to count, when calculating proportions, peaks that don’t have interactions with any TSS!
* Expected for of the answer: “Of the genes that are between 2.5 and 10kb from the nearest TSS, XX % form an interaction with that nearest gene. Of the genes that are more than 10kb away from the nearest TSS, XX % form an interaction with that nearest gene.”


Save your assignment in a R markdown named `assignment.Rmd`, render it, and push the html file to this folder in your github repository
