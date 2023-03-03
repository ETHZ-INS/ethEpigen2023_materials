# Assignment

1. Using AnnotationHub, find and download the following annotations data:       
- The mouse (Mus Musculus) EnsDb object, version 102, genome build GRCm38  
- The mouse genome sequence ( dna_sm ) in TwoBit/2bit format for GRCm38  
- The drosophila melanogaster genome sequence ( dna_sm ) in TwoBit/2bit format for BDGP6  

2. Using the mouse EnsDb, find the following: 
  - How many different ensembl gene IDs and gene symbols are there for protein-coding genes?  
  - Plot the distribution of the (spliced) length of protein-coding transcripts  
      - (tip: this will require you to extract exons of protein-coding transcripts from the database, and split them by transcript, before summing the width of the exons of each transcript)  
    
    
Name your markdown file `assignment.Rmd`, render it, and put it (along with the produced html) in the `week02` folder of your repository, and push!
