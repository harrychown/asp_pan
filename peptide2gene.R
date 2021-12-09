library("biomaRt")
library(rentrez)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
peptide_acc <- readLines("centroid.pep.names")
centroid_id <- readLines("centroid.order")
gene_acc <- c()
for(pep in peptide_acc){
  linked_seq_ids <- entrez_link(dbfrom="protein", id=pep, db="protein")
  protein_link <- linked_seq_ids$links$protein_protein_cdart
  all_recs <- entrez_fetch(db="protein", id=protein_link, rettype="xml")
  gene_tag <- str_match(all_recs, "locus_tag \"(.*?)\"\n")[1,2]
  if(is.na(gene_tag)){
    gene_tag <- str_match(all_recs, "locus-tag \"(.*?)\"\n")[1,2]
    if(is.na(gene_tag)){
      print(all_recs)
      print(pep)
      poop
    }
  }
  gene_acc <- c(gene_acc, gene_tag)
}
df <- data.frame("centroid"=centroid_id,
                 "peptide" = peptide_acc,
                 "gene" = gene_acc)