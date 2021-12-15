library("rstudioapi")   
library(tibble)
setwd(dirname(getActiveDocumentContext()$path))
setwd("/data/OneDrive/asp_pan/post_panoct/scoary_out")
results_list <- list.files(pattern=".results")
centroid_key <- read.csv("/home/harry/centroid.info.clean.key", header = F)

# Add centroid information to results data
add_info <- function(in_data, in_key){
  out_data <- in_data
  
  peptide_col <- c()
  info_col <- c()
  rownames(in_key) <- in_key[,1]
  for(centroid in in_data[,1]){
    info <- in_key[centroid, 2]
    peptide <- in_key[centroid, 3]
    info_col <- c(info_col, info)
    peptide_col <- c(peptide_col, peptide)
  }
  out_data$Non.unique.Gene.name <- info_col
  out_data$Annotation <- peptide_col
  colnames(out_data)[2] <- "peptide"
  return(out_data)
}


# Obtain genes associated with the clades
clade_data <- read.csv(results_list[1])
clade_data <- add_info(clade_data, centroid_key)
significant_clade <- clade_data[clade_data$Bonferroni_p < 0.05,]
cladeA <- significant_clade[significant_clade$Odds_ratio < 1, 1]
cladeB <- significant_clade[significant_clade$Odds_ratio > 1, 1]

write.csv(significant_clade, "~/filtered_scoary_clade.csv")
write(cladeA, "~/cladeA_associated.txt")
write(cladeB, "~/cladeB_associated.txt")

itra_data <- read.csv(results_list[2])
itra_data <- add_info(itra_data, centroid_key)
significant_itra <- itra_data[itra_data$Bonferroni_p < 0.05,]
write.csv(significant_itra, "~/filtered_scoary_itra.csv")

itra_sus <- significant_itra[significant_itra$Odds_ratio < 1, 1]
itra_res <- significant_itra[significant_itra$Odds_ratio > 1, 1]
write(itra_sus, "~/susceptible_associated.txt")
write(itra_res, "~/resistance_associated.txt")

mut_data <- read.csv(results_list[3])
mut_data <- add_info(mut_data, centroid_key)
significant_mut <- mut_data[mut_data$Bonferroni_p < 0.05,]
mutant <- significant_mut[significant_mut$Odds_ratio < 1, 1]
write(mutant, "~/mutant_associated.txt")

wildtype <- significant_mut[significant_mut$Odds_ratio > 1, 1]
write(wildtype, "~/wildtype_associated.txt")
write.csv(significant_mut, "~/filtered_scoary_itra_res.csv")

