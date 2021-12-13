library("rstudioapi")   

setwd(dirname(getActiveDocumentContext()$path))
setwd("/data/OneDrive/asp_pan/post_panoct/scoary_out")
results_list <- list.files(pattern=".results")

# Obtain genes associated with the clades
clade_data <- read.csv(results_list[1])
significant_clade <- clade_data[clade_data$Bonferroni_p < 0.05,]
cladeA <- significant_clade[significant_clade$Odds_ratio < 1, 1]
cladeB <- significant_clade[significant_clade$Odds_ratio > 1, 1]
write.csv(significant_clade, "filtered_scoary_clade.csv")
write(cladeA, "cladeA_associated.txt")
write(cladeB, "cladeB_associated.txt")

itra_data <- read.csv(results_list[2])
significant_itra <- itra_data[itra_data$Bonferroni_p < 0.05,]
write.csv(significant_itra, "filtered_scoary_itra.csv")

itra_sus <- significant_itra[significant_itra$Odds_ratio < 1, 1]
itra_res <- significant_itra[significant_itra$Odds_ratio > 1, 1]
write(itra_sus, "susceptible_associated.txt")
write(itra_res, "resistance_associated.txt")

mut_data <- read.csv(results_list[3])
significant_mut <- mut_data[mut_data$Bonferroni_p < 0.05,]
mutant <- significant_mut[significant_mut$Odds_ratio < 1, 1]
write(mutant, "mutant_associated.txt")

wildtype <- significant_mut[significant_mut$Odds_ratio > 1, 1]
write(wildtype, "wildtype_associated.txt")
write.csv(significant_mut, "filtered_scoary_itra_res.csv")

