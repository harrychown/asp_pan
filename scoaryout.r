library("rstudioapi")   
setwd(dirname(getActiveDocumentContext()$path))
setwd("/home/harry/Documents/imperial/12_panoct/scoary_out")
metadata <- read.csv("clean_meta.csv", row.names =1)
binary <- read.csv("gene_binary.csv", row.names = 1)


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


# resistant_data <- all_scoary_results$trait2_Sitr_Ritr
# mutational_data <- all_scoary_results$trait3_RitrK_RitrU
# 
# itra_res_isolates <- rownames(metadata$itr == 1)
# itra_sus_isolates <- rownames(metadata$itr == 0)
# 
# itra_wt_isolates <- rownames(metadata$AMR == "Wildtype")
# itra_k_isolates <- rownames(metadata$AMR != "Wildtype")
# 
# r_percent <- c()
# s_percent <- c()
# 
# for(centroid in resistant_data$Gene){
#   pav <- binary[centroid, ]
#   pav_id <- which(pav == 1)
#   present_names <- colnames(binary)[pav_id]
#   present_in_resistant <- length(which(present_names %in% itra_res_isolates))
#   percentage_resistant <- present_in_resistant / length(itra_res_isolates) * 100
#   
#   present_in_susceptible <- length(which(present_names %in% itra_sus_isolates))
#   percentage_susceptible <- present_in_susceptible / length(itra_sus_isolates) * 100
#   
#   present_in_wt <- length(which(present_names %in% itra_wt_isolates))
#   percentage_wt <- present_in_wt / length(itra_wt_isolates) * 100
#   
#   present_in_k <- length(which(present_names %in% itra_k_isolates))
#   percentage_k <- present_in_k / length(itra_k_isolates) * 100
#   
#   check_mut <- 0
#   if(centroid %in% mutational_data$Gene){
#     check_mut <- 1
#   }
# }
# 
# 
# 
# 
# 
# 
# results_list <- list.files(pattern=".results")
# all_scoary_results <- list()
# listnames <- c()
# for(i in 1:length(results_list)){
#   name <- strsplit(results_list[i], "_03")[[1]][1]
#   listnames <- c(listnames, name)
#   result_data <- read.csv(results_list[i])
#   significant_data <- result_data[result_data$Bonferroni_p < 0.05,]
#   all_scoary_results[[i]] <- significant_data
# 
# }
# names(all_scoary_results) <- listnames
# 
# # for(i in 1:3){
# #   resistant_data <- all_scoary_results[[i+1]]
# #   mutational_data <- all_scoary_results[[i+2]]
# # }
# 
# resistant_data <- all_scoary_results$trait2_Sitr_Ritr
# mutational_data <- all_scoary_results$trait3_RitrK_RitrU
# 
# total_resistant <- resistant_data[1, 4] + resistant_data[1,6]
# percentage_resistant <- resistant_data[, 4] / total_resistant * 100
# 
# total_susceptible <- resistant_data[1, 5] + resistant_data[1,7]
# percentage_susceptible <- resistant_data[, 5] / total_susceptible * 100
# 
# total_wildtype <- mutational_data[1, 4] + mutational_data[1,6]
# percentage_wildtype <- mutational_data[, 4] / total_wildtype * 100
# 
# total_mutant <- mutational_data[1, 5] + mutational_data[1,7]
# percentage_mutant <- mutational_data[, 5] / total_mutant * 100
# 
# 
# 
# 
# resistant_out <- data.frame("Centroid"=resistant_data[,1],
#                             "Percentage Resistant"=percentage_resistant,
#                             "Percentage Susceptible"=percentage_susceptible)
# 
# 


# itra_data <- read.csv(results_list[2])
# itra_data <- itra_data[itra_data$Bonferroni_p < 0.05,]
# itra_resistant <- itra_data[itra_data$Odds_ratio > 1,]
# itra_susceptible <- itra_data[itra_data$Odds_ratio < 1,]
# 
# itr_amr_data <- read.csv(results_list[3])
# itr_amr_data <- itr_amr_data[itr_amr_data$Bonferroni_p < 0.05,]
# itr_unkown <- itr_amr_data[itr_amr_data$Odds_ratio > 1,]
# itr_known <- itr_amr_data[itr_amr_data$Odds_ratio < 1,]
# 
# itr_unkown[which(! itr_unkown$Gene %in% itra_susceptible$Gene),]
