# PANOCT HISTOGRAMS

library("tidyverse")
library("micropan")
library("ggplot2")
library("VennDiagram")
library("rstudioapi")   
setwd(dirname(getActiveDocumentContext()$path))     
# Read config file for input/output directories
directory_configuration <- readLines("directory.config")
in_dir <- directory_configuration[1]
out_dir <- directory_configuration[2]

setwd(in_dir)
raw_data <- read.table("matchtable.txt")[,-1]
raw_centroids <- readLines("centroids.fasta")
raw_meta <- read.csv("raw_metadata.csv")
essential_centroids <- readLines("essential.centroids")

### ADD HEADER/ROW NAMES BASED ON CENTROIDS ###
# Change row names to centroid index
centroid_names <- c()
for(i in 1:nrow(raw_data)){
  centroid <- paste(c("centroid_", as.character(i)), sep = "", collapse = "")
  centroid_names <- c(centroid_names, centroid)
}
rownames(raw_data) <- centroid_names

# Change column headers to isolate names
isolate_names <- c()
for(c in raw_data[1,]){
  isolate <- str_split(c, "_")[[1]][1]
  isolate_names <- c(isolate_names, isolate)
}
colnames(raw_data) <- isolate_names

## CREATE DATA TABLES ##
# Convert to a binary data format
binary_data <- raw_data
binary_data[binary_data=="----------"]<-0
binary_data[binary_data!="0"]<-1
binary_data <- data.frame(sapply(binary_data, function(x) as.numeric(as.character(x))))
rownames(binary_data) <-  rownames(raw_data)
#write.csv(binary_data, "/home/harry/Documents/imperial/1_11_panoct/gene_binary.csv" )
# Subset the core data
core_data <- binary_data[rowSums(binary_data == 0)==0, , drop = FALSE]
# Get frequency for each peptide
frequency <- rowSums(binary_data)
names(frequency) <- centroid_names
#write.csv(frequency, "/home/harry/Documents/imperial/1_11_panoct/gene_freq.csv" )

## CREATE DATA TO MIMIC ROARY OUTPUT ##
# Replace lacking genes with blank space
data_blanked <- raw_data
data_blanked[data_blanked=="----------"]<-""
db_isolates <- colnames(data_blanked)
mock_data <- data.frame("Gene"=rep(0, nrow(data_blanked)),
                        "Non-unique Gene name"=rep(0, nrow(data_blanked)),
                        "Annotation"=rep(0, nrow(data_blanked)),
                        "No. isolates"=rep(0, nrow(data_blanked)),
                        "No. sequences"=rep(0, nrow(data_blanked)),
                        "Avg sequences per isolate"=rep(0, nrow(data_blanked)),
                        "Genome Fragment"=rep(0, nrow(data_blanked)),
                        "Order within Fragment"=rep(0, nrow(data_blanked)),
                        "Accessory Fragment"=rep(0, nrow(data_blanked)),
                        "Accessory Order with Fragment"=rep(0, nrow(data_blanked)),
                        "QC"=rep(0, nrow(data_blanked)),
                        "Min group size nuc"=rep(0, nrow(data_blanked)),
                        "Max group size nuc"=rep(0, nrow(data_blanked)),
                        "Avg group size nuc"=rep(0, nrow(data_blanked)))
mock_roary <- cbind(mock_data, data_blanked)
mock_roary$Gene <- rownames(data_blanked)
write.csv(mock_roary, paste(c(out_dir,"/mock_roary.csv"), sep="", collapse=""), row.names=F)


## CLEAN METADATA ##
# Remove ">" or "<" strings from mic
# order: itr, vor, pos
clean_meta <- raw_meta
drug_name <- c("itr", "vor", "pos")
drug_col <- 15:17
drug_mic <- c(2, 2, 0.5)
for(i in 1:length(drug_mic)){
  NA_id <- which(clean_meta[,drug_col[i]] == "ND")
  clean_meta[,drug_col[i]] <-  gsub("<","",as.character(clean_meta[,drug_col[i]]))
  clean_meta[,drug_col[i]] <-  gsub(">","",as.character(clean_meta[,drug_col[i]]))
  clean_meta[,drug_col[i]] <-  as.numeric(clean_meta[,drug_col[i]])
  resistant_id <- which(clean_meta[,drug_col[i]] >= drug_mic[i])
  
  new_col <- ncol(clean_meta) + 1
  clean_meta[,new_col] <- rep(0, nrow(clean_meta))
  clean_meta[NA_id,new_col] <- NA
  clean_meta[resistant_id,new_col] <- 1
  
}
colnames(clean_meta) <- c(colnames(raw_meta), drug_name)


write.csv(clean_meta, paste(c(out_dir,"/clean_meta.csv"), sep="", collapse=""), row.names=F)

## CREATE PHENOTYPE DATA FOR SCOARY INPUT ##
# Re-index metadata to match "mock roary"
ordered_meta <- clean_meta[ match(clean_meta$Pangenome.id, colnames(data_blanked)), ]
#strain   trait1    trait2
# Generate data for trait 1: Clade A vs Clade B
# Clade A = 0
# Clade B = 1
trait1_cladeA_cladeB <- ordered_meta$Clade
trait1_cladeA_cladeB[trait1_cladeA_cladeB == "A"] <- 0
trait1_cladeA_cladeB[trait1_cladeA_cladeB == "B"] <- 1

# Generate data for trait 2: Itraconazole susceptible vs resistant
# Itra sus = 0
# Itra res = 1
trait2_Sitr_Ritr <- ordered_meta$itr

# Generate data for trait 3: Itra-resistant known vs unknown
# Known = 0
# Unknown = 1
# Separate the itra-resistant isolates
itra_resistant <- which(ordered_meta$itr == 1)
# Identify the wildtypes to substitute later
itra_R_U <- which(ordered_meta$AMR == "Wildtype" & ordered_meta$itr == 1)
itra_R_K <- which(ordered_meta$AMR != "Wildtype" & ordered_meta$itr == 1)
trait3_RitrK_RitrU <- rep(NA, nrow(ordered_meta))
trait3_RitrK_RitrU[itra_R_U] <-  1
trait3_RitrK_RitrU[itra_R_K] <- 0

# Generate data for trait 4: Voriconazole susceptible vs resistant
# Vor sus = 0
# Vor res = 1
trait4_Svor_Rvor <- ordered_meta$vor

# Generate data for trait 5: Vor-resistant known vs unknown
# Known = 0
# Unknown = 1

# Identify the wildtypes to substitute later
vor_R_U <- which(ordered_meta$AMR == "Wildtype" & ordered_meta$vor == 1)
vor_R_K <- which(ordered_meta$AMR != "Wildtype" & ordered_meta$vor == 1)
trait5_RvorK_RvorU <- rep(NA, nrow(ordered_meta))
trait5_RvorK_RvorU[vor_R_U] <-  1
trait5_RvorK_RvorU[vor_R_K] <-  0

# Generate data for trait 6: posiconazole susceptible vs resistant
# Pos sus = 0
# Pos res = 1
trait6_Spos_Rpos <- ordered_meta$pos

# Generate data for trait 7: Pos-resistant known vs unknown
# Known = 0
# Unknown = 1
# Separate the Pos-resistant isolates
pos_resistant <- which(ordered_meta$pos == 1)
pos_susceptible <- which(ordered_meta$pos == 0)
# Identify the wildtypes to substitute later
pos_R_U <- which(ordered_meta$AMR == "Wildtype" & ordered_meta$pos == 1)
pos_R_K <- which(ordered_meta$AMR != "Wildtype" & ordered_meta$pos == 1)
trait7_RposK_RposU <- rep(NA, nrow(ordered_meta))
trait7_RposK_RposU[pos_R_U] <- 1
trait7_RposK_RposU[pos_R_K] <- 0

# Generate data for trait 8: All-azole susceptible vs resistant
# All sus = 0
# All res = 1
azole_sum <- rowSums(ordered_meta[,21:23])
azole_res <- which(azole_sum == 3)
azole_sus <- which(azole_sum == 0)
trait8_Sall_Rall <- rep(NA, nrow(ordered_meta))
trait8_Sall_Rall[azole_sus] <- 0
trait8_Sall_Rall[azole_res] <- 1

# Generate data for trait 9: All-resistant known vs unknown
# Known = 0
# Unknown = 1
# Identify the wildtypes to substitute later
pos_R_U <- which(ordered_meta$AMR == "Wildtype" & azole_sum == 3)
pos_R_K <- which(ordered_meta$AMR != "Wildtype" & azole_sum == 3)
trait9_RallK_RallU <- rep(NA, nrow(ordered_meta))
trait9_RallK_RallU[pos_R_U] <-  1
trait9_RallK_RallU[pos_R_K] <-  0

# Build the traits file
traits_data <- data.frame("trait1_cladeA_cladeB" = trait1_cladeA_cladeB,
                           "trait2_itra_sensitive_resistant" = trait2_Sitr_Ritr,
                           "trait3_itraresistant_known_unknown" = trait3_RitrK_RitrU,
                           "trait4_vor_sensitive_resistant" = trait4_Svor_Rvor,
                           "trait5_vorresistant_known_unknown" = trait5_RvorK_RvorU,
                           "trait6_pos_sensitive_resistant" = trait6_Spos_Rpos,
                           "trait7_posresistant_known_unknown" = trait7_RposK_RposU,
                           "trait8_all_sensitive_resistant" = trait8_Sall_Rall,
                           "trait9_allresistant_known_unknown" = trait9_RallK_RallU)
traits_data <- rownames(ordered_meta$Pangenome.id)
write.csv(traits_data, paste(c(out_dir,"/traits_data.csv"), sep="", collapse=""), row.names=T)

## CHECK ESSENTIAL GENE PRESENCE ##
essential_freq <- frequency[essential_centroids]
essential_and_core <- essential_freq[which(essential_freq == 218)]
essential_noncore <- essential_freq[which(essential_freq < 218)]

## COUNT NUMBER OF UNIQUE GENES PER STRAIN ##
strain_specific_data <- binary_data[which(frequency==1),]
strain_specific_freq <- colSums(strain_specific_data)
high_specific_freq <- strain_specific_freq[which(strain_specific_freq > 50)]
# write.csv(strain_specific_freq, "/home/harry/Documents/imperial/1_11_panoct/strain_specifc_freq.csv" )

## EXTRACT C30-SPECIFIC GENES ###
c30_specific_data <- strain_specific_data[which(strain_specific_data$C30 == 1),]
c30_centroids <- rownames(c30_specific_data)
c30_raw <- raw_data[c30_centroids, "C30"]

c30_num <- c()
for(i in c30_raw){
  transcript_id <- str_split(i, "_")[[1]][2]
  transcript_num <- as.numeric(transcript_id)
  c30_num <- c(c30_num, transcript_num)
}
c30_data <- as.data.frame(cbind(c30_raw, c30_num))
c30_data$c30_num <- as.numeric(c30_data$c30_num)
tiff(file = "/home/harry/Documents/imperial/1_11_panoct/c30.tiff",
     unit = "in",
     res=600,
     height = 5,
     width = 7)
p<-ggplot(data = c30_data, aes(x=c30_num)) +
  geom_histogram(color="black", fill="white", binwidth=50) +
  xlab("Gene chronology") + ylab("Frequency")

p

dev.off()


## EXTRACT C184-SPECIFIC GENES ###
c184_specific_data <- strain_specific_data[which(strain_specific_data$C184 == 1),]
c184_centroids <- rownames(c184_specific_data)
c184_raw <- raw_data[c184_centroids, "C184"]

c184_num <- c()
for(i in c184_raw){
  transcript_id <- str_split(i, "_")[[1]][2]
  transcript_num <- as.numeric(transcript_id)
  c184_num <- c(c184_num, transcript_num)
}
c184_data <- as.data.frame(cbind(c184_raw, c184_num))
c184_data$c184_num <- as.numeric(c184_data$c184_num)
tiff(file = "/home/harry/Documents/imperial/1_11_panoct/c184.tiff",
     unit = "in",
     res=600,
     height = 5,
     width = 7)
p<-ggplot(data = c184_data, aes(x=c184_num)) +
  geom_histogram(color="black", fill="white", binwidth=50) +
  xlab("Gene chronology") + ylab("Frequency")

p

dev.off()
## GENERATE PANGENOME LINE GRAPH (HEAP) ##

# Generate rarefaction data
line_data <- rarefaction(t(binary_data), n.perm = 1000)
heap_data <- heaps(t(binary_data), n.perm = 1000)
heap=signif(heap_data[2],2)
# Identify the median values across the permutations
perm_data <- line_data[2:nrow(line_data),2:ncol(line_data)]
median_perm <- round(apply(perm_data, 1, median, na.rm=T))
median_data <- as.data.frame(cbind(1:length(median_perm), median_perm))

# Basic line plot with points
tiff(file = "/home/harry/Documents/imperial/1_11_panoct/curve.tiff",
     unit = "in",
     res=600,
     height = 5,
     width = 10)

p <- ggplot(data=median_data, aes(x=V1, y=median_perm)) +
  geom_line() +
  xlim(0, nrow(median_data)) +
  ylim(0, max(median_perm)) +
  xlab("Number of genomes") +
  ylab("Number of gene clusters") +
  annotate("text", x=nrow(median_data)*0.66, y=10, label= paste(c("gamma=", heap), sep = "", collapse ="")) 
print(p)

dev.off()


## GENERATE OVERVIEW HISTOGRAM ##
total <- max(frequency)
highest_frequency <- sort(freq_table,decreasing=TRUE)[1]
freq_table <- table(frequency)
frequency_df <- as.data.frame(frequency)
tiff(file = "/home/harry/Documents/imperial/1_11_panoct/histogram.tiff",
     unit = "in",
     res=600,
     height = 5,
     width = 7)
# Change colors
p<-ggplot(frequency_df, aes(x=frequency)) + 
  geom_histogram(color="black", fill="white", binwidth=1) +
  xlab("Number of genomes") + ylab("Number of gene clusters")
p

dev.off()

## GENERATE SPECIFIC GENE HISTOGRAM ##
strain_specific_frequency_df <- as.data.frame(strain_specific_freq)
strain_specific_frequency_df <- as.data.frame(cbind(rownames(strain_specific_frequency_df), as.numeric(strain_specific_frequency_df$strain_specific_freq)))
tiff(file = "/home/harry/Documents/imperial/1_11_panoct/specific_histogram.tiff",
     unit = "in",
     res=600,
     height = 5,
     width = 7)
# Change colors
strain_specific_frequency_df$V2 <- as.numeric(strain_specific_frequency_df$V2)
p<-ggplot(data = strain_specific_frequency_df, aes(x=V2)) +
  geom_histogram(color="black", fill="white", binwidth=1) +
  labs(x = "Number of strain-specific genes", y="Number of strains")
p 

dev.off()

## SUBSET BASED ON METADATA ##
# Make function to subset things that are  not in a vector
`%nin%` = Negate(`%in%`)
# What unique known resistant mutations are there?
amr_types <- unique(clean_meta$AMR)
t_binary <- t(binary_data)
clade_A_data <- t_binary[clean_meta$Pangenome.id[which(clean_meta$Clade == "A")],]
clade_B_data <- t_binary[clean_meta$Pangenome.id[which(clean_meta$Clade == "B")],]

A_names <- names(frequency)[which(colSums(clade_A_data) == nrow(clade_A_data))]
B_names <- names(frequency)[which(colSums(clade_B_data) == nrow(clade_B_data))]

A_specific  <- A_names[which(A_names %nin% B_names)]
B_specific  <- B_names[which(B_names %nin% A_names)]

write(A_specific, "/home/harry/Documents/imperial/1_11_panoct/clade_A_specific.txt")
write(B_specific, "/home/harry/Documents/imperial/1_11_panoct/clade_B_specific.txt" )

# Itraconazole resistance focus
resistant_data <- t_binary[clean_meta$Pangenome.id[which(clean_meta$itr == 1)],]
susceptible_data <- t_binary[clean_meta$Pangenome.id[which(clean_meta$itr == 0)],]

r_names <- names(frequency)[which(colSums(resistant_data) == nrow(resistant_data))]
s_names <- names(frequency)[which(colSums(susceptible_data) == nrow(susceptible_data))]

r_specific  <- r_names[which(r_names %nin% s_names)]
s_specific  <- s_names[which(s_names %nin% r_names)]

write(r_specific, "/home/harry/Documents/imperial/1_11_panoct/itra_res_specific.txt")
write(s_specific, "/home/harry/Documents/imperial/1_11_panoct/itra_sus_specific.txt" )



# Itra resistance but unknown mechanism
known_r_data <- t_binary[clean_meta$Pangenome.id[which(clean_meta$itr == 1 & clean_meta$AMR != "Wildtype")],]
unknown_r_data <- t_binary[clean_meta$Pangenome.id[which(clean_meta$itr == 1 & clean_meta$AMR == "Wildtype")],]

known_r_specific <- names(frequency)[which(colSums(known_r_data) == nrow(known_r_data))]
unknown_r_specific <- names(frequency)[which(colSums(unknown_r_data) == nrow(unknown_r_data))]

known_specific  <- known_r_specific[which(known_r_specific %nin% unknown_r_specific)]
unknown_specific  <- unknown_r_specific[which(unknown_r_specific %nin% known_r_specific)]

write(known_specific, "/home/harry/Documents/imperial/1_11_panoct/itra_known_amr_specific.txt")
write(unknown_specific, "/home/harry/Documents/imperial/1_11_panoct/itra_unknown_amr_specific.txt" )

