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

fullpath <- function(suffix){
  out <- paste(c(out_dir, "/", suffix), sep="", collapse="")
  return(out)
}

setwd(in_dir)
raw_data <- read.table("matchtable.txt")[,-1]
raw_centroids <- readLines("centroids.fasta")
raw_meta <- read.csv("raw_metadata.csv")
#essential_centroids <- readLines("essential.centroids")

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
if(file.exists(fullpath("gene_binary.csv"))){
  binary_data <- read.csv(fullpath("gene_binary.csv"), row.names = 1)
  binary_data = as.data.frame(sapply(binary_data, as.numeric))
}else{
  binary_data <- raw_data
  binary_data[binary_data=="----------"]<-0
  binary_data[binary_data!="0"]<-1
  binary_data <- data.frame(sapply(binary_data, function(x) as.numeric(as.character(x))))
  rownames(binary_data) <-  rownames(raw_data)
  write.csv(binary_data, fullpath("gene_binary.csv"))
}
# Subset the core data
core_data <- binary_data[rowSums(binary_data == 0)==0, , drop = FALSE]
# Get frequency for each peptide
frequency <- rowSums(binary_data)
names(frequency) <- centroid_names
if(!file.exists(fullpath("gene_freq.csv"))){
  write.csv(frequency, fullpath("gene_freq.csv"))
}


## CREATE DATA TO MIMIC ROARY OUTPUT ##
# Replace lacking genes with blank space
if(!file.exists(fullpath("mock_roary.csv"))){
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
  mock_roary$No..isolates <- frequency
  write.csv(mock_roary, fullpath("mock_roary.csv"), row.names=F)
}


## CLEAN METADATA ##
# Remove ">" or "<" strings from mic
# order: itr, vor, pos
if(!file.exists(fullpath("clean_meta.csv"))){
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
  
  
  write.csv(clean_meta, fullpath("clean_meta.csv"), row.names=F)
} else{
  clean_meta <- read.csv(fullpath("clean_meta.csv"))
}

## CREATE PHENOTYPE DATA FOR SCOARY INPUT ##
if(file.exists(fullpath("traits_data.csv"))){
  # Re-index metadata to match "mock roary"
  data_blanked <- raw_data
  data_blanked[data_blanked=="----------"]<-""
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
  traits_data <- data.frame(trait1_cladeA_cladeB, trait2_Sitr_Ritr,trait3_RitrK_RitrU,
                             trait4_Svor_Rvor, trait5_RvorK_RvorU, trait6_Spos_Rpos,
                             trait7_RposK_RposU, trait8_Sall_Rall, trait9_RallK_RallU)
  
  rownames(traits_data) <- ordered_meta[,1]
  write.csv(traits_data, fullpath("traits_data.csv"), row.names=T)
}


## CHECK ESSENTIAL GENE PRESENCE ##
essential_freq <- frequency[essential_centroids]


## COUNT NUMBER OF UNIQUE GENES PER STRAIN ##
strain_specific_data <- binary_data[which(frequency==1),]
strain_specific_freq <- colSums(strain_specific_data)
high_specific_freq <- strain_specific_freq[which(strain_specific_freq > 50)]
if(!file.exists(fullpath("strain_specifc_freq.csv"))){
  write.csv(strain_specific_freq, fullpath("strain_specifc_freq.csv"))
}


write(names(which(strain_specific_freq == 0)), fullpath("zero_specific.txt"))

## GENERATE PANGENOME LINE GRAPH (HEAP) ##
if(!file.exists(fullpath("curve.tiff"))){
  # Generate rarefaction data
line_data <- rarefaction(t(binary_data), n.perm = 1000)
heap_data <- heaps(t(binary_data), n.perm = 1000)
heap=signif(heap_data[2],2)
# Identify the median values across the permutations
perm_data <- line_data[2:nrow(line_data),2:ncol(line_data)]
# Get median data
median_perm <- round(apply(perm_data, 1, median, na.rm=T))
median_data <- as.data.frame(cbind(1:length(median_perm), median_perm))
# Convert to long
perm_data$X <- rownames(perm_data)
long_data <- melt(setDT(perm_data), id.vars = c("X"), variable.name = "perm")
tiff(file = fullpath("curve.tiff"),
     unit = "in",
     res=600,
     height = 5,
     width = 10)
# Plot raw data
plot(long_data$X, long_data$value, 
     xlim=c(0, 250), ylim=c(9000, 14000), 
     xaxt='n', yaxt='n', 
     xlab = "Number of genomes",
     ylab = "Number of gene clusters",
     col=rgb(0.4,0.4,0.8,0.6),pch=16 , cex=1.3) 
# Now, define a custom axis
axis(side = 1, at=seq(from = 0, to = 250, by = 50))
axis(side = 2, at=seq(from = 9000, to = 14000, by = 1000))
# Add median line
lines(median_data$V1, median_data$median_perm, col=2, lwd=2 )  
text(200, 9000 , labels=(paste(c( "\u03b3=", heap), sep ="", collapse="")))


dev.off()

}


## GENERATE OVERVIEW HISTOGRAM ##
if(!file.exists(fullpath("histogram.tiff"))){
  total <- max(frequency)
#highest_frequency <- sort(freq_table,decreasing=TRUE)[1]
freq_table <- table(frequency)
frequency_df <- as.data.frame(frequency)
tiff(file = fullpath("histogram.tiff"),
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
}


## GENERATE SPECIFIC GENE HISTOGRAM ##
strain_specific_frequency_df <- as.data.frame(strain_specific_freq)
strain_specific_frequency_df <- as.data.frame(cbind(rownames(strain_specific_frequency_df), as.numeric(strain_specific_frequency_df$strain_specific_freq)))
if(!file.exists(fullpath("specific_histogram.tiff"))){
  tiff(file = fullpath("specific_histogram.tiff"),
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
}

