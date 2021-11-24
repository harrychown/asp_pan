# PANOCT HISTOGRAMS
library("tidyverse")
library("micropan")
library("ggplot2")
library("VennDiagram")
# Read config file for input/output directories
directory_configuration <- readLines("directory.config")
in_dir <- directory_configuration[1]
out_dir <- directory_configuration[2]

setwd(in_dir)
raw_data <- read.table("matchtable.txt")[,-1]
raw_centroids <- readLines("centroids.fasta")
raw_meta <- read.csv("raw_metadata.csv")

essential_centroids <- readLines("/home/harry/Documents/imperial/essential.centroids")
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


write.csv(clean_meta, "/home/harry/Documents/imperial/1_11_panoct/clean_meta.csv", row.names=F)

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

