library(DESeq2)
library(tximport)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files <- list.files(pattern="*.sf")
names(files) <- lapply(files, function(x) paste(strsplit(x, "_")[[1]][3:4], sep="", collapse="_"))

gene_id <- read.table(files[1], header = T)[,1]
tx2gene <- data.frame("GENEID"=gene_id, "TXNAME"=gene_id)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

condition_order <- unlist(lapply(files, function(x) paste(strsplit(x, "_")[[1]][3], sep="", collapse="_")))
sampleTable <- data.frame(condition = factor(condition_order))
rownames(sampleTable) <- colnames(txi.salmon$counts)

ddsSalmon <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
dds <- DESeq(ddsSalmon)

results(dds)