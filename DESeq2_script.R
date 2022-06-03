################ DESeq2 ################
## Script by Lucas Miguel de Carvalho ##
######## UNICAMP - 2020 ################
########################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)

#Read counts table
data = read.table("Counts.matrix",header=T, row.names=1)
head(data)


col_ordering = c(1,2,3,7,8,9) #selecting columns  
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix) >= 2,] #filtering
head(rnaseqMatrix)

#Factor names. Equals to number of columns selected (line 18)
conditions = data.frame(conditions=factor(c(rep("Control", 3), rep("Treated", 3))))
rownames(conditions) = colnames(rnaseqMatrix)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)

dds = DESeq(ddsFullCountTable)
#constrast: Treated over Control means that a positive log2FC indicates more
#            reads in Treated. Check the ratio between the the baseMeanB/baseMeanA
contrast=c("conditions","Treated","Control")
res = results(dds, contrast)

baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Control"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Treated"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Control", sampleB="Treated", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

head(res)

#plot read counts boxplot for a specific gene
library(ggplot2)

#get normalized counts
g1_counts <- log2(counts(dds['YML100W',], normalized = TRUE))
m <- list(counts = as.numeric(g1_counts), 
          group = as.factor(conditions$conditions))

df <- data.frame(matrix(unlist(m), nrow=length(m), byrow=TRUE))
df <- as.data.frame(t(df))

colnames(df) <- c("normalized_counts","groups")
df$groups <- c(rep("Control",3),rep("Treated",3))
df$groups <- factor(df$groups, levels=unique(df$groups))

ggplot(df, aes(groups, normalized_counts)) + 
  geom_boxplot(aes(fill=groups)) + geom_jitter(width = 0.1) +
  ylab("Expression(log2 normalized counts)") 

# Volcano plot
pdf("DESeq2_Volcano.pdf")
plot(res$log2FoldChange, -1*log10(res$padj), col=ifelse(res$padj<=0.05, "red", "black"),xlab="logCounts", ylab="logFC", title="Volcano plot", pch=20)
dev.off()

#### Differential expression results ####
# Lets filter Differentially expressed genes (DEGs)
# Filtering conditions:  padj <= 0.05 AND |log2FC| > 1.5
# Remember: |log2FC| > 1.5  is equals to log2FC > 1.5 or log2FC < -1.5

library(dplyr)

#res_filtered <- res %>% select(c("log2FoldChange","padj"))
res_filtered <- res %>% filter(padj <= 0.05 & 
                                 (log2FoldChange > 1.5 | log2FoldChange < -1.5)
                               )
res_filtered

write.table(res_filtered,file="deseq2_results_DEGs.txt")