#install.packages("WGCNA")
library(WGCNA)
#options(stringsAsFactors = FALSE)

setwd("C:/Users/Lucas/WGCNA/")

expression = read.table("expression.txt", header=TRUE,row.names = 1)

#Select some columns
Exp=t(expression[,c(1:17)])

gsg=goodSamplesGenes(Exp, verbose = 3)
gsg$allOK

#if is false:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(Exp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(Exp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  Exp = Exp[gsg$goodSamples, gsg$goodGenes]
}

sampleTree=hclust(dist(Exp), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


######Construction and modules detection ####

powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(Exp, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence (IN)"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.82,col="blue")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity (IN)"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#selected manually
power_value = 12

#matriz de adjacencia
adjacency = adjacency(Exp, power = power_value, type="signed")

#TOM matrix
TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity (LS)",
     labels = FALSE, hang = 0.04)

minModuleSize = 30

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)

table(dynamicColors)

sizeGrWindow(8,6)
png(filename = "Gene_dendrogram_and_module_colors.png", width = 972, height = 600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors (IN)")
dev.off()

#merge
MEList = moduleEigengenes(Exp, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.10
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(Exp, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)

png(filename = "Gene_dendrogram_and_module_colors_merged.png", width = 972, height = 600)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main="Cluster dendrogram")
dev.off()

table(merge$colors)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#save.image("WGCNA_no_traits.RData")

#### Traits #####

##metabolites data##
traits=read.table("traits_conditions.txt",sep="\t", header = TRUE, row.names = 1)
colnames(traits) <- c("1G","2G") #colnames 

prodcolours=numbers2colors(traits, signed = TRUE)

plotDendroAndColors(sampleTree, prodcolours,
                    groupLabels = names(traits),
                    main= "Sample dendogram and tissues heatmap")

##### Association #####    
ngenes=ncol(Exp)
nsamples=nrow(Exp)

MEs0 = moduleEigengenes(Exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleProdCor = cor(MEs, traits, use = "p");
moduleProdPvalue = corPvalueStudent(moduleProdCor, nsamples);

sizeGrWindow(15,6)
textMatrix = paste(signif(moduleProdCor, 2), " (",signif(moduleProdPvalue, 1), ")", sep = "")
print(textMatrix)

write.table(textMatrix,"Matrix_traits.txt")

dim(textMatrix) = dim(moduleProdCor)

png(filename = "Module-sample_relationships.png", width = 8, height = 14,units='in',res=600)
labeledHeatmap(Matrix = moduleProdCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 1,cex.lab.y = 0.9,cex.lab.x = 0.9,
               zlim = c(-1,1),
               main = paste("Module-sample relationships (IN)"))
dev.off()


##### This is for modular membership scores #####
datKME = signedKME(Exp, MEs, outputColumnName = "kME")
write.csv(datKME, file="datKME.csv")

##### This is for gene significance scores #####
#select this by your colname
data_1G = as.data.frame(traits$`1G`)
GS.1G = as.numeric(cor(Exp, data_1G, use = "p"))
write.csv(GS.1G, file="GS.1G.csv")

data_2G = as.data.frame(traits$`2G`)
GS.2G = as.numeric(cor(Exp, data_2G, use = "p"))
write.csv(GS.2G, file="GS.2G.csv")

###genes-module##
probes=colnames(Exp)
write.table(as.data.frame(cbind(probes, mergedColors)),
            file = 'gene_modules_1G_2G.tsv', col.names = c("Gene_ID", "module"),
            row.names = FALSE, sep = "\t", quote = FALSE)



##### CYTOSCAPE ######

table(mergedColors)

#Only key modules in heatmap (r>0.8 e p<=0.05)

modules_total = c("blue","yellow","lightgreen","orange",
                  "cyan","darkred","darkturquoise","purple",
                  "brownm","turquoise","magenta","black")

probes=colnames(Exp)

if(FALSE){
 
  for (modules in modules_total){
    inModule=is.finite(match(moduleColors, modules))
    modProbes=probes[inModule]
   
    modTOM=TOM[inModule, inModule]
    dimnames(modTOM)=list(modProbes, modProbes)
   
    t=(0.5+0.5*0.8)^power_value
    cyt=exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("cytoscape/CytoscapeInput-edges-", paste(modules,collapse = "-"), ".txt", sep = ""),
                                 nodeFile = paste("cytoscape/CytoscapeInput-nodes-", paste(modules,collapse = "-"), ".txt", sep = ""),
                                 weighted = TRUE,
                                 threshold = t,
                                 nodeNames=modProbes,
                                 nodeAttr=moduleColors[inModule])
  }
   
}


png("sampleTree.png",width=5,height=6,units='in',res=600)

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()


png(filename = "Gene dendrogram and module colors.png", width = 5, height = 8,units='in',res=600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


#save.image("WGCNA_traits.RData")
