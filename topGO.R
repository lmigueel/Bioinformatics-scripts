# How to run 'topGO' enrichment analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("topGO")

# topGO.R
# eg. run:
# Rscript topGO.R annotations.txt interestinggenes.txt output_file

# You can change the 'args[i]' by the file path in your computer
args <- commandArgs(trailingOnly = TRUE)
universeFile = args[1]
interestingGenesFile = args[2]
output_file = args[3]

# set the output file
sink(output_file)

# load topGO
library("topGO")

# read in the 'gene universe' file
geneID2GO <- readMappings(file = universeFile)
geneUniverse <- names(geneID2GO)

# read in the genes of interest
genesOfInterest <- read.table(interestingGenesFile,header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# build the GOdata object in topGO
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gen
e2GO = geneID2GO)
myGOdata

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")


# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.001)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

# print out the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes

# print a graph (to a pdf file) with the top 'numsignif' results:
output_file2 = paste(output_file,"Topgo", sep="_")
#printGraph(myGOdata, resultTopgo, firstSigNodes = numsignif, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)

# print out the genes that are annotated with the significantly enriched GO terms:
myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms)){
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
sink() 