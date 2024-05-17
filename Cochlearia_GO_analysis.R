# Install packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")

# Install TopGO
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")

# load topGO
library("topGO")

# Open my files

# This function converts the file to the gene2GO format that topGO requires
go <- readMappings(file="/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_GO_universe_restrictive.tsv", 
                   sep = "\t", IDsep = ",")

# All the genes in my scan (i.e. all the genes in regions from the bpm file)
genesInScan <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_in_scan_clean.txt", 
                        header = FALSE)

# My genes of interest (i.e. Fst top 1%)
top1perc <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top_1perc_minus_inversion.txt", 
                     header = FALSE)

# Create a factor where the indicies are the gene names and the factors are 1/0 for "candidate gene"/"scan gene but not candidate gene"

# First convert top1perc from a dataframe to a vector
top1perc_vector <- as.vector(top1perc[,'V1'])

# Then convert genesInScan to a vector
genesInScan_vector <- as.vector(genesInScan[,'V1'])

# Then make the factor
allGenes_factor <- factor(as.integer(genesInScan_vector %in% top1perc_vector))
names(allGenes_factor) = genesInScan_vector

# Website used for help: https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
# Make the topGOdata object
go_data_BP <- new("topGOdata", 
               ontology="BP", 
               allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan
               annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"
               gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms
               nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene

go_data_CC <- new("topGOdata", 
               ontology="CC", 
               allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan
               annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"
               gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms
               nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene

go_data_MF <- new("topGOdata", 
               ontology="MF", 
               allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan
               annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"
               gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms
               nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene


# Run the GO analysis classic fisher
res_BP_classicFisher <- runTest(go_data_BP, statistic = "fisher")
res_CC_classicFisher <- runTest(go_data_CC, statistic = "fisher")
res_MF_classicFisher <- runTest(go_data_MF, statistic = "fisher")

# Run the GO analysis conserative elim fisher
resBP_elimFisher <- runTest(go_data_BP, algorithm = "elim", statistic = "fisher")
resCC_elimFisher <- runTest(go_data_CC, algorithm = "elim", statistic = "fisher")
resMF_elimFisher <- runTest(go_data_MF, algorithm = "elim", statistic = "fisher")

# Check it
res_BP_classicFisher
res_CC_classicFisher
res_MF_classicFisher

# Check it
resBP_elimFisher
resCC_elimFisher
resMF_elimFisher

# Show a table classic
GenTable(go_data_BP, classicFisher=res_BP_classicFisher, topNodes=51) # 17 is the number of sig at 0.01, 51 at 0.05
GenTable(go_data_CC, classicFisher=res_CC_classicFisher, topNodes=13) # 4 is the number of sig at 0.01, 13 at 0.05
GenTable(go_data_MF, classicFisher=res_MF_classicFisher, topNodes=12) # 2 is the number of sig at 0.01, 12 at 0.05

# Show a table elim
GenTable(go_data_BP, elimFisher=resBP_elimFisher, topNodes=113) # 30 is the number of sig at 0.01, 113 at 0.05
GenTable(go_data_CC, elimFisher=resCC_elimFisher, topNodes=32) # 3 is the number of sig at 0.01, 32 at 0.05
GenTable(go_data_MF, elimFisher=resMF_elimFisher, topNodes=36) # 6 is the number of sig at 0.01, 36

# Show combined table for 0.01 sig cutoff
allResBP <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,
                   elimFisher = resBP_elimFisher,
                   orderBy = "elimFisher",
                   topNodes = 30,
                   numChar = 250)

allResCC <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,
                     elimFisher = resCC_elimFisher,
                     orderBy = "elimFisher",
                     topNodes = 3,
                     numChar = 250)

allResMF <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,
                     elimFisher = resMF_elimFisher,
                     orderBy = "elimFisher",
                     topNodes = 6,
                     numChar = 250)

# Show combined table for 0.05 sig cutoff
allResBP_0.05 <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,
                     elimFisher = resBP_elimFisher,
                     orderBy = "elimFisher",
                     topNodes = 113,
                     numChar = 250)

allResCC_0.05 <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,
                     elimFisher = resCC_elimFisher,
                     orderBy = "elimFisher",
                     topNodes = 32,
                     numChar = 250)

allResMF_0.05 <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,
                     elimFisher = resMF_elimFisher,
                     orderBy = "elimFisher",
                     topNodes = 36,
                     numChar = 250)

# since the elim gives me better results and it was the method used in ...
# ... Madjas paper I will go ahead with this method and not classic

# I think one of the previous problems might be the MTC, the top GO manual ...
# ... states that "in many cases adjusted p-values might be misleading".

# Write the combined tables to files for 0.01 cutoff
write.csv(allResBP, "allResBP.csv")
write.csv(allResCC, "allResCC.csv")
write.csv(allResMF, "allResMF.csv")

# Write the combined tables to files for 0.05 cutoff
write.csv(allResBP_0.05, "allResBP.csv")
write.csv(allResCC_0.05, "allResCC.csv")
write.csv(allResMF_0.05, "allResMF.csv")

# Add the genes to these tables:
# Example from https://support.bioconductor.org/p/29775/

# creates a factor (1's and 0's for true and false)
universe = factor(as.integer(genesInScan_vector %in% top1perc_vector))

# gives the factor the names of the genes
names(universe) = genesInScan_vector

# Creates the opposite list of go, i.e. GOtoGenes not GenesToGO
reverse_go_BP = genesInTerm(go_data_BP)
reverse_go_CC = genesInTerm(go_data_CC)
reverse_go_MF = genesInTerm(go_data_MF)

# strips away the genes from reverse_go that are not in top1perc
go_with_only_1perc_BP = lapply(reverse_go_BP,function(x) x[x %in% top1perc_vector] )
go_with_only_1perc_CC = lapply(reverse_go_CC,function(x) x[x %in% top1perc_vector] )
go_with_only_1perc_MF = lapply(reverse_go_MF,function(x) x[x %in% top1perc_vector] )

# Can now retrieve from go_with_only_1perc using the GO terms as the index ...
# .. and get back only genes in my 1perc list

# Need a vector that is the GO terms in my results, for the 0.01
sig_GOs_BP <- allResBP$GO.ID
sig_GOs_CC <- allResCC$GO.ID
sig_GOs_MF <- allResMF$GO.ID

# Need a vector that is the GO terms in my results, for the 0.05
sig_GOs_BP_0.05 <- allResBP_0.05$GO.ID
sig_GOs_CC_0.05 <- allResCC_0.05$GO.ID
sig_GOs_MF_0.05 <- allResMF_0.05$GO.ID

# get vector to add to the dataframe, for the 0.01
sig_genes_BP <- lapply(sig_GOs_BP, function(x) go_with_only_1perc_BP[x])
sig_genes_CC <- lapply(sig_GOs_CC, function(x) go_with_only_1perc_CC[x])
sig_genes_MF <- lapply(sig_GOs_MF, function(x) go_with_only_1perc_MF[x])

# get vector to add to the dataframe, for the 0.05
sig_genes_BP_0.05 <- lapply(sig_GOs_BP_0.05, function(x) go_with_only_1perc_BP[x])
sig_genes_CC_0.05 <- lapply(sig_GOs_CC_0.05, function(x) go_with_only_1perc_CC[x])
sig_genes_MF_0.05 <- lapply(sig_GOs_MF_0.05, function(x) go_with_only_1perc_MF[x])

# Convert the vectors in the list to single strings, for the 0.01
sig_genes_for_tabel_BP <- c()
for(GO in sig_genes_BP){
  temp_gene <- paste(unlist(GO), collapse = ',')
  temp_vec <- c(sig_genes_for_tabel_BP,temp_gene)
  sig_genes_for_tabel_BP <- temp_vec
}

sig_genes_for_tabel_CC <- c()
for(GO in sig_genes_CC){
  temp_gene <- paste(unlist(GO), collapse = ',')
  temp_vec <- c(sig_genes_for_tabel_CC,temp_gene)
  sig_genes_for_tabel_CC <- temp_vec
}

sig_genes_for_tabel_MF <- c()
for(GO in sig_genes_MF){
  temp_gene <- paste(unlist(GO), collapse = ',')
  temp_vec <- c(sig_genes_for_tabel_MF,temp_gene)
  sig_genes_for_tabel_MF <- temp_vec
}

# Convert the vectors in the list to single strings, for the 0.05
sig_genes_for_tabel_BP_0.05 <- c()
for(GO in sig_genes_BP_0.05){
  temp_gene <- paste(unlist(GO), collapse = ',')
  temp_vec <- c(sig_genes_for_tabel_BP_0.05,temp_gene)
  sig_genes_for_tabel_BP_0.05 <- temp_vec
}

sig_genes_for_tabel_CC_0.05 <- c()
for(GO in sig_genes_CC_0.05){
  temp_gene <- paste(unlist(GO), collapse = ',')
  temp_vec <- c(sig_genes_for_tabel_CC_0.05,temp_gene)
  sig_genes_for_tabel_CC_0.05 <- temp_vec
}

sig_genes_for_tabel_MF_0.05 <- c()
for(GO in sig_genes_MF_0.05){
  temp_gene <- paste(unlist(GO), collapse = ',')
  temp_vec <- c(sig_genes_for_tabel_MF_0.05,temp_gene)
  sig_genes_for_tabel_MF_0.05 <- temp_vec
}

# Add the genes to the dataframe, for the 0.01
allResBP_genes <-(cbind(allResBP, Genes=sig_genes_for_tabel_BP))
allResCC_genes <-(cbind(allResCC, Genes=sig_genes_for_tabel_CC))
allResMF_genes <-(cbind(allResMF, Genes=sig_genes_for_tabel_MF))

# Add the genes to the dataframe, for the 0.05
allResBP_genes_0.05 <-(cbind(allResBP_0.05, Genes=sig_genes_for_tabel_BP_0.05))
allResCC_genes_0.05 <-(cbind(allResCC_0.05, Genes=sig_genes_for_tabel_CC_0.05))
allResMF_genes_0.05 <-(cbind(allResMF_0.05, Genes=sig_genes_for_tabel_MF_0.05))

# Write to file for the 0.01
write.csv(allResBP_genes, "allResBP_genes.csv")
write.csv(allResCC_genes, "allResCC_genes.csv")
write.csv(allResMF_genes, "allResMF_genes.csv")

# Write to file for the 0.05
write.csv(allResBP_genes_0.05, "allResBP_genes_0.05.csv")
write.csv(allResCC_genes_0.05, "allResCC_genes_0.05.csv")
write.csv(allResMF_genes_0.05, "allResMF_genes_0.05.csv")

# Try out a few more ways to vizualise the GOs

# Make higher level interface with runTest
resultBP <- runTest(go_data_BP, algorithm = "elim", statistic = "fisher")
resultCC <- runTest(go_data_CC, algorithm = "elim", statistic = "fisher")
resultMF <- runTest(go_data_MF, algorithm = "elim", statistic = "fisher")

# Histogram of pvalues
pvalBP <- score(resultBP)
hist(pvalBP, 50, xlab = "p-values")

pvalCC <- score(resultCC)
hist(pvalCC, 50, xlab = "p-values")

pvalMF <- score(resultMF)
hist(pvalMF, 50, xlab = "p-values")

# Visulise GO structure for the 0.01
# BP has 30 significant GOs at 0.01
showSigOfNodes(go_data_BP, score(resultBP), firstSigNodes = 30, useInfo = 'all')
printGraph(go_data_BP, resultBP, firstSigNodes = 30, fn.prefix = "BP", useInfo = "all", pdfSW = TRUE)

# CC has 3 significant GOs at 0.01
showSigOfNodes(go_data_CC, score(resultCC), firstSigNodes = 3, useInfo = 'all')
printGraph(go_data_CC, resultCC, firstSigNodes = 3, fn.prefix = "CC", useInfo = "all", pdfSW = TRUE)

# MF has 6 significant GOs at 0.01
showSigOfNodes(go_data_MF, score(resultMF), firstSigNodes = 6, useInfo = 'all')
printGraph(go_data_MF, resultMF, firstSigNodes = 6, fn.prefix = "MF", useInfo = "all", pdfSW = TRUE)

# Visulise GO structure for the 0.05
# BP has 113 significant GOs at 0.05
showSigOfNodes(go_data_BP, score(resultBP), firstSigNodes = 113, useInfo = 'all')
printGraph(go_data_BP, resultBP, firstSigNodes = 113, fn.prefix = "BP", useInfo = "all", pdfSW = TRUE)

# CC has 32 significant GOs at 0.05
showSigOfNodes(go_data_CC, score(resultCC), firstSigNodes = 32, useInfo = 'all')
printGraph(go_data_CC, resultCC, firstSigNodes = 32, fn.prefix = "CC", useInfo = "all", pdfSW = TRUE)

# MF has 36 significant GOs at 0.05
showSigOfNodes(go_data_MF, score(resultMF), firstSigNodes = 36, useInfo = 'all')
printGraph(go_data_MF, resultMF, firstSigNodes = 36, fn.prefix = "MF", useInfo = "all", pdfSW = TRUE)
