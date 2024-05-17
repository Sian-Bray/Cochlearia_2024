# load topGO
library("topGO")

# Open my files

# This function converts the file to the gene2GO format that topGO requires
go <- readMappings(file="/Users/sian_bray/Dropbox/Bray/R_Projects/Permutation_Test_Tests/gene_GO.tsv", 
                   sep = "	", IDsep = ",")

# All the genes in my scan (i.e. all the genes in regions from the bpm file)
genesInScan <- read.csv("/Users/sian_bray/Dropbox/Bray/R_Projects/Permutation_Test_Tests/all_genes.csv", 
                        header = FALSE)

# My genes of interest (i.e. Fst top 1%)
top1perc <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/GO_Permutation/GO_Permutation_Test/random_genes_species2_round10.txt", 
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
#res_BP_classicFisher <- runTest(go_data_BP, statistic = "fisher")
#res_CC_classicFisher <- runTest(go_data_CC, statistic = "fisher")
#res_MF_classicFisher <- runTest(go_data_MF, statistic = "fisher")

# Run the GO analysis conserative elim fisher
resBP_elimFisher <- runTest(go_data_BP, algorithm = "elim", statistic = "fisher")
resCC_elimFisher <- runTest(go_data_CC, algorithm = "elim", statistic = "fisher")
resMF_elimFisher <- runTest(go_data_MF, algorithm = "elim", statistic = "fisher")

###########################

# Extract p-values from the results
#p_values_BP_classicFisher <- score(res_BP_classicFisher)
#p_values_CC_classicFisher <- score(res_CC_classicFisher)
#p_values_MF_classicFisher <- score(res_MF_classicFisher)

p_values_BP_elimFisher <- score(resBP_elimFisher)
p_values_CC_elimFisher <- score(resCC_elimFisher)
p_values_MF_elimFisher <- score(resMF_elimFisher)

# Define the significance threshold
threshold <- 0.05

# Count the number of significant genes
#significant_genes_BP_classicFisher <- sum(p_values_BP_classicFisher < threshold)
#significant_genes_CC_classicFisher <- sum(p_values_CC_classicFisher < threshold)
#significant_genes_MF_classicFisher <- sum(p_values_MF_classicFisher < threshold)

significant_genes_BP_elimFisher <- sum(p_values_BP_elimFisher < threshold)
significant_genes_CC_elimFisher <- sum(p_values_CC_elimFisher < threshold)
significant_genes_MF_elimFisher <- sum(p_values_MF_elimFisher < threshold)

###########################

# Find the number of significant nodes at 0.05 cutoff
####

# Show combined table for 0.05 sig cutoff
if (significant_genes_BP_elimFisher == 0) {
  allResBP <- data.frame()
} else if (significant_genes_BP_elimFisher == 1) {
  allResBP <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,
                       elimFisher = resBP_elimFisher,
                       orderBy = "elimFisher",
                       topNodes = 2,
                       numChar = 250)
  allResBP <- allResBP[-nrow(allResBP), ]
} else if (significant_genes_BP_elimFisher > 1) {
  allResBP <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,
                       elimFisher = resBP_elimFisher,
                       orderBy = "elimFisher",
                       topNodes = significant_genes_BP_elimFisher,
                       numChar = 250)
}

if (significant_genes_CC_elimFisher == 0) {
  allResCC <- data.frame()
} else if (significant_genes_CC_elimFisher == 1) {
  allResCC <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,
                       elimFisher = resCC_elimFisher,
                       orderBy = "elimFisher",
                       topNodes = 2,
                       numChar = 250)
  allResCC <- allResCC[-nrow(allResCC), ]
} else if (significant_genes_CC_elimFisher > 1) {
  allResCC <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,
                       elimFisher = resCC_elimFisher,
                       orderBy = "elimFisher",
                       topNodes = significant_genes_CC_elimFisher,
                       numChar = 250)
}

if (significant_genes_MF_elimFisher == 0) {
  allResMF <- data.frame()
} else if (significant_genes_MF_elimFisher == 1) {
  allResMF <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,
                       elimFisher = resMF_elimFisher,
                       orderBy = "elimFisher",
                       topNodes = 2,
                       numChar = 250)
  allResMF <- allResMF[-nrow(allResMF), ]
} else if (significant_genes_MF_elimFisher > 1) {
  allResMF <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,
                       elimFisher = resMF_elimFisher,
                       orderBy = "elimFisher",
                       topNodes = significant_genes_MF_elimFisher,
                       numChar = 250)
}

# since the elim gives me better results and it was the method used in ...
# ... Madjas paper I will go ahead with this method and not classic

# I think one of the previous problems might be the MTC, the top GO manual ...
# ... states that "in many cases adjusted p-values might be misleading".

# Write the combined tables to files for 0.05 cutoff
write.csv(allResBP, "/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/GO_Permutation/GO_Permutation_Test/allResBP_species_2_10.csv")
write.csv(allResCC, "/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/GO_Permutation/GO_Permutation_Test/allResCC_species_2_10.csv")
write.csv(allResMF, "/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/GO_Permutation/GO_Permutation_Test/allResMF_species_2_10.csv")
