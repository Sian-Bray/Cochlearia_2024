# Working!!!!!
# install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::install("topGO")
# Clostridioides difficile

library(topGO)
library("biomaRt")
library("data.table")

geneList <-fread("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/02_Visualisation_Folder/UK_1perc.txt", header=FALSE)


all_genes <-fread("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/bpm_genes.txt", header=FALSE)


geneID2GO <- readMappings(file = "/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_GO_universe.tsv", sep = "\t", IDsep = ",")
geneID2GO <- readMappings(file = "/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/mini_test_universe.txt", sep = "\t", IDsep = ",")

## #collect gene names from biomart
## mart <- biomaRt::useMart(biomart = "plants_mart",
##                          dataset = "athaliana_eg_gene",
##                          host = 'plants.ensembl.org')
## 
## # Get ensembl gene ids and GO terms
## GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
##                                      "go_id"), mart = mart)
## 
## #Remove blank entries
## GTOGO <- GTOGO[GTOGO$go_id != '',]
## 
## # convert from table format to list format
## geneID2GO <- by(GTOGO$go_id,
##                 GTOGO$ensembl_gene_id,
##                 function(x) as.character(x))
## 
## # Generate gene list for the GO object
## all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
## attempt<-unlist(geneList, recursive = TRUE, use.names = FALSE)
## int.genes <- factor(as.integer(all.genes %in% attempt))
## names(int.genes) = all.genes

# make the GO object
go.obj <- new("topGOdata", ontology='BP',
                allGenes = all_genes, #originally int.genes
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO,
                geneSel = geneList
                 )

# run topGO and make a file of results
resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 1000)
sigRes<-subset(allRes, allRes$elimFisher <=0.05)
write.table(file="/Users/sian_bray/Dropbox/Salt/6_GO_Analysis/GO_Analysis_New/bio_genes_genic_topgo.txt",sigRes,sep="\t",row.names=F)

# make a file that matches genes to GO terms
allGO = genesInTerm(go.obj)
sigGO <- subset(allGO, names(allGO) %in% sigRes$GO.ID)
options(max.print=6100)
capture.output(sigGO, file = "/Users/sian_bray/Dropbox/Salt/6_GO_Analysis/GO_Analysis_New/sig_GO.txt")








### Remake and re-run for mol ###

# make the GO object
go.obj <- new("topGOdata", ontology='MF'
                 , allGenes = int.genes
                 , annot = annFUN.gene2GO
                 , gene2GO = geneID2GO
                 )

# run topGO and make a file of results
resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 1000)
sigRes<-subset(allRes, allRes$elimFisher <=0.05)
write.table(file="/Users/sian_bray/Dropbox/Salt/6_GO_Analysis/GO_Analysis_New/mol_genes_genic_topgo.txt",sigRes,sep="\t",row.names=F)

# make a file that matches genes to GO terms
sigID <- sigRes["GO.ID"]
sigGenes<-subset(GTOGO, GTOGO$go_id %in% sigID$GO.ID)
write.table(file="/Users/sian_bray/Dropbox/Salt/6_GO_Analysis/GO_Analysis_New/mol_sig_genes.tsv",sigGenes,sep="\t",row.names=F)

### Remake and re-run for cel ###

# make the GO object
go.obj <- new("topGOdata", ontology='CC'
                 , allGenes = int.genes
                 , annot = annFUN.gene2GO
                 , gene2GO = geneID2GO
                 )

# run topGO and make a file of results
resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 1000)
sigRes<-subset(allRes, allRes$elimFisher <=0.05)
write.table(file="/Users/sian_bray/Dropbox/Salt/6_GO_Analysis/GO_Analysis_New/cel_genes_genic_topgo.txt",sigRes,sep="\t",row.names=F)

# make a file that matches genes to GO terms
sigID <- sigRes["GO.ID"]
sigGenes<-subset(GTOGO, GTOGO$go_id %in% sigID$GO.ID)
write.table(file="/Users/sian_bray/Dropbox/Salt/6_GO_Analysis/GO_Analysis_New/cel_sig_genes.tsv",sigGenes,sep="\t",row.names=F)







