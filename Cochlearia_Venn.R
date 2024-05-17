# Install packages
install.packages("ggVennDiagram")
install.packages("sf")
install.packages('veccompare')

# Load libraries
library("ggVennDiagram")
library("sf")
library("ggplot2")
library("veccompare")

# Open the files
# Gene lists
coch_genes <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/thalianaIDs_top_1perc_minus_inversion.txt",
                       header = FALSE)
thal_genes <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top1perc_Fst_thaliana.txt",
                       header = FALSE)
amar_genes <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top1perc_Fst_amara.txt",
                       header = FALSE)

# GO lists
coch_GO <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/GO_terms_cochlearia_0.05.txt",
                    header = FALSE)
thal_GO <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/GO_terms_thaliana_0.05.txt",
                    header = FALSE)
amar_GO <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/GO_terms_amara_0.05.txt",
                    header = FALSE)

# Orthogroups lists
coch_ortho <-read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/orthogroups_in_top1perc_cochlearia_NR.txt",
                      header = FALSE)
thal_ortho <-read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/orthogroups_in_top1perc_arenosa_NR.txt",
                      header = FALSE)
amar_ortho <-read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/orthogroups_in_top1perc_amara_NR.txt",
                      header = FALSE)

# Convert to vector
coch_genes_vector <- coch_genes$V1
thal_genes_vector <- thal_genes$V1
amar_genes_vector <- amar_genes$V1

coch_GO_vector <- coch_GO$V1
thal_GO_vector <- thal_GO$V1
amar_GO_vector <- amar_GO$V1

coch_ortho_vector <- coch_ortho$V1
thal_ortho_vector <- thal_ortho$V1
amar_ortho_vector <- amar_ortho$V1

# Make a list and name headers
genes_list <- list("Cochlearia" = coch_genes_vector, "Arabidopsis" = thal_genes_vector, "Cardamine" = amar_genes_vector)
GO_list <- list("Cochlearia" = coch_GO_vector, "Arabidopsis" = thal_GO_vector, "Cardamine" = amar_GO_vector)
ortho_list <- list("Cochlearia" = coch_ortho_vector, "Arabidopsis" = thal_ortho_vector, "Cardamine" = amar_ortho_vector)

# Make Venn diagrams

ggVennDiagram(genes_list,
              label_alpha = 0) +
              scale_fill_gradient(low="white",high = "green") +
              scale_color_manual(values = c("black", "black","black")) +
              ggtitle("Gene Overlaps") +
              theme(plot.title = element_text(hjust = 0.5))

ggVennDiagram(GO_list,
              label_alpha = 0) +
              scale_fill_gradient(low="white",high = "green") +
              scale_color_manual(values = c("black", "black","black"))+
              ggtitle("GO Overlaps") +
              theme(plot.title = element_text(hjust = 0.5))

ggVennDiagram(ortho_list,
              label_alpha = 0) +
              scale_fill_gradient(low="white",high = "green") +
              scale_color_manual(values = c("black", "black","black"))+
              ggtitle("Orthogroup Overlaps") +
              theme(plot.title = element_text(hjust = 0.5))

# What genes are in the overlaps?
gene_comparison <- compare.vectors(genes_list)
gene_comparison_MD <-compare.vectors.and.return.text.analysis.of.overlap(genes_list)

# Write gene output to files
sink("gene_comparison_MD.txt")
cat(gene_comparison_MD)
sink()

# What GOs are in the overlaps?
GO_comparison <- compare.vectors(GO_list)
GO_comparison_MD <-compare.vectors.and.return.text.analysis.of.overlap(GO_list)

# Write gene output to files
sink("GO_comparison_MD.txt")
cat(GO_comparison_MD)
sink()

# What orthogroups are in the overlaps?
ortho_comparison <- compare.vectors(ortho_list)
ortho_comparison_MD <-compare.vectors.and.return.text.analysis.of.overlap(ortho_list)

# Write orthogroup output to files
sink("orthogroup_comparison_MD.txt")
cat(ortho_comparison_MD)
sink()




# Testing commands
# https://www.datanovia.com/en/blog/beautiful-ggplot-venn-diagram-with-r/
set.seed(20190708)
genes <- paste("gene",1:1000,sep="")

# x is a list, the names should be the species (A-D here) and the list contents should be the genes/GO terms
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)

# Default plot
ggVennDiagram(x)

ggVennDiagram(
  x, label_alpha = 0,
  category.names = c("Stage 1","Stage 2","Stage 3", "Stage4")
) +
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")
