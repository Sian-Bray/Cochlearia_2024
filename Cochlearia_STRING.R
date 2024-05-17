# Install the packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")

install.packages("GGally")
install.packages("network")
install.packages("sna")
install.packages("igraph")

# Load modules
library(STRINGdb)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library("topGO")
library("dplyr")

# List all methods
STRINGdb$methods()
# see documentation for get_graph
STRINGdb$help("get_graph")

# Make the reference database, contains all thaliana proteins and interations
# 3702 is the species number for Arabidopsis thaliana
string_db <- STRINGdb$new( version="11.5", 
                           species=3702, 
                           score_threshold=200,
                           input_directory="")

# Import the data, this is simply a list of the genes, one gene per line
# AT1G01010, not AT1G01010.1
data <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/bothIDs_top_1perc_minus_inversion.txt",
                 header=FALSE)

data_3spp <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_to_three_species.csv",
                      header=TRUE)

# Import data for the network firgures
# List of gene IDs to species
genes2species <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_to_three_species_nr.csv",
            header=TRUE)

# Data frame of Cochlearia GO terms to genes
GO2genes <- read.csv("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/GO2gene_for_STRING.csv",
                     header=TRUE)

# Mapps STRING IDs to IDs in my list
# STRING IDs add a species number before and a alternate splice number after
# i.e. AT3G09520 becomes 3702.AT3G09520.1
# V1 is the automatic column name when header=FALSE above
data_STRING <- string_db$map(data, "V2", removeUnmappedRows = TRUE )

data_STRING_3spp <- string_db$map(data_3spp, "Gene", removeUnmappedRows = TRUE )

# Take just the STRING IDs
hits <- data_STRING$STRING_id

hits_3spp <- data_STRING_3spp$STRING_id

# plot the STRING interaction network
# Name the pdf to be plotted
pdf("STRING_Cochlearia_genes.pdf")
# What to plot
# can use network_flavor to specify the flavor of the network
# Options are :"evidence", "confidence" or "actions". default is "evidence"
string_db$plot_network( hits , network_flavor="confidence")
# Stop plotting
dev.off()

pdf("STRING_three_species_genes.pdf")
string_db$plot_network( hits_3spp , network_flavor="confidence")
dev.off()

# Output the hits as a csv file
# What to output
out_hits<-string_db$get_interactions(hits)
# Write the csv
write.csv(out_hits,"STRING_Cochlearia_genes.csv", row.names = FALSE)

out_hits_3spp<-string_db$get_interactions(hits_3spp)
write.csv(out_hits_3spp,"STRING_Cochlearia_genes.csv", row.names = FALSE)

# Get clusters
# In the Cochlearia data there are 5 major clusters (=>16 genes),
# 5 minor clusters (3-4 genes), and 4 pairs (total 14)
clustersList <- string_db$get_clusters(data_STRING)

clustersList_3spp <- string_db$get_clusters(data_STRING_3spp)


# plot clusters - Cocohlearia
# par = a vector of the form c(nr, nc). Subsequent figures will be drawn in an nr-by-nc array on the device by columns (mfcol), or rows (mfrow), respectively
par(mfrow=c(2,2))

# plot cluster 1
pdf("STRING_Cochlearia_cluster_1.pdf")
string_db$plot_network(clustersList[[1]])
dev.off()

# plot cluster 2
pdf("STRING_Cochlearia_cluster_2.pdf")
string_db$plot_network(clustersList[[2]])
dev.off()

# plot cluster 3
pdf("STRING_Cochlearia_cluster_3.pdf")
string_db$plot_network(clustersList[[3]])
dev.off()

# plot cluster 4
pdf("STRING_Cochlearia_cluster_4.pdf")
string_db$plot_network(clustersList[[4]])
dev.off()

# plot cluster 5
pdf("STRING_Cochlearia_cluster_5.pdf")
string_db$plot_network(clustersList[[5]])
dev.off()

# plot cluster 6
pdf("STRING_Cochlearia_cluster_6.pdf")
string_db$plot_network(clustersList[[6]])
dev.off()

# plot cluster 7
pdf("STRING_Cochlearia_cluster_7.pdf")
string_db$plot_network(clustersList[[7]])
dev.off()

# plot cluster 8
pdf("STRING_Cochlearia_cluster_8.pdf")
string_db$plot_network(clustersList[[8]])
dev.off()

# plot cluster 9
pdf("STRING_Cochlearia_cluster_9.pdf")
string_db$plot_network(clustersList[[9]])
dev.off()

# plot cluster 10
pdf("STRING_Cochlearia_cluster_10.pdf")
string_db$plot_network(clustersList[[10]])
dev.off()

# plot cluster 11
pdf("STRING_Cochlearia_cluster_11.pdf")
string_db$plot_network(clustersList[[11]])
dev.off()

# plot cluster 12
pdf("STRING_Cochlearia_cluster_12.pdf")
string_db$plot_network(clustersList[[12]])
dev.off()

# plot cluster 13
pdf("STRING_Cochlearia_cluster_13.pdf")
string_db$plot_network(clustersList[[13]])
dev.off()

# plot cluster 14
pdf("STRING_Cochlearia_cluster_14.pdf")
string_db$plot_network(clustersList[[14]])
dev.off()


# Plot clusters for the three species
# 4 clusters >= 29, 5 clusters == 3, 4 clusters == 2
# Total 13 clusters

# par = a vector of the form c(nr, nc). Subsequent figures will be drawn in an nr-by-nc array on the device by columns (mfcol), or rows (mfrow), respectively
par(mfrow=c(2,2))

# plot cluster 1
pdf("STRING_3species_cluster_1.pdf")
string_db$plot_network(clustersList_3spp[[1]])
dev.off()

# plot cluster 2
pdf("STRING_3species_cluster_2.pdf")
string_db$plot_network(clustersList_3spp[[2]])
dev.off()

# plot cluster 3
pdf("STRING_3species_cluster_3.pdf")
string_db$plot_network(clustersList_3spp[[3]])
dev.off()

# plot cluster 4
pdf("STRING_3species_cluster_4.pdf")
string_db$plot_network(clustersList_3spp[[4]])
dev.off()

# plot cluster 5
pdf("STRING_3species_cluster_5.pdf")
string_db$plot_network(clustersList_3spp[[5]])
dev.off()

# plot cluster 6
pdf("STRING_3species_cluster_6.pdf")
string_db$plot_network(clustersList_3spp[[6]])
dev.off()

# plot cluster 7
pdf("STRING_3species_cluster_7.pdf")
string_db$plot_network(clustersList_3spp[[7]])
dev.off()

# plot cluster 8
pdf("STRING_3species_cluster_8.pdf")
string_db$plot_network(clustersList_3spp[[8]])
dev.off()

# plot cluster 9
pdf("STRING_3species_cluster_9.pdf")
string_db$plot_network(clustersList_3spp[[9]])
dev.off()

# plot cluster 10
pdf("STRING_3species_cluster_10.pdf")
string_db$plot_network(clustersList_3spp[[10]])
dev.off()

# plot cluster 11
pdf("STRING_3species_cluster_11.pdf")
string_db$plot_network(clustersList_3spp[[11]])
dev.off()

# plot cluster 12
pdf("STRING_3species_cluster_12.pdf")
string_db$plot_network(clustersList_3spp[[12]])
dev.off()

# plot cluster 13
pdf("STRING_3species_cluster_13.pdf")
string_db$plot_network(clustersList_3spp[[13]])
dev.off()


# I want to better visualise the nodes, including colour by GO term and shape by original species
# https://briatte.github.io/ggnet/

# Make a dataframe of 2 columns - one interacting gene in each
out_hits_2cols <- out_hits[1:2]
# Remove duplicate columns as network cannot handle them
out_hits_nodups <- out_hits_2cols[!duplicated(out_hits_2cols), ]
# Make the network
net <- network(out_hits_nodups)
# plot the network simply
ggnet2(net)

# Extract outhits for one node only
# vector of genes in the cluster
cluster5 <- clustersList[[5]]
# cluster genes that are in the from column of the interactions
cluster5_from <- out_hits_nodups[out_hits_nodups$from %in% cluster5,]
# cluster genes that are in the to column of the interactions
cluster5_to <- out_hits_nodups[out_hits_nodups$to %in% cluster5,]
# Merge "to" and "from" into one dataframs
cluster5_net <- merge(cluster5_from, cluster5_to)
# simple plot of the network
ggnet2(cluster5_net, label = TRUE)



# For the three species

# input data at the top
# Make a dataframe of 2 columns - one interacting gene in each
out_hits_2cols_3spp <- out_hits_3spp[1:2]
# Remove duplicate columns as network cannot handle them
out_hits_nodups_3spp <- out_hits_2cols_3spp[!duplicated(out_hits_2cols_3spp), ]
# Make the network
net_3spp <- network(out_hits_nodups_3spp)
# plot the network simply
ggnet2(net_3spp)

# Extract outhits for one node only
# vector of genes in the cluster
cluster4_3spp <- clustersList_3spp[[4]]
# cluster genes that are in the from column of the interactions
cluster4_from_3spp <- out_hits_nodups_3spp[out_hits_nodups_3spp$from %in% cluster4_3spp,]
# cluster genes that are in the to column of the interactions
cluster4_to_3spp <- out_hits_nodups_3spp[out_hits_nodups_3spp$to %in% cluster4_3spp,]
# Merge "to" and "from" into one dataframs
cluster4_merge_3spp <- merge(cluster4_from_3spp, cluster4_to_3spp)
# Remove extra bits of the label, i.e. 3702.AT1G34810.1 > AT1G34810
# Remove the last two characters
cluster4_merge_3spp$from <- substr(cluster4_merge_3spp$from,1,nchar(cluster4_merge_3spp$from)-1)
cluster4_merge_3spp$to <- substr(cluster4_merge_3spp$to,1,nchar(cluster4_merge_3spp$to)-1)
# Remove the first 5 chararcters
cluster4_merge_3spp$from <- substr(cluster4_merge_3spp$from,6,nchar(cluster4_merge_3spp$from)-1)
cluster4_merge_3spp$to <- substr(cluster4_merge_3spp$to,6,nchar(cluster4_merge_3spp$to)-1)
# Make network object
cluster4_net_3spp <- network(cluster4_merge_3spp)
# Get a vector of network nodes
cluster4_nodes_3spp <- network.vertex.names(cluster4_net_3spp)

# Create a vector of matching shapes for species
# Need to remake the gene list with the Venn output for overlaps - do this manually
# I want to make a vector of species from genes2species that match cluster4_nodes_3spp

# No idea why the below does not work! Returns the whole dataframe...
# cluster4_nodes_shape_3spp <- genes2species[genes2species$Gene %in% cluster4_nodes_3spp, ]
# cluster4_nodes_shape_3spp <- genes2species[is.element(genes2species$Gene, cluster4_nodes_3spp), ]

# Returns data frame of values only where Gene column data is in the vector cluster4_nodes_3spp
cluster4_nodes_shape_3spp <- filter(genes2species, Gene %in% cluster4_nodes_3spp)

# Order by the vector
cluster4_nodes_shape_ordered_3spp <- cluster4_nodes_shape_3spp[match(cluster4_nodes_3spp, cluster4_nodes_shape_3spp$Gene),]

# Get a vector of the species column
cluster4_nodes_shape_vector_3spp <- cluster4_nodes_shape_ordered_3spp$Species

# Create a vector of weather the genes are in a specific GO or not
# Here the test GO will be GO:0016070 because it has 113 members

# Make a Vector of GO terms
GOs <- GO2genes$GO.ID
GOs_NR <- unique(GOs)

# Add column to 
# Vector to dataframe, GOs as other colums, 1 or 0
cluster4_GO_df_3spp <- data.frame(cluster4_nodes_3spp)

# Vector of genes in GO:0016070
GO_genes_temp <- GO2genes %>% filter(GO.ID == "GO:0016070")
GO_term_genes <- GO_genes_temp$"Genes"

# This creates a data frame
# GO_term_genes <- GO_genes_temp["Genes"]
# And this creates a 'list'
# GO_term_genes <-as.vector(GO_genes_temp["Genes"])

# Get the thaliana IDs for these genes!
# THIS ONLY WORKS IF GO_TERM_GENES IS A CHARACTER CLASS, NOT A DATA FRAME OR LIST (SHOULD HAVE BEEN A VECTOR) WTF?!?!?!?!
genes_in_GO <- filter(data, V1 %in% GO_term_genes)

# Order by the vector
genes_in_GO_ordered <- genes_in_GO[match(cluster4_nodes_3spp, genes_in_GO$V2),]

# Make a TRUE/FALSE vector in the correct order
GO_colour_vector <- genes_in_GO_ordered$V2 %in% cluster4_nodes_3spp


# simple plot of the network with shapes as species
ggnet2(cluster4_net_3spp,
       label = TRUE, label.size = 2,
       shape = cluster4_nodes_shape_vector_3spp,
       shape.palette = c("Cochlearia" = 15, "Amara" = 16, "Thaliana" = 17, "Cochlearia_Thaliana" = 2, "Thaliana_Amara" = 5, "Amara_Cochlearia" = 0, "All three" = 18),
       color = GO_colour_vector,
       color.palette = c("TRUE" = "pink", "FALSE" = "grey"))

# Maybe add a line to size also by species, i.e. the three distinct species can be different shapes, then the rest can be different sizes/larger shapes

# simple plot of the network with shapes colours as species
ggnet2(cluster4_net_3spp,
       label = TRUE, label.size = 2,
       color = cluster4_nodes_shape_vector_3spp,
       color.palette = c("Cochlearia" = "red", "Amara" = "blue", "Thaliana" = "yellow", "Cochlearia_Thaliana" = "orange", "Thaliana_Amara" = "green", "Amara_Cochlearia" = "purple", "All three" = 'grey'))

# 296



# Massively Multiplayer

# Make a plot for each GO/cluster combination
node_no <- 0

for(node in clustersList_3spp){
  node_no <- node_no + 1

  # cluster genes that are in the from column of the interactions
  node_from_3spp <- out_hits_nodups_3spp[out_hits_nodups_3spp$from %in% node,]
  # cluster genes that are in the to column of the interactions
  node_to_3spp <- out_hits_nodups_3spp[out_hits_nodups_3spp$to %in% node,]
  
  # Merge "to" and "from" into one dataframs
  node_merge_3spp <- merge(node_from_3spp, node_to_3spp)
  
  # Remove the last two characters
  node_merge_3spp$from <- substr(node_merge_3spp$from,1,nchar(node_merge_3spp$from)-1)
  node_merge_3spp$to <- substr(node_merge_3spp$to,1,nchar(node_merge_3spp$to)-1)
  
  # Remove the first 5 chararcters
  node_merge_3spp$from <- substr(node_merge_3spp$from,6,nchar(node_merge_3spp$from)-1)
  node_merge_3spp$to <- substr(node_merge_3spp$to,6,nchar(node_merge_3spp$to)-1)
  
  # Make network object
  node_net_3spp <- network(node_merge_3spp)
  
  # Get a vector of network nodes
  node_nodes_3spp <- network.vertex.names(node_net_3spp)
  
  # Returns data frame of values only where Gene column data is in the vector cluster4_nodes_3spp
  node_nodes_shape_3spp <- filter(genes2species, Gene %in% node_nodes_3spp)
  
  # Order by the vector
  node_nodes_shape_ordered_3spp <- node_nodes_shape_3spp[match(node_nodes_3spp, node_nodes_shape_3spp$Gene),]
  
  # Get a vector of the species column
  node_nodes_shape_vector_3spp <- node_nodes_shape_ordered_3spp$Species
  
  for(GO_term in GOs_NR){
    
    # Vector of genes in GO:0016070
    current_GO_genes_temp <- GO2genes %>% filter(GO.ID == GO_term)
    current_GO_term_genes <- current_GO_genes_temp$"Genes"
    
    # Get the thaliana IDs for these genes!
    current_genes_in_GO <- filter(data, V1 %in% current_GO_term_genes)
    
    # Order by the vector
    current_genes_in_GO_ordered <- current_genes_in_GO[match(node_nodes_3spp, current_genes_in_GO$V2),]
    
    # Make a TRUE/FALSE vector in the correct order
    current_GO_colour_vector <- current_genes_in_GO_ordered$V2 %in% node_nodes_3spp
    
    # pdf name
    pdf_name <- c(as.character(node_no), '_', as.character(GO_term), '.pdf')
    pdf_name_string <- paste(pdf_name, collapse='')
    
    # simple plot of the network with shapes as species
    pdf(pdf_name_string)
    ggnet2(node_net_3spp,
           label = TRUE, label.size = 2,
           shape = node_nodes_shape_vector_3spp,
           shape.palette = c("Cochlearia" = 15, "Amara" = 16, "Thaliana" = 17, "Cochlearia_Thaliana" = 2, "Thaliana_Amara" = 5, "Amara_Cochlearia" = 0, "All three" = 18, "NA" = 8),
           color = current_GO_colour_vector,
           color.palette = c("TRUE" = "pink", "FALSE" = "grey"))
    dev.off()
    
  }
}




# test stuff
# make fake net
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# I can use this bit by createing a vector of genes in an OG or set of OGs 
# Must be able to refine to include multiple OGs, probably with a vector the same size as the network
net %v% "phono" = ifelse(letters[1:10] %in% c("a", "e", "i"), "vowel", "consonant")
ggnet2(net, color = "phono", label = TRUE)
# Try colouring with a vector
col_vec = c('red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue')
ggnet2(net, color = col_vec, label = TRUE)
# It works, order has to be the same as in the net
# How to get a vector of node names out that is in the right order?
# see current names with
names_vector <- network.vertex.names(net)






