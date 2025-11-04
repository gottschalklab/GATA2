setwd("/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/AM_Results/TF_target_analysis")

# libraries 
library("qgraph")

# input
net_info_f = "GATA2_DEG_Targets_TFinteractors_annotation.csv"

# functions:
# create an edge list from from gene to gene weight matrix
make_edgelist <- function(weights) {
  # create edgelist
  edgelist <- which(weights != 0, arr.ind = TRUE)
  edgelist[,2] <- dim(weights)[1] + edgelist[,2]
  return(edgelist)
}

# Load network information
net_info <- read.csv(net_info_f, header = TRUE, row.names = 1, check.names = FALSE)

# Drop certain TFs (columns)
net_info <- net_info[, !(names(net_info) %in% c("STAT3", "CEBPA", "CEBPB", "NR1H3", "FLI1", "KLF4", "IRF8", "FOXO1", "MECOM"))]

# Filter genes (rows) where column 2 == 1
# Assume column 2 is "fibrosis", or grab it generically:
col2_name <- names(net_info)[3]
net_info <- net_info[net_info[[col2_name]] == 1, ]

# Now define early and late matrices
early <- t(net_info[, 2])  # column 2 is now filtered, still named "fibrosis"?
late <- t(net_info[, 4:ncol(net_info)])  # rest of TFs

# Filter late
late <- late[rowSums(late) != 0, colSums(late) != 0]

# Find shared sites and genes
common_cols <- colnames(late)

# Subset full edge matrix using common indices
edge_df <- t(net_info[, 4:ncol(net_info)])
edge_df <- edge_df[, common_cols]



# network setup
# create edgelist
edgelist <- make_edgelist(edge_df)

# network formatting 
color_vec <- ifelse(
  net_info[common_cols, "metabolism"] == 1, "hotpink", "gray"
)

node_colors <- c(rep(c('white'),each=dim(edge_df)[1]),color_vec)

node_shapes <- c(rep(c('diamond'),each=dim(edge_df)[1]),rep(c('circle'),each=dim(edge_df)[2]))

node_sizes <- c(rep(c(2.5),each=dim(edge_df)[1]),rep(c(1),each=dim(edge_df)[2]))

node_labels <- c(rownames(edge_df),rep(c(""),each=dim(edge_df)[2]))

edge_color <- '#dbd9d9'


# plot network
pdf("20250325_top6_TF_targets_and_GATA2_DEGs_MetabolismPathway_network.pdf")
qgraph(edgelist, layout = "spring", directed = FALSE, labels = node_labels, color = node_colors, 
       shape = node_shapes, node.width = node_sizes, node.height = node_sizes, edge.color = edge_color, 
       label.scale = FALSE, label.cex = 0.5)
dev.off()
