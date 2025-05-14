
library(clusterProfiler)
library(org.Dm.eg.db)  
library(enrichplot)
library(ggplot2)

# Get significant genes
sig_gene_ids <- rownames(sig_genes_matched)

# Convert FlyBase IDs to ENTREZ IDs
entrez_ids <- bitr(sig_gene_ids, 
                   fromType = "FLYBASE", 
                   toType = c("ENTREZID", "SYMBOL"),  # Get both ENTREZ ID and gene symbol
                   OrgDb = org.Dm.eg.db)

# Keep only genes that were successfully mapped
entrez_list <- entrez_ids$ENTREZID

# Create a named vector for gene labeling
gene_labels <- entrez_ids$SYMBOL
names(gene_labels) <- entrez_ids$ENTREZID

# Replace NA or empty symbols with the original ENTREZID
empty_symbols <- is.na(gene_labels) | gene_labels == ""
gene_labels[empty_symbols] <- names(gene_labels)[empty_symbols]

# Run GO analysis for Biological Process
go_bp <- enrichGO(gene = entrez_list,
                  OrgDb = org.Dm.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

# Run KEGG analysis
kegg_result <- enrichKEGG(gene = entrez_list,
                          organism = 'dme',
                          keyType = 'ncbi-geneid',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")

# Function to create a network plot with explicit node labels
plot_network_with_labels <- function(enrichment_result, title, gene_labels) {
  # Skip if no significant results
  if (nrow(enrichment_result@result) == 0) {
    cat("No significant results for", title, "\n")
    return(NULL)
  }
  
  # Create the network plot with both category and gene labels
  p <- cnetplot(enrichment_result,
                showCategory = min(5, nrow(enrichment_result@result)),
                categorySize = "pvalue",
                nodeLabel = "all",  # Show both category and gene labels
                cex_label_category = 1.0,
                cex_label_gene = 0.8,
                node_label_size = 5,
                label_format = function(id) {
                  # For genes, use the gene symbol
                  if (id %in% names(gene_labels)) {
                    return(gene_labels[id])
                  }
                  # For categories, keep original ID
                  return(id)
                })
  
  # Add a title
  p <- p + ggtitle(title)
  
  return(p)
}

# Alternative approach: create a network plot with editable nodes
plot_network_alt <- function(enrichment_result, title, gene_labels) {
  # Skip if no significant results
  if (nrow(enrichment_result@result) == 0) {
    cat("No significant results for", title, "\n")
    return(NULL)
  }
  
  # Get the data for the plot
  n <- min(5, nrow(enrichment_result@result))
  p <- cnetplot(enrichment_result, 
                showCategory = n, 
                foldChange = NULL, 
                node_label = "none")  # Initially no labels
  
  # Extract the plot data
  plot_data <- ggplot_build(p)
  
  # Get the node coordinates
  nodes <- plot_data$data[[1]]
  
  # Separate category and gene nodes
  category_nodes <- nodes[1:n,]
  gene_nodes <- nodes[(n+1):nrow(nodes),]
  
  # Map gene IDs to labels
  gene_ids <- rownames(gene_nodes)
  gene_text <- sapply(gene_ids, function(id) {
    if (id %in% names(gene_labels)) {
      return(gene_labels[id])
    }
    return(id)
  })
  
  # Add text labels manually
  p <- p + 
    geom_text(data = category_nodes, 
              aes(x = x, y = y, label = group), 
              size = 4, color = "black") +
    geom_text(data = gene_nodes,
              aes(x = x, y = y, label = gene_text),
              size = 3, color = "black") +
    ggtitle(title)
  
  return(p)
}

# Third approach: direct edgeR plot with explicit node labels
direct_network_plot <- function(enrichment_result, title, gene_labels) {
  if (nrow(enrichment_result@result) == 0) {
    return(NULL)
  }
  
  # Get top categories
  n <- min(5, nrow(enrichment_result@result))
  categories <- enrichment_result@result$ID[1:n]
  category_names <- enrichment_result@result$Description[1:n]
  
  # Set up node data
  nodes <- data.frame()
  edges <- data.frame()
  
  # Add category nodes
  for (i in 1:length(categories)) {
    cat_id <- categories[i]
    cat_name <- category_names[i]
    
    # Add category node
    cat_node <- data.frame(id = cat_id, label = cat_name, type = "category")
    nodes <- rbind(nodes, cat_node)
    
    # Get genes for this category
    genes <- enrichment_result@geneSets[[cat_id]]
    
    # Add gene nodes and edges
    for (gene in genes) {
      # Add gene node if not already added
      if (!gene %in% nodes$id) {
        gene_label <- ifelse(gene %in% names(gene_labels), gene_labels[gene], gene)
        gene_node <- data.frame(id = gene, label = gene_label, type = "gene")
        nodes <- rbind(nodes, gene_node)
      }
      
      # Add edge
      edge <- data.frame(from = cat_id, to = gene)
      edges <- rbind(edges, edge)
    }
  }
  
  # Create the network using igraph
  library(igraph)
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # Set node colors
  V(g)$color <- ifelse(V(g)$type == "category", "tomato3", "steelblue")
  
  # Set node sizes
  V(g)$size <- ifelse(V(g)$type == "category", 15, 10)
  
  # Plot the network
  plot(g, 
       layout = layout_with_fr(g),
       vertex.label = V(g)$label,
       vertex.label.cex = ifelse(V(g)$type == "category", 1.0, 0.8),
       vertex.label.color = "black",
       edge.width = 1,
       main = title)
  
  return(g)
}

# Generate the plots

# Biological Process network
cat("Creating Biological Process network plot...\n")
bp_network <- plot_network_with_labels(go_bp, "GO Biological Process Network", gene_labels)
if (!is.null(bp_network)) {
  print(bp_network)
}

# Alternate approach for BP
cat("Creating alternate BP network plot...\n")
bp_network_alt <- plot_network_alt(go_bp, "GO Biological Process Network (Alt)", gene_labels)
if (!is.null(bp_network_alt)) {
  print(bp_network_alt)
}

# Direct igraph approach for BP
cat("Creating direct BP network plot...\n")
bp_network_direct <- direct_network_plot(go_bp, "GO Biological Process Network (Direct)", gene_labels)

# KEGG pathway network
cat("Creating KEGG pathway network plot...\n")
kegg_network <- plot_network_with_labels(kegg_result, "KEGG Pathway Network", gene_labels)
if (!is.null(kegg_network)) {
  print(kegg_network)
}

# Alternate approach for KEGG
cat("Creating alternate KEGG network plot...\n")
kegg_network_alt <- plot_network_alt(kegg_result, "KEGG Pathway Network (Alt)", gene_labels)
if (!is.null(kegg_network_alt)) {
  print(kegg_network_alt)
}

# Direct igraph approach for KEGG
cat("Creating direct KEGG network plot...\n")
kegg_network_direct <- direct_network_plot(kegg_result, "KEGG Pathway Network (Direct)", gene_labels)

