# Add a function to map FlyBase IDs to gene names
# If you have a mapping file, you can load it here
# For demonstration, we'll create a function that could use packages like biomaRt or org.Dm.eg.db

map_flybase_to_gene_names <- function(flybase_ids) {
  
library(AnnotationDbi)
library(org.Dm.eg.db)
gene_names <- mapIds(org.Dm.eg.db, 
                      keys = flybase_ids,
                      column = "SYMBOL", 
                      keytype = "FLYBASE",
                      multiVals = "first")
  
  # For demonstration purposes, we'll just create dummy gene names
  # In your actual code, replace this with one of the methods above
  gene_names <- vapply(flybase_ids, function(id) {
    # Extract just the numeric part from FBgn identifiers and add "gene_" prefix
    if (grepl("^FBgn", id)) {
      num_part <- sub("^FBgn0*", "", id)
      return(paste0("gene_", num_part))
    } else {
      return(id) # Return original ID if it doesn't match FBgn pattern
    }
  }, character(1))
  
  return(gene_names)
}

# Create a mapping from FlyBase IDs to gene names for all genes in the results
all_gene_ids <- rownames(res_matched)
gene_name_mapping <- map_flybase_to_gene_names(all_gene_ids)
names(gene_name_mapping) <- all_gene_ids

# Add gene names to results tables
res_matched$gene_name <- gene_name_mapping[rownames(res_matched)]
res_time$gene_name <- gene_name_mapping[rownames(res_time)]# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(VennDiagram)

# Assuming metadata and gene_count are already loaded in your environment
# Filter metadata for first three timepoints (time 1, 2, and 3)
metadata_filtered <- metadata %>%
  filter(time %in% c(1, 2, 3))

# Get the sample IDs for the filtered metadata
samples_to_keep <- metadata_filtered$sample_id

# Filter gene count matrix to keep only the samples in the filtered metadata
gene_count_filtered <- gene_count %>%
  select(gene_id, all_of(samples_to_keep))

# Make sure gene_id is a character vector (not a factor) before using as rownames
gene_count_filtered$gene_id <- as.character(gene_count_filtered$gene_id)

# Check for duplicate gene_ids and handle them if present
if(any(duplicated(gene_count_filtered$gene_id))) {
  cat("Warning: Duplicate gene_ids found. Adding suffix to make them unique.\n")
  # Create a frequency table of gene_ids
  gene_freq <- table(gene_count_filtered$gene_id)
  # Find duplicated gene_ids
  dup_genes <- names(gene_freq[gene_freq > 1])
  
  # For each duplicated gene_id, add a suffix
  for(gene in dup_genes) {
    # Find all rows with this gene_id
    idx <- which(gene_count_filtered$gene_id == gene)
    # Add suffix to these gene_ids
    gene_count_filtered$gene_id[idx] <- paste0(gene, "_", seq_along(idx))
  }
}

# Convert to matrix format for DESeq2
gene_count_matrix <- gene_count_filtered %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# Create DESeq2 dataset
# Using time as continuous variable and including matched status
dds <- DESeqDataSetFromMatrix(
  countData = gene_count_matrix,
  colData = metadata_filtered,
  design = ~ matched + time 
)

# Filter low count genes
keep <- rowSums(counts(dds) >= 10) >= min(3, ncol(dds))
dds <- dds[keep,]

# Set reference levels for categorical variables
# Keep time as numeric (continuous variable)
dds$time <- as.numeric(dds$time)
dds$group <- factor(dds$group)
dds$matched <- factor(dds$matched)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
# Time as continuous variable
res_time <- results(dds, name="time")
summary(res_time)

# Matched vs unmatched comparison
res_matched <- results(dds, contrast = c("matched", "TRUE", "FALSE"))
summary(res_matched)

# Extract significant genes
sig_genes_time <- subset(res_time, padj < 0.05)
sig_genes_matched <- subset(res_matched, padj < 0.05)

# Print results with gene names instead of just FlyBase IDs
cat("Top significant genes with time as continuous variable:\n")
sig_time_with_names <- data.frame(
  FlyBase_ID = rownames(sig_genes_time),
  Gene_Name = gene_name_mapping[rownames(sig_genes_time)],
  log2FoldChange = sig_genes_time$log2FoldChange,
  padj = sig_genes_time$padj,
  stringsAsFactors = FALSE
)
head(sig_time_with_names[order(sig_time_with_names$padj),], 10)

cat("\nTop significant genes for matched vs unmatched:\n")
sig_matched_with_names <- data.frame(
  FlyBase_ID = rownames(sig_genes_matched),
  Gene_Name = gene_name_mapping[rownames(sig_genes_matched)],
  log2FoldChange = sig_genes_matched$log2FoldChange,
  padj = sig_genes_matched$padj,
  stringsAsFactors = FALSE
)
head(sig_matched_with_names[order(sig_matched_with_names$padj),], 10)


# Get normalized data for visualization
vsd <- vst(dds, blind = FALSE)

# Get top differentially expressed genes
res_time_ordered <- res_time[order(res_time$padj),]
top_time_genes <- row.names(res_time_ordered)[1:50]

res_matched_ordered <- res_matched[order(res_matched$padj),]
top_matched_genes <- row.names(res_matched_ordered)[1:50]

#
# VISUALIZATION SECTION
#

## Make sure gene names are assigned to the res_matched object
if (!exists("gene_name_mapping")) {
  # Create a simple mapping if it doesn't exist
  gene_name_mapping <- setNames(
    paste0("CG", sub("^FBgn0*", "", rownames(res_matched))),
    rownames(res_matched)
  )
}

# Add gene names column to res_matched
res_matched$gene_name <- gene_name_mapping[rownames(res_matched)]

# Create a data frame with all required columns
volcano_data <- data.frame(
  gene_name = res_matched$gene_name,
  log2FoldChange = res_matched$log2FoldChange,
  padj = res_matched$padj
)

# Select top genes to label
top_genes <- rownames(res_matched)[order(res_matched$padj)[1:20]]
top_gene_names <- gene_name_mapping[top_genes]

# Generate the volcano plot
EnhancedVolcano(res_matched,
                lab = res_matched$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0,
                title = 'Matched vs Unmatched',
                subtitle = paste0('Significant genes: ', nrow(sig_genes_matched)),
                selectLab = top_gene_names,
                legendLabels = c('Not significant', 'Log2FC', 'p-value', 'p-value & Log2FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                maxoverlapsConnectors = 20,
                boxedLabels = TRUE,
                colConnectors = 'black')
# 2. Heatmaps

# Annotation dataframe for heatmaps
heatmap_anno <- metadata_filtered %>% 
  select(time, group, matched) %>%
  mutate(time = as.factor(time)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

# Matched vs unmatched genes heatmap
top_matched_gene_counts <- assay(vsd)[top_matched_genes,]
# Rename rows with gene names instead of FlyBase IDs
rownames(top_matched_gene_counts) <- gene_name_mapping[rownames(top_matched_gene_counts)]

pheatmap(
  top_matched_gene_counts,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = heatmap_anno,
  main = "Top 50 Matched vs Unmatched Differentially Expressed Genes",
  scale = "row"
)

# 3. MA plot for matched vs unmatched
plotMA(res_matched, main = "MA Plot: Matched vs Unmatched")

# 4. Visualize top matched vs unmatched genes with boxplots
# Get normalized counts for top matched genes
top_matched_genes <- rownames(res_matched_ordered)[1:10]
top_matched_gene_names <- gene_name_mapping[top_matched_genes]
top_matched_norm_counts <- counts(dds, normalized=TRUE)[top_matched_genes,]

# Create a dataframe with gene counts
top_matched_df <- as.data.frame(t(top_matched_norm_counts))
colnames(top_matched_df) <- top_matched_gene_names  # Use gene names for column names
top_matched_df$sample_id <- rownames(top_matched_df)

# Merge with metadata
top_matched_df <- merge(top_matched_df, 
                        metadata_filtered[, c("sample_id", "time", "group", "matched")], 
                        by="sample_id")

# Reshape for plotting
top_matched_long <- pivot_longer(top_matched_df, 
                                 cols = -c(sample_id, time, group, matched), 
                                 names_to = "gene", 
                                 values_to = "normalized_count")

# Plot matched vs unmatched for top 10 genes
ggplot(top_matched_long, 
       aes(x = matched, y = normalized_count, fill = matched)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  labs(title = "Top 10 Differentially Expressed Genes: Matched vs Unmatched",
       x = "Matched Status", 
       y = "Normalized Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








# Analyze upregulated and downregulated genes in matched vs unmatched comparison
# and perform functional enrichment analysis

# First, separate significant genes into upregulated and downregulated
# Remember: positive log2FoldChange means upregulated in matched compared to unmatched
upregulated_in_matched <- subset(res_matched, padj < 0.05 & log2FoldChange > 0)
downregulated_in_matched <- subset(res_matched, padj < 0.05 & log2FoldChange < 0)

# Which means upregulated in unmatched compared to matched
upregulated_in_unmatched <- downregulated_in_matched

# Count the number of genes in each category
cat("Number of genes upregulated in matched samples:", nrow(upregulated_in_matched), "\n")
cat("Number of genes upregulated in unmatched samples:", nrow(upregulated_in_unmatched), "\n")

# Create data frames with gene names for better readability
upregulated_in_matched_df <- data.frame(
  FlyBase_ID = rownames(upregulated_in_matched),
  Gene_Name = gene_name_mapping[rownames(upregulated_in_matched)],
  log2FoldChange = upregulated_in_matched$log2FoldChange,
  padj = upregulated_in_matched$padj,
  stringsAsFactors = FALSE
)

upregulated_in_unmatched_df <- data.frame(
  FlyBase_ID = rownames(upregulated_in_unmatched),
  Gene_Name = gene_name_mapping[rownames(upregulated_in_unmatched)],
  log2FoldChange = upregulated_in_unmatched$log2FoldChange,
  padj = upregulated_in_unmatched$padj,
  stringsAsFactors = FALSE
)

# Sort by fold change magnitude
upregulated_in_matched_df <- upregulated_in_matched_df[order(-upregulated_in_matched_df$log2FoldChange),]
upregulated_in_unmatched_df <- upregulated_in_unmatched_df[order(upregulated_in_unmatched_df$log2FoldChange),]


## ID Conversion, GO and KEGG Analysis, and Gene Interaction Networks for Drosophila

# Load required libraries
library(clusterProfiler)
library(org.Dm.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(igraph)
library(STRINGdb)
library(tidyverse)

#######################################
## 1. FlyBase ID Conversion Functions
#######################################

# Convert FlyBase IDs to different ID types
convert_flybase_ids <- function(flybase_ids, to_type) {
  # Valid to_type options: "ENTREZID", "SYMBOL", "GENENAME", "ENSEMBL", "REFSEQ", "UNIPROT"
  
  # Remove any NA or empty values
  flybase_ids <- flybase_ids[!is.na(flybase_ids) & flybase_ids != ""]
  
  # Perform the conversion
  converted_ids <- mapIds(org.Dm.eg.db,
                          keys = flybase_ids,
                          column = to_type,
                          keytype = "FLYBASE",
                          multiVals = "first")
  
  # Create a result dataframe
  result_df <- data.frame(
    FlyBase_ID = flybase_ids,
    stringsAsFactors = FALSE
  )
  
  # Add the converted IDs
  result_df[[to_type]] <- converted_ids
  
  return(result_df)
}

# Create a comprehensive mapping table with multiple ID types
create_id_mapping_table <- function(flybase_ids) {
  # Convert to Entrez IDs
  entrez_mapping <- convert_flybase_ids(flybase_ids, "ENTREZID")
  
  # Convert to gene symbols
  symbol_mapping <- convert_flybase_ids(flybase_ids, "SYMBOL")
  
  # Convert to gene names
  genename_mapping <- convert_flybase_ids(flybase_ids, "GENENAME")
  
  # Merge the results
  mapping_table <- entrez_mapping %>%
    left_join(symbol_mapping, by = "FlyBase_ID") %>%
    left_join(genename_mapping, by = "FlyBase_ID")
  
  return(mapping_table)
}

#######################################
## 2. GO Analysis Functions
#######################################

# Perform GO enrichment analysis with visualization
perform_go_analysis <- function(gene_list, gene_universe = NULL, ont = "BP", title = "") {
  # Convert FlyBase IDs to Entrez IDs
  entrez_ids <- mapIds(org.Dm.eg.db,
                       keys = gene_list,
                       column = "ENTREZID",
                       keytype = "FLYBASE",
                       multiVals = "first")
  
  # Remove NAs
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # If universe is provided, convert to Entrez IDs
  if (!is.null(gene_universe)) {
    universe_entrez <- mapIds(org.Dm.eg.db,
                              keys = gene_universe,
                              column = "ENTREZID",
                              keytype = "FLYBASE",
                              multiVals = "first")
    
    universe_entrez <- universe_entrez[!is.na(universe_entrez)]
  } else {
    universe_entrez <- NULL
  }
  
  #Perform GO enrichment analysis with visualization
  # Perform GO enrichment analysis with visualization
  perform_go_analysis <- function(gene_list, gene_universe = NULL, ont = "BP", title = "") {
    # Convert FlyBase IDs to Entrez IDs
    entrez_ids <- mapIds(org.Dm.eg.db,
                         keys = gene_list,
                         column = "ENTREZID",
                         keytype = "FLYBASE",
                         multiVals = "first")
    
    # Remove NAs
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    # If universe is provided, convert to Entrez IDs
    if (!is.null(gene_universe)) {
      universe_entrez <- mapIds(org.Dm.eg.db,
                                keys = gene_universe,
                                column = "ENTREZID",
                                keytype = "FLYBASE",
                                multiVals = "first")
      
      universe_entrez <- universe_entrez[!is.na(universe_entrez)]
    } else {
      universe_entrez <- NULL
    }
    
    # Perform GO enrichment analysis
    go_result <- enrichGO(gene = entrez_ids,
                          universe = universe_entrez,
                          OrgDb = org.Dm.eg.db,
                          ont = ont,
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)
    
    # Return early if no enrichment found
    if (is.null(go_result) || nrow(go_result@result) == 0) {
      message(paste("No significant GO", ont, "terms found for", title))
      return(NULL)
    }
    
    # Create a dotplot
    dot_plot <- dotplot(go_result, 
                        showCategory = min(20, nrow(go_result@result)),
                        title = paste("GO", ont, "Enrichment:", title))
    
    # Display the dotplot
    print(dot_plot)
    
    # Create an enrichment map - with error handling
    if (nrow(go_result@result) >= 5) {
      # Try to create an enrichment map, but handle errors
      tryCatch({
        # Set showCategory to a smaller number to avoid issues with term similarity calculation
        enrichment_map <- emapplot(go_result, 
                                   showCategory = min(15, nrow(go_result@result)),
                                   color = "p.adjust")
        
        # Display the enrichment map with a separate title
        print(enrichment_map + ggtitle(paste("GO", ont, "Network:", title)))
      }, error = function(e) {
        message(paste("Note: Could not create enrichment map plot due to error:", e$message))
        message("This is often due to issues with term similarity calculations or term counts.")
        message("Trying alternative approach with fewer categories...")
        
        # Try with even fewer categories as an alternative
        tryCatch({
          if (nrow(go_result@result) >= 3) {
            alt_map <- emapplot(pairwise_termsim(go_result), 
                                showCategory = min(8, nrow(go_result@result)))
            print(alt_map + ggtitle(paste("GO", ont, "Network (Alt):", title)))
          }
        }, error = function(e2) {
          message("Alternative approach also failed. Skipping enrichment map.")
        })
      })
      
      # Create a category network plot - with error handling
      tryCatch({
        category_net <- cnetplot(go_result,
                                 categorySize = "pvalue",
                                 foldChange = NULL,
                                 showCategory = min(8, nrow(go_result@result)))
        
        # Display the category network
        print(category_net)
      }, error = function(e) {
        message(paste("Note: Could not create category network plot due to error:", e$message))
        message("This is often due to issues with term-gene connections or gene mappings.")
        
        # Try with even fewer categories as an alternative
        tryCatch({
          if (nrow(go_result@result) >= 3) {
            message("Trying with fewer categories...")
            alt_net <- cnetplot(go_result,
                                categorySize = "pvalue",
                                showCategory = 3)
            print(alt_net)
          }
        }, error = function(e2) {
          message("Alternative approach also failed. Skipping category network plot.")
        })
      })
    }
    
    # Return the enrichment result
    return(go_result)
  }
}
  #######################################
  ## 3. KEGG Analysis Functions
  #######################################
  
  # Helper function to process and visualize KEGG results
  process_kegg_results <- function(kegg_result, title) {
    # Make it readable
    kegg_result <- setReadable(kegg_result, OrgDb = org.Dm.eg.db, keyType = "ENTREZID")
    
    # Create a dotplot
    dot_plot <- dotplot(kegg_result, 
                        showCategory = min(20, nrow(kegg_result@result)),
                        title = paste("KEGG Pathway Enrichment:", title))
    
    # Display the dotplot
    print(dot_plot)
    
    # Create a pathway-gene network plot
    if (nrow(kegg_result@result) >= 3) {
      tryCatch({
        pathway_net <- cnetplot(kegg_result,
                                categorySize = "pvalue",
                                foldChange = NULL)
        
        # Display the pathway network
        print(pathway_net)
      }, error = function(e) {
        message(paste("Could not create pathway network plot:", e$message))
      })
    }
    
    # Print the top pathways
    cat("\nTop KEGG pathways:\n")
    print(kegg_result@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "Count")])
    
    return(kegg_result)
  }
  
  # Perform KEGG pathway analysis with visualization
  perform_kegg_analysis <- function(gene_list, gene_universe = NULL, title = "") {
    # Convert FlyBase IDs to Entrez IDs first
    entrez_ids <- mapIds(org.Dm.eg.db,
                         keys = gene_list,
                         column = "ENTREZID",
                         keytype = "FLYBASE",
                         multiVals = "first")
    
    # Remove NAs
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    # If universe is provided, convert to Entrez IDs
    if (!is.null(gene_universe)) {
      universe_entrez <- mapIds(org.Dm.eg.db,
                                keys = gene_universe,
                                column = "ENTREZID",
                                keytype = "FLYBASE",
                                multiVals = "first")
      
      universe_entrez <- universe_entrez[!is.na(universe_entrez)]
    } else {
      universe_entrez <- NULL
    }
    
    # If we don't have any genes, return early
    if (length(entrez_ids) == 0) {
      message("No genes could be mapped to Entrez IDs for KEGG analysis")
      return(NULL)
    }
    
    # Print how many genes we're using
    message(paste("Performing KEGG analysis with", length(entrez_ids), "mapped genes"))
    
    # Try multiple approaches for KEGG analysis
    
    # First attempt: Standard KEGG with Entrez IDs
    tryCatch({
      message("Attempting KEGG analysis with Entrez IDs...")
      kegg_resut}

#######################################
## 4. Gene Interaction Network Functions
#######################################

# Generate a STRING gene interaction network
generate_string_network <- function(gene_list, title = "", species = 7227) {
  # Create a new STRINGdb object
  string_db <- STRINGdb$new(version = "11.5", species = species, score_threshold = 400)
  
  # Convert FlyBase IDs to gene symbols
  symbols <- mapIds(org.Dm.eg.db,
                    keys = gene_list,
                    column = "SYMBOL",
                    keytype = "FLYBASE",
                    multiVals = "first")
  
  # Remove NAs
  symbols <- symbols[!is.na(symbols)]
  
  # Map gene symbols to STRING IDs
  mapped_genes <- string_db$map(data.frame(gene = symbols), "gene", removeUnmappedRows = TRUE)
  
  # If no genes mapped, return NULL
  if (nrow(mapped_genes) == 0) {
    message("No genes could be mapped to STRING IDs")
    return(NULL)
  }
  
  # Get interactions between the mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # If no interactions found, return NULL
  if (nrow(interactions) == 0) {
    message("No interactions found between the provided genes")
    return(NULL)
  }
  
  # Create an igraph object from the interactions
  network <- graph_from_data_frame(interactions[, c("from", "to", "combined_score")], directed = FALSE)
  
  # Add gene names as vertex labels
  V(network)$name <- string_db$get_aliases(V(network)$name)$alias
  
  # Plot the network
  plot(network, 
       layout = layout_with_fr(network),
       vertex.size = 8,
       vertex.label.cex = 0.8,
       vertex.label.dist = 1.5,
       vertex.label.color = "black",
       vertex.color = "lightblue",
       edge.width = E(network)$combined_score / 500,
       main = paste("STRING Gene Interaction Network:", title))
  
  return(network)
}

# Alternative method using ggraph for prettier network visualization
generate_ggraph_network <- function(gene_list, title = "", species = 7227) {
  # Requires additional packages
  if (!requireNamespace("ggraph", quietly = TRUE))
    install.packages("ggraph")
  
  library(ggraph)
  
  # Create a new STRINGdb object
  string_db <- STRINGdb$new(version = "11.5", species = species, score_threshold = 400)
  
  # Convert FlyBase IDs to gene symbols
  symbols <- mapIds(org.Dm.eg.db,
                    keys = gene_list,
                    column = "SYMBOL",
                    keytype = "FLYBASE",
                    multiVals = "first")
  
  # Remove NAs
  symbols <- symbols[!is.na(symbols)]
  
  # Map gene symbols to STRING IDs
  mapped_genes <- string_db$map(data.frame(gene = symbols), "gene", removeUnmappedRows = TRUE)
  
  # If no genes mapped, return NULL
  if (nrow(mapped_genes) == 0) {
    message("No genes could be mapped to STRING IDs")
    return(NULL)
  }
  
  # Get interactions between the mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # If no interactions found, return NULL
  if (nrow(interactions) == 0) {
    message("No interactions found between the provided genes")
    return(NULL)
  }
  
  # Create a node data frame with gene names
  nodes <- data.frame(
    id = unique(c(interactions$from, interactions$to)),
    stringsAsFactors = FALSE
  )
  
  nodes$label <- string_db$get_aliases(nodes$id)$alias
  
  # Create an igraph object
  network <- graph_from_data_frame(interactions, directed = FALSE, vertices = nodes)
  
  # Set edge weights based on combined score
  E(network)$weight <- E(network)$combined_score / 1000
  
  # Plot with ggraph
  ggraph(network, layout = "fr") +
    geom_edge_link(aes(width = weight, alpha = weight)) +
    geom_node_point(size = 5, color = "lightblue") +
    geom_node_text(aes(label = label), repel = TRUE) +
    theme_graph() +
    labs(title = paste("Gene Interaction Network:", title)) +
    theme(legend.position = "none")
}

#######################################
## 5. Run the Analysis for your Data
#######################################

# Let's assume you've already run the DESeq2 analysis and have:
# upregulated_in_matched_df and upregulated_in_unmatched_df

# 1. Create ID mapping tables
cat("\n\n--- Creating ID Mapping Tables ---\n")
matched_id_mapping <- create_id_mapping_table(upregulated_in_matched_df$FlyBase_ID)
unmatched_id_mapping <- create_id_mapping_table(upregulated_in_unmatched_df$FlyBase_ID)

# 2. Perform GO analysis for matched upregulated genes
cat("\n\n--- GO Analysis for Genes Upregulated in Matched Samples ---\n")
matched_go_bp <- perform_go_analysis(
  upregulated_in_matched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  ont = "BP",
  title = "Upregulated in Matched"
)

matched_go_mf <- perform_go_analysis(
  upregulated_in_matched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  ont = "MF",
  title = "Upregulated in Matched"
)

matched_go_cc <- perform_go_analysis(
  upregulated_in_matched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  ont = "CC",
  title = "Upregulated in Matched"
)

# 3. Perform GO analysis for unmatched upregulated genes
cat("\n\n--- GO Analysis for Genes Upregulated in Unmatched Samples ---\n")
unmatched_go_bp <- perform_go_analysis(
  upregulated_in_unmatched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  ont = "BP",
  title = "Upregulated in Unmatched"
)

unmatched_go_mf <- perform_go_analysis(
  upregulated_in_unmatched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  ont = "MF",
  title = "Upregulated in Unmatched"
)

unmatched_go_cc <- perform_go_analysis(
  upregulated_in_unmatched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  ont = "CC",
  title = "Upregulated in Unmatched"
)

# 4. Perform KEGG pathway analysis for matched upregulated genes
cat("\n\n--- KEGG Pathway Analysis for Genes Upregulated in Matched Samples ---\n")
matched_kegg <- perform_kegg_analysis(
  upregulated_in_matched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  title = "Upregulated in Matched"
)

# 5. Perform KEGG pathway analysis for unmatched upregulated genes
cat("\n\n--- KEGG Pathway Analysis for Genes Upregulated in Unmatched Samples ---\n")
unmatched_kegg <- perform_kegg_analysis(
  upregulated_in_unmatched_df$FlyBase_ID,
  gene_universe = rownames(res_matched),
  title = "Upregulated in Unmatched"
)

# 6. Create gene interaction networks for top genes
cat("\n\n--- Gene Interaction Networks ---\n")

# For matched: Take top 30 genes by p-value
top_matched_genes <- upregulated_in_matched_df$FlyBase_ID[
  order(upregulated_in_matched_df$padj)
][1:min(30, nrow(upregulated_in_matched_df))]

cat("\nGenerating interaction network for top genes upregulated in matched samples...\n")
matched_network <- generate_string_network(
  top_matched_genes,
  title = "Top Genes Upregulated in Matched"
)

# For unmatched: Take top 30 genes by p-value
top_unmatched_genes <- upregulated_in_unmatched_df$FlyBase_ID[
  order(upregulated_in_unmatched_df$padj)
][1:min(30, nrow(upregulated_in_unmatched_df))]

cat("\nGenerating interaction network for top genes upregulated in unmatched samples...\n")
unmatched_network <- generate_string_network(
  top_unmatched_genes,
  title = "Top Genes Upregulated in Unmatched"
)

# 7. Create ggraph-based pretty network visualizations
cat("\n\n--- Enhanced Gene Interaction Network Visualizations ---\n")

cat("\nGenerating ggraph network for top genes upregulated in matched samples...\n")
generate_ggraph_network(
  top_matched_genes,
  title = "Top Genes Upregulated in Matched"
)

cat("\nGenerating ggraph network for top genes upregulated in unmatched samples...\n")
generate_ggraph_network(
  top_unmatched_genes,
  title = "Top Genes Upregulated in Unmatched"
)

# 8. Create combined enhanced dotplots for GO terms with gene ratio
cat("\n\n--- Enhanced GO Term Dotplots ---\n")

# Function to create enhanced dotplots
create_enhanced_dotplot <- function(go_result, title) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    message(paste("No significant GO terms found for", title))
    return(NULL)
  }
  
  # Get top 15 terms
  top_terms <- go_result@result %>%
    arrange(p.adjust) %>%
    head(15)
  
  # Create the plot
  ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, -p.adjust))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = paste("GO Term Enrichment:", title),
      x = "Gene Ratio",
      y = "GO Term",
      size = "Gene Count",
      color = "Adjusted p-value"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}

# Create enhanced dotplots for matched samples
if (!is.null(matched_go_bp)) {
  print(create_enhanced_dotplot(matched_go_bp, "BP - Upregulated in Matched"))
}

if (!is.null(matched_go_mf)) {
  print(create_enhanced_dotplot(matched_go_mf, "MF - Upregulated in Matched"))
}

# Create enhanced dotplots for unmatched samples
if (!is.null(unmatched_go_bp)) {
  print(create_enhanced_dotplot(unmatched_go_bp, "BP - Upregulated in Unmatched"))
}

if (!is.null(unmatched_go_mf)) {
  print(create_enhanced_dotplot(unmatched_go_mf, "MF - Upregulated in Unmatched"))
}

# 9. Compare GO terms between matched and unmatched
cat("\n\n--- GO Term Comparison Between Conditions ---\n")

# Function to compare GO terms between two conditions
compare_go_terms <- function(go_result1, go_result2, title1, title2) {
  if (is.null(go_result1) || is.null(go_result2)) {
    message("Cannot compare GO terms: At least one result set is empty")
    return(NULL)
  }
  
  # Extract results
  res1 <- go_result1@result %>% 
    arrange(p.adjust) %>% 
    head(30)
  
  res2 <- go_result2@result %>% 
    arrange(p.adjust) %>% 
    head(30)
  
  # Combine the data
  res1$condition <- title1
  res2$condition <- title2
  
  combined <- rbind(res1, res2)
  
  # Find shared GO terms
  shared_terms <- intersect(res1$ID, res2$ID)
  
  cat("Number of shared GO terms between conditions:", length(shared_terms), "\n")
  
  if (length(shared_terms) > 0) {
    cat("Shared GO terms:\n")
    shared_data <- combined %>% 
      filter(ID %in% shared_terms) %>%
      select(ID, Description, condition, p.adjust) %>%
      arrange(ID, condition)
    
    print(shared_data)
    
    # Create a plot for shared terms
    shared_wide <- shared_data %>%
      select(Description, condition, p.adjust) %>%
      pivot_wider(names_from = condition, values_from = p.adjust
#