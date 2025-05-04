library(ImpulseDE2)
library(writexl)

# Subset metadata for T1-T3 (as in your original analysis)
metadata_early <- metadata[metadata$time <= 3, ]

# Ensure counts are a matrix and match metadata order
gene_count_matrix <- as.matrix(gene_count[, metadata_early$sample_id])  # Align columns

# Prefilter: Keep genes with ≥10 counts in ≥3 samples
keep_genes <- rowSums(gene_count_matrix >= 10) >= 3
gene_count_filtered <- gene_count_matrix[keep_genes, ]

# Case-control annotation (case = mismatched, control = matched)
dfAnnotation <- data.frame(
  Sample = metadata_early$sample_id,
  Time = as.numeric(metadata_early$time),  # Must be numeric for ImpulseDE2
  Condition = ifelse(metadata_early$matched, "control", "case"),  # TRUE=matched=control
  stringsAsFactors = FALSE
)

# Verify sample alignment
stopifnot(all(colnames(gene_count_filtered) == dfAnnotation$Sample))

# Run ImpulseDE2 (case-control mode)
results <- runImpulseDE2(
  matCountData = gene_count_filtered,
  dfAnnotation = dfAnnotation,
  boolCaseCtrl = TRUE,  # Compare case (FA) vs. control (FF/AA)
  scaNProc = 4,         # Parallel processing
  scaQThres = 0.05      # FDR threshold
)

# Extract significant genes (padj < 0.05, |log2FC| > 0.5 if needed)
sig_genes <- results$dfImpulseDE2Results[results$dfImpulseDE2Results$padj < 0.05, ]

# Save results
write_xlsx(
  list(DEGs = as.data.frame(sig_genes)),
  path = "ImpulseDE2_mitonuclear_mismatch_T1-T3.xlsx"
)


# Load required libraries
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(RColorBrewer)
library(org.Dm.eg.db)  # Drosophila-specific annotation package
library(AnnotationDbi)
library(readxl)      # For reading your Excel results
library(writexl)     # For writing results



# If you already have it in your R environment, you can use it directly
results_df <- ImpulseDE2_mitonuclear_mismatch_T1_T3

# Add regulation information
results_df <- results_df %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Upregulated in Mismatched (FA)",
    padj < 0.05 & log2FoldChange < 0 ~ "Upregulated in Matched (FF+AA)",
    TRUE ~ "Not Significant"
  ))

# Get the gene symbols for up and down regulated genes
up_mismatched <- results_df %>% 
  filter(regulation == "Upregulated in Mismatched (FA)") %>%
  pull(symbol)

up_matched <- results_df %>% 
  filter(regulation == "Upregulated in Matched (FF+AA)") %>%
  pull(symbol)

# Print summary counts
cat("Total DEGs:", nrow(results_df %>% filter(padj < 0.05)), "\n")
cat("Upregulated in Mismatched (FA):", length(up_mismatched), "\n")
cat("Upregulated in Matched (FF+AA):", length(up_matched), "\n")

# First, check if symbol column exists and contains valid symbols
if(!"symbol" %in% colnames(results_df)) {
  stop("Symbol column not found in the results dataframe. Check your data.")
}

# Convert symbols to Entrez IDs for GO analysis
entrez_up_mismatched <- mapIds(
  org.Dm.eg.db, 
  keys = up_mismatched,
  column = "ENTREZID", 
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_up_matched <- mapIds(
  org.Dm.eg.db, 
  keys = up_matched,
  column = "ENTREZID", 
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove any NA values
entrez_up_mismatched <- entrez_up_mismatched[!is.na(entrez_up_mismatched)]
entrez_up_matched <- entrez_up_matched[!is.na(entrez_up_matched)]

# Print gene ID conversion statistics
cat("\nGene ID conversion statistics:\n")
cat("Upregulated in Mismatched - Original:", length(up_mismatched), 
    "Mapped to Entrez:", length(entrez_up_mismatched), 
    "Success rate:", round(length(entrez_up_mismatched)/length(up_mismatched)*100, 1), "%\n")

cat("Upregulated in Matched - Original:", length(up_matched), 
    "Mapped to Entrez:", length(entrez_up_matched), 
    "Success rate:", round(length(entrez_up_matched)/length(up_matched)*100, 1), "%\n")

# Create a list for compareCluster
gene_list <- list(
  "Upregulated in Mismatched (FA)" = entrez_up_mismatched,
  "Upregulated in Matched (FF+AA)" = entrez_up_matched
)

# Run compareCluster to perform GO enrichment analysis
compare_result <- compareCluster(
  gene_list,
  fun = "enrichGO",
  OrgDb = org.Dm.eg.db,
  ont = "BP",        # Biological Process
  pAdjustMethod = "BH",  # Benjamini-Hochberg adjustment
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE    # Convert gene IDs to symbols for readability
)

# Define key terms for filtering the biological processes of interest
key_processes <- c(
  # Immune response
  "immune", "inflammatory", "defense", "leukocyte", "cytokine",
  
  # Body morphogenesis
  "morphogenesis", "development", "pattern formation", "body",
  
  # Autophagy
  "autophagy", "lysosome", "vacuole", "phagosome",
  
  # DNA damage
  "DNA damage", "DNA repair", "genome stability",
  
  # Oxidative stress
  "oxidative stress", "reactive oxygen", "antioxidant", "ROS", "redox",
  
  # Amino acid biosynthesis
  "amino acid", "protein synthesis", "translation", "protein folding",
  
  # Mitochondrial expression
  "mitochondri", "respiration",
  
  # Carbohydrate metabolism
  "carbohydrate", "glucose", "glycolysis", "gluconeogenesis", "glycogen", "sugar",
  
  # OXPHOS
  "oxidative phosphorylation", "electron transport", "respiratory chain", "ETC",
  
  # TCA cycle
  "TCA", "tricarboxylic", "citrate", "Krebs", "citric acid cycle"
)

# Filter the enrichment results to focus on these processes
enrichment_results <- compare_result@compareClusterResult

# Create a function to check if a term contains any of our keywords
contains_key_process <- function(term) {
  any(sapply(key_processes, function(process) {
    grepl(process, term, ignore.case = TRUE)
  }))
}

# Filter the results
filtered_results <- enrichment_results %>%
  filter(sapply(Description, contains_key_process))

# Add a new column to categorize each term
filtered_results$Pathway <- sapply(filtered_results$Description, function(term) {
  if (grepl("immune|inflammatory|defense|leukocyte|cytokine", term, ignore.case = TRUE))
    return("Immune Response")
  else if (grepl("morphogenesis|development|pattern formation|body", term, ignore.case = TRUE))
    return("Body Morphogenesis")
  else if (grepl("autophagy|lysosome|vacuole|phagosome", term, ignore.case = TRUE))
    return("Autophagy")
  else if (grepl("DNA damage|DNA repair|genome stability", term, ignore.case = TRUE))
    return("DNA Damage")
  else if (grepl("oxidative stress|reactive oxygen|antioxidant|ROS|redox", term, ignore.case = TRUE))
    return("Oxidative Stress")
  else if (grepl("amino acid|protein synthesis|translation|protein folding", term, ignore.case = TRUE))
    return("Amino Acid Biosynthesis")
  else if (grepl("mitochondri|respiration", term, ignore.case = TRUE))
    return("Mitochondrial Expression")
  else if (grepl("carbohydrate|glucose|glycolysis|gluconeogenesis|glycogen|sugar", term, ignore.case = TRUE))
    return("Carbohydrate Metabolism")
  else if (grepl("oxidative phosphorylation|electron transport|respiratory chain|ETC", term, ignore.case = TRUE))
    return("OXPHOS")
  else if (grepl("TCA|tricarboxylic|citrate|Krebs|citric acid cycle", term, ignore.case = TRUE))
    return("TCA Cycle")
  else
    return("Other")
})

# If we have results, create the plots
if (nrow(filtered_results) > 0) {
  # Define the pathway order for visualization
  pathway_order <- c(
    "Immune Response", 
    "Body Morphogenesis", 
    "Autophagy", 
    "DNA Damage", 
    "Oxidative Stress", 
    "Amino Acid Biosynthesis", 
    "Mitochondrial Expression", 
    "Carbohydrate Metabolism", 
    "OXPHOS", 
    "TCA Cycle",
    "Other"
  )
  
  # Get only the pathways present in our data
  pathways_in_data <- unique(filtered_results$Pathway)
  sorted_pathways <- pathway_order[pathway_order %in% pathways_in_data]
  
  # Order factor levels
  filtered_results$Pathway <- factor(filtered_results$Pathway, levels = sorted_pathways)
  
  # Sort terms by pathway and significance
  filtered_results <- filtered_results %>%
    arrange(Pathway, p.adjust)
  
  # Update the compareCluster result with our filtered and categorized data
  compare_result@compareClusterResult <- filtered_results
  
  # Print a summary of filtered terms
  cat("\nFiltered GO terms by pathway:\n")
  pathway_summary <- filtered_results %>%
    group_by(Pathway, Cluster) %>%
    summarize(
      Terms = n(),
      TopTerm = Description[which.min(p.adjust)],
      BestPadj = min(p.adjust),
      .groups = "drop"
    )
  
  print(pathway_summary)
  
  # Create a dotplot with facets by pathway
  dotplot_facet <- dotplot(
    compare_result,
    showCategory = 50,
    title = "GO Enrichment: Biological Processes in Matched vs. Mismatched",
    font.size = 10
  ) +
    theme_bw() +
    facet_grid(Pathway ~ ., scales = "free_y", space = "free_y") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 8)
    ) +
    scale_color_gradientn(
      colors = rev(brewer.pal(9, "RdYlBu")),
      name = "p.adjust"
    )
  
  # Save the faceted dotplot
  ggsave(
    "GO_enrichment_faceted_by_pathway.png",
    dotplot_facet,
    width = 12,
    height = max(10, length(unique(filtered_results$Description)) * 0.25), # Dynamic height
    dpi = 300
  )
  
  # Create a simplified dotplot without facets
  dotplot_simple <- dotplot(
    compare_result,
    showCategory = 30,
    title = "GO Enrichment: Matched vs. Mismatched",
    font.size = 10
  ) +
    theme_bw()
  
  ggsave(
    "GO_enrichment_simple.png",
    dotplot_simple,
    width = 10,
    height = 12,
    dpi = 300
  )
  
  # Create an enrichment map to show relationships
  set.seed(123) # For reproducible layout
  emap <- emapplot(
    pairwise_termsim(compare_result),
    showCategory = 25, 
    color = "p.adjust",
    layout = "kk" # Force-directed layout
  )
  
  ggsave(
    "GO_enrichment_map.png",
    emap,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  # Create separate analyses for each condition to get more detailed information
  
  # Upregulated in Mismatched
  go_up_mismatched <- enrichGO(
    gene = entrez_up_mismatched,
    OrgDb = org.Dm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # Upregulated in Matched
  go_up_matched <- enrichGO(
    gene = entrez_up_matched,
    OrgDb = org.Dm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # Filter these results for our pathways of interest
  go_up_mismatched_results <- as.data.frame(go_up_mismatched)
  go_up_matched_results <- as.data.frame(go_up_matched)
  
  filtered_mismatched <- go_up_mismatched_results %>%
    filter(sapply(Description, contains_key_process))
  
  filtered_matched <- go_up_matched_results %>%
    filter(sapply(Description, contains_key_process))
  
  # Add pathway categories
  filtered_mismatched$Pathway <- sapply(filtered_mismatched$Description, function(term) {
    if (grepl("immune|inflammatory|defense|leukocyte|cytokine", term, ignore.case = TRUE))
      return("Immune Response")
    else if (grepl("morphogenesis|development|pattern formation|body", term, ignore.case = TRUE))
      return("Body Morphogenesis")
    else if (grepl("autophagy|lysosome|vacuole|phagosome", term, ignore.case = TRUE))
      return("Autophagy")
    else if (grepl("DNA damage|DNA repair|genome stability", term, ignore.case = TRUE))
      return("DNA Damage")
    else if (grepl("oxidative stress|reactive oxygen|antioxidant|ROS|redox", term, ignore.case = TRUE))
      return("Oxidative Stress")
    else if (grepl("amino acid|protein synthesis|translation|protein folding", term, ignore.case = TRUE))
      return("Amino Acid Biosynthesis")
    else if (grepl("mitochondri|respiration", term, ignore.case = TRUE))
      return("Mitochondrial Expression")
    else if (grepl("carbohydrate|glucose|glycolysis|gluconeogenesis|glycogen|sugar", term, ignore.case = TRUE))
      return("Carbohydrate Metabolism")
    else if (grepl("oxidative phosphorylation|electron transport|respiratory chain|ETC", term, ignore.case = TRUE))
      return("OXPHOS")
    else if (grepl("TCA|tricarboxylic|citrate|Krebs|citric acid cycle", term, ignore.case = TRUE))
      return("TCA Cycle")
    else
      return("Other")
  })
  
  filtered_matched$Pathway <- sapply(filtered_matched$Description, function(term) {
    if (grepl("immune|inflammatory|defense|leukocyte|cytokine", term, ignore.case = TRUE))
      return("Immune Response")
    else if (grepl("morphogenesis|development|pattern formation|body", term, ignore.case = TRUE))
      return("Body Morphogenesis")
    else if (grepl("autophagy|lysosome|vacuole|phagosome", term, ignore.case = TRUE))
      return("Autophagy")
    else if (grepl("DNA damage|DNA repair|genome stability", term, ignore.case = TRUE))
      return("DNA Damage")
    else if (grepl("oxidative stress|reactive oxygen|antioxidant|ROS|redox", term, ignore.case = TRUE))
      return("Oxidative Stress")
    else if (grepl("amino acid|protein synthesis|translation|protein folding", term, ignore.case = TRUE))
      return("Amino Acid Biosynthesis")
    else if (grepl("mitochondri|respiration", term, ignore.case = TRUE))
      return("Mitochondrial Expression")
    else if (grepl("carbohydrate|glucose|glycolysis|gluconeogenesis|glycogen|sugar", term, ignore.case = TRUE))
      return("Carbohydrate Metabolism")
    else if (grepl("oxidative phosphorylation|electron transport|respiratory chain|ETC", term, ignore.case = TRUE))
      return("OXPHOS")
    else if (grepl("TCA|tricarboxylic|citrate|Krebs|citric acid cycle", term, ignore.case = TRUE))
      return("TCA Cycle")
    else
      return("Other")
  })
  
  # Save detailed results to Excel
  write_xlsx(
    list(
      "All_Filtered_GO_Terms" = filtered_results,
      "Upregulated_in_Mismatched" = filtered_mismatched,
      "Upregulated_in_Matched" = filtered_matched,
      "Pathway_Summary" = pathway_summary
    ),
    path = "GO_enrichment_analysis_results.xlsx"
  )
  
  # Create barplots for top terms in each condition
  if (nrow(filtered_mismatched) > 0) {
    # For mismatched
    filtered_mismatched <- filtered_mismatched %>%
      arrange(Pathway, p.adjust)
    
    # Get top 3 terms per pathway
    top_mismatched <- filtered_mismatched %>%
      group_by(Pathway) %>%
      slice_head(n = 3) %>%
      ungroup()
    
    # Create barplot
    if (nrow(top_mismatched) > 0) {
      barplot_mismatched <- ggplot(top_mismatched, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Pathway)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_bw() +
        labs(
          title = "Top GO Terms Upregulated in Mismatched (FA)",
          x = NULL,
          y = "-log10(p.adjust)"
        ) +
        scale_fill_brewer(palette = "Set2") +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 9)
        )
      
      ggsave(
        "GO_top_terms_mismatched.png",
        barplot_mismatched,
        width = 10,
        height = max(6, nrow(top_mismatched) * 0.3),
        dpi = 300
      )
    }
  }
  
  if (nrow(filtered_matched) > 0) {
    # For matched
    filtered_matched <- filtered_matched %>%
      arrange(Pathway, p.adjust)
    
    # Get top 3 terms per pathway
    top_matched <- filtered_matched %>%
      group_by(Pathway) %>%
      slice_head(n = 3) %>%
      ungroup()
    
    # Create barplot
    if (nrow(top_matched) > 0) {
      barplot_matched <- ggplot(top_matched, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Pathway)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_bw() +
        labs(
          title = "Top GO Terms Upregulated in Matched (FF+AA)",
          x = NULL,
          y = "-log10(p.adjust)"
        ) +
        scale_fill_brewer(palette = "Set3") +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 9)
        )
      
      ggsave(
        "GO_top_terms_matched.png",
        barplot_matched,
        width = 10,
        height = max(6, nrow(top_matched) * 0.3),
        dpi = 300
      )
    }
  }
  
  # Create a heatmap visualization
  # First, reshape the data to create a heatmap
  top_terms_heatmap <- filtered_results %>%
    group_by(Pathway) %>%
    slice_min(order_by = p.adjust, n = 3) %>%  # Top 3 terms per pathway
    ungroup()
  
  # Create the heatmap
  heatmap_plot <- ggplot(top_terms_heatmap, 
                         aes(x = Cluster, y = Description, fill = -log10(p.adjust))) +
    geom_tile(color = "white") +
    facet_grid(Pathway ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors = c("#FFFFFF", "#FFF7BC", "#FD8D3C", "#E31A1C", "#800026"),
      name = "-log10(p.adjust)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(
      title = "Top GO Terms by Pathway",
      subtitle = "Comparing Matched vs. Mismatched",
      x = NULL,
      y = NULL
    ) +
    geom_text(aes(label = Count), size = 3)
  
  ggsave(
    "GO_heatmap_by_pathway.png",
    heatmap_plot,
    width = 10,
    height = max(8, nrow(top_terms_heatmap) * 0.4),
    dpi = 300
  )
  
} else {
  cat("\nNo terms found that match your biological processes of interest.\n")
  cat("Try broadening your key process terms or check your enrichment results.\n")
}

