# Analysis comparing timepoint 5 vs timepoint 3 in FA group
# Including differential expression, enrichment analysis, and visualization

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(RColorBrewer)
library(ggrepel)

#
# Part 1: Filter data and perform differential expression analysis
#

# Filter metadata to include only FA group, timepoints 3 and 5
timepoint_comparison <- function() {
  # Check if we have the necessary data
  if (!exists("metadata") || !exists("gene_count")) {
    cat("Required data (metadata, gene_count) not found.\n")
    return(NULL)
  }
  
  # Filter metadata
  metadata_filtered <- metadata %>%
    filter(group == "FA" & time %in% c(3, 5))
  
  # Check if we have enough samples
  if (nrow(metadata_filtered) < 2) {
    cat("Not enough FA samples at timepoints 3 and 5.\n")
    return(NULL)
  }
  
  cat("Filtered data includes", nrow(metadata_filtered), "samples.\n")
  
  # Get the sample IDs for the filtered metadata
  samples_to_keep <- metadata_filtered$sample_id
  
  # Filter gene count matrix to keep only the samples in the filtered metadata
  gene_count_filtered <- gene_count %>%
    select(gene_id, all_of(samples_to_keep))
  
  # Make sure gene_id is a character vector (not a factor) before using as rownames
  gene_count_filtered$gene_id <- as.character(gene_count_filtered$gene_id)
  
  # Check for duplicate gene_ids and handle them if present
  if (any(duplicated(gene_count_filtered$gene_id))) {
    cat("Warning: Duplicate gene_ids found. Adding suffix to make them unique.\n")
    # Create a frequency table of gene_ids
    gene_freq <- table(gene_count_filtered$gene_id)
    # Find duplicated gene_ids
    dup_genes <- names(gene_freq[gene_freq > 1])
    
    # For each duplicated gene_id, add a suffix
    for (gene in dup_genes) {
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
  
  # Convert timepoint to factor for DESeq2
  metadata_filtered$time <- factor(metadata_filtered$time)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = gene_count_matrix,
    colData = metadata_filtered,
    design = ~ time
  )
  
  # Filter low count genes
  keep <- rowSums(counts(dds) >= 10) >= min(3, ncol(dds))
  dds <- dds[keep,]
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Get results for timepoint 5 vs 3
  res <- results(dds, contrast = c("time", "5", "3"))
  
  # Add gene symbols to results
  gene_symbols <- mapIds(org.Dm.eg.db, 
                         keys = rownames(res),
                         column = "SYMBOL", 
                         keytype = "FLYBASE",
                         multiVals = "first")
  
  # For missing gene symbols, use the gene ID
  gene_symbols[is.na(gene_symbols)] <- rownames(res)[is.na(gene_symbols)]
  
  # Add gene symbols to results
  res$symbol <- gene_symbols
  
  # Return results
  return(list(dds = dds, res = res, metadata = metadata_filtered))
}

# Run timepoint comparison
tp_results <- timepoint_comparison()

if (!is.null(tp_results)) {
  # Extract results and DESeq2 object
  dds <- tp_results$dds
  res <- tp_results$res
  
  # Summary of results
  cat("\nSummary of differential expression results (FA timepoint 5 vs 3):\n")
  print(summary(res))
  
  # Get significant genes
  padj_cutoff <- 0.05
  sig_genes <- subset(res, padj < padj_cutoff)
  
  cat("\nFound", nrow(sig_genes), "significantly differentially expressed genes (padj <", padj_cutoff, ")\n")
  
  if (nrow(sig_genes) > 0) {
    # Print top significant genes
    cat("\nTop significantly differentially expressed genes:\n")
    sig_genes_ordered <- sig_genes[order(sig_genes$padj), ]
    top_genes <- head(sig_genes_ordered, 20)
    print(data.frame(
      gene = rownames(top_genes),
      symbol = top_genes$symbol,
      log2FoldChange = top_genes$log2FoldChange,
      padj = top_genes$padj
    ))
    
    #
    # Part 2: Enrichment Analysis
    #
    
    # Function to perform enrichment analysis
    run_enrichment <- function(res, padj_cutoff = 0.05, log2FC_cutoff = 0) {
      # Get significant genes
      sig_genes <- subset(res, padj < padj_cutoff & abs(log2FoldChange) > log2FC_cutoff)
      
      if (nrow(sig_genes) == 0) {
        cat("No significant genes found for enrichment analysis.\n")
        return(NULL)
      }
      
      # Get gene IDs
      gene_ids <- rownames(sig_genes)
      
      # Convert FlyBase IDs to ENTREZ IDs for functional enrichment analysis
      id_mappings <- bitr(gene_ids, 
                          fromType = "FLYBASE", 
                          toType = c("ENTREZID", "SYMBOL"),
                          OrgDb = org.Dm.eg.db)
      
      # Check if mapping was successful
      if (nrow(id_mappings) == 0) {
        cat("No genes could be mapped to ENTREZ IDs for enrichment analysis.\n")
        return(NULL)
      }
      
      # Create named vector for fold changes
      gene_list <- sig_genes$log2FoldChange
      names(gene_list) <- rownames(sig_genes)
      
      # Keep only mapped genes
      mapped_gene_ids <- id_mappings$FLYBASE
      mapped_entrez_ids <- id_mappings$ENTREZID
      
      # Perform GO enrichment analysis
      go_bp <- enrichGO(gene = mapped_entrez_ids,
                        OrgDb = org.Dm.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
      
      # Perform KEGG pathway enrichment
      kegg_result <- enrichKEGG(gene = mapped_entrez_ids,
                                organism = 'dme',
                                keyType = 'ncbi-geneid',
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")
      
      # Return results
      return(list(
        go_bp = go_bp,
        kegg = kegg_result,
        id_mappings = id_mappings,
        gene_list = gene_list
      ))
    }
    
    # Run enrichment analysis
    enrichment_results <- run_enrichment(res, padj_cutoff = 0.05)
    
    # Check if enrichment analysis was successful
    if (!is.null(enrichment_results)) {
      # Extract enrichment results
      go_bp <- enrichment_results$go_bp
      kegg_result <- enrichment_results$kegg
      
      # Print GO enrichment results
      cat("\nGO Biological Process enrichment results:\n")
      if (nrow(go_bp@result) > 0) {
        print(head(go_bp@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
      } else {
        cat("No significant GO terms found.\n")
      }
      
      # Print KEGG enrichment results
      cat("\nKEGG pathway enrichment results:\n")
      if (nrow(kegg_result@result) > 0) {
        print(head(kegg_result@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
      } else {
        cat("No significant KEGG pathways found.\n")
      }
      
      #
      # Part 3: Visualizations
      #
      
      # Volcano plot
      volcano_plot <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = padj < 0.05), size = 1.5, alpha = 0.7) +
        scale_color_manual(values = c("grey", "red"), name = "Significant", labels = c("No", "Yes")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "Volcano Plot - FA Timepoint 5 vs 3",
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom"
        )
      
      # Add gene symbols for top genes
      top_genes_volcano <- rbind(
        head(res[order(res$padj), ], 10),  # Top significant 
        head(res[order(res$log2FoldChange, decreasing = TRUE), ], 5),  # Top upregulated
        head(res[order(res$log2FoldChange), ], 5)  # Top downregulated
      )
      
      top_genes_df <- as.data.frame(top_genes_volcano)
      top_genes_df$symbol <- top_genes_volcano$symbol
      
      volcano_plot <- volcano_plot +
        geom_text_repel(
          data = top_genes_df,
          aes(x = log2FoldChange, y = -log10(padj), label = symbol),
          max.overlaps = 20,
          size = 3,
          box.padding = 0.5,
          point.padding = 0.3,
          force = 3
        )
      
      print(volcano_plot)
      
      # GO enrichment visualization
      if (nrow(go_bp@result) > 0) {
        # Dotplot
        dotplot_go <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment")
        print(dotplot_go)
        
        # Enrichment map
        if (nrow(go_bp@result) >= 3) {
          ego_map <- pairwise_termsim(go_bp)
          emap_plot <- emapplot(ego_map, showCategory = min(30, nrow(go_bp@result)))
          print(emap_plot)
        }
        
        # Category-gene network plot
        cnet_plot <- cnetplot(go_bp, 
                              showCategory = 5,
                              categorySize = "pvalue", 
                              node_label = "all")
        print(cnet_plot)
      }
      
      # KEGG pathway visualization
      if (nrow(kegg_result@result) > 0) {
        # Dotplot
        dotplot_kegg <- dotplot(kegg_result, showCategory = 15, title = "KEGG Pathway Enrichment")
        print(dotplot_kegg)
        
        # Category-gene network plot
        cnet_plot_kegg <- cnetplot(kegg_result, 
                                   showCategory = 5,
                                   categorySize = "pvalue", 
                                   node_label = "all")
        print(cnet_plot_kegg)
      }
      
      # Heatmap of top differentially expressed genes
      # Get normalized counts
      vsd <- vst(dds, blind = FALSE)
      
      # Select top differentially expressed genes
      top_de_genes <- rownames(sig_genes_ordered)[1:min(50, nrow(sig_genes_ordered))]
      
      # Extract counts for these genes
      top_gene_counts <- assay(vsd)[top_de_genes, ]
      
      # Add gene symbols as row names
      rownames(top_gene_counts) <- sig_genes_ordered$symbol[1:length(top_de_genes)]
      
      # Create annotation for heatmap
      anno_col <- data.frame(
        Timepoint = factor(tp_results$metadata$time),
        row.names = colnames(top_gene_counts)
      )
      
      # Generate heatmap
      pheatmap(
        top_gene_counts,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = anno_col,
        scale = "row",
        main = "Top Differentially Expressed Genes - FA T5 vs T3",
        fontsize_row = 8,
        fontsize_col = 8
      )
      
      # Respiration genes analysis
      # Define list of respiration genes
      respiration_genes <- c(
        # Electron Transport Chain Complex I (NADH:ubiquinone oxidoreductase)
        "ND75", "ND42", "ND51", "ND23", "ND24", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
        "NDUFS2", "NDUFS3", "NDUFA8", "NDUFV1", 
        
        # Complex II (Succinate dehydrogenase)
        "SdhA", "SdhB", "SdhC", "SdhD", 
        
        # Complex III (Ubiquinol-cytochrome c reductase)
        "UQCR-C1", "UQCR-C2", "UQCR-Q", "Ucrh", "Cyc1", 
        
        # Complex IV (Cytochrome c oxidase)
        "COX1", "COX2", "COX3", "COX4", "COX5A", "COX5B", "COX6A", "COX6B", "COX7A", "COX7C", "COX8",
        
        # Complex V (ATP synthase)
        "ATPsynα", "ATPsynalpha", "ATPsyn-alpha", "ATP5A", 
        "ATPsynβ", "ATPsynbeta", "ATPsyn-beta", "ATP5B",
        "ATPsynγ", "ATPsyngamma", "ATPsyn-gamma", "ATP5C1",
        "ATPsynδ", "ATPsyndelta", "ATPsyn-delta", "ATP5D",
        "ATPsynε", "ATPsynepsilon", "ATPsyn-epsilon", "ATP5E",
        "ATPsynB", "ATPsynC", "ATPsynD", "ATPsynE", "ATPsynF", "ATPsynG", "ATPsynO",
        
        # TCA cycle genes
        "Idh", "Mdh1", "Mdh2", "SdhA", "SdhB", "Acon", "Fum", "l(1)G0255", "Pdha", "Pdhb", "Dlat", "Dld",
        
        # Mitochondrial DNA maintenance 
        "tam", "mtSSB", "Pol32", "DNApol-γ35", "MtTFB2", "MtTFB1", "TFAM", "polG", "polG2", "polGA", "polGB",
        
        # Key oxygen sensors and respiration regulators
        "sima", "tango", "Hph", "Alas", "Pdp1", "dNRF2", "cbt", "HIF1A",
        
        # Mitochondrial transporters
        "sesB", "Ant2", "Mfn", "Marf", "Drp1", "Fis1", "Opa1",
        
        # Other important respiratory genes
        "Cyp", "mRpL1", "mRpL2", "mRpL3", "mRpL4", "mRpS5", "mRpS6", "sesB", "bou", "Pdsw",
        "mt:ND1", "mt:ND2", "mt:ND3", "mt:ND4", "mt:ND5", "mt:ND6", "mt:Cyt-b", "mt:CoI", "mt:CoII", "mt:CoIII",
        "mt:ATPase6", "mt:ATPase8"
      )
      
      # Check which respiration genes are differentially expressed
      resp_genes_in_results <- res$symbol %in% respiration_genes
      resp_gene_results <- res[resp_genes_in_results, ]
      
      # Filter for significant respiration genes
      sig_resp_genes <- subset(resp_gene_results, padj < 0.05)
      
      if (nrow(sig_resp_genes) > 0) {
        cat("\nSignificantly differentially expressed respiration genes:\n")
        sig_resp_df <- data.frame(
          gene_id = rownames(sig_resp_genes),
          symbol = sig_resp_genes$symbol,
          log2FoldChange = sig_resp_genes$log2FoldChange,
          padj = sig_resp_genes$padj,
          stringsAsFactors = FALSE
        )
        print(sig_resp_df[order(sig_resp_df$padj), ])
        
        # Create a specialized volcano plot for respiration genes
        resp_volcano <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
          # Add background of all genes
          geom_point(color = "grey85", size = 1, alpha = 0.5) +
          # Add respiration genes with special highlighting
          geom_point(data = as.data.frame(resp_gene_results), 
                     aes(color = padj < 0.05, size = padj < 0.05)) +
          scale_color_manual(values = c("blue", "red"), name = "Significant", labels = c("No", "Yes")) +
          scale_size_manual(values = c(2, 3), name = "Significant", labels = c("No", "Yes")) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
          labs(title = "Respiration Genes - FA Timepoint 5 vs 3",
               x = "Log2 Fold Change",
               y = "-Log10 Adjusted P-value") +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "bottom"
          )
        
        # Add gene labels for respiration genes
        resp_volcano <- resp_volcano +
          geom_text_repel(
            data = as.data.frame(resp_gene_results),
            aes(x = log2FoldChange, y = -log10(padj), label = symbol),
            max.overlaps = 50,
            size = 3,
            box.padding = 0.5,
            point.padding = 0.3,
            force = 3
          )
        
        print(resp_volcano)
      } else {
        cat("\nNo significantly differentially expressed respiration genes found.\n")
      }
    } else {
      cat("Enrichment analysis did not produce results.\n")
    }
  } else {
    cat("No significantly differentially expressed genes found.\n")
  }
} else {
  cat("Analysis could not be completed.\n")
}

# Detailed Analysis of Upregulation in T5 vs T3 in FA Group
# This code extends your existing analysis with more focused direction-specific analysis

# Load required libraries if not already loaded
if (!requireNamespace("DESeq2", quietly = TRUE)) library(DESeq2)
if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
if (!requireNamespace("pheatmap", quietly = TRUE)) library(pheatmap)
if (!requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) library(clusterProfiler)
if (!requireNamespace("enrichplot", quietly = TRUE)) library(enrichplot)
if (!requireNamespace("org.Dm.eg.db", quietly = TRUE)) library(org.Dm.eg.db)
if (!requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
if (!requireNamespace("RColorBrewer", quietly = TRUE)) library(RColorBrewer)
if (!requireNamespace("VennDiagram", quietly = TRUE)) library(VennDiagram)
if (!requireNamespace("UpSetR", quietly = TRUE)) library(UpSetR)

# Assuming tp_results object is already available from your previous code
# If not, run your timepoint_comparison() function first to get tp_results

# PART 1: Separate and analyze up and down regulated genes
#----------------------------------------------------------------------------

# Function to analyze DEGs by direction and create detailed reports
analyze_degs_by_direction <- function(res, padj_cutoff = 0.05, log2FC_cutoff = 0.5) {
  # Ensure results are in data frame format
  res_df <- as.data.frame(res)
  
  # Add a column for regulation direction
  res_df$regulation <- "Not significant"
  res_df$regulation[res_df$padj < padj_cutoff & res_df$log2FoldChange > log2FC_cutoff] <- "Upregulated in T5"
  res_df$regulation[res_df$padj < padj_cutoff & res_df$log2FoldChange < -log2FC_cutoff] <- "Upregulated in T3"
  
  # Convert to factor with specific order for plotting
  res_df$regulation <- factor(res_df$regulation, 
                              levels = c("Upregulated in T5", "Not significant", "Upregulated in T3"))
  
  # Count genes in each category
  up_t5 <- sum(res_df$regulation == "Upregulated in T5", na.rm = TRUE)
  up_t3 <- sum(res_df$regulation == "Upregulated in T3", na.rm = TRUE)
  not_sig <- sum(res_df$regulation == "Not significant", na.rm = TRUE)
  
  cat("Classification of differentially expressed genes:\n")
  cat("- Upregulated in T5 (log2FC >", log2FC_cutoff, "and padj <", padj_cutoff, "):", up_t5, "genes\n")
  cat("- Upregulated in T3 (log2FC < -", log2FC_cutoff, "and padj <", padj_cutoff, "):", up_t3, "genes\n")
  cat("- Not significantly differentially expressed:", not_sig, "genes\n")
  
  # Filter for up and down regulated genes
  up_t5_genes <- rownames(res_df[res_df$regulation == "Upregulated in T5", ])
  up_t3_genes <- rownames(res_df[res_df$regulation == "Upregulated in T3", ])
  
  # Get the gene symbols
  up_t5_symbols <- res_df$symbol[res_df$regulation == "Upregulated in T5"]
  up_t3_symbols <- res_df$symbol[res_df$regulation == "Upregulated in T3"]
  
  # Return the results
  list(
    complete_results = res_df,
    up_t5_genes = up_t5_genes,
    up_t3_genes = up_t3_genes,
    up_t5_symbols = up_t5_symbols,
    up_t3_symbols = up_t3_symbols
  )
}

# Run the analysis
direction_analysis <- analyze_degs_by_direction(tp_results$res, padj_cutoff = 0.05, log2FC_cutoff = 0.5)

# Get the dataframe with complete results
res_with_direction <- direction_analysis$complete_results

# PART 2: Enhanced Visualizations for Direction Analysis
#----------------------------------------------------------------------------

# Create enhanced volcano plot with direction colors
enhanced_volcano <- ggplot(res_with_direction, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated in T5" = "red", 
                                "Not significant" = "grey", 
                                "Upregulated in T3" = "blue")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.5, 0, 0.5), linetype = "dashed", color = c("grey70", "grey50", "grey70")) +
  labs(title = "Enhanced Volcano Plot - FA Timepoint 5 vs 3",
       subtitle = paste0("Red: Upregulated in T5, Blue: Upregulated in T3\n",
                         "padj < 0.05, |log2FC| > 0.5"),
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

# Top genes to label (combining both directions)
top_t5_genes <- head(res_with_direction[res_with_direction$regulation == "Upregulated in T5", ], 15)
top_t3_genes <- head(res_with_direction[res_with_direction$regulation == "Upregulated in T3", ], 15)
top_genes_combined <- rbind(top_t5_genes, top_t3_genes)

# Add gene labels to volcano plot
enhanced_volcano_labeled <- enhanced_volcano +
  geom_text_repel(
    data = top_genes_combined,
    aes(label = symbol, color = regulation),
    max.overlaps = 30,
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 3
  )

print(enhanced_volcano_labeled)

# MA Plot with direction colors
ma_plot <- ggplot(res_with_direction, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = regulation), size = 1.2, alpha = 0.8) +
  scale_x_log10() +
  scale_color_manual(values = c("Upregulated in T5" = "red", 
                                "Not significant" = "grey", 
                                "Upregulated in T3" = "blue")) +
  geom_hline(yintercept = c(-0.5, 0, 0.5), linetype = "dashed", color = c("grey70", "grey50", "grey70")) +
  labs(title = "MA Plot - FA Timepoint 5 vs 3",
       subtitle = "Mean expression vs Log2 Fold Change",
       x = "Mean Expression (log10 scale)",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

print(ma_plot)

# PART 3: Direction-Specific Enrichment Analysis
#----------------------------------------------------------------------------

# Function to perform enrichment analysis for a gene list
direction_enrichment <- function(gene_ids, gene_symbols, direction_name) {
  cat("\n\n======================================================\n")
  cat("Enrichment Analysis for", direction_name, "\n")
  cat("======================================================\n")
  
  cat("Number of genes:", length(gene_ids), "\n\n")
  
  # Print top genes
  cat("Top genes in this category:\n")
  top_genes_df <- data.frame(
    FlyBase_ID = gene_ids[1:min(20, length(gene_ids))],
    Symbol = gene_symbols[1:min(20, length(gene_ids))]
  )
  print(top_genes_df)
  
  # If no genes, return early
  if (length(gene_ids) == 0) {
    cat("No genes available for enrichment analysis.\n")
    return(NULL)
  }
  
  # Convert FlyBase IDs to ENTREZ IDs for functional enrichment analysis
  id_mappings <- tryCatch({
    bitr(gene_ids, 
         fromType = "FLYBASE", 
         toType = c("ENTREZID", "SYMBOL"),
         OrgDb = org.Dm.eg.db)
  }, error = function(e) {
    cat("Error in ID mapping:", e$message, "\n")
    return(NULL)
  })
  
  # Check if mapping was successful
  if (is.null(id_mappings) || nrow(id_mappings) == 0) {
    cat("No genes could be mapped to ENTREZ IDs for enrichment analysis.\n")
    return(NULL)
  }
  
  cat("Successfully mapped", nrow(id_mappings), "out of", length(gene_ids), "genes to ENTREZ IDs.\n\n")
  
  # Get ENTREZ IDs
  entrez_ids <- id_mappings$ENTREZID
  
  # Perform GO enrichment analysis - Biological Process
  go_bp <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Dm.eg.db,
             keyType = "ENTREZID",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.1)
  }, error = function(e) {
    cat("Error in GO BP enrichment:", e$message, "\n")
    return(NULL)
  })
  
  # Perform GO enrichment analysis - Molecular Function
  go_mf <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Dm.eg.db,
             keyType = "ENTREZID",
             ont = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.1)
  }, error = function(e) {
    cat("Error in GO MF enrichment:", e$message, "\n")
    return(NULL)
  })
  
  # Perform GO enrichment analysis - Cellular Component
  go_cc <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Dm.eg.db,
             keyType = "ENTREZID",
             ont = "CC",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.1)
  }, error = function(e) {
    cat("Error in GO CC enrichment:", e$message, "\n")
    return(NULL)
  })
  
  # Perform KEGG pathway enrichment
  kegg_result <- tryCatch({
    enrichKEGG(gene = entrez_ids,
               organism = 'dme',
               keyType = 'ncbi-geneid',
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
  }, error = function(e) {
    cat("Error in KEGG enrichment:", e$message, "\n")
    return(NULL)
  })
  
  # Print GO BP results
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    cat("\nGO Biological Process enrichment results:\n")
    print(head(go_bp@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
    
    # Visualize top GO terms
    dotplot(go_bp, title = paste0("GO BP Enrichment - ", direction_name), showCategory = 15)
  } else {
    cat("\nNo significant GO Biological Process terms found.\n")
  }
  
  # Print GO MF results
  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    cat("\nGO Molecular Function enrichment results:\n")
    print(head(go_mf@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
    
    # Visualize top GO terms
    dotplot(go_mf, title = paste0("GO MF Enrichment - ", direction_name), showCategory = 15)
  } else {
    cat("\nNo significant GO Molecular Function terms found.\n")
  }
  
  # Print GO CC results
  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    cat("\nGO Cellular Component enrichment results:\n")
    print(head(go_cc@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
    
    # Visualize top GO terms
    dotplot(go_cc, title = paste0("GO CC Enrichment - ", direction_name), showCategory = 15)
  } else {
    cat("\nNo significant GO Cellular Component terms found.\n")
  }
  
  # Print KEGG results
  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    cat("\nKEGG pathway enrichment results:\n")
    print(head(kegg_result@result[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
    
    # Visualize top KEGG pathways
    dotplot(kegg_result, title = paste0("KEGG Pathway Enrichment - ", direction_name), showCategory = 15)
  } else {
    cat("\nNo significant KEGG pathways found.\n")
  }
  
  # Return results
  return(list(
    go_bp = go_bp,
    go_mf = go_mf,
    go_cc = go_cc,
    kegg = kegg_result,
    mapped_genes = id_mappings
  ))
}

# Run enrichment for upregulated in T5
up_t5_enrichment <- direction_enrichment(
  direction_analysis$up_t5_genes, 
  direction_analysis$up_t5_symbols,
  "Upregulated in T5"
)
dotplot(up_t5_enrichment$go_bp)


# Run enrichment for upregulated in T3
up_t3_enrichment <- direction_enrichment(
  direction_analysis$up_t3_genes, 
  direction_analysis$up_t3_symbols,
  "Upregulated in T3"
)





# PART 4: Specialized Analysis for Metabolism and Respiration
#----------------------------------------------------------------------------

# Define broader lists of functional gene categories
metabolism_genes <- c(
  # Respiration genes from your original code
  "ND75", "ND42", "ND51", "ND23", "ND24", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
  "NDUFS2", "NDUFS3", "NDUFA8", "NDUFV1", 
  "SdhA", "SdhB", "SdhC", "SdhD", 
  "UQCR-C1", "UQCR-C2", "UQCR-Q", "Ucrh", "Cyc1", 
  "COX1", "COX2", "COX3", "COX4", "COX5A", "COX5B", "COX6A", "COX6B", "COX7A", "COX7C", "COX8",
  "ATPsynα", "ATPsynalpha", "ATPsyn-alpha", "ATP5A", 
  "ATPsynβ", "ATPsynbeta", "ATPsyn-beta", "ATP5B",
  "ATPsynγ", "ATPsyngamma", "ATPsyn-gamma", "ATP5C1",
  "ATPsynδ", "ATPsyndelta", "ATPsyn-delta", "ATP5D",
  "ATPsynε", "ATPsynepsilon", "ATPsyn-epsilon", "ATP5E",
  "ATPsynB", "ATPsynC", "ATPsynD", "ATPsynE", "ATPsynF", "ATPsynG", "ATPsynO",
  "Idh", "Mdh1", "Mdh2", "SdhA", "SdhB", "Acon", "Fum", "l(1)G0255", "Pdha", "Pdhb", "Dlat", "Dld",
  
  # Add broader metabolism categories
  # Glycolysis
  "Hex-A", "Hex-C", "Pfk", "Ald", "Tpi", "Gapdh1", "Gapdh2", "Pgk", "Pglym78", "Eno", "PyK",
  
  # Fatty acid metabolism
  "ACC", "FASN1", "FASN2", "ACSL", "CPT1", "CPT2", "ACOX1", "ACADL", "ACADM", "ACADS", 
  "ECHS1", "HADHA", "HADHB", "ACAA1", "ACAA2",
  
  # Pentose phosphate pathway
  "G6pd", "Pgd", "Tkt", "Tal", "Rpe", "Rpi",
  
  # Amino acid metabolism
  "Got1", "Got2", "Gdh", "Ast", "Gs1", "Gs2", "ProDH", "ArgK", "Argl", "Oat",
  
  # Lipid metabolism
  "Lsd-1", "Lsd-2", "bmm", "Lsd", "Lip4", "mdy", "ATGL", "Lpin", "Dgat1", "Dgat2",
  
  # Gluconeogenesis 
  "Pepck", "Fbp", "G6pc",
  
  # Insulin/TOR signaling
  "InR", "chico", "Pi3K92E", "Akt1", "foxo", "Thor", "S6k", "Tor", "Tsc1", "Tsc2", "Rheb",
  
  # Mitochondrial dynamics
  "Marf", "Opa1", "Drp1", "Fis1", "Mfn", "Mfn2", "MitoPLD"
)

oxidative_stress_genes <- c(
  # Antioxidant enzymes
  "Sod1", "Sod2", "Cat", "Trx-2", "TrxR-1", "TrxR-2", "Prx2540-1", "Prx2540-2", "Prx3", "Prx5", "Prx6005",
  "Gpx1", "Gst", "GstD1", "GstS1", "GstE1", "GstT1", "GstZ1", "GstO1", "Gr", "Trx",
  
  # Transcription factors involved in oxidative stress response
  "Cnc", "Keap1", "Hr96", "dFOXO", "Mtf-1", "Hsf", "Nrf2", "ATF-2",
  
  # Metabolite sensors and stress response
  "Thor", "JNK", "hep", "bsk", "p38a", "p38b", "Atf3", "Gadd45", "DHR96", "Hsp22", "Hsp23", "Hsp26", "Hsp27",
  "Hsp60", "Hsp68", "Hsp70", "Hsp83", "GSTe", "GSTd", "GCLC", "GCLM",
  
  # DNA damage and repair
  "p53", "rad50", "rad51", "lig3", "mus201", "mus209", "mus308", "mei-41", "lok", "mnk",
  
  # Apoptosis and cell death
  "rpr", "hid", "grim", "debcl", "buffy", "dark", "dronc", "drICE", "Diap1", "Diap2",
  
  # Autophagy
  "Atg1", "Atg5", "Atg6", "Atg7", "Atg8a", "Atg8b", "Atg9", "Atg12", "Atg13", "Atg18",
  
  # Others
  "pink1", "parkin", "BNIP3", "DHR38", "EcR"
)

# Function to analyze specific gene categories
analyze_gene_category <- function(res_df, gene_list, category_name) {
  cat("\n\n======================================================\n")
  cat("Analysis of", category_name, "Genes\n")
  cat("======================================================\n")
  
  # Find genes in the category that are in results
  category_genes_in_results <- which(res_df$symbol %in% gene_list)
  
  if (length(category_genes_in_results) == 0) {
    cat("No", category_name, "genes found in the results.\n")
    return(NULL)
  }
  
  cat("Found", length(category_genes_in_results), category_name, "genes in the results.\n\n")
  
  # Extract results for these genes
  category_results <- res_df[category_genes_in_results, ]
  
  # Count by regulation direction
  up_t5 <- sum(category_results$regulation == "Upregulated in T5", na.rm = TRUE)
  up_t3 <- sum(category_results$regulation == "Upregulated in T3", na.rm = TRUE)
  not_sig <- sum(category_results$regulation == "Not significant", na.rm = TRUE)
  
  cat("Direction of regulation:\n")
  cat("- Upregulated in T5:", up_t5, "genes\n")
  cat("- Upregulated in T3:", up_t3, "genes\n")
  cat("- Not significantly changed:", not_sig, "genes\n\n")
  
  # Filter for significant genes in this category
  sig_category_genes <- category_results[category_results$padj < 0.05, ]
  
  if (nrow(sig_category_genes) > 0) {
    cat("Significantly differentially expressed", category_name, "genes:\n")
    sig_df <- data.frame(
      gene_id = rownames(sig_category_genes),
      symbol = sig_category_genes$symbol,
      log2FoldChange = sig_category_genes$log2FoldChange,
      padj = sig_category_genes$padj,
      regulation = sig_category_genes$regulation,
      stringsAsFactors = FALSE
    )
    print(sig_df[order(sig_df$padj), ])
    
    # Create a specialized plot for this category
    category_plot <- ggplot(category_results, aes(x = log2FoldChange, y = -log10(padj))) +
      # Add background of all genes in grey
      geom_point(data = res_df, color = "grey90", size = 1, alpha = 0.3) +
      # Add category genes with regulation colors
      geom_point(aes(color = regulation, size = padj < 0.05)) +
      scale_color_manual(values = c("Upregulated in T5" = "red", 
                                    "Not significant" = "grey", 
                                    "Upregulated in T3" = "blue")) +
      scale_size_manual(values = c(2, 3), guide = "none") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = c(-0.5, 0, 0.5), linetype = "dashed", color = c("grey70", "grey50", "grey70")) +
      labs(title = paste0(category_name, " Genes - FA Timepoint 5 vs 3"),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    # Add gene labels
    category_plot <- category_plot +
      geom_text_repel(
        data = sig_category_genes,
        aes(label = symbol, color = regulation),
        max.overlaps = 30,
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        force = 3
      )
    
    print(category_plot)
    
    # Return results
    return(list(
      all_genes = category_results,
      significant_genes = sig_category_genes
    ))
  } else {
    cat("No significantly differentially expressed", category_name, "genes found.\n")
    return(list(all_genes = category_results))
  }
}

# Analyze metabolism genes
metabolism_analysis <- analyze_gene_category(res_with_direction, metabolism_genes, "Metabolism")

# Analyze oxidative stress genes
oxidative_analysis <- analyze_gene_category(res_with_direction, oxidative_stress_genes, "Oxidative Stress")

# PART 5: Normalized Expression Heatmaps by Direction
#----------------------------------------------------------------------------

# Get normalized counts from the DESeq object
vsd <- vst(tp_results$dds, blind = FALSE)

# Function to create a heatmap for a set of genes
create_gene_heatmap <- function(gene_ids, gene_symbols, title) {
  if (length(gene_ids) == 0 || length(gene_symbols) == 0) {
    cat("No genes available for heatmap:", title, "\n")
    return(NULL)
  }
  
  # Subset gene IDs that are in the data
  gene_ids_in_data <- gene_ids[gene_ids %in% rownames(vsd)]
  
  if (length(gene_ids_in_data) == 0) {
    cat("None of the provided gene IDs are found in the data for heatmap:", title, "\n")
    return(NULL)
  }
  
  # Limit to the top 50 genes if there are more
  if (length(gene_ids_in_data) > 50) {
    # Sort by adjusted p-value
    ordered_idx <- order(res_with_direction[gene_ids_in_data, "padj"])
    gene_ids_in_data <- gene_ids_in_data[ordered_idx[1:50]]
    gene_symbols <- gene_symbols[ordered_idx[1:50]]
  }
  
  # Extract counts for these genes
  gene_counts <- assay(vsd)[gene_ids_in_data, ]
  
  # Use gene symbols as row names if available
  if (length(gene_symbols) >= nrow(gene_counts)) {
    rownames(gene_counts) <- gene_symbols[1:nrow(gene_counts)]
  }
  
  # Create annotation for heatmap
  anno_col <- data.frame(
    Timepoint = factor(tp_results$metadata$time),
    row.names = colnames(gene_counts)
  )
  
  # Generate heatmap
  pheatmap(
    gene_counts,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = anno_col,
    scale = "row",
    main = title,
    fontsize_row = 8,
    fontsize_col = 8,
    cutree_rows = 2,
    cutree_cols = 2
  )
}

# Create heatmaps for up and down regulated genes
create_gene_heatmap(
  direction_analysis$up_t5_genes, 
  direction_analysis$up_t5_symbols,
  "Genes Upregulated in T5 vs T3"
)

create_gene_heatmap(
  direction_analysis$up_t3_genes,
  direction_analysis$up_t3_symbols,
  "Genes Upregulated in T3 vs T5"
)

# PART 6: Summary Table of Top Genes by Direction
#----------------------------------------------------------------------------

# Function to create a detailed table of top genes by direction
create_detailed_gene_table <- function(res_df, direction, n = 25) {
  # Filter for the specific direction
  genes <- res_df[res_df$regulation == direction, ]
  
  # Sort by fold change magnitude if there are genes in this direction
  if (nrow(genes) > 0) {
    if (direction == "Upregulated in T5") {
      genes <- genes[order(genes$log2FoldChange, decreasing = TRUE), ]
    } else if (direction == "Upregulated in T3") {
      genes <- genes[order(genes$log2FoldChange), ] # Already negative values
    }
    
    # Take top N genes
    top_genes <- head(genes, n)
    
    # Create table
    gene_table <- data.frame(
      Rank = 1:nrow(top_genes),
      FlyBase_ID = rownames(top_genes),
      Symbol = top_genes$symbol,
      Log2FC = round(top_genes$log2FoldChange, 3),
      Padj = signif(top_genes$padj, 3),
      BaseMean = round(top_genes$baseMean, 1)
    )
    
    # Add direction to title
    cat("\n\nTop", n, "genes", direction, ":\n")
    print(gene_table)
    
    return(gene_table)
  } else {
    cat("\n\nNo genes found for direction:", direction, "\n")
    return(NULL)
  }
}
# PART 7: Gene Expression Patterns Analysis 
#----------------------------------------------------------------------------

# Extract z-scores for visualization of expression patterns
extract_z_scores <- function() {
  # Get normalized counts
  norm_counts <- assay(vsd)
  
  # Combine upregulated genes from both directions
  genes_of_interest <- unique(c(
    direction_analysis$up_t5_genes[1:min(50, length(direction_analysis$up_t5_genes))],
    direction_analysis$up_t3_genes[1:min(50, length(direction_analysis$up_t3_genes))]
  ))
  
  # Subset counts for these genes
  if (length(genes_of_interest) > 0) {
    # Keep only genes that are in the expression data
    genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(norm_counts)]
    if (length(genes_of_interest) == 0) {
      return(NULL)
    }
    
    gene_counts <- norm_counts[genes_of_interest, ]
    
    # Calculate z-scores for each gene
    gene_means <- rowMeans(gene_counts)
    gene_sds <- apply(gene_counts, 1, sd)
    
    # Avoid division by zero
    gene_sds[gene_sds == 0] <- 1
    
    z_scores <- sweep(gene_counts, 1, gene_means, "-")
    z_scores <- sweep(z_scores, 1, gene_sds, "/")
    
    # Get regulation direction for each gene
    gene_direction <- res_with_direction[genes_of_interest, "regulation"]
    
    # Create a data frame for plotting
    z_df <- data.frame(
      Gene = rep(rownames(z_scores), each = ncol(z_scores)),
      Sample = rep(colnames(z_scores), times = nrow(z_scores)),
      Z_score = as.vector(z_scores),
      Direction = rep(gene_direction, each = ncol(z_scores))
    )
    
    # Add timepoint information
    sample_timepoint <- tp_results$metadata$time
    names(sample_timepoint) <- tp_results$metadata$sample_id
    z_df$Timepoint <- factor(sample_timepoint[z_df$Sample])
    
    return(z_df)
  } else {
    return(NULL)
  }
}

# Extract z-scores
z_scores_data <- extract_z_scores()

if (!is.null(z_scores_data)) {
  # Create boxplot of z-scores by timepoint and regulation direction
  z_boxplot <- ggplot(z_scores_data, aes(x = Timepoint, y = Z_score, fill = Direction)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~Direction, scales = "free_y") +
    scale_fill_manual(values = c("Upregulated in T5" = "red", 
                                 "Not significant" = "grey", 
                                 "Upregulated in T3" = "blue")) +
    labs(title = "Expression Patterns of Differentially Expressed Genes",
         subtitle = "Z-scores by timepoint and regulation direction",
         x = "Timepoint",
         y = "Z-score (normalized expression)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  print(z_boxplot)
}

# PART 8: Pathway Impact Analysis
#----------------------------------------------------------------------------

# Function to analyze pathway impact using GSEA-based approach
pathway_impact_analysis <- function() {
  # First create a named vector of log2 fold changes
  gene_list <- res_with_direction$log2FoldChange
  names(gene_list) <- rownames(res_with_direction)
  
  # Sort the list in decreasing order
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert IDs from FlyBase to ENTREZ
  id_mappings <- bitr(names(gene_list), 
                      fromType = "FLYBASE", 
                      toType = c("ENTREZID", "SYMBOL"),
                      OrgDb = org.Dm.eg.db)
  
  # Create a new gene list with ENTREZ ids
  entrez_gene_list <- gene_list[id_mappings$FLYBASE]
  names(entrez_gene_list) <- id_mappings$ENTREZID
  
  cat("\nRunning Gene Set Enrichment Analysis (GSEA)...\n")
  
  # Perform GSEA with GO
  gsea_go_bp <- gseGO(geneList = entrez_gene_list,
                      OrgDb = org.Dm.eg.db,
                      ont = "BP",
                      keyType = "ENTREZID",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")
  
  # Perform GSEA with KEGG
  gsea_kegg <- gseKEGG(geneList = entrez_gene_list,
                       organism = 'dme',
                       keyType = 'ncbi-geneid',
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
  
  # Print results
  cat("\nGSEA results for GO Biological Process:\n")
  if (!is.null(gsea_go_bp) && nrow(gsea_go_bp@result) > 0) {
    # Extract top enriched GO terms
    top_go <- head(gsea_go_bp@result, 15)
    print(top_go[, c("ID", "Description", "enrichmentScore", "NES", "pvalue", "p.adjust")])
    
    # Visualize GSEA results
    dotplot(gsea_go_bp, showCategory = 15, title = "GO BP GSEA Results")
    
    # Create enrichment plot
    tryCatch({
      gseaplot2(gsea_go_bp, geneSetID = 1:min(3, nrow(gsea_go_bp@result)), 
                title = "GSEA Enrichment Plot for Top GO BP Pathways")
    }, error = function(e) {
      cat("Could not generate GSEA plot:", e$message, "\n")
    })
  } else {
    cat("No significant GO Biological Process terms found in GSEA.\n")
  }
  
  cat("\nGSEA results for KEGG pathways:\n")
  if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    # Extract top enriched KEGG pathways
    top_kegg <- head(gsea_kegg@result, 15)
    print(top_kegg[, c("ID", "Description", "enrichmentScore", "NES", "pvalue", "p.adjust")])
    
    # Visualize GSEA results
    dotplot(gsea_kegg, showCategory = 15, title = "KEGG GSEA Results")
    
    # Create enrichment plot
    tryCatch({
      gseaplot2(gsea_kegg, geneSetID = 1:min(3, nrow(gsea_kegg@result)), 
                title = "GSEA Enrichment Plot for Top KEGG Pathways")
    }, error = function(e) {
      cat("Could not generate GSEA plot:", e$message, "\n")
    })
  } else {
    cat("No significant KEGG pathways found in GSEA.\n")
  }
  
  return(list(gsea_go_bp = gsea_go_bp, gsea_kegg = gsea_kegg))
}

# Run pathway impact analysis
gsea_results <- pathway_impact_analysis()

# PART 9: Focused Analysis of Biological Processes of Interest
#----------------------------------------------------------------------------

# Define key biological processes that might be relevant to your study
# These can be adjusted based on your research question
key_processes <- list(
  "Respiration" = c(
    "mt:ND1", "mt:ND2", "mt:ND3", "mt:ND4", "mt:ND5", "mt:ND6", "mt:Cyt-b", "mt:CoI", "mt:CoII", "mt:CoIII",
    "mt:ATPase6", "mt:ATPase8", "ND75", "ND42", "ND51", "ND23", "ND24",
    "ATPsynα", "ATPsynβ", "ATPsynγ", "ATPsynδ", "ATPsynε", "Cyc1", "Cyt-c-p", "Cyt-c-d"
  ),
  
  "Hypoxia_Response" = c(
    "sima", "tango", "Hph", "HIF1A", "Alas", "Pdp1", "dNRF2", "cbt"
  ),
  
  "Stress_Response" = c(
    "Sod1", "Sod2", "Cat", "Trx-2", "TrxR-1", "TrxR-2", "Hsp22", "Hsp23", "Hsp26", "Hsp27",
    "Hsp60", "Hsp68", "Hsp70", "Hsp83", "Atf3", "Gadd45", "p53", "JNK"
  ))

# Function to analyze specific biological processes
analyze_biological_process <- function(res_df, process_list, process_name) {
  cat("\n\n======================================================\n")
  cat("Analysis of", process_name, "Process\n")
  cat("======================================================\n")
  
  # Get genes in this process category
  genes_in_process <- process_list[[process_name]]
  
  # Find process genes in the results
  process_genes_in_results <- res_df$symbol %in% genes_in_process
  
  if (sum(process_genes_in_results) == 0) {
    cat("No", process_name, "genes found in the results.\n")
    return(NULL)
  }
  
  # Extract results for these genes
  process_results <- res_df[process_genes_in_results, ]
  
  cat("Found", nrow(process_results), process_name, "genes in the results.\n\n")
  
  # Count by regulation direction
  up_t5 <- sum(process_results$regulation == "Upregulated in T5", na.rm = TRUE)
  up_t3 <- sum(process_results$regulation == "Upregulated in T3", na.rm = TRUE)
  not_sig <- sum(process_results$regulation == "Not significant", na.rm = TRUE)
  
  cat("Direction of regulation:\n")
  cat("- Upregulated in T5:", up_t5, "genes\n")
  cat("- Upregulated in T3:", up_t3, "genes\n")
  cat("- Not significantly changed:", not_sig, "genes\n\n")
  
  # Filter for significant genes in this process
  sig_process_genes <- process_results[process_results$padj < 0.05, ]
  
  if (nrow(sig_process_genes) > 0) {
    cat("Significantly differentially expressed", process_name, "genes:\n")
    sig_df <- data.frame(
      gene_id = rownames(sig_process_genes),
      symbol = sig_process_genes$symbol,
      log2FoldChange = sig_process_genes$log2FoldChange,
      padj = sig_process_genes$padj,
      regulation = sig_process_genes$regulation,
      stringsAsFactors = FALSE
    )
    print(sig_df[order(sig_df$padj), ])
    
    # Create a specialized plot for this process
    process_plot <- ggplot(process_results, aes(x = log2FoldChange, y = -log10(padj))) +
      # Add background of all genes in grey
      geom_point(data = res_df, color = "grey90", size = 1, alpha = 0.3) +
      # Add process genes with regulation colors
      geom_point(aes(color = regulation, size = padj < 0.05)) +
      scale_color_manual(values = c("Upregulated in T5" = "red", 
                                    "Not significant" = "grey", 
                                    "Upregulated in T3" = "blue")) +
      scale_size_manual(values = c(2, 3), guide = "none") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = c(-0.5, 0, 0.5), linetype = "dashed", color = c("grey70", "grey50", "grey70")) +
      labs(title = paste0(process_name, " Genes - FA Timepoint 5 vs 3"),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    # Add gene labels
    process_plot <- process_plot +
      geom_text_repel(
        data = sig_process_genes,
        aes(label = symbol, color = regulation),
        max.overlaps = 30,
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        force = 3
      )
    
    print(process_plot)
    
    return(list(
      all_genes = process_results,
      significant_genes = sig_process_genes
    ))
  } else {
    cat("No significantly differentially expressed", process_name, "genes found.\n")
    return(list(all_genes = process_results))
  }
}

# Analyze each biological process
process_results <- list()
for (process in names(key_processes)) {
  process_results[[process]] <- analyze_biological_process(
    res_with_direction, key_processes, process
  )
}

# PART 10: Compare T5 Upregulated vs T3 Upregulated Genes
#----------------------------------------------------------------------------

# Create a summary comparison of T5 vs T3 upregulated genes
gene_comparison_summary <- function() {
  cat("\n\n======================================================\n")
  cat("Comparison Summary: T5 Upregulated vs T3 Upregulated Genes\n")
  cat("======================================================\n")
  
  # Get the sets of genes
  up_t5_genes <- direction_analysis$up_t5_genes
  up_t3_genes <- direction_analysis$up_t3_genes
  
  # Basic summary statistics
  cat("Number of genes upregulated in T5:", length(up_t5_genes), "\n")
  cat("Number of genes upregulated in T3:", length(up_t3_genes), "\n\n")
  
  # Get the mapped gene symbols for each set
  up_t5_symbols <- direction_analysis$up_t5_symbols
  up_t3_symbols <- direction_analysis$up_t3_symbols
  
  # Print top genes in each direction with fold changes
  cat("Top 10 genes upregulated in T5:\n")
  if (length(up_t5_genes) > 0) {
    top_t5 <- data.frame(
      FlyBase_ID = up_t5_genes[1:min(10, length(up_t5_genes))],
      Symbol = up_t5_symbols[1:min(10, length(up_t5_symbols))],
      Log2FC = res_with_direction[up_t5_genes[1:min(10, length(up_t5_genes))], "log2FoldChange"],
      Padj = res_with_direction[up_t5_genes[1:min(10, length(up_t5_genes))], "padj"]
    )
    print(top_t5)
  }
  
  cat("\nTop 10 genes upregulated in T3:\n")
  if (length(up_t3_genes) > 0) {
    top_t3 <- data.frame(
      FlyBase_ID = up_t3_genes[1:min(10, length(up_t3_genes))],
      Symbol = up_t3_symbols[1:min(10, length(up_t3_symbols))],
      Log2FC = res_with_direction[up_t3_genes[1:min(10, length(up_t3_genes))], "log2FoldChange"],
      Padj = res_with_direction[up_t3_genes[1:min(10, length(up_t3_genes))], "padj"]
    )
    print(top_t3)
  }
  
  # Compare gene lists from T5 and T3 for distinct biological functions
  if (!is.null(up_t5_enrichment$go_bp) && !is.null(up_t3_enrichment$go_bp)) {
    cat("\nDistinct biological processes in T5 upregulated genes:\n")
    if (nrow(up_t5_enrichment$go_bp@result) > 0) {
      print(head(up_t5_enrichment$go_bp@result[, c("ID", "Description", "GeneRatio", "p.adjust")], 5))
    } else {
      cat("No significant GO terms found.\n")
    }
    
    cat("\nDistinct biological processes in T3 upregulated genes:\n")
    if (nrow(up_t3_enrichment$go_bp@result) > 0) {
      print(head(up_t3_enrichment$go_bp@result[, c("ID", "Description", "GeneRatio", "p.adjust")], 5))
    } else {
      cat("No significant GO terms found.\n")
    }
  }
  
  # Create a comparison barplot
  comparison_data <- data.frame(
    Category = c(
      rep("Metabolism", 2),
      rep("TCA Cycle", 2),
      rep("Oxidative Stress", 2),
      rep("Signaling", 2)
    ),
    Direction = rep(c("Upregulated in T5", "Upregulated in T3"), 4),
    Count = c(
      # Metabolism counts
      sum(direction_analysis$up_t5_symbols %in% metabolism_genes),
      sum(direction_analysis$up_t3_symbols %in% metabolism_genes),
      
      # TCA cycle counts (subset of metabolism genes)
      sum(direction_analysis$up_t5_symbols %in% c("Idh", "Mdh1", "Mdh2", "SdhA", "SdhB", "Acon", "Fum")),
      sum(direction_analysis$up_t3_symbols %in% c("Idh", "Mdh1", "Mdh2", "SdhA", "SdhB", "Acon", "Fum")),
      
      # Oxidative stress counts
      sum(direction_analysis$up_t5_symbols %in% oxidative_stress_genes),
      sum(direction_analysis$up_t3_symbols %in% oxidative_stress_genes),
      
      # Signaling counts
      sum(direction_analysis$up_t5_symbols %in% key_processes$Insulin_Signaling),
      sum(direction_analysis$up_t3_symbols %in% key_processes$Insulin_Signaling)
    )
  )
  
  # Create comparison barplot
  comparison_plot <- ggplot(comparison_data, aes(x = Category, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Upregulated in T5" = "red", "Upregulated in T3" = "blue")) +
    labs(title = "Comparison of Gene Categories Between T5 and T3 Upregulation",
         x = "Functional Category",
         y = "Number of Genes") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(comparison_plot)
}

# Run the comparison summary
gene_comparison_summary()

