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
 keep<-rowSums(counts(dds) >= 10) >= min(3, ncol(dds))
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
      
     
