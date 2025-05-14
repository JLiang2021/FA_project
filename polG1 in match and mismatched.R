library(ggplot2)
library(org.Dm.eg.db)
library(dplyr)

# Function to find a gene by name or partial name in results
find_gene <- function(gene_pattern, results_obj, full_results = TRUE) {
  pattern <- tolower(gene_pattern)
  if ("gene_name" %in% colnames(results_obj)) {
    matches <- grep(pattern, tolower(results_obj$gene_name), value = FALSE)
    if (length(matches) > 0) {
      if (full_results) {
        return(results_obj[matches, ])
      } else {
        return(rownames(results_obj)[matches])
      }
    }
  }
  

  matches <- grep(pattern, tolower(rownames(results_obj)), value = FALSE)
  if (length(matches) > 0) {
    if (full_results) {
      return(results_obj[matches, ])
    } else {
      return(rownames(results_obj)[matches])
    }
  }
  
  if (gene_pattern %in% rownames(results_obj)) {
    if (full_results) {
      return(results_obj[gene_pattern, , drop = FALSE])
    } else {
      return(gene_pattern)
    }
  }
  
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    tryCatch({
      fb_id <- AnnotationDbi::select(org.Dm.eg.db, 
                                     keys = gene_pattern,
                                     columns = "FLYBASE", 
                                     keytype = "SYMBOL")
      if (nrow(fb_id) > 0) {
        for (id in fb_id$FLYBASE) {
          if (id %in% rownames(results_obj)) {
            if (full_results) {
              return(results_obj[id, , drop = FALSE])
            } else {
              return(id)
            }
          }
        }
      }
    }, error = function(e) {
    })
  }
  cat("Gene", gene_pattern, "not found in results.\n")
  return(NULL)
}


# Try different possible names/IDs for polG1
polg_results <- find_gene("polG", res_matched)
if (is.null(polg_results)) {
  polg_results <- find_gene("polG", res_matched)
}
if (is.null(polg_results)) {
  polg_results <- find_gene("tam", res_matched)
}
if (is.null(polg_results)) {
  polg_results <- find_gene("tamas", res_matched)
}
if (is.null(polg_results)) {
  polg_results <- find_gene("FBgn0004406", res_matched)
}

# Print the results if found
if (!is.null(polG_results)) {
  cat("\nDifferential expression results for polG (tamas):\n")
  print(polg_results)
  
  # Extract key statistics
  log2FC <- polg_results$log2FoldChange
  padj <- polg_results$padj
  
  cat("\nLog2 fold change:", log2FC, "\n")
  cat("Adjusted p-value:", padj, "\n")
  cat("Significant:", ifelse(padj < 0.05, "YES", "NO"), "\n")
} else {
  cat("\npolG1 gene not found in the results. Check if it's filtered out or using a different identifier.\n")
}

# Get normalized counts for visualization
if (!is.null(polg1_results) && exists("dds") && exists("metadata_filtered")) {
  # Get the FlyBase ID
  polg1_id <- rownames(polg1_results)[1]
  
  # Extract normalized counts
  normalized_counts <- counts(dds, normalized=TRUE)[polg1_id, , drop=FALSE]
  
  # Create a data frame for plotting
  count_data <- data.frame(
    sample_id = colnames(normalized_counts),
    count = as.numeric(normalized_counts),
    stringsAsFactors = FALSE
  )
  
  # Merge with metadata to get matched status and other info
  plot_data <- merge(count_data, metadata_filtered, by = "sample_id")
  
  # Create a boxplot of polG1 expression by matched status
  p1 <- ggplot(plot_data, aes(x = matched, y = count, fill = matched)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    labs(title = "polG1 Expression: Matched vs Unmatched",
         x = "Matched Status",
         y = "Normalized Expression Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Print the plot
  print(p1)
  
  # Create a plot showing expression across time and groups
  p2 <- ggplot(plot_data, aes(x = factor(time), y = count, color = group, shape = matched)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    stat_summary(aes(group = interaction(group, matched)), 
                 fun = mean, geom = "line", position = position_dodge(width = 0.3)) +
    labs(title = "polG1 Expression Over Time by Group and Matched Status",
         x = "Time Point",
         y = "Normalized Expression Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Print the plot
  print(p2)
  
  # Check if polG1 appears in any significant GO terms or KEGG pathways
  if (exists("go_bp") && exists("entrez_ids")) {
    # First, find the Entrez ID for polG1
    polg1_entrez <- entrez_ids$ENTREZID[entrez_ids$FLYBASE == polg1_id]
    
    if (length(polg1_entrez) > 0) {
      # Check if this gene appears in any significant GO terms
      cat("\nChecking if polG1 appears in significant GO Biological Process terms...\n")
      
      found_in_terms <- FALSE
      
      if (nrow(go_bp@result) > 0) {
        for (i in 1:nrow(go_bp@result)) {
          term_id <- go_bp@result$ID[i]
          term_name <- go_bp@result$Description[i]
          genes <- go_bp@geneSets[[term_id]]
          
          if (polg1_entrez %in% genes) {
            cat("Found in GO term:", term_name, "(", term_id, ")\n")
            cat("  - Adjusted p-value:", go_bp@result$p.adjust[i], "\n")
            found_in_terms <- TRUE
          }
        }
      }
      
      if (!found_in_terms) {
        cat("polG1 not found in any significant GO Biological Process terms.\n")
      }
      
      # Check KEGG pathways
      if (exists("kegg_result")) {
        cat("\nChecking if polG1 appears in significant KEGG pathways...\n")
        
        found_in_pathways <- FALSE
        
        if (nrow(kegg_result@result) > 0) {
          for (i in 1:nrow(kegg_result@result)) {
            pathway_id <- kegg_result@result$ID[i]
            pathway_name <- kegg_result@result$Description[i]
            genes <- kegg_result@geneSets[[pathway_id]]
            
            if (polg1_entrez %in% genes) {
              cat("Found in KEGG pathway:", pathway_name, "(", pathway_id, ")\n")
              cat("  - Adjusted p-value:", kegg_result@result$p.adjust[i], "\n")
              found_in_pathways <- TRUE
            }
          }
        }
        
        if (!found_in_pathways) {
          cat("polG1 not found in any significant KEGG pathways.\n")
        }
      }
    } else {
      cat("\nCouldn't find Entrez ID for polG1, so cannot check pathway involvement.\n")
    }
  }
} else {
  cat("\nCannot generate expression plots as required data is missing.\n")
}

#mitochondria genes 
# Analysis of respiration genes comparing mismatched to matched samples

# Gene Expression Fold Change Plot with Gene Symbols
# X-axis: FF/AA, Y-axis: FA/AA

# Load required libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggrepel)
library(org.Dm.eg.db)
library(AnnotationDbi)

# Improved function to map FlyBase IDs to gene symbols
map_ids_to_symbols <- function(gene_ids) {
  # Create a vector to store the mapped symbols
  symbols <- rep(NA, length(gene_ids))
  names(symbols) <- gene_ids
  
  # First try direct mapping using org.Dm.eg.db
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    # Check if gene_ids look like FlyBase IDs
    flybase_pattern <- "^FB"
    if (length(gene_ids) > 0 && any(grepl(flybase_pattern, gene_ids))) {
      # Filter IDs that match the FlyBase pattern
      fb_ids <- gene_ids[grepl(flybase_pattern, gene_ids)]
      
      # Map these IDs to symbols
      mapped_symbols <- AnnotationDbi::mapIds(org.Dm.eg.db, 
                                              keys = fb_ids,
                                              column = "SYMBOL", 
                                              keytype = "FLYBASE",
                                              multiVals = "first")
      
      # Update the symbols vector
      symbols[names(mapped_symbols)] <- mapped_symbols
    }
  }
  
  # For any remaining unmapped IDs, extract a symbolic name from the ID itself
  na_indices <- is.na(symbols)
  if (any(na_indices)) {
    unmapped_ids <- names(symbols)[na_indices]
    for (i in seq_along(unmapped_ids)) {
      id <- unmapped_ids[i]
      # Try to extract a meaningful part from the ID
      if (grepl("^FB", id)) {
        # For FlyBase IDs, use the numeric part
        symbols[id] <- sub("^FBgn0*", "g", id)
      } else if (grepl("^mt:", id)) {
        # For mitochondrial genes, use the name part
        symbols[id] <- sub("^mt:", "mt", id)
      } else {
        # For other IDs, use as is
        symbols[id] <- id
      }
    }
  }
  
  return(symbols)
}

#
# Get fold changes for the correct comparisons
#

get_correct_fold_changes <- function() {
  # Check if we have the necessary DESeq2 object
  if (!exists("dds")) {
    cat("DESeq2 object 'dds' not found. Cannot calculate fold changes.\n")
    return(NULL)
  }
  
  # Check if group variable exists in colData
  if (!"group" %in% colnames(colData(dds))) {
    cat("Group variable not found in DESeq2 object.\n")
    return(NULL)
  }
  
  # Extract groups
  groups <- unique(as.character(colData(dds)$group))
  
  # Check if we have the needed groups
  if (!all(c("FA", "AA", "FF") %in% groups)) {
    cat("Not all required groups (FA, AA, FF) found in data.\n")
    cat("Available groups:", paste(groups, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Create result objects for the needed comparisons
  try({
    res_FF_AA <- results(dds, contrast = c("group", "FF", "AA"))
    res_FA_AA <- results(dds, contrast = c("group", "FA", "AA"))
    
    # Extract log2 fold changes
    fc_FF_AA <- res_FF_AA$log2FoldChange
    fc_FA_AA <- res_FA_AA$log2FoldChange
    
    # Create data frame
    fc_data <- data.frame(
      gene_id = rownames(res_FF_AA),
      FC_FF_AA = fc_FF_AA,
      FC_FA_AA = fc_FA_AA,
      padj_FF_AA = res_FF_AA$padj,
      padj_FA_AA = res_FA_AA$padj,
      stringsAsFactors = FALSE
    )
    
    # Add significance columns
    fc_data$sig_FF_AA <- fc_data$padj_FF_AA < 0.05
    fc_data$sig_FA_AA <- fc_data$padj_FA_AA < 0.05
    fc_data$sig_both <- fc_data$sig_FF_AA & fc_data$sig_FA_AA
    
    # Remove rows with NA values
    fc_data <- fc_data[complete.cases(fc_data[, c("FC_FF_AA", "FC_FA_AA")]), ]
    
    return(fc_data)
  }, silent = TRUE)
}

# Get fold changes
fc_data <- get_correct_fold_changes()

  gene_symbols <- map_ids_to_symbols(fc_data$gene_id)
  fc_data$gene_symbol <- gene_symbols
  
  # Function to label a subset of genes in the plot
  label_genes <- function(df, n_top = 10, method = "both") {
    if (method == "both") {
      # Label genes significant in both comparisons
      sig_both <- df[df$sig_both, ]
      # Sort by the sum of absolute fold changes
      sig_both$total_fc <- abs(sig_both$FC_FF_AA) + abs(sig_both$FC_FA_AA)
      sig_both <- sig_both[order(sig_both$total_fc, decreasing = TRUE), ]
      # Take top n
      return(head(sig_both, n_top)$gene_id)
    } else if (method == "extreme") {
      # Label genes with extreme fold changes in either direction
      df$extreme_score <- abs(df$FC_FF_AA) + abs(df$FC_FA_AA)
      df <- df[order(df$extreme_score, decreasing = TRUE), ]
      return(head(df, n_top)$gene_id)
    } else if (method == "diagonal") {
      # Label genes furthest from the diagonal (showing different patterns)
      df$diag_dist <- abs(df$FC_FF_AA - df$FC_FA_AA)
      df <- df[order(df$diag_dist, decreasing = TRUE), ]
      return(head(df, n_top)$gene_id)
    }
    return(NULL)
  }
  
  # Define colors for significance
  if ("sig_both" %in% colnames(fc_data)) {
    fc_data$color_group <- "Not significant"
    fc_data$color_group[fc_data$sig_FF_AA & !fc_data$sig_FA_AA] <- "Significant in FF/AA only"
    fc_data$color_group[!fc_data$sig_FF_AA & fc_data$sig_FA_AA] <- "Significant in FA/AA only"
    fc_data$color_group[fc_data$sig_both] <- "Significant in both"
    
    # Set colors
    color_mapping <- c(
      "Not significant" = "gray", 
      "Significant in FF/AA only" = "blue", 
      "Significant in FA/AA only" = "red",
      "Significant in both" = "purple"
    )
    
    # Select genes to label
    genes_to_label <- label_genes(fc_data, 15, "both")
    if (length(genes_to_label) < 10) {
      # If not enough genes significant in both, label genes with extreme fold changes
      genes_to_label <- c(genes_to_label, setdiff(label_genes(fc_data, 15, "extreme"), genes_to_label))
      if (length(genes_to_label) > 15) genes_to_label <- genes_to_label[1:15]
    }
  } else {
    # Define colors based on quadrants
    fc_data$color_group <- "Q1: Up in both"
    fc_data$color_group[fc_data$FC_FF_AA < 0 & fc_data$FC_FA_AA >= 0] <- "Q2: Down in FF/AA, Up in FA/AA"
    fc_data$color_group[fc_data$FC_FF_AA < 0 & fc_data$FC_FA_AA < 0] <- "Q3: Down in both"
    fc_data$color_group[fc_data$FC_FF_AA >= 0 & fc_data$FC_FA_AA < 0] <- "Q4: Up in FF/AA, Down in FA/AA"
    
    # Set colors
    color_mapping <- c(
      "Q1: Up in both" = "purple", 
      "Q2: Down in FF/AA, Up in FA/AA" = "blue", 
      "Q3: Down in both" = "green",
      "Q4: Up in FF/AA, Down in FA/AA" = "red"
    )
    
    # Select genes to label
    genes_to_label <- label_genes(fc_data, 15, "extreme")
  }
  
  # Create a scatter plot
  p1 <- ggplot(fc_data, aes(x = FC_FF_AA, y = FC_FA_AA, color = color_group)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    scale_color_manual(values = color_mapping) +
    labs(title = "Gene Expression Fold Changes",
         subtitle = "Comparing FF/AA (x-axis) vs FA/AA (y-axis)",
         x = "Log2 Fold Change (FF/AA)",
         y = "Log2 Fold Change (FA/AA)",
         color = "Significance") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 10)
    )
  
  # Add labels for selected genes
  label_df <- fc_data[fc_data$gene_id %in% genes_to_label, ]
  
  p1 <- p1 + 
    geom_text_repel(data = label_df, 
                    aes(label = gene_symbol),
                    box.padding = 0.5,
                    point.padding = 0.3,
                    max.overlaps = 30,
                    segment.color = "gray50",
                    min.segment.length = 0,
                    force = 5)
  
  print(p1)
  
  #
  # Define a comprehensive list of respiration genes with proper symbols
  #
  
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
  
  # Get mappings for respiration genes
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    gene_mappings <- AnnotationDbi::select(org.Dm.eg.db, 
                                           keys = respiration_genes,
                                           columns = c("FLYBASE", "SYMBOL"),
                                           keytype = "SYMBOL")
  } else {
    gene_mappings <- data.frame(
      SYMBOL = respiration_genes,
      FLYBASE = NA,
      stringsAsFactors = FALSE
    )
  }
  
  # Find respiration genes in our fold change data
  # Use both gene_id and gene_symbol for matching
  resp_gene_ids <- fc_data$gene_id[fc_data$gene_id %in% gene_mappings$FLYBASE | 
                                     fc_data$gene_symbol %in% gene_mappings$SYMBOL |
                                     fc_data$gene_symbol %in% respiration_genes |
                                     fc_data$gene_id %in% respiration_genes]
  
  if (length(resp_gene_ids) > 0) {
    cat("Found", length(resp_gene_ids), "respiration genes in the fold change data.\n")
    
    # Create dataset for respiration genes
    resp_fc_data <- fc_data[fc_data$gene_id %in% resp_gene_ids, ]
    
    # Highlight downregulated genes in mismatched (FA) in red
    # This means genes with negative FA/AA fold change
    resp_fc_data$downreg <- resp_fc_data$FC_FA_AA < 0
    
    # Create a scatter plot for respiration genes
    p2 <- ggplot(resp_fc_data, aes(x = FC_FF_AA, y = FC_FA_AA, color = downreg)) +
      geom_point(alpha = 0.7, size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                         labels = c("TRUE" = "Downregulated in FA", "FALSE" = "Upregulated in FA")) +
      labs(title = "Respiration Genes Fold Changes",
           subtitle = "Comparing FF/AA (x-axis) vs FA/AA (y-axis)",
           x = "Log2 Fold Change (FF/AA)",
           y = "Log2 Fold Change (FA/AA)",
           color = "Regulation") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 10)
      )
    
    # Add labels for all respiration genes
    p2 <- p2 + 
      geom_text_repel(data = resp_fc_data, 
                      aes(label = gene_symbol),
                      box.padding = 0.5,
                      point.padding = 0.3,
                      max.overlaps = 50,
                      segment.color = "gray50",
                      min.segment.length = 0,
                      force = 3)
    
    print(p2)
  } else {
    cat("No respiration genes found in the fold change data.\n")
  }
} else {
  cat("Could not calculate fold changes. Please ensure you have the required data.\n")
}

# mitochondria biosynthesis 
# Analysis of mitochondrial biosynthesis gene expression over time
# Examining changes across the three timepoints (similar to polG1 analysis)

# Load required libraries
library(ggplot2)
library(dplyr)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# Define list of mitochondrial biosynthesis genes
# Including genes involved in mitochondrial biogenesis, translation, and assembly
mito_biosynthesis_genes <- c(
  # Mitochondrial DNA replication and maintenance
  "tam", "mtSSB", "Pol32", "DNApol-γ35", "MtTFB2", "MtTFB1", "TFAM", "polG", "polG2", "polGA", "polGB",
  
  # Mitochondrial transcription
  "mtTFB1", "mtTFB2", "TFAM", "mTerf3", "mTerf5", "POLRMT", "TFB1M", "TFB2M",
  
  # Mitochondrial translation and ribosomes
  "mRpL1", "mRpL2", "mRpL3", "mRpL4", "mRpL9", "mRpL10", "mRpL11", "mRpL12", "mRpL13", "mRpL14", 
  "mRpL15", "mRpL16", "mRpL17", "mRpL18", "mRpL19", "mRpL20", "mRpL21", "mRpL22", "mRpL23", "mRpL24",
  "mRpL27", "mRpL28", "mRpL30", "mRpL32", "mRpL33", "mRpL34", "mRpL35", "mRpL36", "mRpL37", "mRpL38",
  "mRpL39", "mRpL40", "mRpL41", "mRpL42", "mRpL43", "mRpL44", "mRpL45", "mRpL46", "mRpL47", "mRpL48",
  "mRpL49", "mRpL50", "mRpL51", "mRpL52", "mRpL53", "mRpL54", "mRpL55",
  "mRpS1", "mRpS2", "mRpS5", "mRpS6", "mRpS7", "mRpS8", "mRpS9", "mRpS10", "mRpS11", "mRpS12", 
  "mRpS14", "mRpS15", "mRpS16", "mRpS17", "mRpS18", "mRpS21", "mRpS22", "mRpS23", "mRpS24", "mRpS25",
  "mRpS26", "mRpS27", "mRpS28", "mRpS29", "mRpS30", "mRpS31", "mRpS33", "mRpS34", "mRpS35",
  
  # Mitochondrial import machinery
  "Tom20", "Tom22", "Tom40", "Tom70", "Tim8", "Tim9", "Tim10", "Tim13", "Tim17", "Tim22", "Tim23", "Tim44",
  "Tim50", "Mia40", "Erv1", "Sam50", "metaxin", "Hsp60", "Hsp10", "mtHsp70", "Mdj1", "Mge1",
  
  # Mitochondrial dynamics
  "Mfn", "Marf", "Drp1", "Fis1", "Opa1", "Pink1", "park", "Miro", "Milton",
  
  # Nuclear-encoded mitochondrial assembly factors
  "Coa3", "Coa6", "Coa7", "Cox10", "Cox11", "Cox14", "Cox15", "Cox16", "Cox17", "Cox18", "Cox19", "Cox20",
  "Surf1", "Sco1", "Sco2", "Oxa1", "Oxa1l",
  
  # Mitochondrial metabolism regulators
  "Pgc-1", "Pgc1", "Nrf1", "Nrf2", "Tfb1m", "Tfb2m", "Ppargc1a", "Ppargc1b", "Esrra", "Gabpa", "Gabpb",
  
  # Mitochondrial quality control
  "Pink1", "park", "PARL", "PGAM5", "USP30", "Atg1", "Atg8a", "Atg8b", "Ref(2)P", "Optineurin", "Parkin"
)

# Function to find mitochondrial biosynthesis genes in the dataset
find_mito_biosynthesis_genes <- function() {
  # Get mappings for genes
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    gene_mappings <- AnnotationDbi::select(org.Dm.eg.db, 
                                           keys = mito_biosynthesis_genes,
                                           columns = c("FLYBASE", "SYMBOL"),
                                           keytype = "SYMBOL")
  } else {
    gene_mappings <- data.frame(
      SYMBOL = mito_biosynthesis_genes,
      FLYBASE = NA,
      stringsAsFactors = FALSE
    )
  }
  
  # Find genes in our expression data
  if (exists("dds")) {
    gene_ids <- rownames(counts(dds))
    
    # Match genes by FlyBase ID or symbol
    matched_ids <- gene_ids[gene_ids %in% gene_mappings$FLYBASE]
    
    # If we have gene symbols in our data, try matching those too
    if (exists("gene_symbols") && is.character(gene_symbols) && length(gene_symbols) == length(gene_ids)) {
      matched_by_symbol <- gene_ids[gene_symbols %in% mito_biosynthesis_genes | 
                                      gene_symbols %in% gene_mappings$SYMBOL]
      matched_ids <- unique(c(matched_ids, matched_by_symbol))
    }
    
    # Return the matched gene IDs
    return(matched_ids)
  } else {
    cat("DESeq2 object 'dds' not found. Cannot retrieve expression data.\n")
    return(NULL)
  }
}

# Get normalized counts for visualization
get_gene_expression_over_time <- function(gene_ids) {
  if (!exists("dds") || !exists("metadata_filtered")) {
    cat("Required data (dds or metadata_filtered) not found.\n")
    return(NULL)
  }
  
  # Extract normalized counts for the genes
  normalized_counts <- counts(dds, normalized=TRUE)[gene_ids, , drop=FALSE]
  
  # Create a data frame for plotting
  count_data <- data.frame(
    sample_id = rep(colnames(normalized_counts), each = length(gene_ids)),
    gene_id = rep(rownames(normalized_counts), times = ncol(normalized_counts)),
    count = as.numeric(t(normalized_counts)),
    stringsAsFactors = FALSE
  )
  
  # Merge with metadata to get time, matched status, and group info
  plot_data <- merge(count_data, metadata_filtered, by = "sample_id")
  
  # Try to add gene symbols if available
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    # Filter for FlyBase IDs
    flybase_pattern <- "^FB"
    fb_ids <- unique(plot_data$gene_id[grepl(flybase_pattern, plot_data$gene_id)])
    
    if (length(fb_ids) > 0) {
      # Map these IDs to symbols
      symbols <- AnnotationDbi::mapIds(org.Dm.eg.db, 
                                       keys = fb_ids,
                                       column = "SYMBOL", 
                                       keytype = "FLYBASE",
                                       multiVals = "first")
      
      # Create a lookup table
      symbol_lookup <- data.frame(
        gene_id = names(symbols),
        gene_symbol = as.character(symbols),
        stringsAsFactors = FALSE
      )
      
      # Merge with plot data
      plot_data <- merge(plot_data, symbol_lookup, by = "gene_id", all.x = TRUE)
    }
  }
  
  # If gene_symbol column doesn't exist or has NAs, create a simplified version from gene_id
  if (!"gene_symbol" %in% colnames(plot_data) || any(is.na(plot_data$gene_symbol))) {
    if (!"gene_symbol" %in% colnames(plot_data)) {
      plot_data$gene_symbol <- NA
    }
    
    # Fill in NAs with simplified gene ID
    na_symbols <- is.na(plot_data$gene_symbol)
    plot_data$gene_symbol[na_symbols] <- gsub("^FBgn0*", "g", plot_data$gene_id[na_symbols])
  }
  
  return(plot_data)
}

# Function to create plots for gene expression across timepoints
create_time_plots <- function(plot_data) {
  # Ensure time is treated as a factor for proper plotting
  plot_data$time <- factor(plot_data$time)
  
  # 1. Overall expression trends by group and matched status
  p1 <- ggplot(plot_data, aes(x = time, y = count, color = group, linetype = matched)) +
    stat_summary(aes(group = interaction(group, matched)), 
                 fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(aes(group = interaction(group, matched)), 
                 fun = mean, geom = "point", size = 3) +
    labs(title = "Mitochondrial Biosynthesis Gene Expression Over Time",
         subtitle = "Average expression across all biosynthesis genes",
         x = "Time Point",
         y = "Normalized Expression Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right"
    )
  
  print(p1)
  
  # 2. Expression by gene for the top genes (most variable)
  # Calculate variance for each gene
  gene_variance <- plot_data %>%
    group_by(gene_symbol) %>%
    summarize(variance = var(count)) %>%
    arrange(desc(variance))
  
  # Select top variable genes
  top_genes <- head(gene_variance$gene_symbol, 12)
  
  # Filter data for top genes
  top_genes_data <- plot_data %>%
    filter(gene_symbol %in% top_genes)
  
  # Create a faceted plot
  p2 <- ggplot(top_genes_data, aes(x = time, y = count, color = group, linetype = matched)) +
    stat_summary(aes(group = interaction(group, matched)), 
                 fun = mean, geom = "line", linewidth = 1) +
    stat_summary(aes(group = interaction(group, matched)), 
                 fun = mean, geom = "point", size = 2) +
    facet_wrap(~gene_symbol, scales = "free_y") +
    labs(title = "Top Variable Mitochondrial Biosynthesis Genes",
         subtitle = "Expression patterns across timepoints by group and matched status",
         x = "Time Point",
         y = "Normalized Expression Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p2)
  

# Function to compare expression fold changes between groups
analyze_fold_changes <- function(mito_genes) {
  if (!exists("dds")) {
    cat("DESeq2 object 'dds' not found. Cannot calculate fold changes.\n")
    return(NULL)
  }
  
  # Define the comparisons we want to make
  comparisons <- list(
    c("group", "FF", "AA"),
    c("group", "FA", "AA"),
    c("group", "FF", "FA")
  )
  
  # Create a list to store results
  results_list <- list()
  
  # Calculate results for each comparison
  for (comp in comparisons) {
    comp_name <- paste(comp[2], "vs", comp[3])
    results_list[[comp_name]] <- results(dds, contrast = comp)
  }
  
  # Extract log2 fold changes for mitochondrial biosynthesis genes
  fc_data <- data.frame(
    gene_id = mito_genes,
    FC_FF_AA = results_list[["FF vs AA"]]$log2FoldChange[mito_genes],
    FC_FA_AA = results_list[["FA vs AA"]]$log2FoldChange[mito_genes],
    FC_FF_FA = results_list[["FF vs FA"]]$log2FoldChange[mito_genes],
    padj_FF_AA = results_list[["FF vs AA"]]$padj[mito_genes],
    padj_FA_AA = results_list[["FA vs AA"]]$padj[mito_genes],
    padj_FF_FA = results_list[["FF vs FA"]]$padj[mito_genes],
    stringsAsFactors = FALSE
  )
  
  # Add significance columns
  fc_data$sig_FF_AA <- fc_data$padj_FF_AA < 0.05
  fc_data$sig_FA_AA <- fc_data$padj_FA_AA < 0.05
  fc_data$sig_FF_FA <- fc_data$padj_FF_FA < 0.05
  
  # Add gene symbols
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
    # Filter for FlyBase IDs
    flybase_pattern <- "^FB"
    fb_ids <- fc_data$gene_id[grepl(flybase_pattern, fc_data$gene_id)]
    
    if (length(fb_ids) > 0) {
      # Map these IDs to symbols
      symbols <- AnnotationDbi::mapIds(org.Dm.eg.db, 
                                       keys = fb_ids,
                                       column = "SYMBOL", 
                                       keytype = "FLYBASE",
                                       multiVals = "first")
      
      # Create a lookup table and merge
      symbol_lookup <- data.frame(
        gene_id = names(symbols),
        gene_symbol = as.character(symbols),
        stringsAsFactors = FALSE
      )
      
      fc_data <- merge(fc_data, symbol_lookup, by = "gene_id", all.x = TRUE)
    }
  }
  
  # If gene_symbol column doesn't exist or has NAs, create from gene_id
  if (!"gene_symbol" %in% colnames(fc_data)) {
    fc_data$gene_symbol <- gsub("^FBgn0*", "g", fc_data$gene_id)
  } else if (any(is.na(fc_data$gene_symbol))) {
    na_symbols <- is.na(fc_data$gene_symbol)
    fc_data$gene_symbol[na_symbols] <- gsub("^FBgn0*", "g", fc_data$gene_id[na_symbols])
  }
  
  # Create plots
  # 1. FF/AA vs FA/AA (mismatched effect)
  p1 <- ggplot(fc_data, aes(x = FC_FF_AA, y = FC_FA_AA)) +
    geom_point(aes(color = sig_FA_AA | sig_FF_AA), size = 3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(title = "Mitochondrial Biosynthesis Gene Expression Fold Changes",
         subtitle = "Comparing FF/AA (x-axis) vs FA/AA (y-axis)",
         x = "Log2 Fold Change (FF/AA)",
         y = "Log2 Fold Change (FA/AA)",
         color = "Significant") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Add labels for significant genes
  sig_genes <- fc_data$sig_FA_AA | fc_data$sig_FF_AA
  sig_data <- fc_data[sig_genes, ]
  
  if (nrow(sig_data) > 0) {
    p1 <- p1 + 
      geom_text_repel(data = sig_data, 
                      aes(label = gene_symbol),
                      box.padding = 0.5,
                      point.padding = 0.3,
                      max.overlaps = 30,
                      segment.color = "gray50")
  }
  
  print(p1)
  
  # 2. Show timepoint-specific effects if time information is available
  if (exists("metadata_filtered") && "time" %in% colnames(metadata_filtered)) {
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)[mito_genes, , drop = FALSE]
    
    # Create a data frame of counts with metadata
    count_data <- data.frame(
      sample_id = colnames(norm_counts),
      t(norm_counts),
      stringsAsFactors = FALSE
    )
    
    # Merge with metadata
    meta_counts <- merge(count_data, metadata_filtered, by = "sample_id")
    
    # Calculate mean expression by group, time, and matched status
    time_effects <- list()
    for (g in unique(meta_counts$group)) {
      for (m in unique(meta_counts$matched)) {
        # Filter data
        subset_data <- meta_counts[meta_counts$group == g & meta_counts$matched == m, ]
        
        # Calculate means for each gene and time point
        means <- aggregate(subset_data[, mito_genes], 
                           by = list(time = subset_data$time), 
                           FUN = mean)
        
        # Store results
        time_effects[[paste(g, m, sep = "_")]] <- means
      }
    }
    
    # Compare timepoint effects between groups
    # For example, compare FF matched vs FA matched over time
    if (all(c("FF_matched", "FA_matched") %in% names(time_effects))) {
      ff_matched <- time_effects[["FF_matched"]]
      fa_matched <- time_effects[["FA_matched"]]
      
      # Ensure the data has the same time points
      common_times <- intersect(ff_matched$time, fa_matched$time)
      
      if (length(common_times) > 0) {
        ff_matched <- ff_matched[ff_matched$time %in% common_times, ]
        fa_matched <- fa_matched[fa_matched$time %in% common_times, ]
        
        # Sort by time
        ff_matched <- ff_matched[order(ff_matched$time), ]
        fa_matched <- fa_matched[order(fa_matched$time), ]
        
        # Calculate log2 fold changes for each gene over time
        fc_over_time <- data.frame(time = ff_matched$time)
        
        for (gene in mito_genes) {
          fc_over_time[[gene]] <- log2(ff_matched[[gene]] / fa_matched[[gene]])
        }
        
        # Reshape for plotting
        fc_long <- reshape2::melt(fc_over_time, id.vars = "time", 
                                  variable.name = "gene_id", value.name = "log2FC")
        
        # Add gene symbols
        if (requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
          flybase_pattern <- "^FB"
          fb_ids <- unique(fc_long$gene_id[grepl(flybase_pattern, fc_long$gene_id)])
          
          if (length(fb_ids) > 0) {
            symbols <- AnnotationDbi::mapIds(org.Dm.eg.db, 
                                             keys = fb_ids,
                                             column = "SYMBOL", 
                                             keytype = "FLYBASE",
                                             multiVals = "first")
            
            symbol_lookup <- data.frame(
              gene_id = names(symbols),
              gene_symbol = as.character(symbols),
              stringsAsFactors = FALSE
            )
            
            fc_long <- merge(fc_long, symbol_lookup, by = "gene_id", all.x = TRUE)
          }
        }
        
        # If gene_symbol column doesn't exist or has NAs, create from gene_id
        if (!"gene_symbol" %in% colnames(fc_long)) {
          fc_long$gene_symbol <- as.character(fc_long$gene_id)
        } else if (any(is.na(fc_long$gene_symbol))) {
          na_symbols <- is.na(fc_long$gene_symbol)
          fc_long$gene_symbol[na_symbols] <- as.character(fc_long$gene_id[na_symbols])
        }
        
        # Identify genes with largest changes
        gene_changes <- aggregate(abs(log2FC) ~ gene_symbol, fc_long, FUN = max)
        top_changing_genes <- head(gene_changes[order(-gene_changes$`abs(log2FC)`), "gene_symbol"], 10)
        
        # Filter for top changing genes
        top_gene_data <- fc_long[fc_long$gene_symbol %in% top_changing_genes, ]
        
        # Plot fold changes over time for top genes
        p2 <- ggplot(top_gene_data, aes(x = factor(time), y = log2FC, group = gene_symbol, color = gene_symbol)) +
          geom_line(linewidth = 1.2) +
          geom_point(size = 3) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
          labs(title = "Temporal Changes in Mitochondrial Biosynthesis Genes",
               subtitle = "Log2 Fold Change (FF/FA) over time for top changing genes",
               x = "Time Point",
               y = "Log2 Fold Change (FF/FA)",
               color = "Gene") +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.position = "right"
          )
        
        print(p2)
      }
    }
  }
  
  return(fc_data)
}

# Main analysis workflow
cat("Starting analysis of mitochondrial biosynthesis genes...\n")

# Find mitochondrial biosynthesis genes in our dataset
mito_genes <- find_mito_biosynthesis_genes()

if (!is.null(mito_genes) && length(mito_genes) > 0) {
  cat("Found", length(mito_genes), "mitochondrial biosynthesis genes in the dataset.\n")
  
  # Get expression data over time
  expression_data <- get_gene_expression_over_time(mito_genes)
  
  if (!is.null(expression_data)) {
    cat("Creating expression plots across timepoints...\n")
    time_plots <- create_time_plots(expression_data)
    
    cat("Analyzing fold changes between groups...\n")
    fc_analysis <- analyze_fold_changes(mito_genes)
  } else {
    cat("Could not retrieve expression data for plotting.\n")
  }
} else {
  cat("No mitochondrial biosynthesis genes found in the dataset.\n")
}

cat("Analysis complete.\n")

