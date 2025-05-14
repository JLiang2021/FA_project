polg_id <- rownames(res_matched)[grep("PolG-1", res_matched$gene_symbol)]

if (length(polg_id) > 0) {
  # Extract the DESeq2 results for PolG-1
  polg_result <- res_matched[polg_id, ]
  
  cat("\nPolG-1 expression (proxy for mtDNA copy number):\n")
  cat("Log2 Fold Change:", polg_result$log2FoldChange, "\n")
  cat("p-value:", polg_result$pvalue, "\n")
  cat("Adjusted p-value:", polg_result$padj, "\n")
  
  # Get normalized counts for matched vs mismatched comparison
  if (exists("dds")) {
    polg_counts <- counts(dds, normalized=TRUE)[polg_id,]
    matched_samples <- colData(dds)$matched == TRUE
    unmatched_samples <- colData(dds)$matched == FALSE
    
    # Compare expression using t-test
    polg_matched <- polg_counts[, matched_samples]
    polg_unmatched <- polg_counts[, unmatched_samples]
    
    polg_ttest <- t.test(polg_matched, polg_unmatched)
    
    cat("t-test p-value:", polg_ttest$p.value, "\n")
    cat("Mean expression in matched samples:", mean(polg_matched), "\n")
    cat("Mean expression in mismatched samples:", mean(polg_unmatched), "\n")
  }
} else {
  cat("\nPolG-1 gene not found in the dataset.\n")
}
