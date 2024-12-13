**Authors:** Bhargav Pulugundla and Anjali Rai
**Last Updated:** 11-Dec-2024


Load Dependencies
``` r
library(ggplot2)
```

Set paths to Ribo-seq and RNA-seq data
``` r
ribo_file <- "results_riboseq_zscore"
rna_file <- "results_rnaseq_zscore"

# Load the data
ribo_data <- read.csv(ribo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rna_data <- read.csv(rna_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
```

Filter out -inf values
``` r
# Filter out -Inf values from ribo-seq data
ribo_data_filtered <- subset(ribo_data, 
                             Avg_RPKM_norm != -Inf & Avg_RPKM_dep != -Inf & 
                             !is.na(Avg_RPKM_norm) & !is.na(Avg_RPKM_dep))

# Filter out -Inf values from RNA-seq data
rna_data_filtered <- subset(rna_data, 
                             Avg_RPKM_norm != -Inf & Avg_RPKM_dep != -Inf & 
                             !is.na(Avg_RPKM_norm) & !is.na(Avg_RPKM_dep))
```

``` r
# Subset to retain only relevant columns for ribo-seq data
ribo_data_filtered <- ribo_data_filtered[, c("Gene", "Avg_RPKM_norm", "Avg_RPKM_dep")]
# Subset to retain only relevant columns for RNA-seq data
rna_data_filtered <- rna_data_filtered[, c("Gene", "Avg_RPKM_norm", "Avg_RPKM_dep")]
```

Merge based on common genes
``` r
# Find common genes between ribo-seq and RNA-seq data
common_genes <- intersect(ribo_data_filtered$Gene, rna_data_filtered$Gene)
```

``` r
# Subset both datasets to include only common genes
ribo_common <- ribo_data_filtered[ribo_data_filtered$Gene %in% common_genes, ]
rna_common <- rna_data_filtered[rna_data_filtered$Gene %in% common_genes, ]
```

``` r
# Merge datasets on Gene column to align rows
merged_data <- merge(ribo_common, rna_common, by = "Gene", suffixes = c("_ribo", "_rna"))
```

Compute Log2 Ribosomal Load for both conditions
``` r
# Calculate ribosomal load for normal and depleted conditions
merged_data$ribo_load_norm <- merged_data$Avg_RPKM_norm_ribo / merged_data$Avg_RPKM_norm_rna
merged_data$ribo_load_dep <- merged_data$Avg_RPKM_dep_ribo / merged_data$Avg_RPKM_dep_rna
```

``` r
# Log2 transform the ribosomal load values
merged_data$ribo_load_norm <- log2(merged_data$ribo_load_norm)
merged_data$ribo_load_dep <- log2(merged_data$ribo_load_dep)
```

``` r
# Filter out rows with non-finite ribosomal load values
merged_data <- merged_data[is.finite(merged_data$ribo_load_norm) & is.finite(merged_data$ribo_load_dep), ]
```


``` r
# Calculate fold-change as the difference (ribo_load_dep - ribo_load_norm)
merged_data$fold_change <- merged_data$ribo_load_dep - merged_data$ribo_load_norm
```


``` r
# Tag genes based on fold-change thresholds
merged_data$tag <- ifelse(merged_data$fold_change > 1.5, "Upregulated",
                                 ifelse(merged_data$fold_change < -1.5, "Downregulated", "Not Significant"))
```

Colour code genes based on fold-change thresholds for both conditions
``` r
# Assign colors based on tags
merged_data$color <- ifelse(merged_data$tag == "Upregulated", "red",
                     ifelse(merged_data$tag == "Downregulated", "blue", "black"))
```

Count upregulated and downregulated genes
``` r
# Count upregulated genes (Log2 fold-change > 1.5)
upregulated_count <- sum(merged_data$color == "red")

# Count downregulated genes (Log2 fold-change < -1.5)
downregulated_count <- sum(merged_data$color == "blue")
```


Create plot
``` r
plot <- ggplot(merged_data, 
       aes(x = ribo_load_norm, 
           y = ribo_load_dep, 
           color = color)) +
  geom_point(alpha = 0.6) +
  scale_color_identity() +  # Use colors as defined in the data
    scale_x_continuous(breaks = seq(floor(min(merged_data$ribo_load_norm)),
                                  ceiling(max(merged_data$ribo_load_norm)), 
                                  by = 2),
                       limits = c(-11, ceiling(max(merged_data$ribo_load_norm)))
    ) +
  scale_y_continuous(breaks = seq(floor(min(merged_data$ribo_load_dep)), 
                                  ceiling(max(merged_data$ribo_load_dep)), 
                                  by = 2)) +
  labs(
    title = expression("Ribosomal Load"^"ORF"),
    x = "log2(ribosomal load) 0 mM Met",
    y = "log2(ribosomal load) 0.5 mM Met"
  ) +
annotate("text", 
           x = ceiling(min(merged_data$ribo_load_norm, na.rm = TRUE)),
           y = max(merged_data$ribo_load_dep, na.rm = TRUE),
           label = paste(upregulated_count, "↑"), 
           hjust = 0, vjust = 1, size = 5, color = "red", fontface = "bold") +
  annotate("text", 
           x = ceiling(max(merged_data$ribo_load_norm, na.rm = TRUE)),
           y = min(merged_data$ribo_load_dep, na.rm = TRUE),  
           label = paste(downregulated_count, "↓"), 
           hjust = 1, vjust = -0.2, size = 5, color = "blue", fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
```


``` r
# Display the plot
print(plot)
```

``` r
# Save the plot as a PNG file
ggsave("dge_riboload.png", width = 8, height = 6, dpi = 300)
```
