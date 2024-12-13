**Authors:** Bhargav Pulugundla and Anjali Rai
**Last Updated:** 11-Dec-2024


Load Dependencies
``` r
library(ggplot2)
```

Set path to file and read the data
``` r
file <- "results_rnaseq_zscore"
data <- read.csv(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
```

Filter out -inf values
``` r
filtered_data <- subset(data, 
                        Log2_RPKM_avg_norm != -Inf & 
                        Log2_RPKM_avg_dep != -Inf)
```

Colour code tag
``` r
# Assign colors based on the Tag column
filtered_data$color <- ifelse(
  filtered_data$Tag == "upregulated", "red",
  ifelse(filtered_data$Tag == "downregulated", "blue", "black")
)
```

Count upregulated and downregulated genes
``` r
upregulated_count <- sum(filtered_data$Tag == "upregulated")
downregulated_count <- sum(filtered_data$Tag == "downregulated")
```

Plot
``` r
plot <- ggplot(filtered_data, 
       aes(x = Log2_RPKM_avg_norm, 
           y = Log2_RPKM_avg_dep, 
           color = color)) +
  geom_point(alpha = 0.6) +
  scale_color_identity() +  # Use colors as defined in the data
  scale_x_continuous(breaks = seq(floor(min(filtered_data$Log2_RPKM_avg_norm)),
                                  ceiling(max(filtered_data$Log2_RPKM_avg_norm)), 
                                  by = 2)) +
  scale_y_continuous(breaks = seq(floor(min(filtered_data$Log2_RPKM_avg_dep)), 
                                  ceiling(max(filtered_data$Log2_RPKM_avg_dep)), 
                                  by = 2)) +
  labs(
    title = expression("RNA-Seq"^"ORF (Z-score)"),
    x = "log2(RPKM) 0 mM Met",
    y = "log2(RPKM) 0.5 mM Met"
  ) +
  annotate("text", 
           x = ceiling(max(filtered_data$Log2_RPKM_avg_norm, na.rm = TRUE)) - 25,
           y = max(filtered_data$Log2_RPKM_avg_dep, na.rm = TRUE) - 1,
           label = paste(upregulated_count, "↑"), 
           hjust = 0, vjust = 1, size = 5, color = "red", fontface = "bold") +
  annotate("text", 
           x = ceiling(max(filtered_data$Log2_RPKM_avg_norm, na.rm = TRUE)) - 1,
           y = max(filtered_data$Log2_RPKM_avg_dep, na.rm = TRUE) - 25,  
           label = paste(downregulated_count, "↓"), 
           hjust = 1, vjust = -0.2, size = 5, color = "blue", fontface = "bold") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# Display the plot
print(plot)
```

``` r
# Save the plot as a PNG file
ggsave("dge_rnaseq_zscore.eps", width = 8, height = 6, dpi = 300)
```

