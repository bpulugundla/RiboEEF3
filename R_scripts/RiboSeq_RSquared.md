**Author:** Anjali Rai
**Last Updated:** 10-Dec-2024


``` r
# Install packages
# install.packages("BiocManager")
# BiocManager::install("Rsubread")
# BiocManager::install("rtracklayer")
```


``` r
# Load libraries
library(Rsubread)  # For feature counting
library(ggplot2)   # For plotting
library(rtracklayer)  # For annotation loading
```

# Prepare BAM Files
``` r
# Path to BAM files
bam_file_1 <- "ERR2660266.bam"
bam_file_2 <- "ERR2660267.bam"

# Path to the reference annotation GTF file 
annotation_file <- "Genome.gtf"
```

# Perform Feature Counting
Counts how many reads from each BAM file overlap with the features in the annotation file.

``` r
# Count reads for replicate 1
counts_rep1 <- featureCounts(bam_file_1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                             GTF.featureType = "CDS", GTF.attrType = "gene_id", 
                             allowMultiOverlap = TRUE)
# Count reads for replicate 2
counts_rep2 <- featureCounts(bam_file_2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                             GTF.featureType = "CDS", GTF.attrType = "gene_id", 
                             allowMultiOverlap = TRUE)
```

# Extract the Raw Counts
``` r
# Extract raw counts for both replicates
orf_counts_1 <- counts_rep1$counts
orf_counts_2 <- counts_rep2$counts
```

# Calculate RPM 
``` r
total_counts_rep1 <- sum(orf_counts_1)  # Total reads in replicate 1
total_counts_rep2 <- sum(orf_counts_2)  # Total reads in replicate 2
```

``` r
# Function to calculate RPM (Reads Per Million)
calc_rpm <- function(counts, total_counts) {
  return((counts / total_counts) * 1e6)  # Normalise counts to RPM
}

# Calculate RPM values for both replicates
rpm_rep1 <- calc_rpm(orf_counts_1, total_counts_rep1)
rpm_rep2 <- calc_rpm(orf_counts_2, total_counts_rep2)
```

# Create Data Frame for Plotting
``` r
# Create a data frame with the RPM values
data <- data.frame(rpm_rep1, rpm_rep2)
```

``` r
# Rename columns
colnames(data) <- c("Replicate1_RPM", "Replicate2_RPM")

# Check the first few rows of the data frame
print(head(data))
```


# Filter and Fit Linear Model
``` r
# Remove rows where the RPM values are 0, as log(0) is undefined
filtered_data <- data[data$Replicate1_RPM > 0 & data$Replicate2_RPM > 0, ]

# Fit a linear model on log-transformed RPM values from both replicates
lm_model <- lm(log10(Replicate2_RPM) ~ log10(Replicate1_RPM), data = filtered_data)
```

``` r
# Calculate R-squared value
r_squared <- summary(lm_model)$r.squared
```

# Plot
``` r
# Generate the scatter plot with R-squared value
ggplot(filtered_data, aes(x = Replicate1_RPM, y = Replicate2_RPM)) +
  geom_point(alpha = 0.5) +  # Add points with transparency
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear regression line
  theme_bw() +  # Set the background to white
  scale_x_log10(
    limits = c(1e-2, 1e6),  # X-axis range
    breaks = scales::trans_breaks("log10", function(x) 10^x),  # Log10 breaks
    labels = scales::trans_format("log10", scales::math_format(10^.x))  # Log10 labels
  ) +
  scale_y_log10(
    limits = c(1e-2, 1e6),  # Y-axis range
    breaks = scales::trans_breaks("log10", function(x) 10^x),  # Log10 breaks
    labels = scales::trans_format("log10", scales::math_format(10^.x))  # Log10 labels
  ) +
  labs(
    x = expression("1st Replicate, RPM "),  # X-axis label
    y = expression("2nd Replicate, RPM "),  # Y-axis label
    title = "ORF coverage: Ribo-Seq (0.5 mM methionine)"
  ) +
  # Annotate the plot with R-squared value
  annotate(
    "text", 
    x = 0.1, y = 1e4,  # Position of the text annotation
    label = bquote(R^2 == .(signif(r_squared, 3))),  # Format R-squared to 3 significant figures
    parse = FALSE,
    size = 5, color = "blue"
  )
```

``` r
# Save the plot as a PNG file
ggsave("Dep_Ribo_Seq.png", width = 8, height = 6, dpi = 300)
```

