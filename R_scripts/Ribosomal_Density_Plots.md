**Author:** Bhargav Pulugundla  
**Last Updated:** 13-Dec-2024

``` r
library(ggplot2)
library(dplyr) 
``


``` r
# Load data
ribo_normal <- read.csv("density_results/ERR2660262_RiboSeq_normalized_pos_reads.csv")
ribo_depleted <- read.csv("density_results/ERR2660266_RiboSeq_normalized_pos_reads.csv")
rna_normal <- read.csv("density_results/ERR2660264_RNASeq_normalized_pos_reads.csv")
rna_depleted <- read.csv("density_results/ERR2660269_RNASeq_normalized_pos_reads.csv")
```


``` r
# Add source labels
ribo_normal$Source <- "RiboSeq"
ribo_depleted$Source <- "RiboSeq"
rna_normal$Source <- "RNASeq"
rna_depleted$Source <- "RNASeq"

# Combine data into one data frame
ribo_normal$Variant <- "Normal"
ribo_depleted$Variant <- "Depleted"
rna_normal$Variant <- "Normal"
rna_depleted$Variant <- "Depleted"
```


``` r
# Combine all data into one dataframe
all_data <- rbind(ribo_normal, ribo_depleted, rna_normal, rna_depleted)

# Restructure data to calculate ratio of RiboSeq to RNASeq
ribo_data <- all_data %>% filter(Source == "RiboSeq")
rna_data <- all_data %>% filter(Source == "RNASeq")

# Join RiboSeq and RNASeq data by Relative Position and Variant
merged_data <- merge(ribo_data, rna_data, 
                     by = c("Relative.Position", "Variant", "Codon.Type"), 
                     suffixes = c("_RiboSeq", "_RNASeq"))

# Calculate Ratio
merged_data <- merged_data %>% 
  mutate(Ratio = Normalized.Read.Count..RPM._RiboSeq /
                    Normalized.Read.Count..RPM._RNASeq)

# Filter data for start_codon and stop_codon separately
start_codon_data <- merged_data %>% 
  filter(Codon.Type == "start_codon" & Relative.Position >= -5 
                        & Relative.Position <= 40)

stop_codon_data <- merged_data %>% 
  filter(Codon.Type == "stop_codon" & Relative.Position >= -40 
                        & Relative.Position <= 5)
```


``` r
# Create a line plot for start_codon
start_plot <- ggplot(start_codon_data, aes(x = Relative.Position, 
y = Ratio, color = Variant)) +
  geom_line(size = 1) +
  labs(
    title = "1st Replicate Start Codon",
    x = "Relative Position",
    y = "Ribosomal Load (RiboSeq RPM/RNASeq RPM)",
    color = "Variant"
  ) +
  scale_color_manual(values = c("Normal" = "black", "Depleted" = "red")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
```


``` r
# Create a line plot for stop_codon
stop_plot <- ggplot(stop_codon_data, aes(x = Relative.Position, 
y = Ratio, color = Variant)) +
  geom_line(size = 1) +
  labs(
    title = "1st Replicate Stop Codon",
    x = "Relative Position",
    y = "Ribosomal Load (RiboSeq RPM/RNASeq RPM)",
    color = "Variant"
  ) +
  scale_color_manual(values = c("Normal" = "black", "Depleted" = "red")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
```


``` r
# Display the plots
print(start_plot)
print(stop_plot)
```


``` r
# Save plots as image files
ggsave("StartCodon_Plot.eps", plot = start_plot, width = 10, height = 6)
ggsave("StopCodon_Plot.eps", plot = stop_plot, width = 10, height = 6)
```





