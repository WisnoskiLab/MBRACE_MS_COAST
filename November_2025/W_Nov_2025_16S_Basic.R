# ---- Core data handling & manipulation ----
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
library(magrittr)
library(parallel)
library(reshape2)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(MicrobiomeStat)
library(meconetcomp)
library(WGCNA)
library(ggClusterNet)
library(ape)
library(picante)
library(Biostrings)

# ---- Differential abundance & compositional analysis ----
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)

# ---- Visualization ----
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(paletteer)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggraph)

# ---- Statistical modeling ----
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(multcompView)
library(dplyr)
library(usethis)
library(nlMS)
library(iCAMP)
library(minpack.lm)
library(Hmisc)
library(Biostrings)
library(meconetcomp)

# ---- Package manager ----
library(BiocManager)


##### Importing files ########

biom = import_biom("Data files/w_nov_2025_table.biom")

metadata = import_qiime_sample_data("Data files/metadata.txt")

tree = read_tree("Data files/w_nov_2025-unrooted_tree.nwk")

rep_fasta = readDNAStringSet("Data files/w_nov_2025_repseqs.fasta", format = "fasta")

w_biom = merge_phyloseq(biom, rep_fasta, metadata, tree)

#rename columns in taxonomy table
colnames(tax_table(w_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

w_nov_2025 <- phyloseq2meco(w_biom)

# removing zero-abundance taxa and make sure sample and taxonomy tables are consistent.
w_nov_2025$tidy_dataset()

# check sequencing depth 

sample_depth <- colSums(w_nov_2025$otu_table)

summary(sample_depth)
min(sample_depth)
quantile(sample_depth, probs = c(0.05, 0.10, 0.25, 0.50, 0.75))

hist(sample_depth, breaks = 30,
     main = "Sequencing depth per sample",
     xlab = "Reads per sample")

tmp <- trans_norm$new(dataset = w_nov_2025)

# rarefaction
w_nov_2025_rarefied <- tmp$norm(method = "rarefy", sample.size = 40000)

w_nov_2025_rarefied$sample_table$Salinity <- factor(
  w_nov_2025_rarefied$sample_table$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

w_nov_2025_rarefied$tidy_dataset()


w_nov_2025_rarefied$cal_abund()
w_nov_2025_rarefied$cal_alphadiv()
w_nov_2025_rarefied$cal_betadiv()


loc_alpha <- trans_alpha$new(dataset = w_nov_2025_rarefied, group = "Location")
loc_alpha$cal_diff(method = "anova", formula = "Location+Salinity")

loc_alpha <- trans_alpha$new(dataset = w_nov_2025_rarefied, group = "Location", by_group = "Salinity")
loc_alpha$cal_diff(method = "wilcox")
loc_alpha$plot_alpha(measure = "Shannon")

alpha = loc_alpha[["data_alpha"]]

write.csv(alpha, file = "Alpha.csv")

alpha_dat = read_csv("Alpha-dat.csv")

alpha_dat$Salinity <- factor(alpha_dat$Salinity,
                             levels = c("Freshwater", "Moderate Salinity", "High Salinity"))

w_nov_2025_rarefied$sample_table <- factor(w_nov_2025_rarefied$sample_table,
                         levels = c("Freshwater", "Moderate Salinity", "High Salinity"))


alpha_shannon <- ggplot(alpha_dat, 
                                  aes(x = factor(Salinity), 
                                      y = Shannon, 
                                      fill = factor(Salinity))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  facet_wrap(~ Location, scales = "free_y") +
  labs(x = "Salinity",
       y = "Shannon",
       fill = "Salinity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_shannon

alpha_chao1 <- ggplot(alpha_dat, 
                        aes(x = factor(Salinity), 
                            y = Chao1, 
                            fill = factor(Salinity))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  facet_wrap(~ Location, scales = "free_y") +
  labs(x = "Salinity",
       y = "Chao1",
       fill = "Salinity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_chao1

ggsave("Output figures/alpha_shannon.pdf", plot = alpha_shannon, width = 6, height = 4, dpi = 1000)
ggsave("Output figures/alpha_chao1.pdf", plot = alpha_chao1, width = 6, height = 4, dpi = 1000)

#Beta diversity

beta_loc <- trans_beta$new(dataset = w_nov_2025_rarefied, group = "Location", measure = "bray")

beta_loc$cal_ordination(method = "PCoA")

nov_2025_loc = beta_loc$plot_ordination(plot_color = "Salinity", plot_shape = "Location", plot_type = c("point", "ellipse"), point_size = 2) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Increase x-axis text size
        axis.text.y = element_text(size = 18), # Increase y-axis text size
        axis.title.x = element_text(size = 20), # Increase x-axis label size
        axis.title.y = element_text(size = 20), # Increase y-axis label size
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

nov_2025_loc

ggsave("Output figures/nov_2025_loc.pdf", plot = nov_2025_loc, device = "pdf", width = 7, 
       height = 5, units = "in", dpi = 1000)


top15_palette <- c("#290AD8", "#264DFF", 
                   "#3FA0FF", "#AAF7FF", 
                   "#B2FFB2", "#FFFFBF", 
                   "#FFE099", "#FFAD72", 
                   "#F76D5E", "#D82632", 
                   "#A50021", '#F3B79C',
                   '#D64D72', '#841859',
                   '#312A56', '#6B6100',
                   '#004F78', '#0096B5',
                   "#427e83", "#686a47")

Abundance_L <- trans_abund$new(dataset = w_nov_2025_rarefied, taxrank = "Phylum", ntaxa = 15)


box_abundance_loc = Abundance_L$plot_box(group = "Location",  xtext_angle = 30)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

box_abundance_sal = Abundance_L$plot_box(group = "Salinity",  xtext_angle = 30)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

ggsave("Output figures/box_abundance_loc.pdf", plot = box_abundance_loc, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)

ggsave("Output figures/box_abundance_sal.pdf", plot = box_abundance_sal, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)

stack_loc <- trans_abund$new(dataset = w_nov_2025_rarefied, taxrank = "Phylum", ntaxa = 15, groupmean = "Location")

stack_sal <- trans_abund$new(dataset = w_nov_2025_rarefied, taxrank = "Phylum", ntaxa = 15, groupmean = "Salinity")

stack_loc_plot <- stack_loc$plot_bar(others_color = "grey70", legend_text_italic = FALSE, color_values = top15_palette)  +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border



stack_sal_plot <- stack_sal$plot_bar(others_color = "grey70", legend_text_italic = FALSE, color_values = top15_palette)  +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

ggsave("Output figures/stack_loc_plot.pdf", plot = stack_loc_plot, device = "pdf", width = 8, 
       height = 7, units = "in", dpi = 1000)

ggsave("Output figures/stack_sal_plot.pdf", plot = stack_sal_plot, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)
