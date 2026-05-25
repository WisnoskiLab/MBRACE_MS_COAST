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
library(patchwork)

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

loc_alpha <- trans_alpha$new(dataset = w_nov_2025_rarefied, group = "Location")
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

alpha_shannon_loc <- ggplot(alpha_dat, 
                        aes(x = Location, 
                            y = Shannon, 
                            fill = factor(Location))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  labs(x = "Location",
       y = "Shannon",
       fill = "Location") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_shannon_loc

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
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_chao1

alpha_chao1_loc <- ggplot(alpha_dat, 
                      aes(x = Location, 
                          y = Chao1, 
                          fill = factor(Location))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2)  +
  labs(x = "Location",
       y = "Chao1",
       fill = "Location") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

ggsave("Output figures/alpha_shannon.pdf", plot = alpha_shannon, width = 6, height = 4, dpi = 1000)
ggsave("Output figures/alpha_chao1.pdf", plot = alpha_chao1, width = 6, height = 4, dpi = 1000)
ggsave("Output figures/alpha_chao1_loc.pdf", plot = alpha_chao1_loc, width = 6, height = 4, dpi = 1000)
ggsave("Output figures/alpha_shannon_loc.pdf", plot = alpha_shannon_loc, width = 6, height = 4, dpi = 1000)

alpha_chao1_loc

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



library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

make_stamp_plot <- function(dataset,
                            tax_rank = "Phylum",
                            group_var = "Location",
                            groups = NULL,
                            top_n = 10) {
  
  otu <- as.data.frame(dataset$otu_table)
  tax <- as.data.frame(dataset$tax_table)
  meta <- as.data.frame(dataset$sample_table)
  
  meta$SampleID <- rownames(meta)
  meta[[group_var]] <- trimws(as.character(meta[[group_var]]))
  
  if (!tax_rank %in% colnames(tax)) {
    stop("Taxonomy rank not found. Available columns: ",
         paste(colnames(tax), collapse = ", "))
  }
  
  if (!group_var %in% colnames(meta)) {
    stop("Grouping variable not found. Available metadata columns: ",
         paste(colnames(meta), collapse = ", "))
  }
  
  if (is.null(groups)) {
    groups <- unique(meta[[group_var]])
  }
  
  groups <- trimws(as.character(groups))
  
  if (length(groups) != 2) {
    stop("This plot requires exactly two groups. Current groups are: ",
         paste(groups, collapse = ", "))
  }
  
  meta_sub <- meta[meta[[group_var]] %in% groups, , drop = FALSE]
  
  cat("\nSamples per location:\n")
  print(table(meta_sub[[group_var]]))
  
  if (length(unique(meta_sub[[group_var]])) != 2) {
    stop("Only one Location was retained after filtering. Check exact spelling in groups.")
  }
  
  # Make sure OTU table has taxa as rows and samples as columns
  if (all(meta_sub$SampleID %in% colnames(otu))) {
    otu_sub <- otu[, meta_sub$SampleID, drop = FALSE]
  } else if (all(meta_sub$SampleID %in% rownames(otu))) {
    otu <- as.data.frame(t(otu))
    otu_sub <- otu[, meta_sub$SampleID, drop = FALSE]
  } else {
    stop("Sample IDs do not match between otu_table and sample_table.")
  }
  
  # Relative abundance percentage
  otu_rel <- sweep(otu_sub, 2, colSums(otu_sub), FUN = "/") * 100
  
  # Match taxonomy
  tax_sub <- tax[rownames(otu_rel), , drop = FALSE]
  
  taxon_vec <- as.character(tax_sub[[tax_rank]])
  taxon_vec[is.na(taxon_vec) | taxon_vec == "" | taxon_vec == " "] <- "Unclassified"
  
  # Aggregate ASVs/OTUs to selected taxonomic rank
  tax_abund <- rowsum(as.matrix(otu_rel), group = taxon_vec)
  
  # Create long table using base R
  long_df <- expand.grid(
    Taxon = rownames(tax_abund),
    SampleID = colnames(tax_abund),
    stringsAsFactors = FALSE
  )
  
  long_df$Abundance <- as.vector(tax_abund)
  
  long_df <- merge(
    long_df,
    meta_sub[, c("SampleID", group_var)],
    by = "SampleID",
    all.x = TRUE
  )
  
  colnames(long_df)[colnames(long_df) == group_var] <- "Group"
  long_df$Group <- factor(long_df$Group, levels = groups)
  
  long_df <- long_df[!is.na(long_df$Group), ]
  
  cat("\nLong table columns:\n")
  print(colnames(long_df))
  
  cat("\nSamples in final long table:\n")
  print(table(long_df$Group))
  
  # Select top taxa
  mean_taxa <- aggregate(
    Abundance ~ Taxon,
    data = long_df,
    FUN = mean
  )
  
  mean_taxa <- mean_taxa[order(mean_taxa$Abundance, decreasing = TRUE), ]
  top_taxa <- head(mean_taxa$Taxon, top_n)
  
  long_df <- long_df[long_df$Taxon %in% top_taxa, ]
  
  # Mean abundance by group
  mean_df <- aggregate(
    Abundance ~ Taxon + Group,
    data = long_df,
    FUN = mean
  )
  
  colnames(mean_df)[colnames(mean_df) == "Abundance"] <- "mean_abund"
  
  # Statistics
  stat_list <- lapply(unique(long_df$Taxon), function(tx) {
    
    dat <- long_df[long_df$Taxon == tx, ]
    
    g1 <- dat$Abundance[dat$Group == groups[1]]
    g2 <- dat$Abundance[dat$Group == groups[2]]
    
    pval <- wilcox.test(g1, g2)$p.value
    
    data.frame(
      Taxon = tx,
      mean_group1 = mean(g1, na.rm = TRUE),
      mean_group2 = mean(g2, na.rm = TRUE),
      diff = mean(g1, na.rm = TRUE) - mean(g2, na.rm = TRUE),
      p_value = pval
    )
  })
  
  stat_df <- do.call(rbind, stat_list)
  
  stat_df$sig <- ifelse(stat_df$p_value < 0.001, "***",
                        ifelse(stat_df$p_value < 0.01, "**",
                               ifelse(stat_df$p_value < 0.05, "*", "")))
  
  # Order taxa
  tax_order <- aggregate(
    mean_abund ~ Taxon,
    data = mean_df,
    FUN = sum
  )
  
  tax_order <- tax_order[order(tax_order$mean_abund), "Taxon"]
  
  mean_df$Taxon <- factor(mean_df$Taxon, levels = tax_order)
  stat_df$Taxon <- factor(stat_df$Taxon, levels = tax_order)
  
  # Plot A: mean relative abundance
  p1 <- ggplot(mean_df, aes(x = mean_abund, y = Taxon, fill = Group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    labs(x = "Proportions (%)", y = NULL, fill = "Location") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, colour = "black"),
      axis.title.x = element_text(size = 14),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
    )
  
  # Plot B: difference
  p2 <- ggplot(stat_df, aes(x = diff, y = Taxon)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.8) +
    geom_text(aes(label = sig), nudge_x = 0.3, color = "red", size = 5) +
    labs(
      x = paste0("Difference between proportions (%)"),
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14, colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
    )
  
  # Plot C: p-values
  p3 <- ggplot(stat_df, aes(y = Taxon, x = 1, label = signif(p_value, 3))) +
    geom_text(size = 3.5) +
    labs(x = NULL, y = NULL, title = "P-value") +
    theme_void() +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  final_plot <- p1 + p2 + p3 +
    plot_layout(widths = c(1.5, 1.1, 0.45))
  
  return(list(
    plot = final_plot,
    stats = stat_df,
    abundance_table = mean_df,
    long_table = long_df
  ))
}

unique(w_nov_2025_rarefied$sample_table$Location)
table(w_nov_2025_rarefied$sample_table$Location)

phylum_location <- make_stamp_plot( w_nov_2025_rarefied,
  tax_rank = "Phylum",
  group_var = "Location",
  groups = c("Biloxi Bay", "Pascagoula Bay"),
  top_n = 15
)

loc_stat_abund = phylum_location$plot

ggsave("Output figures/loc_stat_abund.pdf", plot = loc_stat_abund, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)

