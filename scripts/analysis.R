################################################################################
###################### INSTALL AND LOAD REQUIRED PACKAGES ######################
################################################################################

# Define a vector with all required packages
packages <- c("ggplot2", "vegan", "indicspecies", "gridExtra", "tidyverse",
              "patchwork", "scales", "tibble", "cowplot", "ggrepel")

# Install only the packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))

# Install Bioconductor manager if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install phyloseq and microbiome packages from Bioconductor
BiocManager::install("phyloseq")
BiocManager::install("microbiome")

################################################################################
################################ LOAD PACKAGES  ################################
################################################################################

library(phyloseq)      # For microbial community analysis
library(ggplot2)       # For creating plots
library(vegan)         # For ecological analysis
library(microbiome)    # For microbiome-specific utilities
library(indicspecies)  # For indicator species analysis
library(gridExtra)     # For combining plots
library(tidyverse)     # For data manipulation and visualization
library(stringr)       # For string operations
library(dplyr)         # Data manipulation
library(tidyr)         # Data reshaping
library(patchwork)     # Plot composition
library(scales)        # Scales for ggplot2
library(tibble)        # Enhanced data frames
library(forcats)       # Factor handling
library(readr)         # Reading files
library(cowplot)       # Another plot combining package
library(ggrepel)       # For better label placement in ggplot

################################################################################
############################### LOAD INPUT DATA  ###############################
################################################################################

# Read OTU (Operational Taxonomic Unit) count table
otu <- read.table("results/OTUs.97.rep.count_table", header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE)

# Read taxonomy assignment
taxonomy_raw <- read.table("results/OTUs.97.cons.taxonomy", sep = "\t", header = FALSE,
                           col.names = c("OTU", "Size", "Taxonomy"),
                           stringsAsFactors = FALSE, fill = TRUE)

# Read sample metadata
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t", row.names = 1)

# Split taxonomy into hierarchical levels
tax_split <- str_split_fixed(taxonomy_raw$Taxonomy, ";", 7)
colnames(tax_split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_split) <- taxonomy_raw$OTU
tax_mat <- as.matrix(tax_split)  # Convert to matrix for use in phyloseq

# Convert OTU table to matrix
otu_mat <- as.matrix(otu)

# Create phyloseq object
physeq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(metadata)
)

# Create directory for figures
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

################################################################################
############# FIGURE 1: ALPHA DIVERSITY (Chao1, Shannon, Simpson)  #############
################################################################################

alpha_df <- estimate_richness(physeq, measures = c("Chao1", "Shannon", "Simpson"))
alpha_df$Primer <- sample_data(physeq)$Primer  # Add primer information

# Normality tests (Shapiro-Wilk) for each alpha diversity metric
shapiro.test(alpha_df$Chao1)   # Normality test for Chao1
shapiro.test(alpha_df$Shannon) # Normality test for Shannon
shapiro.test(alpha_df$Simpson) # Normality test for Simpson

# Hypothesis testing between primers
# ANOVA for metrics that follow normal distribution
anova_chao1 <- aov(Chao1 ~ Primer, data = alpha_df)
summary(anova_chao1)

anova_shannon <- aov(Shannon ~ Primer, data = alpha_df)
summary(anova_shannon)

anova_simpson <- aov(Simpson ~ Primer, data = alpha_df)
summary(anova_simpson)


# Chao1 Richness
p1 <- ggplot(alpha_df, aes(x = Primer, y = Chao1, fill = Primer)) +
  geom_boxplot() +
  ggtitle("Chao1 Richness") +
  scale_fill_manual(values = c("V3V4" = "darkslategray4", "V4V5" = "darkorange"))

# Shannon Diversity
p2 <- ggplot(alpha_df, aes(x = Primer, y = Shannon, fill = Primer)) +
  geom_boxplot() +
  ggtitle("Shannon Diversity") +
  scale_fill_manual(values = c("V3V4" = "darkslategray4", "V4V5" = "darkorange"))

# Simpson Evenness
p3 <- ggplot(alpha_df, aes(x = Primer, y = Simpson, fill = Primer)) +
  geom_boxplot() +
  ggtitle("Simpson Evenness") +
  scale_fill_manual(values = c("V3V4" = "darkslategray4", "V4V5" = "darkorange"))

ggsave(file.path(fig_dir, "Figura1_Alpha_Diversity.pdf"),
       plot = plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3), width = 12, height = 4)

ggsave(file.path(fig_dir, "Figura1_Alpha_Diversity.png"), 
       plot = plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3), width = 12, height = 4, dpi = 300)

################################################################################
########## FIGURE 2: BETA DIVERSITY (PCoA with Bray-Curtis distance)  ##########
################################################################################

ord <- ordinate(physeq, method = "PCoA", distance = "bray")

# Calculate Bray-Curtis distance matrix
dist_bray <- phyloseq::distance(physeq, method = "bray")

# Extract sample metadata as data.frame
sample_df <- as(sample_data(physeq), "data.frame")

# Check homogeneity of group dispersions
beta_disp <- betadisper(dist_bray, sample_df$Primer)
permutest(beta_disp)

# PERMANOVA to test for significant differences between groups (Primers)
permanova_result <- adonis2(dist_bray ~ Primer, data = sample_df, permutations = 999)
print(permanova_result)


pcoa_plot <- plot_ordination(physeq, ord, color = "Primer") +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("V3V4" = "darkslategray4", "V4V5" = "darkorange"))

ggsave(file.path(fig_dir,"Figura2_PCoA.pdf"), pcoa_plot, width = 6, height = 5)
ggsave(file.path(fig_dir,"Figura2_PCoA.png"), pcoa_plot, width = 6, height = 5, dpi = 300)

################################################################################
########## FIGURES 3–5: TAXONOMIC COMPOSITION AND RELATIVE ABUNDANCE  ##########
################################################################################

#### Generate Taxonomic Plots of Microbial Composition

# This function generates comparative bar and pie plots of ASV (Amplicon Sequence Variants)
# counts and relative abundances across primer groups. It highlights taxonomic differences
# by phylum or class and marks taxa with a ≥10x difference in relative abundance.

#' @param phyobj A `phyloseq` object containing OTU table, taxonomy, and metadata.
#' @param kingdom_filter Optional. A character vector to filter by taxonomic kingdom 
#'        (e.g., `"Bacteria"` or `"Archaea"`). If `NULL`, no filtering is applied.
#' @param tax_rank The taxonomic rank to summarize at. Options: `"Phylum"`, `"Class"`, etc.
#' @param filename The output file name (must include `.pdf` or `.png` extension).
#' @param top_cutoff Minimum relative abundance threshold for inclusion in the plot.
#' @param max_taxa Maximum number of taxa (e.g., phyla or classes) to show individually.
#'        Remaining taxa will be lumped into `"Other"`.
#' @return Saves the plot as a file (PDF or PNG) and prints a message upon success.
#' @examples
#' make_fig(physeq, kingdom_filter = "Bacteria", tax_rank = "Phylum",
#'          filename = "Figura4_Bacteria.pdf", top_cutoff = 0.01, max_taxa = 15)
#'
make_fig <- function(phyobj,
                     kingdom_filter = NULL,
                     tax_rank = "Phylum",
                     filename = "Figura.pdf",
                     top_cutoff = 0.01,
                     max_taxa = 15,
                     fig_dir = ".") {
  # 1. Extract relevant tables from phyloseq object
  tax_tab <- as.data.frame(tax_table(phyobj))
  otu_tab <- as.data.frame(otu_table(phyobj))
  sample_tab <- as.data.frame(sample_data(phyobj))
  
  # 2. Optionally filter by kingdom
  if (!is.null(kingdom_filter)) {
    keep_otus <- rownames(tax_tab)[tax_tab$Kingdom %in% kingdom_filter]
    tax_tab <- tax_tab[keep_otus, ]
    otu_tab <- otu_tab[keep_otus, , drop = FALSE]
  }
  
  # 3. Reshape OTU table into long format
  otu_long <- otu_tab %>%
    rownames_to_column("OTU") %>%
    pivot_longer(-OTU, names_to = "Sample", values_to = "Count") %>%
    filter(Count > 0) %>%
    left_join(sample_tab %>% rownames_to_column("Sample"), by = "Sample") %>%
    left_join(tax_tab %>% rownames_to_column("OTU") %>%
                select(OTU, Rank = all_of(tax_rank)), by = "OTU") %>%
    filter(!is.na(Rank), Rank != "", !grepl("unclassified", Rank, ignore.case = TRUE)) %>%
    mutate(Rank = fct_lump(fct_infreq(Rank), n = max_taxa)) %>%
    filter(!is.na(Primer))
  
  # 4. Create custom color palette for taxa
  all_ranks <- fct_drop(otu_long$Rank) %>% fct_infreq() %>% fct_drop() %>% levels()
  palette <- RColorBrewer::brewer.pal(max(length(all_ranks), 3), "Set3")
  phylum_colors <- setNames(palette[seq_along(all_ranks)], all_ranks)
  
  # 5. Count ASVs per primer group and taxonomic rank
  asv_counts <- otu_long %>%
    distinct(Primer, OTU, Rank) %>%
    group_by(Primer, Rank) %>%
    summarise(ASVs = n(), .groups = "drop")
  
  # 6. Calculate relative proportions for 10x difference test
  total_asvs <- asv_counts %>%
    group_by(Primer) %>%
    summarise(Total = sum(ASVs))
  
  asv_counts <- asv_counts %>%
    left_join(total_asvs, by = "Primer") %>%
    mutate(prop = ASVs / Total)
  
  # 7. Identify taxa with ≥10-fold difference in relative abundance
  prop_wide <- asv_counts %>%
    select(Primer, Rank, prop) %>%
    pivot_wider(names_from = Primer, values_from = prop, values_fill = 1e-6) %>%
    mutate(fold = pmax(V3V4 / V4V5, V4V5 / V3V4),
           mark = ifelse(fold >= 10, "†", ""))
  
  label_map <- prop_wide %>%
    mutate(RankLabel = paste0(Rank, mark)) %>%
    select(Rank, RankLabel)
  
  # 8. Update OTU table with modified rank labels
  otu_long <- otu_long %>%
    left_join(label_map, by = "Rank")
  
  asv_counts <- asv_counts %>%
    left_join(label_map, by = "Rank") %>%
    mutate(ASVs = ifelse(Primer == "V3V4", -ASVs, ASVs))
  
  # 9. Create bar plot: ASV counts per taxonomic group
  p_bar <- ggplot(asv_counts, aes(x = ASVs, y = fct_rev(RankLabel), fill = RankLabel)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = abs(ASVs)),
              hjust = ifelse(asv_counts$Primer == "V3V4", 1.1, -0.1), size = 3) +
    scale_x_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
    scale_fill_manual(values = phylum_colors[label_map$Rank]) +
    labs(x = "ASVs per Group", y = NULL) +
    facet_wrap(~Primer, ncol = 2, scales = "free_x") +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle("A")
  
  # 10. Calculate relative abundances per group and taxon
  rel_abund <- otu_long %>%
    group_by(Sample) %>%
    mutate(RelAbund = Count / sum(Count)) %>%
    group_by(Primer, RankLabel) %>%
    summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
    group_by(Primer) %>%
    mutate(prop = RelAbund / sum(RelAbund),
           ypos = cumsum(prop) - 0.5 * prop,
           label = NA)
  
  # 11. Define pie chart theme
  pie_theme <- theme_void() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # 12. Pie chart for V3V4
  p_pie_v3v4 <- ggplot(rel_abund %>% filter(Primer == "V3V4"),
                       aes(x = "", y = prop, fill = RankLabel)) +
    geom_col(color = "white", width = 1) +
    coord_polar("y") +
    ggrepel::geom_text_repel(aes(y = ypos, label = label),
                             size = 3, show.legend = FALSE,
                             nudge_x = 1, segment.color = "grey50",
                             direction = "y", max.overlaps = 20) +
    scale_fill_manual(values = phylum_colors[label_map$Rank]) +
    ggtitle("B") +
    pie_theme
  
  # 13. Pie chart for V4V5
  p_pie_v4v5 <- ggplot(rel_abund %>% filter(Primer == "V4V5"),
                       aes(x = "", y = prop, fill = RankLabel)) +
    geom_col(color = "white", width = 1) +
    coord_polar("y") +
    ggrepel::geom_text_repel(aes(y = ypos, label = label),
                             size = 3, show.legend = FALSE,
                             nudge_x = 1, segment.color = "grey50",
                             direction = "y", max.overlaps = 20) +
    scale_fill_manual(values = phylum_colors[label_map$Rank]) +
    ggtitle("C") +
    pie_theme
  
  # 14. Combine bar and pie plots into one layout
  final_fig <- p_bar / (p_pie_v3v4 | p_pie_v4v5) +
    plot_layout(heights = c(1, 1.2), guides = "collect") &
    theme(legend.position = "bottom")
  
  # 15. Save the figure
  # Full path for the PDF file
  pdf_path <- file.path(fig_dir, filename)
  
  # Construct the PNG path by replacing the .pdf extension with .png
  png_path <- sub("\\.pdf$", ".png", pdf_path)
  
  # Save the figure as PDF
  ggsave(pdf_path, final_fig, width = 10, height = 9)
  
  # Save the figure as PNG with 300 dpi resolution
  ggsave(png_path, final_fig, width = 10, height = 9, dpi = 300)
  
  message("Figure saved as ", pdf_path, " and ", png_path)
  
}

# FIGURE 3: Archaea only
tryCatch({
  make_fig(physeq,
         kingdom_filter = "Archaea",
         tax_rank = "Phylum",
         filename = "Figura3_Archaea.pdf",
         top_cutoff = 0.001,
         fig_dir = fig_dir)
}, error = function(e) {
  message("Error in Figura3_Archaea: ", e$message)
})

# FIGURE 4: Bacteria only
tryCatch({
  make_fig(physeq,
         kingdom_filter = "Bacteria",
         tax_rank = "Phylum",
         filename = "Figura4_Bacteria.pdf",
         top_cutoff = 0.01,
         max_taxa = 15,
         fig_dir = fig_dir)
}, error = function(e) {
  message("Error in Figura4_Bacteria: ", e$message)
})

# FIGURE 5: Proteobacteria + Campylobacterota at class level
phy_protcampy <- prune_taxa(tax_table(physeq)[, "Phylum"] %in%
                              c("Proteobacteria", "Campylobacterota"), physeq)
tryCatch({
  make_fig(phy_protcampy,
         tax_rank = "Class",
         filename = "Figura5_ProteoCampy.pdf",
         top_cutoff = 0.01,
         max_taxa = 10,
         fig_dir = fig_dir)
}, error = function(e) {
  message("Error in Figura5_ProteoCampy: ", e$message)
})

################################################################################
#################### FIGURE 6: INDICATOR SPECIES ANALYSIS   ####################
################################################################################

# Transpose OTU table if needed
otu_mat <- as(otu_table(physeq), "matrix")
if (taxa_are_rows(physeq)) {
  otu_mat <- t(otu_mat)
}

group_vector <- sample_data(physeq)$Primer

# Run indicator species analysis
indval <- multipatt(otu_mat, group_vector, func = "IndVal.g", control = how(nperm=999))

# Get significant OTUs (p <= 0.05)
indval_df <- as.data.frame(indval$sign)
indval_df$OTU <- rownames(indval_df)
indval_sig <- indval_df %>% filter(p.value <= 0.05)

# Extract primer labels (requires correction of the index column manually if needed)
group_names <- colnames(indval$str)
indval_sig <- indval_sig %>%
  mutate(
    Primer = group_names[index] 
  ) %>%
  select(OTU, Primer, stat, p.value)

# Merge with taxonomy (class and phylum)
tax_subset <- as.data.frame(tax_mat[indval_sig$OTU, c("Class", "Phylum"), drop = FALSE])
bubble_df <- indval_sig %>%
  left_join(tibble(OTU = rownames(tax_subset), Class = tax_subset$Class, Phylum = tax_subset$Phylum), by = "OTU")

# Get relative abundance per OTU
otu_table <- read.delim("results/OTUs.97.rep.count_table", row.names = 1, check.names = FALSE)
otu_table <- otu_table[, !colnames(otu_table) %in% "total"]
otu_rel <- sweep(otu_table, 2, colSums(otu_table), FUN = "/")
otu_mean_rel_abund <- rowMeans(otu_rel) * 100
rel_abund_df <- tibble(OTU = rownames(otu_table), rel_abund = otu_mean_rel_abund)

# Combine abundance data
bubble_df <- bubble_df %>%
  left_join(rel_abund_df, by = "OTU") %>%
  filter(!is.na(Class), !is.na(Phylum))

# Order factor levels
bubble_df$Class <- fct_rev(fct_infreq(bubble_df$Class))
bubble_df$Phylum <- fct_infreq(bubble_df$Phylum)

# Labels for primer plots
primer_labels <- c("V3V4" = "A", "V4V5" = "B")
plots <- list()

for (primer in names(primer_labels)) {
  df_primer <- bubble_df %>% filter(Primer == primer)
  if (nrow(df_primer) == 0) next
  
  label <- primer_labels[primer]
  
  p <- ggplot(df_primer, aes(x = stat, y = Class, size = rel_abund, color = Phylum)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(
      title = label,
      x = "Correlation Coefficient",
      y = "Class",
      size = "% Relative Abundance",
      color = "Phylum"
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 16, hjust = -0.2)
    )
  
  plots[[primer]] <- p
}

# Combine both panels vertically
figura6 <- plot_grid(plots[["V3V4"]], plots[["V4V5"]], nrow = 2, align = "v")

# Save output
ggsave(file.path(fig_dir,"Figura6_Indicadoras.pdf"), plot = figura6, width = 10, height = 14)
ggsave(file.path(fig_dir,"Figura6_Indicadoras.png"), plot = figura6, width = 10, height = 14, dpi = 300)

################################################################################
################# FIGURE 7: ZETAPROTEOBACTERIA zOTUs ANALYSIS  #################
################################################################################

# Align OTU and taxonomy matrices
common_otus <- intersect(rownames(otu_mat), rownames(tax_mat))
otu_mat_aligned <- otu_mat[common_otus, ]
tax_mat_aligned <- tax_mat[common_otus, ]

stopifnot(all(rownames(otu_mat_aligned) == rownames(tax_mat_aligned)))

# Recreate phyloseq object with filtered data
physeq <- phyloseq(
  otu_table(otu_mat_aligned, taxa_are_rows = TRUE),
  tax_table(tax_mat_aligned),
  sample_data(metadata)
)

# Filter for Zetaproteobacteria class
physeq_zeta <- subset_taxa(physeq, Class == "Zetaproteobacteria")

# Melt phyloseq object into long format
zeta_df <- psmelt(physeq_zeta)

# Create barplot
p7 <- ggplot(zeta_df, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Primer, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Figure 7 - Zetaproteobacteria zOTUs Diversity")

# Save output
ggsave(file.path(fig_dir,"Figura7_Zetaproteobacteria.pdf"), p7, width = 10, height = 6)

