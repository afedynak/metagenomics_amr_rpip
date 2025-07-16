### AMR/RPIP Data Visualization Script
# Author: Anastasia Amber Fedynak
# Description: This R script reads RPIP AMR summary data and creates a series of plots for antimicrobial resistance (AMR) genes,
# microbial community profiles, and variant analysis in lake samples from Ontario. Each block is commented with purpose, steps, and outputs.

# --- Load Required Libraries ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(gridExtra)
library(RColorBrewer)
library(forcats)
library(ggbeeswarm)

# --- Set Working Directory and Load Data ---
setwd("/Users/admin/Desktop/AMR_RPIP/")
amr_data <- read_excel("Lake_samples_summary_RPIP_AMR.xlsx", sheet = "AMR (RPIP)")

# --- Filter for high-confidence genes and select important columns ---
high_conf_amr <- amr_data %>%
  filter(Confidence_Interpretation == "high") %>%
  select(Accession, Best_Match, AMR_Name, NT_RPKM, AA_RPKM)

# --- Plot 1: Gene vs Abundance (AA_RPKM) ---
ggplot(high_conf_amr, aes(x = Best_Match, y = AA_RPKM, color = AMR_Name)) +
  geom_point(alpha = 0.7) +
  labs(title = "Gene Identified vs. AA Abundance", x = "Gene Identified", y = "AA_RPKM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Plot 2: Gene Count Distribution ---
ggplot(high_conf_amr, aes(x = AMR_Name, fill = AMR_Name)) +
  geom_bar() +
  labs(title = "Distribution of AMR Genes", x = "AMR Gene", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- AMR Genes by Site ---
site_gene_counts <- amr_data %>%
  group_by(Accession, AMR_Name) %>%
  tally(name = "Gene_Count") %>%
  mutate(Accession = recode(Accession, "Wilcox_Lake" = "Wilcox Lake", "Guelph_Lake" = "Guelph Lake"))

# --- Plot 3: Stacked Barplot of AMR Gene Count by Site ---
jpeg("gene_count_per_site.jpg", height = 5, width = 7, units = "in", res = 300)

ggplot(site_gene_counts, aes(x = Accession, y = Gene_Count, fill = AMR_Name)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Gene Count", fill = "AMR Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

dev.off()

# --- Quality Plots: Coverage, PID, Depth for Each Site ---
# Function to generate dot/bar plots by gene feature
guelph_data <- selected_data %>%
  filter(Accession == "Guelph Lake")

wilcox_data <- selected_data %>%
  filter(Accession == "Wilcox Lake")

wilcox_color <- "#8DD3C7"   # greenish
guelph_color <- "#FB8072"   # reddish


# Create the bar plot
# ------------------- Guelph Lake: Coverage -------------------
pdf("amr_coverage_guelph.pdf", height=5)
ggplot(guelph_data, aes(x = AMR_Name, y = AA_Coverage, fill = Accession)) +
  geom_bar(stat = "identity") +
  labs(
    x = "AMR Gene",
    y = "Coverage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Guelph Lake" = "#1f77b4"))
dev.off()

# ------------------- Wilcox: Coverage -------------------
pdf("amr_coverage_wilcox.pdf", height=5)
ggplot(wilcox_data, aes(x = AMR_Name, y = AA_Coverage, fill = Accession)) +
  geom_bar(stat = "identity") +
  labs(
    x = "AMR Gene",
    y = "Coverage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Wilcox Lake" = "#1f77b4"))
dev.off()

# ------------------- Guelph: AA_PID -------------------
pdf("amr_pid_guelph.pdf", height=5)
ggplot(guelph_data, aes(x = AMR_Name, y = AA_PID, fill = Accession)) +
  geom_bar(stat = "identity") +
  labs(x = "AMR Gene", y = "Percent Identity (AA PID)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("Guelph Lake" = guelph_color))
dev.off()

# ------------------- Wilcox: AA_PID -------------------
pdf("amr_pid_wilcox.pdf", height=5)
ggplot(wilcox_data, aes(x = AMR_Name, y = AA_PID, fill = Accession)) +
  geom_bar(stat = "identity") +
  labs(x = "AMR Gene", y = "Percent Identity (AA PID)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("Wilcox Lake" = wilcox_color))
dev.off()

# ------------------- Guelph: AA_Median_Depth -------------------
pdf("amr_depth_guelph.pdf", height=5)
ggplot(guelph_data, aes(x = AMR_Name, y = AA_Median_Depth, fill = Accession)) +
  geom_bar(stat = "identity") +
  labs(x = "AMR Gene", y = "Median Depth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("Guelph Lake" = guelph_color))
dev.off()

# ------------------- Wilcox: AA_Median_Depth -------------------
pdf("amr_depth_wilcox.pdf", height=5)
ggplot(wilcox_data, aes(x = AMR_Name, y = AA_Median_Depth, fill = Accession)) +
  geom_bar(stat = "identity") +
  labs(x = "AMR Gene", y = "Median Depth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("Wilcox Lake" = wilcox_color))
dev.off()

# Prepare data
coverage_data <- selected_data %>%
  select(Accession, AMR_Name, AA_Coverage) %>%
  filter(!is.na(AA_Coverage))

# Plot
pdf("amr_coverage_violin.pdf", height = 5)
ggplot(coverage_data, aes(x = Accession, y = AA_Coverage, fill = Accession)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  scale_fill_manual(values = c("Guelph Lake" = guelph_color, "Wilcox Lake" = wilcox_color)) +
  labs(x = "Site", y = "AA Coverage") +
  theme_minimal()
dev.off()

pid_data <- selected_data %>%
  select(Accession, AMR_Name, AA_PID) %>%
  filter(!is.na(AA_PID))

# Plot
pdf("amr_pid_boxplot.pdf", height = 5)
ggplot(pid_data, aes(x = Accession, y = AA_PID, fill = Accession)) +
  geom_boxplot(alpha = 0.9, outlier.color = "black") +
  scale_fill_manual(values = c("Guelph Lake" = guelph_color, "Wilcox Lake" = wilcox_color)) +
  labs(x = "Site", y = "Percent Identity (PID)") +
  theme_minimal()
dev.off()

# Prepare data
depth_data <- selected_data %>%
  select(Accession, AMR_Name, AA_Median_Depth) %>%
  filter(!is.na(AA_Median_Depth))

# Plot for Wilcox Lake
wilcox_data <- depth_data %>% filter(Accession == "Wilcox Lake")

pdf("amr_median_depth_wilcox.pdf", height = 5)
ggplot(wilcox_data, aes(x = reorder(AMR_Name, -AA_Median_Depth), y = AA_Median_Depth)) +
  geom_segment(aes(xend = AMR_Name, y = 0, yend = AA_Median_Depth), color = wilcox_color) +
  geom_point(color = wilcox_color, size = 3) +
  coord_flip() +
  labs(x = "AMR Gene", y = "Median Depth") +
  theme_minimal()
dev.off()

# Plot for Guelph Lake
guelph_data <- depth_data %>% filter(Accession == "Guelph Lake")

pdf("amr_median_depth_guelph.pdf", height = 5)
ggplot(guelph_data, aes(x = reorder(AMR_Name, -AA_Median_Depth), y = AA_Median_Depth)) +
  geom_segment(aes(xend = AMR_Name, y = 0, yend = AA_Median_Depth), color = guelph_color) +
  geom_point(color = guelph_color, size = 3) +
  coord_flip() +
  labs(x = "AMR Gene", y = "Median Depth") +
  theme_minimal()
dev.off()


# --- Microorganism Data Analysis ---
microbe_data <- read_excel("Lake_samples_summary_RPIP_AMR.xlsx", sheet = "Microorganisms (RPIP)")

clean_microbes <- microbe_data %>%
  filter(`Passed Cutoffs` == "y") %>%
  mutate(Accession = recode(Accession, "Wilcox_Lake" = "Wilcox Lake", "Guelph_Lake" = "Guelph Lake"),
         Microorganism = fct_recode(`Microorganism Name`, "SARS-CoV-2" = "SARS-CoV-2 (2019-nCoV)"),
         Abundance = as.numeric(Abundance),
         Coverage = as.numeric(Coverage),
         ANI = as.numeric(ANI))

# Plot: Abundance per Microorganism by Site
# Guelph Lake
jpeg("guelph_microorg.jpg", height = 5, width = 7, units = "in", res = 300)

ggplot(filter(clean_microbes, Accession == "Guelph Lake"),
       aes(x = reorder(Microorganism, -Abundance), y = Abundance, fill = `Class Type`)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(x = "Microorganism", y = "Abundance")

dev.off()

# Wilcox Lake
jpeg("wilcox_microorg.jpg", height = 5, width = 7, units = "in", res = 300)

ggplot(filter(clean_microbes, Accession == "Wilcox Lake"),
       aes(x = reorder(Microorganism, -Abundance), y = Abundance, fill = `Class Type`)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Microorganism", y = "Abundance")

dev.off()

# --- SARS-CoV-2 Variant Analysis ---
variant_data <- read_excel("Lake_samples_summary_RPIP_AMR.xlsx", sheet = "Variants (RPIP)") %>%
  rename_with(~str_replace_all(., "\\s+", "_")) %>%
  mutate(Accession = gsub("_", " ", Accession),
         Variant_Type = case_when(
           str_detect(NT_Change, "del") ~ "Deletion",
           str_detect(NT_Change, "ins") ~ "Insertion",
           TRUE ~ "Substitution"
         ))

# Plot: Variant Type per Site
jpeg("variant_type_custom_colors.jpg", height = 5, width = 7, units = "in", res = 300)

ggplot(variant_data, aes(x = Variant_Type, y = Allele_Frequency, color = Variant_Type)) +
  geom_beeswarm(size = 2) +
  facet_wrap(~ Accession) +
  labs(x = "Variant Type", y = "Allele Frequency") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(values = c("#8DD3C7", "#FB8072", "#BC80BD"))  # Green, redish, and purple

dev.off()

# Plot: Variant Annotation Count per Site
variant_annotation <- variant_data %>%
  group_by(Accession, Annotation) %>%
  summarise(Count = n(), .groups = "drop")

jpeg("variant_annotation.jpg", height = 5, width = 7, units = "in", res = 300)
ggplot(variant_annotation, aes(x = Accession, y = Count, fill = Annotation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")
dev.off()

