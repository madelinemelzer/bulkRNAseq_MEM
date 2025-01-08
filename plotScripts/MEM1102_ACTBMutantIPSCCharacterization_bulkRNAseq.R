#Plotting DESeq results
# created on 20250107 by MEM
# last modified on 20250107 by MEM

set.seed(23)
options(future.globals.maxSize = 32 * 1024 * 1024^2)  # 32 GB

library(ggplot2)
library(ggrepel)
library(svglite)

output <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/MEM1102_ACTBMutantIPSCCharacterization/bulkRNAseq/extData/deseq/"
plotDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/MEM1102_ACTBMutantIPSCCharacterization/bulkRNAseq/plots/"
plotDataDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/MEM1102_ACTBMutantIPSCCharacterization/bulkRNAseq/plotData/"

wt_het = read.csv(paste0(output, "DESeqResult_wt_het_withGeneNames.csv"))
wt_homo = read.csv(paste0(output, "DESeqResult_wt_homo_withGeneNames.csv"))

wt_het_subset <- wt_het[order(wt_het$padj)[1:1000], ]  # Top 1000 most significant genes to reduce plot complexity
wt_het_plot = ggplot(wt_het_subset, aes(x = log2FoldChange, 
                                        y = -log10(padj), 
                                        color = ifelse(padj < 0.05 & log2FoldChange > 1, "Upregulated",
                                                       ifelse(padj < 0.05 & log2FoldChange < -1, "Downregulated", "Not Significant")))) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray"),
                     name = "Regulation") +
  
  # Add vertical lines for fold change cutoff
  geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", 
             color = "darkgray") +
  
  # Add horizontal line for significance
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "darkgray") +
  
  # Add labels for significant genes
  geom_text_repel(data = subset(wt_het_subset, padj < 0.05 & abs(log2FoldChange) > 1),
                  aes(label = external_gene_name),
                  max.overlaps = 10,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50") +
  
  theme_classic() +
  labs(
    x = "log2FoldChange",
    y = "-log10(adj. p-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )
wt_het_plot

ggsave(wt_het_plot, file = paste0(plotDir, "wt_het_volcano.svg"), width = 6, height = 6) #20250107
ggsave(wt_het_plot, file = paste0(plotDir, "wt_het_volcano.tif"), width = 6, height = 6) #20250107


wt_homo_subset <- wt_homo[order(wt_homo$padj)[1:1000], ]  # Top 1000 most significant genes to reduce plot complexity
wt_homo_plot = ggplot(wt_homo_subset, aes(x = log2FoldChange, 
                                        y = -log10(padj), 
                                        color = ifelse(padj < 0.05 & log2FoldChange > 1, "Upregulated",
                                                       ifelse(padj < 0.05 & log2FoldChange < -1, "Downregulated", "Not Significant")))) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray"),
                     name = "Regulation") +
  
  # Add vertical lines for fold change cutoff
  geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", 
             color = "darkgray") +
  
  # Add horizontal line for significance
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "darkgray") +
  
  # Add labels for significant genes
  geom_text_repel(data = subset(wt_homo_subset, padj < 0.05 & abs(log2FoldChange) > 1),
                  aes(label = external_gene_name),
                  max.overlaps = 15,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50") +
  
  theme_classic() +
  labs(
    x = "log2FoldChange",
    y = "-log10(adj. p-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )
wt_homo_plot

ggsave(wt_homo_plot, file = paste0(plotDir, "wt_homo_volcano.svg"), width = 6, height = 6) #20250107
ggsave(wt_homo_plot, file = paste0(plotDir, "wt_homo_volcano.tif"), width = 6, height = 6) #20250107

