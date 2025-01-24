#Code to analyze Count Input Matrix using DESeq2 - created by Aurelia on 2024-08-09 and is based on Keerthana's code 
# last modified on 20250123 by MEM
library(Seurat)
library(dplyr)
library(dbplyr)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(readr)


#Main ---------------------------------
#extData <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/MEM1102_ACTBMutantIPSCCharacterization/bulkRNAseq/extData/matrix/"
Scripts<- "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/extractionScripts/bulkRNAseq/DESeqbulk/"
output <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/extractedData/MEM1201_Grimes_Urp/bulkRNAseq_20250123/deseq/"

source(paste0(Scripts, "deseqbulkfunction.R"))

#"GeneID"
mastercounts = as_tibble(read_table(file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/extractedData/MEM1201_Grimes_Urp/bulkRNAseq_20250123/mastercounts3.txt"))

#to get "LLgeneID" column
geneInformation = as_tibble(read_csv(file = "/Users/mem3579/p32655/melzer/MEM1201_Grimes_Urp/refData/v4_3_2geneinfo.csv"))

#adding gene information to mastercounts.
mastercounts_meta <- mastercounts %>%
  left_join(geneInformation, by = c("GeneID" = "LLgeneID"))

#write.csv(mastercounts_meta, file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/transcriptionalAdaptation/MEM1201_Grimes_Urp/data/mastercounts3_meta.csv")


countmatrix = mastercounts_meta %>% select(-Length, -LLchr, -LLstart, -LLend, -LLstrand, -geneName, -Ens99geneIDversion, -EntrezGeneID, -ZFINgeneID, -RefSeqAnnotationNote, -Ens99annotationNote, -'LLgeneSymbol (V4.3.2.gtf)')

#Create directory 
create_dir_if_not_exists(output)

#countmatrix <- read_csv(paste0(extData, "countmatrix.csv"))

#Look at what the column names are
print(colnames(countmatrix))

#  [1] "ensembl_gene"             "ID_HS_017_B1_sot_S84"     "ID_HS_017_B2_sot_S52"     "ID_HS_017_B3_sot_S59"    
# [5] "ID_HS_017_gem1_S51"       "ID_HS_017_gem2_S83"       "ID_HS_017_p2a2_S85"       "ID_HS_025_tram1_S60"     
# [9] "ID_HS_025_tram2_S69"      "ID_HS_025_tram3_S77"      "ID_HS_026_20g_20s1_S62"   "ID_HS_026_20g_20s2_S93"  
# [13] "ID_HS_026_20g_50s_S80"    "ID_HS_026_50g_20s1_S48"   "ID_HS_026_50g_20s2_S49"   "ID_HS_026_50g_50s1_S67"  
# [17] "ID_HS_026_50g_50s2_S76"   "ID_HS_027_12daysS1_S65"   "ID_HS_027_12daysS2_S86"   "ID_HS_027_12daysS3_S96"  
# [21] "ID_HS_027_48hCis1_S56"    "ID_HS_027_48hCis2_S61"    "ID_HS_027_48hCis3_S88"    "ID_HS_027_48hGem1_S66"   
# [25] "ID_HS_027_48hGem2_S73"    "ID_HS_027_48hGem3_S75"    "ID_HS_027_48hSot1_S64"    "ID_HS_027_48hSot2_S78"   
# [29] "ID_HS_027_48hSot3_S95"    "ID_HS_027_8daysG1_S71"    "ID_HS_027_8daysG2_S92"    "ID_HS_027_8daysG3_S54"   
# [33] "ID_HS_027_cis1_S58"       "ID_HS_027_cis2_S68"       "ID_HS_027_cis3_S90"       "ID_SA2-100_20g_20s_S79"  
# [37] "ID_SA2-100_20g_50s_S47"   "ID_SA2-100_50g_50s_S82"   "ID_SA2-101_50g_20s_S81"   "ID_SA2-102_A1_sotGem_S91"
# [41] "ID_SA2-102_A2_sotGem_S70" "ID_SA2-102_C1_sotR_S45"   "ID_SA2-102_C2_sotR_S57"   "ID_SA2-102_C3_sotR_S89"  
# [45] "ID_SA2-102_P_S53"         "ID_SA2-103_A2_sotGem_S63" "ID_SA2-103_P1_S74"        "ID_SA2-104_A1_gemSot_S94"
# [49] "ID_SA2-104_A2_gemSot_S72" "ID_SA2-104_B2_gem_S50"    "ID_SA2-104_C1_gemR_S55"   "ID_SA2-104_C2_gemR_S87"  
# [53] "ID_SA2-104_C3_gemR_S46"  

# Make ensembl_gene as the rownames for later processing 
countmatrixdf <- column_to_rownames(as.data.frame(countmatrix),'GeneID')

# Select the columns you would like to have in your DESeq 
# if you want to delete the columns use -c("colname1", "colname2")
# If you want to include the columns use c("colname1", "colname2")
# (madeline did not use because i want all?) countmatrixdf<- dplyr::select(countmatrixdf,c("ID_HS_017_B1_sot_S84", "ID_HS_017_B2_sot_S52", "ID_HS_017_B3_sot_S59", "ID_HS_017_p2a2_S85","ID_SA2-102_P_S53", "ID_SA2-103_P1_S74"))
# countmatrixdf<- dplyr::select(countmatrixdf,-c("Cisplatin_1_S7", "Cisplatin_2_S8", "Gemcitabine_3_S14", "Gemcitabine_1_S12"))
# cts <- as.matrix.data.frame(countmatrixdf)

#Make an annotation matrix 
# First, let's create the annotation matrix
sample_names <- colnames(countmatrixdf)[colnames(countmatrixdf) != "ensembl_gene"]
# Specify the condition and the replicate 
condition_maps <- list(
  "abc_1_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 1),
  "abc_2_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 2),
  "abc_3_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 3),
  "abc_4_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 4),
  "ABC_5_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 5),
  "ABC_6_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 6),
  "ABC_7_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 7),
  "ABC_8_F.fastqsanger.gz.subread.BAM" = list(condition = "control", replicate = 8),
  
  "URP_1_F.fastqsanger.gz.subread.BAM" = list(condition = "urp1urp2", replicate = 1),
  "URP_2_F.fastqsanger.gz.subread.BAM" = list(condition = "urp1urp2", replicate = 2),
  "URP_3_F.fastqsanger.gz.subread.BAM" = list(condition = "urp1urp2", replicate = 3),
  "URP_4_F.fastqsanger.gz.subread.BAM" = list(condition = "urp1urp2", replicate = 4),
  "URP_5_F.fastqsanger.gz.subread.BAM" = list(condition = "urp1urp2", replicate = 5),
  "URP_6_F.fastqsanger.gz.subread.BAM" = list(condition = "urp1urp2", replicate = 6),
  
  "urp1_1_S12_R1_001.fastq.gz.subread.BAM" = list(condition = "urp1", replicate = 1),
  "urp1_2_S13_R1_001.fastq.gz.subread.BAM" = list(condition = "urp1", replicate = 2),
  "urp1_3_S14_R1_001.fastq.gz.subread.BAM" = list(condition = "urp1", replicate = 3),
  "urp1_4_S15_R1_001.fastq.gz.subread.BAM" = list(condition = "urp1", replicate = 4),
  "urp1_6_S16_R1_001.fastq.gz.subread.BAM" = list(condition = "urp1", replicate = 6),
  
  "urp2_3_S19_R1_001.fastq.gz.subread.BAM" = list(condition = "urp2", replicate = 3),
  "urp2_4_S20_R1_001.fastq.gz.subread.BAM" = list(condition = "urp2", replicate = 4),
  "urp2_5_S21_R1_001.fastq.gz.subread.BAM" = list(condition = "urp2", replicate = 5),
  "urp2_a_S17_R1_001.fastq.gz.subread.BAM" = list(condition = "urp2", replicate = NA),
  "urp2_b_S18_R1_001.fastq.gz.subread.BAM" = list(condition = "urp2", replicate = NA),
  
  "uts_1_F.fastqsanger.gz.subread.BAM" = list(condition = "uts", replicate = 1),
  "uts_2_F.fastqsanger.gz.subread.BAM" = list(condition = "uts", replicate = 2),
  "uts_3_F.fastqsanger.gz.subread.BAM" = list(condition = "uts", replicate = 3),
  "uts_4_F.fastqsanger.gz.subread.BAM" = list(condition = "uts", replicate = 4),
  "uts_5_F.fastqsanger.gz.subread.BAM" = list(condition = "uts", replicate = 5),
  "uts_6_F.fastqsanger.gz.subread.BAM" = list(condition = "uts", replicate = 6)
)
# Use the function to create annotation matrix 
annotation_data <- create_annotation_matrix(sample_names, condition_mapping=condition_maps)

# View the first few rows of the annotation matrix
head(annotation_data)

# Create a factor level so the comparison will be Naive vs Sotorasib 
# annotate_data <- column_to_rownames(as.data.frame(annotate_data),'sample') 
# annotate_data <- annotate_data %>% filter(!(condition == "Gemcitabine" & replicate == 2))
# annotation_data$condition <- factor(annotation_data$condition, levels = c("naive", "Gemcitabine", "Sotorasib", "Cisplatin"))
annotation_data$condition <- factor(annotation_data$condition, levels = c("control", "urp1", "urp2", "urp1urp2", "uts"))
write_csv(annotation_data, paste0(output,"annotationmatrix.csv"))



# Start deseq 
dds <- DESeqDataSetFromMatrix(countData = countmatrixdf,
                              colData = annotation_data,
                              design = ~ condition)


# Calculate the low count threshold based on a quantile (added 20250121)
low_count_threshold <- quantile(rowMeans(counts(dds)), 0.25) # 4.233


#Pre-filtering process-----------------------
smallestGroupSize <- 1
count_threshold <- 4

# Get the count matrix once
# Use base R's rowSums, which is still quite efficient
keep <- rowSums(counts(dds) >= count_threshold) >= smallestGroupSize

# Use logical indexing
dds <- dds[keep,]

##Differential Expression Analysis--------------------------------
# dds_naive <- dds
# dds_naive$condition <- relevel(x = dds_naive$condition, ref = "naive")
# 
# dds_naive <- DESeq(dds_naive)
# res_naive_vs_sotorasib <- results(dds_naive, contrast=c("condition","naive","Sotorasib"))
# res_naive_vs_sotarasib <- data.frame(res_naive_vs_sotorasib)
# res_naive_vs_sotarasib <- res_naive_vs_sotarasib %>% filter(padj < 0.05)

#Naive samples
#contrast_info <- c("condition", "homo", "wt")
padj_threshold <- 0.05

# Relevel the condition factor to indicate that the reference is naive
dds$condition <- relevel(dds$condition, ref = "control")

# Perform DESeq analysis
dds <- DESeq(dds)

# Get results for the specified comparisons
control_vs_urp1 <- results(dds, contrast = c("condition", "urp1", "control"))
control_vs_urp2 <- results(dds, contrast = c("condition", "urp2", "control"))
control_vs_urp1urp2 <- results(dds, contrast = c("condition", "urp1urp2", "control"))
control_vs_uts <- results(dds, contrast = c("condition", "uts", "control"))

# Save RDS files
saveRDS(control_vs_urp1, paste0(output, "DESeqResult__control_urp1.rds")) #20250121
saveRDS(control_vs_urp2, paste0(output, "DESeqResult__control_urp2.rds")) #20250121
saveRDS(control_vs_urp1urp2, paste0(output, "DESeqResult__control_urp1urp2.rds")) #20250121
saveRDS(control_vs_uts, paste0(output, "DESeqResult__control_uts.rds")) #20250121

# Save CSV files
write.csv(control_vs_urp1, paste0(output, "DESeqResult_control_urp1.csv"), row.names = TRUE) #20250123
write.csv(control_vs_urp2, paste0(output, "DESeqResult_control_urp2.csv"), row.names = TRUE) #20250123
write.csv(control_vs_urp1urp2, paste0(output, "DESeqResult_control_urp1urp2.csv"), row.names = TRUE) #20250123
write.csv(control_vs_uts, paste0(output, "DESeqResult_control_uts.csv"), row.names = TRUE) #20250123




#Adding back metadata --------------------------------

# Re-open results for adding back metadata
control_vs_urp1 <- read.csv(paste0(output, "DESeqResult_control_urp1.csv"))
control_vs_urp2 <- read.csv(paste0(output, "DESeqResult_control_urp2.csv"))
control_vs_urp1urp2 <- read.csv(paste0(output, "DESeqResult_control_urp1urp2.csv"))
control_vs_uts <- read.csv(paste0(output, "DESeqResult_control_uts.csv"))

add_metadata_to_deseq_results <- function(deseq_result,comparison_name, output_path) {
  
  # Add metadata
  result_with_meta <- deseq_result %>%
    left_join(
      geneInformation %>% 
        select(LLgeneID, 
               `LLgeneSymbol (V4.3.2.gtf)`,
               LLchr, 
               LLstart, 
               LLend, 
               LLstrand, 
               geneName, 
               Ens99geneIDversion, 
               EntrezGeneID, 
               ZFINgeneID, 
               RefSeqAnnotationNote, 
               Ens99annotationNote), 
      by = c("X" = "LLgeneID"))
  
  # Write CSV
  write.csv(result_with_meta, 
            paste0(output_path, "DESeqResult_", comparison_name, "_withGeneNamesAndMetadata.csv"), 
            row.names = TRUE)
  
  # Return the result in case you want to do further processing
  return(result_with_meta)
}

control_vs_urp1_meta <- add_metadata_to_deseq_results(
  deseq_result = control_vs_urp1, 
  comparison_name = "control_vs_urp1", 
  output_path = output
)

control_vs_urp2_meta <- add_metadata_to_deseq_results(
  deseq_result = control_vs_urp2, 
  comparison_name = "control_vs_urp2", 
  output_path = output
)

control_vs_urp1urp2_meta <- add_metadata_to_deseq_results(
  deseq_result = control_vs_urp1urp2, 
  comparison_name = "control_vs_urp1urp2", 
  output_path = output
)

control_vs_uts_meta <- add_metadata_to_deseq_results(
  deseq_result = control_vs_uts, 
  comparison_name = "control_vs_uts", 
  output_path = output
)






















