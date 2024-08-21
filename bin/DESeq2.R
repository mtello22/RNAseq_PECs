library(DESeq2)
library(data.table)

load(file = "~/GitHub/RNAseq_PECs/data/txi_coldat_object.rds")

DESeqDataset <- DESeqDataSetFromTximport(txi = txi, colData = deseq_coldat, design = ~ condition)

# Filter lowly expressed genes
min_reads <- 1
min_samples <- 12
index <- rowSums(counts(DESeqDataset) >= min_reads) > min_samples
DESeqDataset <- DESeqDataset[index, ]

# Perform DESeq2
DESeqDataset <- DESeq(DESeqDataset)

# resultsNames(DESeqDataset)
# "condition_GEN9_vs_Ctrl" "condition_H2O2_vs_Ctrl" 
# "condition_PNS2_vs_Ctrl" "condition_PTS3_vs_Ctrl"

DEG_results <- results(DESeqDataset, 
                        name="condition_H2O2_vs_Ctrl",
                        pAdjustMethod = "BH", 
                        alpha = 0.05)

DEG_results <- as.data.table(keep.rownames = TRUE, DEG_results)
setnames(DEG_results, "rn", "ENSG")
DEG_results[, Group := "H2O2_Vs_Ctrl"]

####################

deseq_coldat$condition <- relevel(deseq_coldat$condition, ref = "H2O2")
DESeqDataset <- DESeqDataSetFromTximport(txi = txi, colData = deseq_coldat, design = ~ condition)
# Filter lowly expressed genes
DESeqDataset <- DESeqDataset[index, ]

# Perform DESeq2
DESeqDataset <- DESeq(DESeqDataset)

# resultsNames(DESeqDataset)
# "condition_Ctrl_vs_H2O2" "condition_GEN9_vs_H2O2"
# "condition_PNS2_vs_H2O2" "condition_PTS3_vs_H2O2"
DEGs <- results(DESeqDataset, 
                name="condition_GEN9_vs_H2O2",
                pAdjustMethod = "BH", 
                alpha = 0.05)
DEGs <- as.data.table(DEGs, keep.rownames = TRUE)
setnames(DEGs, "rn", "ENSG")
DEGs[, Group := "GEN9_vs_H2O2"]
DEG_results <- rbind(DEG_results, DEGs)

DEGs <- results(DESeqDataset, 
                name="condition_PNS2_vs_H2O2",
                pAdjustMethod = "BH", 
                alpha = 0.05)
DEGs <- as.data.table(DEGs, keep.rownames = TRUE)
setnames(DEGs, "rn", "ENSG")
DEGs[, Group := "PNS2_vs_H2O2"]
DEG_results <- rbind(DEG_results, DEGs)

DEGs <- results(DESeqDataset, 
                name="condition_PTS3_vs_H2O2",
                pAdjustMethod = "BH", 
                alpha = 0.05)
DEGs <- as.data.table(DEGs, keep.rownames = TRUE)
setnames(DEGs, "rn", "ENSG")
DEGs[, Group := "PTS3_vs_H2O2"]
DEG_results <- rbind(DEG_results, DEGs)

fwrite(DEG_results, file = "~/GitHub/RNAseq_PECs/data/DEG_results.tsv", 
       quote = FALSE, append = FALSE, sep = '\t', 
       row.names = FALSE, col.names = TRUE)
