load(file = "~/GitHub/RNAseq_PECs/data/txi_coldat_object.rds")

DESeqDataset <- DESeqDataSetFromTximport(txi = txi, colData = deseq_coldat, design = ~ condition)
DESeqDataset <- DESeq(DESeqDataset)

# resultsNames(DESeqDataset)
# "condition_GEN9_vs_Ctrl" "condition_H2O2_vs_Ctrl" 
# "condition_PNS2_vs_Ctrl" "condition_PTS3_vs_Ctrl"
results <- results(DESeqDataset, 
                   name="condition_GEN9_vs_Ctrl",
                   pAdjustMethod = "BH", 
                   alpha = 0.05)

