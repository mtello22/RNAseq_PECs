EDA
================
Alice
2024-02-22

# Exploratory Data Analysis

Once we have the final count table, we can analyze the distribution of
counts per sample and possible variables that could influence our
results.

## Count distribution

We looked at the density distribution of gene expression values for all
samples. This will help us to identify outlier samples in the gene
expression data.

``` r
# Read count table
# We remove last sample which is not properly identified
count_df <- fread(file.path(data_path, "Salmon_EstCount_ENSG.tsv") )


count_df[1:5,1:5]
```

    ##                  ENSG  Ctrl_1  Ctrl_2  Ctrl_3  GEN9_1
    ##                <char>   <num>   <num>   <num>   <num>
    ## 1: ENSG00000000003.16 286.513  90.259  72.215 104.211
    ## 2:  ENSG00000000005.6   0.000   0.000   0.000   0.000
    ## 3: ENSG00000000419.14 640.854 510.415 587.028 504.972
    ## 4: ENSG00000000457.14 259.760 142.849  94.973 164.791
    ## 5: ENSG00000000460.17 153.036  39.680  30.645  56.586

``` r
count_df %>%
  pivot_longer(!ENSG) %>%
  ggplot(aes(x=value, color = name))+
    geom_density() +
    scale_x_continuous(trans = 'log10', labels = scales::comma)
```

    ## Warning in scale_x_continuous(trans = "log10", labels = scales::comma): log-10
    ## transformation introduced infinite values.

    ## Warning: Removed 254452 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](Alice-Template_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
temp <- count_df %>%
  pivot_longer(!ENSG) %>%
  mutate(condition = str_replace(string = name, 
                                 pattern = "_\\d", 
                                 replacement = "")) %>%
  mutate(condition = case_when(
    condition == "Ctrl" ~ "blue",
    condition == "GEN9" ~ "red",
    condition == "H2O2" ~ "darkgreen",
    condition == "PNS2" ~ "purple",
    condition == "PTS3" ~ "pink")) 

temp %>%
  ggplot(aes(x=value, color = name))+
    geom_density() +
    scale_x_continuous(trans = 'log10',labels = scales::comma)+
  scale_color_manual(values = temp$condition,
                     breaks = temp$name)
```

    ## Warning in scale_x_continuous(trans = "log10", labels = scales::comma): log-10
    ## transformation introduced infinite values.

    ## Warning: Removed 254452 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](Alice-Template_files/figure-gfm/Gene%20expression%20density%20plots-1.png)<!-- -->


    We initially observed two sample in Ctrl and GEN9 that do not completely follow the distribution of counts as the others. Thus we flagged these 2 samples as outlier with a different distribution of counts compared to other samples. The next step is to remove lowly expressed genes and check if the density plots improve.

    ### Filter lowly expressed genes

    To increase data quality we removed lowly expressed genes in two steps:

    1.  Removed genes with no counts across all samples.

    2.  Removed genes with less than two reads in more than half of the samples.

    Following this strategy we kept 14,450 protein coding genes from the initial 20,124. We used a threshold of minimum 2 reads in more than half of the samples for keeping any gene. This value was determined by the average of the first quartile from all samples.


    ``` r
    # We first remove duplicated gene symbols
    exp_mat <- count_df %>%
      column_to_rownames("ENSG")


    # Remove rows with all entries equal zero
    exp_mat <- exp_mat[!(rowSums(exp_mat == 0) == ncol(exp_mat)),]

    # We use the average 1st quartile as threshold (2.1), round down to 2
    # floor(mean(apply(exp_mat, MARGIN = 2, quantile, 0.25)))

    # Select genes to keep
    min_reads <- 5
    min_samples <- 8
    genes_to_keep <- apply(exp_mat >= min_reads,
                           MARGIN = 1, sum) > min_samples

    # Final gene count matrix
    exp_mat <- exp_mat[genes_to_keep,]


    This step made most samples distributions to look mostly similar.

    ## CPM normalization

    Before proceeding to further analysis we need to normalize the gene counts, we will use *counts per million* (CPM) to adjust for library size and transform it to log space for better visualization.


    ``` r
    # CPM normalization and transform to log2

    expr_log2cpm <- cpm(exp_mat, 
                        log = TRUE, 
                        prior.count = 1) %>% 
      data.frame() 

    expr_log2cpm %>%
      rownames_to_column(var = "ENSG") %>%
      pivot_longer(!ENSG) %>%
      dplyr::select(!ENSG) %>%
      ggplot(aes(x=value, color= name))+
      geom_density() 

![](Alice-Template_files/figure-gfm/cpm%20norm-1.png)<!-- -->

### Violin plots

To directly compare sample gene expression distribution without overlap
between density plots we generated the respective violin plots

``` r
temp <- expr_log2cpm %>%
  rownames_to_column(var="ENSG") %>%
  pivot_longer(!ENSG) %>%
  dplyr::select(!ENSG) %>%
  mutate(condition = str_replace(string = name, 
                                 pattern = "_\\d+", 
                                 replacement = "")) 
# Convert condition to factor 
temp <- temp %>%
  mutate(condition = factor(condition, levels = c("Ctrl", "GEN9", "H2O2", "PNS2","PTS3")))

# Plot using ggplot with correct fill mapping
temp %>%
  ggplot(aes(x=name, y=value, fill = condition)) +
  geom_violin()+
  scale_fill_manual(values = c("Ctrl" = "blue", 
                               "GEN9" = "red", 
                               "H2O2" = "darkgreen", 
                               "PNS2" = "purple",
                               "PTS3" = "pink")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
```

![](Alice-Template_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->


    In this way we can easily identify that sample  **Ctrl_1**  and **GEN9_3** has a different distribution, with a higherproportion of genes with high log2CPM values

    ### Sample-sample correlation plot

    Using normalized counts we generated a plot to observe if there is correlation within condition groups which would group the samples accordingly.


    ``` r
    # Annotation

    # Standardize conditon names
    condition <- names(expr_log2cpm) %>% 
      str_replace(pattern = "_\\d+", 
                  replacement = "")

    annot <- data.frame(condition = as.factor(condition),
                        row.names = names(expr_log2cpm))

    annot_colors <- list(condition = c("Ctrl" = "blue", 
                                   "GEN9" = "red", 
                                   "H2O2" = "darkgreen", 
                                   "PNS2" = "purple",
                                   "PTS3" = "pink"))

    # Heatmap            
    expr_log2cpm %>%
      cor() %>%
      pheatmap(annotation_col = annot,
               annotation_row = annot,
               show_rownames = FALSE, 
               annotation_colors = annot_colors, 
               angle_col = 45)

![](Alice-Template_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->


    For the heatmap, one sample from the Ctrl group does not seem to cluster with the other two samples. I wonder if there is a solution to fix it and if it is going to impact the downstream analysis.

    Other than that, I think the samples in most conditions are nicely correlated, meaning that no intragroup variability. Except for GEN9, the correlation seems weaker.

    The next step is to perform a principal component analysis to further investigate this low correlation within condition groups.

### Principal Component Analysis (PCA) analysis

Another way of determining relationship between the samples is through a
PCA analysis, which reduces the dimentionality of our data to a set of
independent variables (principal components) that represent the major
proportion of variability in our data.

``` r
PCs <- prcomp(t(cpm(exp_mat)), center = TRUE, scale = TRUE)
# Scree plot 
fviz_eig(PCs)
```

![](Alice-Template_files/figure-gfm/PCA%20CPM-1.png)<!-- -->

``` r
# Scatter plot
eig_val <- get_eigenvalue(PCs)
PCs <- cbind(annot, PCs$x[,1:10])
PCs$sample_id <- rownames(PCs)

PCs <- PCs %>%
  mutate(color_class = case_when(
    condition == "Ctrl" ~ "blue",
    condition == "GEN9" ~ "red",
    condition == "H2O2" ~ "darkgreen",
    condition == "PNS2" ~ "purple",
    condition == "PTS3" ~ "pink")) 

PCs %>%
  ggplot(aes(x = PC1, y = PC2, 
             color = condition)) + 
  geom_point(aes(size = 8)) + 
  scale_color_manual(values = PCs$color_class, 
                     breaks = PCs$condition) +
  labs(x= paste("PC1 (",round(eig_val$variance.percent[1], 2),"%)", sep = ""),
       y= paste("PC2 (",round(eig_val$variance.percent[2], 2),"%)", sep = ""))+
  guides(size = "none") +
  theme_bw()
```

![](Alice-Template_files/figure-gfm/PCA%20CPM-2.png)<!-- -->

For the PCA plot, I could identify 2 outliers from the Ctrl and GEN9
groups (as what we expected in the gene expression density
plot).Additional metadata is required to evaluate other associations
driving most of the variability in the dataset (i.e. PC1).

Interestingly, although the dots are clustered they do not seem to
differ much on the x-axis but differ more on the y-axis.

Next, I generated another plot based on PC2 and PC3, they looked more
separated from each other

``` r
temp <- exp_mat[, c(1:3,4:6,7:9,10:12,13:15)]

#Top 1000 variable
# index <- which( rownames(temp) %in% names(sort(apply(X = temp, MARGIN = 1, var), decreasing = TRUE)[1:1000]))
# temp <- temp[index,]

PCs <- prcomp(t(cpm(temp)), center = TRUE, scale = TRUE)
# Scree plot 
fviz_eig(PCs)
```

![](Alice-Template_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Scatter plot
eig_val <- get_eigenvalue(PCs)
PCs <- as.data.table(PCs$x[,1:15])
PCs$condition <- as.data.table(annot, keep.rownames = TRUE)[condition %in% c("Ctrl", "GEN9","H2O2","PNS2","PTS3"), "condition"]
PCs$sample_id <- as.data.table(annot, keep.rownames = TRUE)[condition %in% c("Ctrl", "GEN9","H2O2","PNS2","PTS3"), "rn"]

PCs <- PCs %>%
  mutate(color_class = case_when(
    condition == "Ctrl" ~ "blue",
    condition == "GEN9" ~ "red",
    condition == "H2O2" ~ "darkgreen",
    condition == "PNS2" ~ "purple",
    condition == "PTS3" ~ "pink"))

PCs %>%
  ggplot(aes(x = PC2, y = PC3, 
             color = condition)) + 
  geom_point(aes(size = 8)) + 
  scale_color_manual(values = PCs$color_class, 
                     breaks = PCs$condition) +
  labs(x= paste("PC2 (",round(eig_val$variance.percent[2], 2),"%)", sep = ""),
       y= paste("PC3 (",round(eig_val$variance.percent[3], 2),"%)", sep = ""))+
  guides(size = "none") 
```

![](Alice-Template_files/figure-gfm/unnamed-chunk-4-2.png)<!-- --> It is
likely that combined explain the expression variability associated with
the experimental conditions (7% + 19.8%). The first PC could be
associated with the exposure to H2O2, which we expect will induce most
of the expression variability in the dataset.

## EDA conclusions

Initial visualization of gene counts suggested sample **Ctrl_1** and
**GEN0_3** to be outliers due to higher reads. Sample correlation
resulted in clustering by each condition (weaker correaltion in Ctrl and
GEN9 but still \~0.9).

PCA analysis indicated a separation between conditions along the second
principal component. Ctrl and PNS2 have more similarities with each
other compared to other experimental conditions, while GEN9 showed
furthest from the Ctrl condition. On the same line, PTS3 was the least
effective of the treatments since the expression is more similar to the
H2O2 group.

VST counts has not yet been conducted as clustering has shown in the
initial EDA.

## Variance stabilized counts visualization.

We observed a random clustering of the samples in both sample-sample
correlation and PCA analysis. This could be caused because the majority
of genes have a low variance while a handful are highly variable. To
improve visualization we used the `vst` function from `DESeq2`.

``` r
expr_vst <- vst(object = as.matrix(sapply(exp_mat, as.integer))) %>%
  data.frame(row.names = rownames(expr_log2cpm))
```

``` r
expr_vst %>%
  rownames_to_column(var = "ENSG") %>%
  pivot_longer(!ENSG) %>%
  dplyr::select(!ENSG) %>%
    ggplot(aes(x=value + 1, color = name))+
    geom_density() +
    scale_x_continuous(trans = 'log10', labels = scales::comma) +
    labs(title = "Gene expression density plot", 
         subtitle = "Colored by Sample", 
         x = "VST counts + 1", 
         y = "Density",
         color = "Sample") +
    theme_bw()
```

![](Alice-Template_files/figure-gfm/vst%20density-1.png)<!-- -->

``` r
temp <- expr_vst %>%
  rownames_to_column(var = "ENSG") %>%
  pivot_longer(!ENSG) %>%
  mutate(condition = str_replace(string = name, 
                                 pattern = ".\\d+", 
                                 replacement = "")) %>%
  mutate(condition = str_replace(string = condition, 
                                 pattern = "_\\d{6}.\\d+", 
                                 replacement = "")) %>%
  mutate(condition = case_when(
    condition == "Ctrl" ~ "blue",
    condition == "GEN9" ~ "red",
    condition == "H2O2" ~ "darkgreen",
    condition == "PNS2" ~ "purple",
    condition == "PTS3" ~ "pink")) %>% 
  dplyr::select(!ENSG) 

temp %>%
  ggplot(aes(x=value, color= name))+
  geom_density() +
  scale_x_continuous(trans = 'log10', labels = scales::comma) +
  scale_color_manual(values = temp$condition, 
                     breaks = temp$name) 
```

![](Alice-Template_files/figure-gfm/vst%20density-2.png)<!-- -->

``` r
  labs(x = "VST")
```

    ## $x
    ## [1] "VST"
    ## 
    ## attr(,"class")
    ## [1] "labels"

VST count transformation improved the visualization of the gene
expression values for every sample.

### Violin plots

This can also be appreciated using violin plots for representing the VST
counts per sample.

``` r
temp <- expr_vst %>%
  rownames_to_column(var="ENSG") %>%
  pivot_longer(!ENSG) %>%
  mutate(condition = str_replace(string = name, 
                            pattern = "_\\d+", 
                            replacement = "")) %>%
  mutate(color_class = case_when(
    condition == "Ctrl" ~ "blue",
    condition == "GEN9" ~ "red",
    condition == "H2O2" ~ "darkgreen",
    condition == "PNS2" ~ "purple",
    condition == "PTS3" ~ "pink")) %>%
  dplyr::select(!ENSG) 

temp %>%
  ggplot(aes(x=name, y=value, fill = condition)) +
  geom_violin()+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values= c("blue","red","darkgreen", "purple","pink")) 
```

![](Alice-Template_files/figure-gfm/vst%20violin-1.png)<!-- -->

### Sample-sample correlation plot

Next, we calculated the correlation between samples using VST counts.
Unfortunately, sample clustering still looks random to some degree.

``` r
# Annotation
condition <- names(expr_vst) %>% 
  str_replace(pattern = "_\\d+", 
              replacement = "")

annot <- data.frame(condition = as.factor(condition),
                    row.names = names(expr_log2cpm))

annot_colors <- list(condition = c("Ctrl" = "blue", 
                               "GEN9" = "red", 
                               "H2O2" = "darkgreen", 
                               "PNS2" = "purple",
                               "PTS3" = "pink"))

# Heatmap            
expr_vst %>%
  cor() %>%
  pheatmap(annotation_col = annot,
           annotation_row = annot,
           show_rownames = FALSE, 
           annotation_colors = annot_colors,
           angle_col = 45)
```

![](Alice-Template_files/figure-gfm/vst%20sample-sample%20correlation-1.png)<!-- -->

#### Top variable genes

Since we expect subtle changes in gene expression, we can narrow the set
of genes used to calculate the correlation between samples to the 100
variable genes.

``` r
genes_to_keep <- names(sort(apply(exp_mat,MARGIN = 1,FUN = var),
                            decreasing = TRUE)[1:100])
expr_vst_top <- expr_vst %>% 
  rownames_to_column(var = "gene") %>% 
  data.table(key = "gene")

expr_vst_top <- expr_vst_top[genes_to_keep] %>% 
  column_to_rownames(var = "gene")

# Heatmap            
expr_vst_top %>%
  cor() %>%
  pheatmap(annotation_col = annot,
           annotation_row = annot,
           show_rownames = FALSE, 
           annotation_colors = annot_colors,
           angle_col = 45, treeheight_row = 0)
```

![](Alice-Template_files/figure-gfm/vst%20heatmap%20top%20genes-1.png)<!-- -->

Using only the top 100 variable genes improved sample clustering,
although it is still difficult to identify clusters with a single diet
group. It is important to notice that the right-most cluster is composed
primarily of the reference diet (CSAA) and CGA supplemented diet. This
result shows that expression profiles of highly variable genes between
CSAA and CSAA + CGA diets are very similar.

### PCA analysis

Our final EDA exploration is a PCA using VST counts.

``` r
# Remove RSV_14 due to heavy bias on PCA
temp <- expr_vst[, names(expr_vst) != "RSV_14"]

# Annotation
condition <- names(temp) %>% 
  str_replace(pattern = "_\\d+", 
              replacement = "")

annot <- data.frame(condition = as.factor(condition),
                    row.names = names(temp))


# PCA analysis
PCs <- prcomp(t(temp), center = TRUE, scale = TRUE)
# Scree plot 
fviz_eig(PCs)
```

![](Alice-Template_files/figure-gfm/vst%20PCA-1.png)<!-- -->

``` r
# Scatter plot
eig_val <- get_eigenvalue(PCs)
PCs <- cbind(annot, PCs$x[,1:10])
PCs$sample_id <- str_extract(string = rownames(PCs), pattern = "\\d+")


PCs <- PCs %>%
  mutate(color_class = case_when(
    condition == "Ctrl" ~ "blue",
    condition == "GEN9" ~ "red",
    condition == "H2O2" ~ "darkgreen",
    condition == "PNS2" ~ "purple",
    condition == "PTS3" ~ "pink")) 

PCs %>%
  ggplot(aes(x = PC1, y = PC2, color = condition)) + 
  geom_point(aes(size = 8))+
  scale_color_manual(values = PCs$color_class, 
                     breaks = PCs$condition) +
  labs(x= paste("PC1 (",round(eig_val$variance.percent[1], 2),"%)", sep = ""),
       y= paste("PC2 (",round(eig_val$variance.percent[2], 2),"%)", sep = ""))+
  guides(size = "none") 
```

![](Alice-Template_files/figure-gfm/vst%20PCA-2.png)<!-- --> Using VST
counts showed a clear separation along the second PC between CSAA and
CSAA + PTS condition Suggesting that a low
