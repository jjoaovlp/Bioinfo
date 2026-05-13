setwd("C:/Users/João - PC/Desktop/Tópicos Avançados em Bioinformática IV/Nova pasta/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")

library(data.table)
# library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(DESeq2)
library(matrixStats)
library(pheatmap)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(apeglm)
library(ashr)

# Metadados (GEO SOFT)
meta <- read.table(
  "sra.sra.SRP111173.MD/metadados(Sheet1).csv",
  sep = ";",
  fill = TRUE,
  header = TRUE,
  stringsAsFactors = FALSE,
)

# Contagens (CSV de STAR)
counts_raw <- read.table(
  "sra.gene_sums.SRP111173.G029/RNA counts.csv",
  sep = ";",
  fill = TRUE,
  header = TRUE,
  stringsAsFactors = FALSE,
)

# Detecta coluna de genes (primeira coluna) e define rownames
if (!is.numeric(counts_raw[[1]])) {
  rownames(counts_raw) <- counts_raw[[1]]
  counts_raw <- counts_raw[, -1]
}

# Mantem apenas amostras presentes nos metadados
common_samples <- intersect(colnames(counts_raw), meta$código)
counts <- counts_raw[, common_samples, drop = FALSE]
meta <- meta %>% filter(código %in% common_samples)

# Ordena meta para alinhar com as colunas
meta <- meta %>% arrange(match(código, colnames(counts)))

stopifnot(all(meta$código == colnames(counts)))

## Dia 1 - QC e exploracao

#```{r qc-basic}

# Tamanho de biblioteca
lib_size <- colSums(counts)
qc_tbl <- meta %>%
  mutate(lib_size = lib_size, genes_detected = colSums(counts > 0))

p_lib <- ggplot(qc_tbl, aes(x = código, y = lib_size)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Tamanho de biblioteca", y = "Reads", x = "")

p_genes <- ggplot(qc_tbl, aes(x = código, y = genes_detected)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genes detectados", y = "Genes", x = "")

p_lib
p_genes

# ```{r qc-distribution}

# Correlacao entre amostras (log1p)

log_counts <- log1p(counts)
sample_cor <- cor(log_counts, use = "pairwise.complete.obs")

sample_cor <- sample_cor[meta$código, meta$código]

anno_col <- data.frame(
  amostra = factor(meta$amostra),
  tratamento = factor(meta$tratamento),
  replicata.biológica = factor(meta$replicada.biológica),
  row.names = meta$código
)

amostra_levels <- levels(anno_col$amostra)
tratamento_levels <- levels(anno_col$tratamento)
replicada.biológica_levels <- levels(anno_col$replicada.biológica)

amostra_cols <- setNames(colorRampPalette(brewer.pal(3, "Blues"))(length(amostra_levels)), amostra_levels)
tratamento_cols <- setNames(colorRampPalette(brewer.pal(3, "Greens"))(length(tratamento_levels)), tratamento_levels)
rb_cols <- setNames(colorRampPalette(brewer.pal(3, "Set1"))(length(replicada.biológica_levels)), replicada.biológica_levels)

ha <- HeatmapAnnotation(
  df = anno_col,
  col = list(
    amostra = amostra_cols,
    tratamento = tratamento_cols,
    replicada.biológica = rb_cols
  ),
  annotation_name_side = "left"
)

Heatmap(
  sample_cor,
  name = "cor",
  top_annotation = ha,
  col = colorRamp2(c(0.7, 0.85, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Correlacao entre amostras"
)

# ```{r qc-mad}

# Variancia por amostra (VST)

vsd_mad <- vst(DESeqDataSetFromMatrix(counts, meta, ~ 1))
vst_mad_mat <- assay(vsd_mad)
sample_var <- matrixStats::colVars(vst_mad_mat)
mad_tbl <- meta %>% mutate(sample_var = sample_var)

p_mad <- ggplot(mad_tbl, aes(x = código, y = código_var)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Variancia por amostra", y = "Variancia", x = "")

p_mad

# ```{r qc-prepost}

# Normalizacao (VST) para comparar distribuicoes

vsd <- vst(DESeqDataSetFromMatrix(counts, meta, ~ 1))
vst_mat <- assay(vsd)

raw_tbl <- as.data.frame(log1p(counts)) %>%
  pivot_longer(cols = everything(), names_to = "código", values_to = "value") %>%
  mutate(stage = "Pre-normalizacao")

vst_tbl <- as.data.frame(vst_mat) %>%
  pivot_longer(cols = everything(), names_to = "código", values_to = "value") %>%
  mutate(stage = "Pos-normalizacao (VST)")

dist_tbl <- bind_rows(raw_tbl, vst_tbl) %>%
  left_join(meta %>% select(código, tratamento), by = "código")

p_prepost <- ggplot(dist_tbl, aes(x = código, y = value, color = tratamento)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ stage, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Distribuicao pre vs pos normalizacao", x = "Amostras", y = "Expressao")

p_prepost

# ```{r qc-hclust}

# Clustering hierarquico com top 500 genes por MAD

gene_mad <- matrixStats::rowMads(vst_mat)
top_idx <- order(gene_mad, decreasing = TRUE)[1:500]

anno_col <- meta %>%
  dplyr::select(amostra, tratamento, replicada.biológica) %>%
  as.data.frame()
rownames(anno_col) <- meta$código

pheatmap::pheatmap(
  vst_mat[top_idx, ],
  show_rownames = FALSE,
  annotation_col = anno_col,
  clustering_distance_cols = "correlation",
  main = "Clustering hierarquico (Top 500 MAD)"
)

# ```{r qc-pca}

# PCA usando VST

pca <- prcomp(t(vst_mat))

pca_tbl <- as.data.frame(pca$x) %>%
  mutate(código = rownames(pca$x)) %>%
  left_join(meta, by = "código")

pca_tbl <- pca_tbl %>% mutate(amostra = factor(amostra))

p_pca <- ggplot(pca_tbl, aes(x = PC1, y = PC2, color = amostra, shape = tratamento)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA (VST)")

p_pca2 <- ggplot(pca_tbl, aes(x = PC1, y = PC2, color = factor(amostra), shape = tratamento)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA (VST)")

cowplot::plot_grid(p_pca, p_pca2)

# Dia 2

### Setup DESeq2

meta <- meta %>%
  mutate(
    amostra = factor(amostra),
    tratamento = factor(tratamento),
    replicata_biológica = factor(replicada.biológica)
  )

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ tratamento
)

# Filtro basico
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

### Contraste e resultados

res <- lfcShrink(dds, coef = "tratamento_RA_vs_DMSO", type = "apeglm")

res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
head(res_tbl, 10)

# MA plot
plotMA(res, ylim = c(-4, 4))

# Volcano plot
vol_tbl <- res_tbl %>%
  mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "sig", "ns"))

p_volcano <- ggplot(vol_tbl, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey60", "sig" = "red")) +
  theme_bw() +
  labs(
    title = "Volcano plot – RA vs DMSO",
    x = "log2 Fold Change",
    y = "-log10(padj)"
  )
p_volcano

# Volcano SEM shrinkage para diagnóstico
res_sem_shrink <- results(dds, contrast = c("tratamento", "RA", "DMSO"))

vol_tbl_raw <- as.data.frame(res_sem_shrink) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "sig", "ns"))

ggplot(vol_tbl_raw, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey60", "sig" = "red")) +
  theme_bw() +
  labs(
    title = "Volcano plot – RA vs DMSO (sem shrinkage)",
    x = "log2 Fold Change",
    y = "-log10(padj)"
  )

# Resumo para ver quantos genes DEG existem
summary(res_sem_shrink)

# Volcano com type ashr

res <- lfcShrink(dds, 
                 contrast = c("tratamento", "RA", "DMSO"), 
                 type = "ashr")
res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "sig", "ns"))

ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey60", "sig" = "red")) +
  theme_bw() +
  labs(
    title = "Volcano plot – RA vs DMSO",
    x = "log2 Fold Change",
    y = "-log10(padj)"
  )


#Volcano por linhagem

linhagens <- levels(meta$amostra)
res_list <- lapply(linhagens, function(lin) {
  idx <- meta$amostra == lin
  
  # Pular linhagens com menos de 2 amostras por tratamento
  tab <- table(meta$tratamento[idx])
  if (any(tab < 2)) {
    message("Pulando ", lin, " — réplicas insuficientes")
    return(NULL)
  }
  
  message("Processando: ", lin)
  
  dds_l <- DESeqDataSetFromMatrix(
    countData = counts[, idx, drop = FALSE],
    colData   = meta[idx, ],
    design    = ~ tratamento
  )
  dds_l <- dds_l[rowSums(counts(dds_l)) >= 10, ]
  dds_l$tratamento <- relevel(dds_l$tratamento, ref = "DMSO")
  dds_l <- DESeq(dds_l)
  
  res_l <- lfcShrink(dds_l,
                     contrast = c("tratamento", "RA", "DMSO"),
                     type = "ashr")
  
  as.data.frame(res_l) %>%
    rownames_to_column("gene") %>%
    mutate(linhagem = lin)
})

# Remover linhagens puladas
res_all <- bind_rows(Filter(Negate(is.null), res_list))

# Volcano para cada linhagem
plots <- lapply(linhagens, function(lin) {
  df <- res_all %>%
    filter(linhagem == lin) %>%
    mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "sig", "ns"))
  
  if (nrow(df) == 0) return(NULL)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("ns" = "dodgerblue", "sig" = "red")) +
    theme_bw() +
    labs(
      title = paste("Volcano –", lin, "| RA vs DMSO"),
      x = "log2 Fold Change",
      y = "-log10(padj)"
    )
})
plots <- Filter(Negate(is.null), plots)

# Exibir todos
library(patchwork)
wrap_plots(plots, ncol = 2)

# Salvar
ggsave("volcano_por_linhagem.png", 
       wrap_plots(plots, ncol = 2),
       width = 14, height = 18, dpi = 300)