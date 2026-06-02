setwd("C:/Users/João - PC/Desktop/Diretório R/Tópicos Avançados em Bioinformática IV/WGCNA.XPA")

library(data.table)
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
library(ashr)
library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(WGCNA)
library(lionessR)
library(SummarizedExperiment)
library(AnnotationDbi)
library(GSVA)
library(msigdbr)

# Metadados (GEO SOFT)
meta <- read.table(
  "sra.sra.SRP111173.MD/metadados(Sheet1).csv",
  sep = ",",
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

##############################################################################
# ── Mapeamento ENSEMBL → símbolo ──────────────────────────────────────────────

# Remove versão dos rownames do counts (ENSG00000278704.1 → ENSG00000278704)
ensembl_ids_clean <- sub("\\.(\\d+)(_PAR_Y)?$", "", rownames(counts))

# Busca símbolos
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(ensembl_ids_clean),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# Função auxiliar — recebe IDs com ou sem versão
ensembl_to_symbol <- function(ids) {
  ids_clean <- sub("\\.(\\d+)(_PAR_Y)?$", "", ids)
  syms <- map_df$SYMBOL[match(ids_clean, map_df$ENSEMBL)]
  ifelse(is.na(syms) | syms == "", ids_clean, syms)  # já está assim — sem mudança
}

##############################################################################
## Dia 1 - QC e exploracao

# Tamanho de biblioteca
lib_size <- colSums(counts)
qc_tbl <- meta %>%
  mutate(lib_size = lib_size, genes_detected = colSums(counts > 0))

# Ordena por proficiencia + tratamento para agrupar visualmente
qc_tbl <- qc_tbl %>%
  arrange(proficiencia, tratamento, código) %>%
  mutate(código = factor(código, levels = código))

p_lib <- ggplot(qc_tbl, aes(x = código, y = lib_size / 1e6, fill = proficiencia, alpha = tratamento)) +
  geom_col() +
  geom_hline(yintercept = mean(qc_tbl$lib_size / 1e6), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = c("proficiente" = "#2166AC", "deficiente" = "#D62728")) +
  scale_alpha_manual(values = c("DMSO" = 0.7, "RA" = 1.0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title   = "Tamanho de biblioteca",
    y       = "Reads (milhões)",
    x       = NULL,
    fill    = "Proficiência",
    alpha   = "Tratamento"
  )

p_genes <- ggplot(qc_tbl, aes(x = código, y = genes_detected / 1e3, fill = proficiencia, alpha = tratamento)) +
  geom_col() +
  geom_hline(yintercept = mean(qc_tbl$genes_detected / 1e3), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = c("proficiente" = "#2166AC", "deficiente" = "#D62728")) +
  scale_alpha_manual(values = c("DMSO" = 0.7, "RA" = 1.0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom"
  ) +
  labs(
    title = "Genes detectados",
    y     = "Genes (milhares)",
    x     = NULL,
    fill  = "Proficiência",
    alpha = "Tratamento"
  )

p_lib / p_genes

##############################################################################
# Correlacao entre amostras (log1p)
log_counts <- log1p(counts)
sample_cor <- cor(log_counts, use = "pairwise.complete.obs")
sample_cor <- sample_cor[meta$código, meta$código]
# rownames/colnames permanecem como SRR — não substituir

anno_col <- data.frame(
  amostra      = factor(meta$amostra),
  tratamento   = factor(meta$tratamento),
  proficiencia = factor(meta$proficiencia),
  row.names    = meta$código    # usa código SRR direto, sem duplicata
)

amostra_levels    <- levels(anno_col$amostra)
n_amostras        <- length(amostra_levels)
amostra_cols      <- setNames(
  RColorBrewer::brewer.pal(max(n_amostras, 3), "Set2")[1:n_amostras],
  amostra_levels
)
tratamento_cols   <- c("DMSO" = "#A8D5A2", "RA" = "#2D6A4F")
proficiencia_cols <- c("proficiente" = "#2166AC", "deficiente" = "#D62728")

ha <- HeatmapAnnotation(
  amostra      = anno_col$amostra,
  tratamento   = anno_col$tratamento,
  proficiencia = anno_col$proficiencia,
  col = list(
    amostra      = amostra_cols,
    tratamento   = tratamento_cols,
    proficiencia = proficiencia_cols
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9),
  simple_anno_size     = unit(4, "mm")
)

Heatmap(
  sample_cor,
  name              = "cor",
  top_annotation    = ha,
  col               = colorRamp2(c(0.7, 0.85, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
  show_row_names    = TRUE,
  show_column_names = FALSE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = 5.5),
  row_dend_side     = "left",
  row_dend_width    = unit(8, "mm"),
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "pearson",
  column_title    = "Correlação entre amostras",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    title          = "Pearson r",
    direction      = "horizontal",
    legend_width   = unit(3, "cm"),
    title_gp       = gpar(fontsize = 9),
    labels_gp      = gpar(fontsize = 8),
    title_position = "lefttop"
  )
) %>%
  draw(
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    legend_grouping        = "original",
    padding = unit(c(2, 20, 10, 2), "mm")
  )

##############################################################################
# Variancia por amostra (VST)
vsd_mad     <- vst(DESeqDataSetFromMatrix(counts, meta, ~ 1))
vst_mad_mat <- assay(vsd_mad)
sample_var  <- matrixStats::colVars(vst_mad_mat)
mad_tbl     <- meta %>% mutate(sample_var = sample_var)

p_mad <- ggplot(mad_tbl, aes(x = código, y = sample_var, fill = proficiencia)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Variancia por amostra", y = "Variancia", x = "")
p_mad

##############################################################################
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
  left_join(dplyr::select(meta, código, tratamento, proficiencia), by = "código")

p_prepost <- ggplot(dist_tbl, aes(x = código, y = value, color = tratamento, fill = proficiencia)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ stage, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Distribuicao pre vs pos normalizacao", x = "Amostras", y = "Expressao")
p_prepost

##############################################################################
# Clustering hierarquico com top 500 genes por MAD
gene_mad <- matrixStats::rowMads(vst_mat)
top_idx  <- order(gene_mad, decreasing = TRUE)[1:500]

anno_col <- meta %>%
  dplyr::select(amostra, tratamento, proficiencia) %>%
  as.data.frame()
rownames(anno_col) <- meta$código

pheatmap::pheatmap(
  vst_mat[top_idx, ],
  show_rownames = FALSE,
  annotation_col = anno_col,
  clustering_distance_cols = "correlation",
  fontsize = 7,
  fontsize_col = 6,
  main = "Clustering hierarquico (Top 500 MAD)"
)

##############################################################################
# PCA usando VST
pca <- prcomp(t(vst_mat))

pca_tbl <- as.data.frame(pca$x) %>%
  mutate(código = rownames(pca$x)) %>%
  left_join(meta, by = "código") %>%
  mutate(amostra = factor(amostra))

p_pca <- ggplot(pca_tbl, aes(x = PC1, y = PC2, color = amostra, shape = tratamento)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA (VST)")

p_pca2 <- ggplot(pca_tbl, aes(x = PC1, y = PC2, color = factor(amostra), shape = tratamento)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA (VST)")

cowplot::plot_grid(p_pca, p_pca2)

## Dia 2

##############################################################################
### Setup DESeq2

meta <- meta %>%
  mutate(
    amostra      = factor(amostra),
    tratamento   = factor(tratamento),
    proficiencia = factor(proficiencia)
  )

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ tratamento
)
dds$tratamento <- relevel(dds$tratamento, ref = "DMSO")
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

### Contraste global: RA vs DMSO (ashr)

res <- lfcShrink(
  dds,
  contrast = c("tratamento", "RA", "DMSO"),
  type     = "ashr"
)

res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

res_tbl <- res_tbl %>%
  mutate(symbol = ensembl_to_symbol(gene))

head(res_tbl, 10)

##############################################################################
### Volcano global – RA vs DMSO (com ashr)

vol_tbl <- res_tbl %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    color = case_when(
      padj < 0.05 & log2FoldChange >  0 ~ "Positivo (up)",
      padj < 0.05 & log2FoldChange <  0 ~ "Negativo (down)",
      TRUE                               ~ "Não significativo"
    ),
    color = factor(color, levels = c("Positivo (up)", "Negativo (down)", "Não significativo"))
  )

top10_vol <- vol_tbl %>% arrange(padj) %>% slice_head(n = 10)

ggplot(vol_tbl, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c(
    "Positivo (up)"     = "#D62728",
    "Negativo (down)"   = "#1F77B4",
    "Não significativo" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data          = top10_vol,
    aes(label     = symbol),
    size          = 2.8,
    color         = "black",
    fontface      = "bold",
    max.overlaps  = Inf,
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size  = 0.3,
    show.legend   = FALSE
  ) +
  theme_bw() +
  labs(
    title = "Volcano plot – RA vs DMSO (ashr)",
    x     = "log2 Fold Change",
    y     = "-log10(padj)",
    color = NULL
  )

##############################################################################
### GSEA

# Ranking por log2FoldChange do DESeq2 -> Deficiente

meta_group <- meta %>%
  mutate(
    group = paste(proficiencia, tratamento, sep = " | "),
    group = factor(
      group,
      levels = c(
        "proficiente | DMSO",
        "proficiente | RA",
        "deficiente | DMSO",
        "deficiente | RA"
      )
    )
  )

group_levels <- unique(as.character(meta_group$proficiencia))

res_list <- lapply(group_levels, function(g) {
  
  idx      <- meta_group$proficiencia == g
  meta_sub <- meta_group[idx, ] %>%
    mutate(tratamento = droplevels(factor(tratamento)))
  
  dds_g <- DESeqDataSetFromMatrix(
    countData = counts[, idx, drop = FALSE],
    colData   = meta_sub,
    design    = ~ tratamento
  )
  dds_g <- dds_g[rowSums(counts(dds_g)) >= 10, ]
  dds_g$tratamento <- relevel(dds_g$tratamento, ref = "DMSO")
  dds_g <- DESeq(dds_g)
  res_g <- lfcShrink(
    dds_g,
    contrast = c("tratamento", "RA", "DMSO"),
    type     = "ashr"
  )
  
  as.data.frame(res_g) %>%
    rownames_to_column("gene") %>%
    mutate(group = g)
})

res_group_tbl <- bind_rows(res_list)

res_group_tbl <- res_group_tbl %>%
  mutate(symbol = ensembl_to_symbol(gene))

# Top 5 genes por grupo
res_group_tbl %>%
  group_by(group) %>%
  arrange(padj) %>%
  slice_head(n = 5)

rank_tbl <- res_group_tbl %>%
  filter(group == "deficiente", !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
  mutate(gene_clean = sub("\\..*$", "", gene))

rank_entrez <- bitr(
  unique(rank_tbl$gene_clean),
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

rank_vec <- rank_tbl %>%
  inner_join(rank_entrez, by = c("gene_clean" = "ENSEMBL")) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  filter(!is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  { setNames(.$log2FoldChange, .$ENTREZID) }

gsea <- gseGO(
  geneList      = rank_vec,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

gsea_tbl <- as.data.frame(gsea)

top_up   <- gsea_tbl %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = 10)
top_down <- gsea_tbl %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = 10)

gsea_top <- bind_rows(top_up, top_down) %>%
  mutate(Description = factor(Description, levels = Description[order(NES)]))

ggplot(gsea_top, aes(x = NES, y = Description, color = p.adjust, size = setSize)) +
  geom_point() +
  scale_color_viridis_c(option = "C", direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  theme_bw() +
  labs(
    title = "GSEA GO:BP — deficiente | RA vs DMSO",
    x     = "NES (Normalized Enrichment Score)",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene set size"
  )

# Ranking por log2FoldChange do DESeq2 -> Proficiente
rank_tbl <- res_group_tbl %>%
  filter(group == "proficiente", !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
  mutate(gene_clean = sub("\\..*$", "", gene))

rank_entrez <- bitr(
  unique(rank_tbl$gene_clean),
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

rank_vec <- rank_tbl %>%
  inner_join(rank_entrez, by = c("gene_clean" = "ENSEMBL")) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  filter(!is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  { setNames(.$log2FoldChange, .$ENTREZID) }

gsea <- gseGO(
  geneList      = rank_vec,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

gsea_tbl <- as.data.frame(gsea)

top_up   <- gsea_tbl %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = 10)
top_down <- gsea_tbl %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = 10)

gsea_top <- bind_rows(top_up, top_down) %>%
  mutate(Description = factor(Description, levels = Description[order(NES)]))

ggplot(gsea_top, aes(x = NES, y = Description, color = p.adjust, size = setSize)) +
  geom_point() +
  scale_color_viridis_c(option = "C", direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  theme_bw() +
  labs(
    title = "GSEA GO:BP — proficiente | RA vs DMSO",
    x     = "NES (Normalized Enrichment Score)",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene set size"
  )

##############################################################################
### PPI com STRING + igraph

string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 500)

# Mapeia genes significativos para STRING

sig_genes <- res_group_tbl %>%
  filter(padj < 0.05, abs(log2FoldChange) > log2(1.5), group == "deficiente") %>%
  pull(gene)

# Converte para Entrez (IDs estao em ENSEMBL)
sig_genes <- sub("\\..*$", "", sig_genes)
entrez <- bitr(sig_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

mapped <- string_db$map(data.frame(gene = sig_genes), "gene", removeUnmappedRows = TRUE)

ppi <- string_db$get_interactions(mapped$STRING_id)

# Grafo PPI
g <- graph_from_data_frame(ppi[, c("from", "to")], directed = FALSE)

set.seed(123)
coords <- layout_with_kk(g)
node_tbl <- data.frame(
  id = V(g)$name,
  x  = coords[, 1],
  y  = coords[, 2]
)
edge_tbl <- as.data.frame(get.edgelist(g))
colnames(edge_tbl) <- c("from", "to")
edge_tbl <- edge_tbl %>%
  left_join(node_tbl, by = c("from" = "id")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(node_tbl, by = c("to" = "id")) %>%
  rename(x2 = x, y2 = y)

# Grau de conectividade para destacar hubs
# Mapeia STRING ID → símbolo para o grafo
mapped <- mapped %>%
  mutate(symbol = ensembl_to_symbol(gene))
string_id_to_symbol <- setNames(mapped$symbol, mapped$STRING_id)

node_tbl <- data.frame(
  id     = V(g)$name,
  x      = coords[, 1],
  y      = coords[, 2],
  symbol = ifelse(
    !is.na(string_id_to_symbol[V(g)$name]),
    string_id_to_symbol[V(g)$name],
    V(g)$name
  )
)
node_tbl$degree <- degree(g)

edge_tbl <- as.data.frame(get.edgelist(g))
colnames(edge_tbl) <- c("from", "to")
edge_tbl <- edge_tbl %>%
  left_join(node_tbl, by = c("from" = "id")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(node_tbl, by = c("to" = "id")) %>%
  rename(x2 = x, y2 = y)

p_ppi <- ggplot() +
  geom_segment(
    data = edge_tbl,
    aes(x = x1, y = y1, xend = x2, yend = y2),
    alpha = 0.2, color = "grey60"
  ) +
  geom_point(
    data = node_tbl,
    aes(
      x = x,
      y = y,
      size = degree,
      color = degree
    ),
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    data          = node_tbl,
    aes(x = x, y = y, label = symbol),   # <-- corrigido
    size          = 2.5,
    color         = "black",
    max.overlaps  = Inf,
    box.padding   = 0.3,
    segment.color = "grey40",
    segment.size  = 0.2
  ) +
  scale_size_continuous(range = c(1.5, 6)) +
  theme_void() +
  labs(
    title = "PPI (STRING + igraph + ggplot2)",
    size  = "Grau"
  ) +
  theme(
   legend.title = element_text(size = 8),
   legend.text  = element_text(size = 6)
  )

p_ppi

##############################################################################
### WGCNA

# Usa VST para reduzir efeitos de contagem
vsd <- vst(dds)
expr <- t(assay(vsd))

# Escolha do soft-threshold
powers <- 1:10
sft <- pickSoftThreshold(expr, powerVector = powers, verbose = 0)

# Exemplo com power recomendado
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 6

net <- blockwiseModules(
  expr,
  power = soft_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  verbose = 0
)

module_eigengenes <- net$MEs

##############################################################################
# Relacao com fenotipos

trait_tbl <- meta %>%
  mutate(
    amostra      = factor(amostra),
    tratamento   = factor(tratamento,   levels = c("DMSO",         "RA")),       # ref = DMSO → aparece RA
    proficiencia = factor(proficiencia, levels = c("proficiente", "deficiente"))  # ref = proficiente → aparece deficiente
  )

# Matriz de traits (dummy) para correlacao com modulos
trait_mat <- model.matrix(~ tratamento + proficiencia, time, data = trait_tbl)
trait_mat <- trait_mat[, -1, drop = FALSE]

colnames(trait_mat) <- gsub("^tratamento", "", colnames(trait_mat))
colnames(trait_mat) <- gsub("^proficiencia", "", colnames(trait_mat))

rownames(trait_mat) <- meta$código

# Garante mesma ordem de linhas
common_rows  <- intersect(rownames(module_eigengenes), rownames(trait_mat))
ME_aligned   <- module_eigengenes[common_rows, , drop = FALSE]
trait_aligned <- trait_mat[common_rows, , drop = FALSE]

# ── Calcula correlação e p-valor módulo × trait ────────────────────────────────
module_trait_cor <- cor(ME_aligned, trait_aligned, use = "pairwise.complete.obs")

# p-valores via corPvalueStudent() do próprio WGCNA
nSamples        <- nrow(ME_aligned)
module_trait_p  <- corPvalueStudent(module_trait_cor, nSamples)

# Heatmap de correlacao modulo x trait
pheatmap::pheatmap(
  module_trait_cor,
  color           = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50),
  display_numbers = TRUE,
  number_format   = "%.2f",
  number_color    = "black",
  fontsize        = 9,          # labels dos eixos
  fontsize_row    = 7,          # labels dos módulos (eixo direito)
  fontsize_col    = 10,         # labels das traits (eixo inferior)
  fontsize_number = 6,          # números dentro das células
  angle_col       = 45,
  cellwidth       = 40,         # largura fixa por célula
  cellheight      = 8,          # altura fixa por célula — compacta as linhas
  main            = "WGCNA: Módulos vs Traits"
)

# Aplica teto em 10 para não deixar 1-2 módulos dominarem a escala
p_mat_capped <- pmin(-log10(module_trait_p), 10)

pheatmap::pheatmap(
  p_mat_capped,
  color           = colorRampPalette(c("#F7F7F7", "#FDB863", "#B2182B"))(50),
  display_numbers = matrix(
    sprintf("%.1f", -log10(module_trait_p)),
    nrow = nrow(module_trait_p)
  ),
  number_color    = "black",
  fontsize        = 9,
  fontsize_row    = 7,
  fontsize_col    = 10,
  fontsize_number = 6,
  angle_col       = 45,
  cellwidth       = 40,
  cellheight      = 8,
  main            = "WGCNA: -log10(p)"
)

# Scatter de um trait continuo (time) com o modulo mais correlacionado
if ("time" %in% colnames(trait_tbl)) {
  me_names <- colnames(module_eigengenes)
  cor_time <- module_trait_cor[, grepl("^time", colnames(module_trait_cor)), drop = FALSE]
  top_mod <- rownames(cor_time)[which.max(abs(cor_time[, 1]))]
  
  plot_tbl <- data.frame(
    ME = module_eigengenes[, top_mod],
    time = as.numeric(factor(trait_tbl$time))
  )
  
  p_scatter <- ggplot(plot_tbl, aes(x = time, y = ME)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "#B2182B") +
    theme_bw() +
    labs(title = paste("Modulo", top_mod, "vs time"), x = "Time (ordinal)")
  
  p_scatter
}

# Boxplot do modulo mais associado ao fator tratamento
if (any(colnames(module_trait_cor) %in% c("DMSO", "RA"))) {
  trat_cols <- intersect(colnames(module_trait_cor), c("DMSO", "RA"))
  cor_tratamento <- module_trait_cor[, trat_cols, drop = FALSE]
  
  # Módulo com maior correlação absoluta com tratamento
  top_mod_tratamento <- rownames(cor_tratamento)[which.max(abs(cor_tratamento[, 1]))]
  
  box_tbl <- data.frame(
    ME         = module_eigengenes[common_rows, top_mod_tratamento],
    tratamento = trait_tbl$tratamento
  )
  
  y_lims <- box_tbl %>%
    group_by(tratamento) %>%
    summarise(q1 = quantile(ME, 0.25), q3 = quantile(ME, 0.75), iqr = IQR(ME)) %>%
    summarise(ymin = min(q1 - 1.5 * iqr), ymax = max(q3 + 1.5 * iqr))
  
  p_box <- ggplot(box_tbl, aes(x = tratamento, y = ME, fill = tratamento)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.6, color = "grey20") +
    coord_cartesian(ylim = c(y_lims$ymin * 1.2, y_lims$ymax * 1.2)) +
    scale_fill_manual(values = c("DMSO" = "#A8D5A2", "RA" = "#2D6A4F")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      title    = paste("Módulo", top_mod_tratamento, "por tratamento"),
      subtitle = "Outliers fora do intervalo IQR×1.5 omitidos do zoom",
      x        = NULL,
      y        = "Module Eigengene (ME)"
    )
  p_box
}

##############################################################################
# VST a partir do DESeq2

vst_data <- vst(dds, blind = TRUE)
normalized_counts <- assay(vst_data)

# Mantem genes protein-coding
clean_ids <- sub("\\..*$", "", rownames(normalized_counts))
rownames(normalized_counts) <- clean_ids

keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
gene_info <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys,
  columns = c("ENSEMBL", "GENETYPE"),
  keytype = "ENSEMBL"
)

protein_coding_genes <- gene_info$ENSEMBL[gene_info$GENETYPE == "protein-coding"]
protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes)]

valid_genes <- intersect(rownames(normalized_counts), protein_coding_genes)
normalized_counts <- normalized_counts[valid_genes, , drop = FALSE]

# Top 2000 genes mais variaveis
gene_var <- matrixStats::rowVars(normalized_counts)
top_n <- min(2000, length(gene_var))
top_idx <- order(gene_var, decreasing = TRUE)[1:top_n]
expr_top <- normalized_counts[top_idx, , drop = FALSE]
expr_top <- expr_top[, meta$código, drop = FALSE]

se_filtered <- SummarizedExperiment(
  assays = list(counts = expr_top),
  colData = DataFrame(código = colnames(expr_top), row.names = colnames(expr_top)),
  rowData = DataFrame(gene = rownames(expr_top), row.names = rownames(expr_top))
)

lioness_cor_optimized <- function(se) {
  X <- assay(se)
  n <- ncol(X)
  S <- rowSums(X)
  SS <- rowSums(X ^ 2)
  SP <- X %*% t(X)
  
  num_full <- SP - outer(S, S) / n
  den_full <- sqrt((SS - S ^ 2 / n) %o% (SS - S ^ 2 / n))
  C_full <- num_full / den_full
  
  ut_idx <- which(upper.tri(C_full), arr.ind = TRUE)
  genes <- rownames(X)
  reg_names <- genes[ut_idx[, 1]]
  tar_names <- genes[ut_idx[, 2]]
  v_full <- C_full[ut_idx]
  
  out <- matrix(NA_real_, nrow = length(reg_names), ncol = n)
  colnames(out) <- colnames(X)
  
  for (q in seq_len(n)) {
    xq <- X[, q]
    S_q <- S - xq
    SS_q <- SS - xq ^ 2
    SP_q <- SP - xq %*% t(xq)
    
    num_q <- SP_q - outer(S_q, S_q) / (n - 1)
    den_q <- sqrt((SS_q - S_q ^ 2 / (n - 1)) %o% (SS_q - S_q ^ 2 / (n - 1)))
    C_q <- num_q / den_q
    
    out[, q] <- n * (v_full - C_q[ut_idx]) + C_q[ut_idx]
  }
  
  list(W = out, reg = reg_names, tar = tar_names)
}

lioness_out <- lioness_cor_optimized(se_filtered)
W <- lioness_out$W
reg_names <- lioness_out$reg
tar_names <- lioness_out$tar

# Matriz de grau por gene e amostra (threshold mean + 2*sd)
all_genes_in_network <- unique(c(reg_names, tar_names))

deg_list <- lapply(seq_len(ncol(W)), function(i) {
  col_w <- abs(W[, i])
  thr <- mean(col_w, na.rm = TRUE) + (2 * sd(col_w, na.rm = TRUE))
  keep <- !is.na(col_w) & col_w > thr
  
  active_genes <- c(reg_names[keep], tar_names[keep])
  gene_counts <- table(active_genes)
  
  sample_deg <- setNames(integer(length(all_genes_in_network)), all_genes_in_network)
  sample_deg[names(gene_counts)] <- as.integer(gene_counts)
  sample_deg
})

deg_mat <- as.data.frame(do.call(cbind, deg_list))
colnames(deg_mat) <- colnames(W)
rownames(deg_mat) <- all_genes_in_network

deg_mat[1:5, 1:5]

##############################################################################
# DIGs: t-test nos graus por gene (RA vs DMSO)

dig_meta <- meta %>%
  filter(proficiencia == "deficiente")

common_samples <- intersect(colnames(deg_mat), dig_meta$código)
deg_mat_sub <- deg_mat[, common_samples, drop = FALSE]
dig_meta <- dig_meta %>%
  filter(código %in% common_samples) %>%
  arrange(match(código, common_samples))

stopifnot(all(colnames(deg_mat_sub) == dig_meta$código))

sf_idx <- which(dig_meta$tratamento == "RA")
gc_idx <- which(dig_meta$tratamento == "DMSO")

dig_results_list <- lapply(seq_len(nrow(deg_mat_sub)), function(i) {
  gene_name <- rownames(deg_mat_sub)[i]
  deg_vec <- as.numeric(deg_mat_sub[i, ])
  
  group_sf <- deg_vec[sf_idx]
  group_gc <- deg_vec[gc_idx]
  
  if (sd(group_sf) == 0 && sd(group_gc) == 0) return(NULL)
  
  res <- tryCatch({
    t_test <- t.test(group_sf, group_gc, alternative = "two.sided")
    data.frame(
      gene = gene_name,
      mean_deg_sf = mean(group_sf),
      mean_deg_gc = mean(group_gc),
      diff_deg = mean(group_sf) - mean(group_gc),
      p_val = t_test$p.value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) return(NULL))
  
  res
})

dig_results <- bind_rows(dig_results_list) %>%
  mutate(padj = p.adjust(p_val, method = "BH")) %>%
  arrange(p_val)

# Mapeia Ensembl -> simbolo
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(dig_results$gene),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

dig_results <- dig_results %>%
  left_join(gene_map, by = c("gene" = "ENSEMBL")) %>%
  dplyr::rename(gene_symbol = SYMBOL)

DIGs <- dig_results %>% filter(padj < 0.1)

if (nrow(DIGs) == 0) {
  message("Usando p_val < 0.05 (exploratório)")
  DIGs <- dig_results %>% filter(p_val < 0.05)
}

cat("Genes selecionados:", nrow(DIGs), "\n")

top_genes <- DIGs %>% arrange(padj) %>% slice_head(n = 6) %>% pull(gene)

if (length(top_genes) > 0) {
  gene_labels <- dig_results %>%
    dplyr::filter(gene %in% top_genes) %>%
    mutate(gene_label = ifelse(is.na(gene_symbol) | gene_symbol == "",
                               gene, gene_symbol)) %>%
    distinct(gene, gene_label)
  
  plot_tbl <- as.data.frame(t(deg_mat_sub[top_genes, , drop = FALSE])) %>%
    rownames_to_column("código") %>%
    pivot_longer(cols = -código, names_to = "gene", values_to = "degree") %>%
    left_join(
      dig_meta %>% dplyr::select(código, tratamento, proficiencia),
      by = "código"
    ) %>%
    left_join(gene_labels, by = "gene")
  
  p_digs <- ggplot(plot_tbl, aes(x = tratamento, y = degree, fill = tratamento)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.6, color = "grey20") +
    scale_fill_manual(values = c("DMSO" = "#A8D5A2", "RA" = "#2D6A4F")) +
    facet_wrap(~ gene_label, scales = "free_y") +
    theme_bw() +
    labs(
      title    = "Top DIGs — deficiente: RA vs DMSO (FDR < 10%)",
      subtitle = paste(nrow(DIGs), "genes com padj < 0.1"),
      x        = NULL,
      y        = "Grau de conectividade"
    )
  
  p_digs
}

##############################################################################
# ssGSEA com Hallmark e regressao com conectividade media
hallmark_tbl <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_list <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)

gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(normalized_counts),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

expr_sym_tbl <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  left_join(gene_symbols, by = c("gene" = "ENSEMBL")) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

expr_sym_mat <- expr_sym_tbl %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

if (exists("ssgseaParam")) {
  param_names <- names(formals(ssgseaParam))
  if ("geneSets" %in% param_names) {
    ssgsea_params <- ssgseaParam(
      expr = expr_sym_mat,
      geneSets = hallmark_list
    )
  } else {
    ssgsea_params <- ssgseaParam(
      expr = expr_sym_mat,
      gset.idx.list = hallmark_list
    )
  }
  ssgsea_mat <- gsva(ssgsea_params)
} else {
  ssgsea_mat <- gsva(expr_sym_mat, hallmark_list, method = "ssgsea", kcdf = "Gaussian")
}

conn_vec <- colMeans(deg_mat_sub, na.rm = TRUE)
conn_tbl <- data.frame(código = names(conn_vec), connectivity = conn_vec) %>%
  left_join(dig_meta %>% dplyr::select(código, tratamento), by = "código")

ssgsea_sub <- ssgsea_mat[, conn_tbl$código, drop = FALSE]

reg_results <- lapply(rownames(ssgsea_sub), function(gs) {
  df <- conn_tbl
  df$score <- as.numeric(ssgsea_sub[gs, ])
  fit <- lm(connectivity ~ score + tratamento, data = df)
  coef_tbl <- summary(fit)$coefficients
  
  data.frame(
    geneset = gs,
    beta_score = coef_tbl["score", "Estimate"],
    p_score = coef_tbl["score", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
})

reg_tbl <- bind_rows(reg_results) %>%
  mutate(padj = p.adjust(p_score, method = "BH")) %>%
  arrange(padj)

head(reg_tbl, 10)

# Regressao: ssGSEA vs conectividade do gene OGN por amostra
ogn_map <- gene_symbols %>% dplyr::filter(SYMBOL == "OGN") %>% pull(ENSEMBL)
ogn_map <- ogn_map[!is.na(ogn_map)]

if ("OGN" %in% rownames(expr_sym_mat)) {
  
  ogn_conn <- as.numeric(expr_sym_mat["OGN", conn_tbl$código])
  
  top_sets <- reg_tbl %>%
    slice_head(n = 6) %>%
    pull(geneset)
  
  plot_ssgsea_tbl <- data.frame(
    código     = conn_tbl$código,
    tratamento = conn_tbl$tratamento,
    ogn_conn   = ogn_conn
  )
  
  plot_long <- lapply(top_sets, function(gs) {
    data.frame(
      código     = plot_ssgsea_tbl$código,
      tratamento = plot_ssgsea_tbl$tratamento,
      ogn_conn   = plot_ssgsea_tbl$ogn_conn,
      geneset    = gs,
      ssgsea     = as.numeric(ssgsea_sub[gs, plot_ssgsea_tbl$código])
    )
  }) %>% bind_rows()
  
  stats_tbl <- plot_long %>%
    group_by(geneset) %>%
    summarise(
      r2    = summary(lm(ssgsea ~ ogn_conn))$r.squared,
      p_val = summary(lm(ssgsea ~ ogn_conn))$coefficients["ogn_conn", "Pr(>|t|)"],
      .groups = "drop"
    ) %>%
    mutate(label = paste0("R²=", signif(r2, 2), "\np=", signif(p_val, 2)))
  
  p_ssgsea <- ggplot(plot_long, aes(x = ogn_conn, y = ssgsea, color = tratamento)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ geneset, scales = "free_y",
               labeller = labeller(geneset = function(x) {
                 x %>%
                   stringr::str_remove("HALLMARK_") %>%
                   stringr::str_replace_all("_", " ") %>%
                   stringr::str_wrap(width = 20)  # <- quebra em 20 caracteres
               })) +
    geom_text(
      data        = stats_tbl,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = 1.1, size = 3
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8)  # <- reduz fonte dos títulos dos painéis
    ) +
    labs(
      title = "ssGSEA vs expressão de OGN",
      x     = "Expressão normalizada de OGN",
      y     = "ssGSEA score"
    )
  
  print(p_ssgsea)
  
} else {
  message("OGN não encontrado em expr_sym_mat")
}