## ============================================================================
## Modulo 09 -- Anotacao (ChIPseeker)
## ============================================================================
## Descricao:
##   Anota os picos de cada amostra (promotores, exons, introns, regioes
##   distais/enhancer-like, intergenicas, TSS) usando
##   `ChIPseeker::annotatePeak()` contra `TxDb.Hsapiens.UCSC.hg38.knownGene`
##   (todos os 4 datasets sao humanos apos a substituicao registrada em
##   CLAUDE.md S9) e `org.Hs.eg.db` para mapear SYMBOL/ENTREZID/ENSEMBL.
##
## Entradas:
##   Dados/Peaks/<amostra>_peaks.*Peak   (Modulo 07)
##
## Saidas:
##   Arquivos/annotation/<amostra>_annotation.csv   (tabela completa por pico)
##   Arquivos/annotation/annotation_summary.csv     (% por categoria genomica, por amostra)
##   Figuras/annotation/<amostra>_annobar.png
##   Logs/09_annotation.log
##
## Dependencias:
##   00_setup.R (import_peaks())
##   ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db (Bioconductor)
##
## Funcoes definidas neste modulo:
##   annotate_sample_peaks(), summarize_annotation(), run_module_09()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("ChIPseeker")
install_if_missing("TxDb.Hsapiens.UCSC.hg38.knownGene")
install_if_missing("org.Hs.eg.db")
suppressMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

ANNOTATION_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "annotation")
ANNOTATION_FIG_DIR <- file.path(PROJECT_DIRS$figuras, "annotation")

TXDB_HG38 <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

## --- anotacao por amostra ------------------------------------------------------

#' Anota os picos de uma amostra com ChIPseeker::annotatePeak() e salva a
#' tabela completa (uma linha por pico, com categoria genomica e gene mais
#' proximo) em Arquivos/annotation/.
annotate_sample_peaks <- function(peak_file, sample_id, tss_region = c(-3000, 3000)) {
  peaks_gr <- import_peaks(peak_file)
  log_message("09_annotation", sprintf("Anotando %d pico(s) da amostra '%s'.", length(peaks_gr), sample_id))

  peak_anno <- annotatePeak(
    peaks_gr, tssRegion = tss_region, TxDb = TXDB_HG38,
    annoDb = "org.Hs.eg.db", verbose = FALSE
  )

  ensure_dir(ANNOTATION_ARQ_DIR)
  out_file <- file.path(ANNOTATION_ARQ_DIR, sprintf("%s_annotation.csv", sample_id))
  write.csv(as.data.frame(peak_anno), out_file, row.names = FALSE)
  log_message("09_annotation", sprintf("Anotacao de '%s' salva em '%s'.", sample_id, out_file))

  peak_anno
}

## --- resumo por categoria genomica -----------------------------------------------

#' Extrai o percentual de picos por categoria genomica (Promoter, Exon,
#' Intron, Downstream, Distal Intergenic -- usada aqui como proxy de
#' enhancer/regiao distal) de um objeto annotatePeak.
summarize_annotation <- function(peak_anno, sample_id) {
  stats <- as.data.frame(peak_anno@annoStat)
  stats$sample_id <- sample_id
  stats
}

#' Salva o grafico de barras da distribuicao genomica de uma amostra (PNG,
#' 300 dpi) em Figuras/annotation/.
save_annotation_plot <- function(peak_anno, sample_id, output_dir = ANNOTATION_FIG_DIR) {
  ensure_dir(output_dir)
  install_if_missing("ggplot2")
  p <- plotAnnoBar(peak_anno) + ggplot2::ggtitle(sprintf("Distribuicao genomica -- %s", sample_id))
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_annobar.png", sample_id)),
                   p, width = 8, height = 4, dpi = 300)
  invisible(TRUE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 09 para uma tabela de amostras (colunas: sample_id,
#' peak_file), salvando a tabela de anotacao completa, o resumo por
#' categoria genomica e a figura de barras de cada amostra.
run_module_09 <- function(samples_df, tss_region = c(-3000, 3000)) {
  log_message("09_annotation", "Iniciando Modulo 09 -- Anotacao.")
  ensure_dir(ANNOTATION_ARQ_DIR)

  summaries <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    peak_anno <- annotate_sample_peaks(row$peak_file, row$sample_id, tss_region = tss_region)
    save_annotation_plot(peak_anno, row$sample_id)
    summarize_annotation(peak_anno, row$sample_id)
  })
  summary_df <- do.call(rbind, summaries)
  out_file <- file.path(ANNOTATION_ARQ_DIR, "annotation_summary.csv")
  write.csv(summary_df, out_file, row.names = FALSE)
  log_message("09_annotation", sprintf("Resumo de anotacao salvo em '%s'.", out_file))
  log_message("09_annotation", "Modulo 09 concluido.")
  invisible(summary_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_09(samples_df) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R).
