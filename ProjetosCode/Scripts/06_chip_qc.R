## ============================================================================
## Modulo 06 -- QC especifico de ChIP (ChIPQC)
## ============================================================================
## Descricao:
##   Roda controle de qualidade especifico de ChIP-seq por amostra e em lote,
##   usando o pacote `ChIPQC` (substitui deepTools -- ver CLAUDE.md S5.2/S9:
##   sem dependencia de PATH, mesma finalidade). Produz metricas e figuras de
##   fingerprint (SSD), fragment size (cross-coverage), coverage e correlacao
##   entre amostras, alem de PCA quando ha mais de uma amostra.
##
##   Este modulo roda ANTES do peak calling (Modulo 07); por isso FRiP fica
##   NA ate as amostras serem re-analisadas com peaks (ver run_chipqc_sample()
##   com o argumento `peaks`, opcional). As demais metricas (SSD/fingerprint,
##   fragment length, coverage, correlacao) nao dependem de peaks.
##
## Entradas:
##   Dados/BAM/<amostra>.filtered.bam(.bai)   (Modulo 05)
##   Dados/Peaks/<amostra>_peaks.*            (opcional, Modulo 07, para FRiP)
##
## Saidas:
##   Figuras/chip_qc/<amostra>_fragmentsize.png
##   Figuras/chip_qc/<amostra>_fingerprint.png
##   Figuras/chip_qc/correlation_heatmap.png   (so com >=2 amostras)
##   Figuras/chip_qc/pca.png                   (so com >=2 amostras)
##   Arquivos/chip_qc/chipqc_metrics.csv
##   Logs/06_chip_qc.log
##
## Dependencias:
##   00_setup.R
##   05_filtering.R (download_encode_blacklist(), reaproveitada para o
##     ponto "Post_Blacklist" do plotSSD())
##   ChIPQC (Bioconductor)
##
## Funcoes definidas neste modulo:
##   run_chipqc_sample(), save_sample_qc_plots(), run_chipqc_batch(),
##   save_batch_qc_plots(), run_module_06()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
source(here::here("Scripts", "05_filtering.R"))
install_if_missing("ChIPQC")
suppressMessages(library(ChIPQC))

CHIPQC_FIG_DIR <- file.path(PROJECT_DIRS$figuras, "chip_qc")
CHIPQC_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "chip_qc")

## --- QC por amostra ------------------------------------------------------------

#' Roda ChIPQCsample() para um BAM. `peaks` e `blacklist_gr` sao opcionais
#' (NULL antes do Modulo 07 rodar) -- sem peaks, FRiP fica NA; sem
#' blacklist, ChIPQC nao filtra regioes conhecidas de artefato na propria
#' analise de QC (o BAM ja deveria ter sido filtrado no Modulo 05).
run_chipqc_sample <- function(bam_file, sample_id, peaks = NULL, blacklist_gr = NULL) {
  validate_file_exists(bam_file, "BAM filtrado")
  log_message("06_chip_qc", sprintf("Rodando ChIPQCsample para '%s'.", sample_id))
  ChIPQCsample(
    reads = bam_file, peaks = peaks, annotation = NULL,
    blacklist = blacklist_gr, runCrossCor = TRUE, verboseT = FALSE
  )
}

#' Salva as figuras de fragment size (plotCC) e fingerprint/SSD (plotSSD)
#' de uma amostra em Figuras/chip_qc/, em PNG a 300 dpi.
save_sample_qc_plots <- function(qc_sample, sample_id, output_dir = CHIPQC_FIG_DIR) {
  ensure_dir(output_dir)
  install_if_missing("ggplot2")

  cc_plot <- plotCC(qc_sample) +
    ggplot2::ggtitle(sprintf("Fragment size / cross-coverage -- %s", sample_id))
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_fragmentsize.png", sample_id)),
                   cc_plot, width = 7, height = 5, dpi = 300)

  ssd_plot <- plotSSD(qc_sample) +
    ggplot2::ggtitle(sprintf("Fingerprint (SSD) -- %s", sample_id))
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_fingerprint.png", sample_id)),
                   ssd_plot, width = 7, height = 5, dpi = 300)

  log_message("06_chip_qc", sprintf("Figuras de QC salvas para '%s' em '%s'.", sample_id, output_dir))
  invisible(TRUE)
}

## --- QC em lote (correlacao, PCA) ------------------------------------------------

#' Monta a folha de amostras no formato esperado por ChIPQC() (estilo
#' DiffBind) a partir de samples_df (colunas: sample_id, bam, factor,
#' condition, replicate, peaks [opcional]) e roda ChIPQC() para o lote.
#'
#' Forca execucao serial (BiocParallel::SerialParam()) antes de chamar
#' ChIPQC(): os workers SNOW do BiocParallel (paralelismo default no
#' Windows) rodam em processos-filho que nao carregam automaticamente o
#' namespace de GenomeInfoDb, e ChIPQC() usa `seqlevels<-` internamente --
#' isso derruba a execucao inteira com "nao foi possivel encontrar a funcao
#' 'seqlevels<-'" assim que o lote de amostras e' processado em paralelo
#' (visto na pratica com 19 amostras). Rodar serial evita o problema (mais
#' lento, mas roda no processo principal que ja tem o namespace carregado).
run_chipqc_batch <- function(samples_df, blacklist_gr = NULL) {
  install_if_missing("BiocParallel")
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  sample_sheet <- data.frame(
    SampleID = samples_df$sample_id,
    Factor = samples_df$factor,
    Condition = samples_df$condition,
    Replicate = samples_df$replicate,
    bamReads = samples_df$bam,
    Peaks = if ("peaks" %in% names(samples_df)) samples_df$peaks else NA,
    stringsAsFactors = FALSE
  )
  log_message("06_chip_qc", sprintf("Rodando ChIPQC em lote para %d amostra(s).", nrow(sample_sheet)))
  ChIPQC(sample_sheet, annotation = NULL, blacklist = blacklist_gr)
}

#' Salva heatmap de correlacao e PCA entre amostras (so faz sentido com 2+
#' amostras) em Figuras/chip_qc/.
save_batch_qc_plots <- function(chipqc_exp, output_dir = CHIPQC_FIG_DIR) {
  ensure_dir(output_dir)
  install_if_missing("ggplot2")

  corr_plot <- plotCorHeatmap(chipqc_exp)
  ggplot2::ggsave(file.path(output_dir, "correlation_heatmap.png"),
                   corr_plot, width = 7, height = 6, dpi = 300)

  pca_plot <- plotPrincomp(chipqc_exp)
  ggplot2::ggsave(file.path(output_dir, "pca.png"), pca_plot, width = 7, height = 6, dpi = 300)

  log_message("06_chip_qc", sprintf("Figuras de correlacao/PCA salvas em '%s'.", output_dir))
  invisible(TRUE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 06 completo: QC por amostra (fragment size, fingerprint)
#' e, se houver 2+ amostras, correlacao e PCA em lote. Salva as metricas
#' numericas (FRiP, SSD, fragment length estimado) em CSV.
#'
#' Se `blacklist_gr` nao for informado, carrega a blacklist ENCODE do
#' `genome_build` (mesma usada no Modulo 05) para que plotSSD() produza o
#' ponto "Post_Blacklist" alem do "Pre_Blacklist" -- sem isso, o ChIPQC nao
#' tem contra o que comparar e o grafico de fingerprint fica com metade da
#' legenda sem dado (o Modulo 05 ja filtra os BAMs, entao o efeito aqui e'
#' mais para completar o grafico do que para remover reads novos).
run_module_06 <- function(samples_df, blacklist_gr = NULL, genome_build = "hg38") {
  log_message("06_chip_qc", "Iniciando Modulo 06 -- QC especifico de ChIP.")
  ensure_dir(CHIPQC_ARQ_DIR)
  if (is.null(blacklist_gr)) {
    blacklist_gr <- download_encode_blacklist(genome_build)
  }

  qc_samples <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    peaks <- if ("peaks" %in% names(row) && !is.na(row$peaks)) row$peaks else NULL
    qc <- run_chipqc_sample(row$bam, row$sample_id, peaks = peaks, blacklist_gr = blacklist_gr)
    save_sample_qc_plots(qc, row$sample_id)
    qc
  })
  names(qc_samples) <- samples_df$sample_id

  metrics_df <- do.call(rbind, lapply(seq_along(qc_samples), function(i) {
    qc <- qc_samples[[i]]
    data.frame(
      sample_id = names(qc_samples)[i],
      frip = tryCatch(frip(qc), error = function(e) NA_real_),
      ssd = tryCatch(ChIPQC::QCmetrics(qc)[["SSD"]], error = function(e) NA_real_),
      fragment_length = tryCatch(fragmentlength(qc), error = function(e) NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  out_file <- file.path(CHIPQC_ARQ_DIR, "chipqc_metrics.csv")
  write.csv(metrics_df, out_file, row.names = FALSE)
  log_message("06_chip_qc", sprintf("Metricas de ChIP-QC salvas em '%s'.", out_file))

  if (nrow(samples_df) >= 2) {
    chipqc_exp <- run_chipqc_batch(samples_df, blacklist_gr = blacklist_gr)
    save_batch_qc_plots(chipqc_exp)
  } else {
    log_message("06_chip_qc",
      "Apenas 1 amostra -- pulando correlacao/PCA em lote (precisa de 2+).",
      level = "WARN")
  }

  log_message("06_chip_qc", "Modulo 06 concluido.")
  invisible(metrics_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_06(samples_df) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R). `samples_df` precisa das colunas sample_id, bam,
## factor, condition, replicate e, opcionalmente, peaks.
