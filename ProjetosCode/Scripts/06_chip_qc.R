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

## Rotulos manuais para amostras sem entrada no metadata padronizado (ex.
## input ENCODE do ELK1, que nao tem GSM proprio -- ver CLAUDE.md S4 2026-07-18).
MANUAL_SAMPLE_LABELS <- c(
  ELK1input_ENCFF002ECM = "ELK1 Input (ENCODE, rep1)",
  ELK1input_ENCFF002ECL = "ELK1 Input (ENCODE, rep2)",
  ELK1input_pooled       = "ELK1 Input (ENCODE, pooled)"
)

#' Resolve um rotulo legivel "Proteina Genotipo" (ex. "XPC WT") para um
#' sample_id (GSM), lendo Protein/Genotype de chipseq_metadata.csv e, se nao
#' encontrado la (amostras fora do escopo original, ex. IRF9/H3K4me3 -- ver
#' CLAUDE.md S4 2026-07-18), de chipseq_metadata_filtered_out.csv. Amostras
#' sinteticas (sem GSM, ex. inputs ENCODE combinados) usam
#' MANUAL_SAMPLE_LABELS. Se nada for encontrado, devolve NA (chamador decide
#' o fallback).
sample_label <- function(sample_id) {
  if (sample_id %in% names(MANUAL_SAMPLE_LABELS)) {
    return(unname(MANUAL_SAMPLE_LABELS[sample_id]))
  }
  for (f in c("chipseq_metadata.csv", "chipseq_metadata_filtered_out.csv")) {
    path <- file.path(PROJECT_DIRS$metadata, f)
    if (!file.exists(path)) next
    meta <- read.csv(path, stringsAsFactors = FALSE)
    row <- meta[meta$GSM == sample_id, ]
    if (nrow(row) > 0) {
      return(trimws(paste(row$Protein[1], row$Genotype[1])))
    }
  }
  NA_character_
}

#' Monta o titulo de figura "<sample_id> (<rotulo>)" ou so "<sample_id>" se
#' o rotulo nao for encontrado (nunca falha por amostra desconhecida).
sample_title <- function(sample_id) {
  label <- sample_label(sample_id)
  if (is.na(label) || !nzchar(label)) sample_id else sprintf("%s (%s)", sample_id, label)
}

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

  title_suffix <- sample_title(sample_id)
  cc_plot <- plotCC(qc_sample) +
    ggplot2::ggtitle(sprintf("Fragment size / cross-coverage -- %s", title_suffix))
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_fragmentsize.png", sample_id)),
                   cc_plot, width = 7, height = 5, dpi = 300)

  ssd_plot <- plotSSD(qc_sample) +
    ggplot2::ggtitle(sprintf("Fingerprint (SSD) -- %s", title_suffix))
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
  chipqc_exp <- ChIPQC(sample_sheet, annotation = NULL, blacklist = blacklist_gr)
  ## Salva o objeto ChIPQCexperiment em disco -- computa-lo custa horas
  ## (recalcula ChIPQCsample de todas as amostras), entao guardar o RDS
  ## permite regenerar os graficos de correlacao/PCA (save_batch_qc_plots())
  ## depois sem recomputar tudo.
  ensure_dir(CHIPQC_ARQ_DIR)
  saveRDS(chipqc_exp, file.path(CHIPQC_ARQ_DIR, "chipqc_experiment.rds"))
  chipqc_exp
}

#' Salva heatmap de correlacao e PCA entre amostras (so faz sentido com 2+
#' amostras) em Figuras/chip_qc/.
#'
#' plotCorHeatmap()/plotPrincomp() do ChIPQC desenham em BASE GRAPHICS (nao
#' devolvem objeto ggplot), entao precisam ser salvos abrindo um dispositivo
#' png() e fechando com dev.off() -- usar ggplot2::ggsave() aqui gera imagem
#' EM BRANCO (ggsave captura o dispositivo ggplot vazio, nao o base). Bug
#' visto na pratica: correlation_heatmap.png e pca.png sairam identicos e
#' totalmente brancos (2026-07-17).
save_batch_qc_plots <- function(chipqc_exp, output_dir = CHIPQC_FIG_DIR) {
  ensure_dir(output_dir)

  save_base_plot <- function(file, plot_expr) {
    grDevices::png(file, width = 7, height = 6, units = "in", res = 300)
    on.exit(grDevices::dev.off(), add = TRUE)
    tryCatch(plot_expr, error = function(e) {
      log_message("06_chip_qc",
        sprintf("Falha ao gerar '%s': %s", basename(file), conditionMessage(e)),
        level = "WARN")
    })
  }

  save_base_plot(file.path(output_dir, "correlation_heatmap.png"), plotCorHeatmap(chipqc_exp))
  save_base_plot(file.path(output_dir, "pca.png"), plotPrincomp(chipqc_exp))

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
