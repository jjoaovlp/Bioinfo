## ============================================================================
## Modulo 21 -- Validacao Cientifica e Tecnica
## ============================================================================
## Descricao:
##   Bateria de checagens defensivas que qualquer modulo pode chamar antes
##   de prosseguir para uma etapa que dependa da integridade dos dados
##   anteriores: metadata (GSM duplicado, especie/genoma inconsistente,
##   replicas insuficientes), alinhamento (taxa de alinhamento, BAM sem
##   indice), picos (arquivo vazio, FRiP fora da faixa esperada), GRanges
##   (ranges invalidos, orientacao +/- invalida). Cada check devolve
##   status PASS/WARN/FAIL; `run_module_21()` agrega tudo num relatorio e
##   **interrompe a execucao** se qualquer check critico (FAIL) ocorrer --
##   nenhuma etapa prossegue silenciosamente apos uma falha que comprometa a
##   validade cientifica da analise.
##
## Entradas:
##   Dados/Metadata/chipseq_metadata.csv          (Modulo 02)
##   Arquivos/alignment/alignment_stats.csv        (Modulo 04, opcional)
##   Dados/Peaks/*.{narrowPeak,broadPeak}           (Modulo 07, opcional)
##   Arquivos/chip_qc/chipqc_metrics.csv           (Modulo 06, opcional)
##
## Saidas:
##   Arquivos/validation/validation_report.csv
##   Logs/21_validation.log
##
## Dependencias:
##   00_setup.R
##
## Funcoes definidas neste modulo:
##   check_no_duplicate_samples(), check_organism_consistency(),
##   check_genome_consistency(), check_replicate_counts(),
##   check_bam_has_index(), check_alignment_rate(), check_peaks_nonempty(),
##   check_frip(), check_granges_validity(), check_strand_orientation(),
##   run_module_21()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))

VALIDATION_DIR <- file.path(PROJECT_DIRS$arquivos, "validation")

#' Constroi um resultado de check padronizado.
check_result <- function(name, status = c("PASS", "WARN", "FAIL"), message = "") {
  status <- match.arg(status)
  list(check = name, status = status, message = message)
}

## --- checks de metadata --------------------------------------------------------

#' FAIL se houver GSM duplicado no metadata padronizado.
check_no_duplicate_samples <- function(metadata_df) {
  dups <- unique(metadata_df$GSM[duplicated(metadata_df$GSM)])
  if (length(dups) > 0) {
    return(check_result("no_duplicate_samples", "FAIL",
      sprintf("GSM duplicado(s): %s", paste(dups, collapse = ", "))))
  }
  check_result("no_duplicate_samples", "PASS")
}

#' FAIL se o metadata combinar mais de uma especie (nunca deve acontecer
#' apos a substituicao de dataset registrada em CLAUDE.md S9 -- funciona
#' como rede de seguranca contra regressao).
check_organism_consistency <- function(metadata_df, organism_col = "organism") {
  if (!organism_col %in% names(metadata_df)) {
    return(check_result("organism_consistency", "WARN",
      sprintf("Coluna '%s' ausente -- nao foi possivel checar.", organism_col)))
  }
  organisms <- unique(na.omit(metadata_df[[organism_col]]))
  if (length(organisms) > 1) {
    return(check_result("organism_consistency", "FAIL",
      sprintf("Mais de uma especie no metadata: %s. Nunca misturar especies (CLAUDE.md S9).",
              paste(organisms, collapse = ", "))))
  }
  check_result("organism_consistency", "PASS")
}

#' FAIL se `genome_build` tiver valor fora de {hg19, hg38} -- nunca assumir
#' montagem sem chain de liftOver cadastrada (Modulo 12).
check_genome_consistency <- function(samples_df, genome_col = "genome_build") {
  if (!genome_col %in% names(samples_df)) {
    return(check_result("genome_consistency", "WARN",
      sprintf("Coluna '%s' ausente -- nao foi possivel checar.", genome_col)))
  }
  invalid <- setdiff(unique(na.omit(samples_df[[genome_col]])), c("hg19", "hg38"))
  if (length(invalid) > 0) {
    return(check_result("genome_consistency", "FAIL",
      sprintf("genome_build invalido(s): %s (so hg19/hg38 tem chain cadastrada).",
              paste(invalid, collapse = ", "))))
  }
  check_result("genome_consistency", "PASS")
}

#' WARN (nao FAIL -- pode ser intencional, ex. ELK1/STAT1 occupancy-only)
#' se alguma combinacao Protein+Genotype tiver menos de `min_replicates`.
check_replicate_counts <- function(metadata_df, min_replicates = 2) {
  counts <- aggregate(GSM ~ Protein + Genotype, data = metadata_df, FUN = length)
  low <- counts[counts$GSM < min_replicates, ]
  if (nrow(low) > 0) {
    combos <- sprintf("%s/%s (n=%d)", low$Protein, low$Genotype, low$GSM)
    return(check_result("replicate_counts", "WARN",
      sprintf("Combinacoes com menos de %d replica(s): %s", min_replicates, paste(combos, collapse = "; "))))
  }
  check_result("replicate_counts", "PASS")
}

## --- checks de alinhamento --------------------------------------------------------

#' FAIL se um BAM nao tiver o indice .bai correspondente.
check_bam_has_index <- function(bam_file) {
  if (!file.exists(paste0(bam_file, ".bai"))) {
    return(check_result("bam_has_index", "FAIL",
      sprintf("BAM sem indice: '%s'.", bam_file)))
  }
  check_result("bam_has_index", "PASS")
}

#' FAIL se a taxa de alinhamento de alguma amostra estiver abaixo de
#' `min_rate` (%) -- taxa muito baixa compromete a validade dos picos.
check_alignment_rate <- function(stats_df, min_rate = 50) {
  low <- stats_df[!is.na(stats_df$mapping_rate_pct) & stats_df$mapping_rate_pct < min_rate, ]
  if (nrow(low) > 0) {
    return(check_result("alignment_rate", "FAIL",
      sprintf("Amostra(s) com taxa de alinhamento < %d%%: %s",
              min_rate, paste(low$sample_id, collapse = ", "))))
  }
  check_result("alignment_rate", "PASS")
}

## --- checks de picos/GRanges --------------------------------------------------------

#' FAIL se um arquivo de picos estiver vazio.
check_peaks_nonempty <- function(peak_file) {
  validate_file_exists(peak_file, "arquivo de picos")
  n_lines <- length(readLines(peak_file))
  if (n_lines == 0) {
    return(check_result("peaks_nonempty", "FAIL", sprintf("Arquivo de picos vazio: '%s'.", peak_file)))
  }
  check_result("peaks_nonempty", "PASS")
}

#' WARN se o FRiP estiver fora da faixa tipicamente aceita para ChIP-seq de
#' fator de transcricao/histona (ENCODE sugere >= 1%; nao interrompe porque
#' o limiar exato depende do alvo e do desenho experimental).
check_frip <- function(frip_value, min_frip = 0.01) {
  if (is.na(frip_value)) {
    return(check_result("frip", "WARN", "FRiP nao disponivel (picos ainda nao calculados?)."))
  }
  if (frip_value < min_frip) {
    return(check_result("frip", "WARN",
      sprintf("FRiP = %.4f abaixo do minimo sugerido (%.2f).", frip_value, min_frip)))
  }
  check_result("frip", "PASS")
}

#' FAIL se um GRanges tiver ranges invalidos (NA/largura <= 0) ou estiver
#' vazio -- reaproveita a mesma logica de validate_granges_consistency()
#' (Modulo 11), aqui exposta como check generico do Modulo 21.
check_granges_validity <- function(gr, label) {
  tryCatch({
    validate_granges_consistency(gr, label)
    check_result("granges_validity", "PASS")
  }, error = function(e) {
    check_result("granges_validity", "FAIL", conditionMessage(e))
  })
}

#' FAIL se algum valor de strand for diferente de "+", "-" ou "*" (erro de
#' orientacao).
check_strand_orientation <- function(gr) {
  install_if_missing("GenomicRanges")
  invalid <- setdiff(unique(as.character(GenomicRanges::strand(gr))), c("+", "-", "*"))
  if (length(invalid) > 0) {
    return(check_result("strand_orientation", "FAIL",
      sprintf("Valor(es) de strand invalido(s): %s.", paste(invalid, collapse = ", "))))
  }
  check_result("strand_orientation", "PASS")
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa a bateria de checks aplicaveis aos artefatos fornecidos (todos os
#' argumentos sao opcionais -- so roda os checks cujo insumo esteja
#' disponivel), salva o relatorio completo em CSV e **interrompe a
#' execucao** (stop()) se qualquer check tiver status FAIL, listando todos
#' os problemas encontrados na mensagem de erro.
run_module_21 <- function(metadata_df = NULL, alignment_stats_df = NULL, bam_files = NULL,
                            peak_files = NULL, frip_values = NULL) {
  log_message("21_validation", "Iniciando Modulo 21 -- Validacao Cientifica e Tecnica.")
  ensure_dir(VALIDATION_DIR)
  results <- list()

  if (!is.null(metadata_df)) {
    results <- c(results, list(
      check_no_duplicate_samples(metadata_df),
      check_organism_consistency(metadata_df),
      check_genome_consistency(metadata_df),
      check_replicate_counts(metadata_df)
    ))
  }
  if (!is.null(bam_files)) {
    results <- c(results, lapply(bam_files, check_bam_has_index))
  }
  if (!is.null(alignment_stats_df)) {
    results <- c(results, list(check_alignment_rate(alignment_stats_df)))
  }
  if (!is.null(peak_files)) {
    results <- c(results, lapply(peak_files, check_peaks_nonempty))
  }
  if (!is.null(frip_values)) {
    results <- c(results, lapply(frip_values, check_frip))
  }

  report_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
  out_file <- file.path(VALIDATION_DIR, "validation_report.csv")
  write.csv(report_df, out_file, row.names = FALSE)

  for (i in seq_len(nrow(report_df))) {
    if (report_df$status[i] == "WARN") {
      log_message("21_validation", sprintf("%s: %s", report_df$check[i], report_df$message[i]), level = "WARN")
    }
  }

  failures <- report_df[report_df$status == "FAIL", ]
  if (nrow(failures) > 0) {
    msg <- paste(sprintf("[%s] %s", failures$check, failures$message), collapse = "\n")
    log_message("21_validation", sprintf("Validacao falhou com %d problema(s) critico(s).", nrow(failures)),
                level = "ERROR")
    stop(sprintf(
      "Validacao cientifica/tecnica falhou -- execucao interrompida:\n%s", msg
    ), call. = FALSE)
  }

  log_message("21_validation", sprintf("Validacao concluida sem falhas criticas (relatorio em '%s').", out_file))
  log_message("21_validation", "Modulo 21 concluido.")
  invisible(report_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_21(...) explicitamente, a qualquer momento do pipeline
## (tipicamente apos cada modulo critico, e sempre antes do relatorio final
## no Modulo 22), passando os artefatos disponiveis naquele ponto.
