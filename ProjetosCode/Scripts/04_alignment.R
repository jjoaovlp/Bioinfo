## ============================================================================
## Modulo 04 -- Alinhamento (Bowtie2 + samtools)
## ============================================================================
## Descricao:
##   Alinha os FASTQ filtrados/brutos contra um indice Bowtie2 do genoma de
##   referencia, converte para BAM, ordena por coordenada, indexa e registra
##   estatisticas de alinhamento (samtools flagstat: total de reads, reads
##   mapeados, taxa de alinhamento). Cada etapa e' um passo discreto com
##   arquivo intermediario proprio (em vez de pipes de shell encadeados) para
##   manter o modulo portavel entre Windows/Linux/macOS.
##
## Entradas:
##   Dados/FASTQ/<SRR>_{1,2}.fastq.gz     (Modulo 01/03)
##   Indice Bowtie2 do genoma de referencia (ver build_bowtie2_index())
##
## Saidas:
##   Dados/BAM/<amostra>.sorted.bam(.bai)
##   Arquivos/alignment/<amostra>_flagstat.txt
##   Arquivos/alignment/alignment_stats.csv   (resumo de todas as amostras)
##   Logs/04_alignment.log
##
## Dependencias:
##   00_setup.R
##   Bowtie2 e samtools no PATH (verificados via check_external_tool())
##
## Funcoes definidas neste modulo:
##   build_bowtie2_index(), align_bowtie2(), sort_and_index_bam(),
##   parse_flagstat(), align_sample(), run_module_04()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))

ALIGNMENT_DIR <- file.path(PROJECT_DIRS$arquivos, "alignment")

## --- indice do genoma de referencia ------------------------------------------------

#' Constroi um indice Bowtie2 a partir de um FASTA de referencia, caso ainda
#' nao exista em index_prefix (nunca reconstroi um indice ja presente).
build_bowtie2_index <- function(reference_fasta, index_prefix) {
  if (!check_external_tool("bowtie2-build")) {
    stop("Validacao falhou: 'bowtie2-build' ausente do PATH. Execucao interrompida.",
         call. = FALSE)
  }
  validate_file_exists(reference_fasta, "FASTA de referencia")
  if (file.exists(paste0(index_prefix, ".1.bt2"))) {
    log_message("04_alignment", sprintf("Indice Bowtie2 ja existe em '%s', reaproveitando.", index_prefix))
    return(invisible(index_prefix))
  }
  ensure_dir(dirname(index_prefix))
  log_message("04_alignment", sprintf("Construindo indice Bowtie2 em '%s'.", index_prefix))
  status <- system2("bowtie2-build", args = c(reference_fasta, index_prefix))
  if (status != 0) {
    stop(sprintf("Falha ao construir indice Bowtie2 (codigo %d). Execucao interrompida.", status),
         call. = FALSE)
  }
  invisible(index_prefix)
}

## --- alinhamento ------------------------------------------------------------------

#' Alinha um par (ou single-end) de FASTQ contra `index_prefix` com Bowtie2,
#' devolvendo o caminho do SAM gerado. `fastq2` = NULL para single-end.
align_bowtie2 <- function(fastq1, fastq2 = NULL, index_prefix, output_sam, threads = 4) {
  if (!check_external_tool("bowtie2")) {
    stop("Validacao falhou: 'bowtie2' ausente do PATH. Execucao interrompida.", call. = FALSE)
  }
  validate_file_exists(fastq1, "FASTQ (mate 1 ou single-end)")
  if (!is.null(fastq2)) validate_file_exists(fastq2, "FASTQ (mate 2)")
  if (!file.exists(paste0(index_prefix, ".1.bt2"))) {
    stop(sprintf(
      "Validacao falhou: indice Bowtie2 nao encontrado em '%s'. Rode build_bowtie2_index() primeiro.",
      index_prefix
    ), call. = FALSE)
  }
  ensure_dir(dirname(output_sam))
  args <- c("-p", threads, "-x", index_prefix)
  args <- if (!is.null(fastq2)) c(args, "-1", fastq1, "-2", fastq2) else c(args, "-U", fastq1)
  args <- c(args, "-S", output_sam)
  log_message("04_alignment", sprintf("bowtie2 -> '%s'", output_sam))
  status <- system2("bowtie2", args = args)
  if (status != 0) {
    stop(sprintf("Falha no bowtie2 (codigo %d) para '%s'. Execucao interrompida.", status, output_sam),
         call. = FALSE)
  }
  invisible(output_sam)
}

## --- conversao, ordenacao e indexacao -------------------------------------------------

#' Converte SAM -> BAM, ordena por coordenada e indexa. Devolve o caminho do
#' BAM ordenado final.
sort_and_index_bam <- function(sam_file, output_bam, threads = 4) {
  if (!check_external_tool("samtools")) {
    stop("Validacao falhou: 'samtools' ausente do PATH. Execucao interrompida.", call. = FALSE)
  }
  validate_file_exists(sam_file, "arquivo SAM")
  ensure_dir(dirname(output_bam))
  unsorted_bam <- tempfile(fileext = ".bam")

  log_message("04_alignment", sprintf("samtools view (SAM -> BAM) '%s'", sam_file))
  st1 <- system2("samtools", args = c("view", "-@", threads, "-bS", "-o", unsorted_bam, sam_file))
  if (st1 != 0) stop(sprintf("Falha em 'samtools view' (codigo %d).", st1), call. = FALSE)

  log_message("04_alignment", sprintf("samtools sort -> '%s'", output_bam))
  st2 <- system2("samtools", args = c("sort", "-@", threads, "-o", output_bam, unsorted_bam))
  if (st2 != 0) stop(sprintf("Falha em 'samtools sort' (codigo %d).", st2), call. = FALSE)
  unlink(unsorted_bam)

  log_message("04_alignment", sprintf("samtools index '%s'", output_bam))
  st3 <- system2("samtools", args = c("index", output_bam))
  if (st3 != 0) stop(sprintf("Falha em 'samtools index' (codigo %d).", st3), call. = FALSE)

  validate_file_exists(output_bam, "BAM ordenado")
  validate_file_exists(paste0(output_bam, ".bai"), "indice do BAM")
  invisible(output_bam)
}

## --- estatisticas de alinhamento -----------------------------------------------------

#' Roda `samtools flagstat` sobre um BAM e devolve uma lista com total de
#' reads, reads mapeados e taxa de alinhamento (%), alem de gravar a saida
#' bruta em Arquivos/alignment/.
parse_flagstat <- function(bam_file, sample_id) {
  validate_file_exists(bam_file, "BAM")
  ensure_dir(ALIGNMENT_DIR)
  raw_output <- system2("samtools", args = c("flagstat", bam_file), stdout = TRUE)
  out_file <- file.path(ALIGNMENT_DIR, sprintf("%s_flagstat.txt", sample_id))
  writeLines(raw_output, out_file)

  total_line <- grep("in total", raw_output, value = TRUE)
  mapped_line <- grep(" mapped \\(", raw_output, value = TRUE)
  total_reads <- as.numeric(sub("^([0-9]+).*", "\\1", total_line[1]))
  mapped_reads <- as.numeric(sub("^([0-9]+).*", "\\1", mapped_line[1]))
  mapping_rate <- if (length(total_reads) == 1 && total_reads > 0) {
    round(100 * mapped_reads / total_reads, 2)
  } else {
    NA_real_
  }

  list(sample_id = sample_id, total_reads = total_reads,
       mapped_reads = mapped_reads, mapping_rate_pct = mapping_rate)
}

## --- orquestracao por amostra ---------------------------------------------------------

#' Alinha uma amostra do inicio ao fim: bowtie2 -> BAM ordenado+indexado ->
#' flagstat. Devolve uma lista com o caminho do BAM final e as estatisticas.
align_sample <- function(sample_id, fastq1, fastq2 = NULL, index_prefix, threads = 4) {
  log_message("04_alignment", sprintf("Alinhando amostra '%s'.", sample_id))
  sam_file <- file.path(tempdir(), sprintf("%s.sam", sample_id))
  bam_file <- file.path(PROJECT_DIRS$bam, sprintf("%s.sorted.bam", sample_id))

  align_bowtie2(fastq1, fastq2, index_prefix, sam_file, threads = threads)
  sort_and_index_bam(sam_file, bam_file, threads = threads)
  unlink(sam_file)
  stats <- parse_flagstat(bam_file, sample_id)

  log_message("04_alignment", sprintf(
    "Amostra '%s': %s reads totais, taxa de alinhamento %s%%.",
    sample_id, stats$total_reads, stats$mapping_rate_pct
  ))
  c(list(bam = bam_file), stats)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 04 para uma tabela de amostras (colunas esperadas:
#' sample_id, fastq1, fastq2 [opcional]) contra um mesmo index_prefix, e
#' salva o resumo de estatisticas de alinhamento em CSV.
run_module_04 <- function(samples_df, index_prefix, threads = 4) {
  log_message("04_alignment", "Iniciando Modulo 04 -- Alinhamento.")
  ensure_dir(ALIGNMENT_DIR)
  results <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    fastq2 <- if ("fastq2" %in% names(row) && !is.na(row$fastq2)) row$fastq2 else NULL
    align_sample(row$sample_id, row$fastq1, fastq2, index_prefix, threads = threads)
  })
  stats_df <- do.call(rbind, lapply(results, function(r) {
    data.frame(sample_id = r$sample_id, total_reads = r$total_reads,
               mapped_reads = r$mapped_reads, mapping_rate_pct = r$mapping_rate_pct,
               bam = r$bam, stringsAsFactors = FALSE)
  }))
  out_file <- file.path(ALIGNMENT_DIR, "alignment_stats.csv")
  write.csv(stats_df, out_file, row.names = FALSE)
  log_message("04_alignment", sprintf("Estatisticas de alinhamento salvas em '%s'.", out_file))
  log_message("04_alignment", "Modulo 04 concluido.")
  invisible(stats_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_04(samples_df, index_prefix) explicitamente (interativamente ou
## a partir de 22_master_pipeline.R).
