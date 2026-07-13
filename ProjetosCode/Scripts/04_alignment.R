## ============================================================================
## Modulo 04 -- Alinhamento (Rbowtie2 + Rsamtools)
## ============================================================================
## Descricao:
##   Alinha os FASTQ contra um indice Bowtie2 do genoma de referencia usando
##   o pacote `Rbowtie2` (mesmo algoritmo Bowtie2 -- ver CLAUDE.md S5.2/S9 --
##   empacotado como pacote R, sem precisar configurar PATH manualmente),
##   converte/ordena/indexa o BAM com `Rsamtools`, e registra estatisticas de
##   alinhamento (total de reads e taxa de alinhamento, extraidas da saida do
##   proprio bowtie2 e de Rsamtools::countBam()).
##
## Entradas:
##   Dados/FASTQ/<SRR>_{1,2}.fastq.gz     (Modulo 01/03)
##   Indice Bowtie2 do genoma de referencia (ver build_bowtie2_index())
##
## Saidas:
##   Dados/BAM/<amostra>.sorted.bam(.bai)
##   Arquivos/alignment/alignment_stats.csv   (resumo de todas as amostras)
##   Logs/04_alignment.log
##
## Dependencias:
##   00_setup.R
##   Rbowtie2, Rsamtools (Bioconductor)
##
## Observacao (2026-07-13): Rbowtie2::bowtie2_build()/bowtie2_samtools()
## chamam um script Python (`python3`) internamente -- no Windows isso exige
## um `python3.exe` real no PATH (o instalador oficial do python.org so cria
## `python.exe`; foi criada uma copia `python3.exe` ao lado dela nesta
## maquina). Ver CLAUDE.md S5.2 para o registro dessa dependencia.
##
## Funcoes definidas neste modulo:
##   build_bowtie2_index(), align_sample(), parse_alignment_log(),
##   run_module_04()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("Rbowtie2")
install_if_missing("Rsamtools")
suppressMessages({
  library(Rbowtie2)
  library(Rsamtools)
})

ALIGNMENT_DIR <- file.path(PROJECT_DIRS$arquivos, "alignment")

## --- indice do genoma de referencia ------------------------------------------------

#' Constroi um indice Bowtie2 a partir de um FASTA de referencia, caso ainda
#' nao exista em index_prefix (nunca reconstroi um indice ja presente).
build_bowtie2_index <- function(reference_fasta, index_prefix, threads = 4) {
  validate_file_exists(reference_fasta, "FASTA de referencia")
  if (file.exists(paste0(index_prefix, ".1.bt2"))) {
    log_message("04_alignment", sprintf("Indice Bowtie2 ja existe em '%s', reaproveitando.", index_prefix))
    return(invisible(index_prefix))
  }
  ensure_dir(dirname(index_prefix))
  log_message("04_alignment", sprintf("Construindo indice Bowtie2 em '%s'.", index_prefix))
  bowtie2_build(references = reference_fasta, bt2Index = index_prefix,
                sprintf("--threads %d", threads), overwrite = TRUE)
  if (!file.exists(paste0(index_prefix, ".1.bt2")) && !file.exists(paste0(index_prefix, ".1.bt2l"))) {
    stop(sprintf(
      "Validacao falhou: indice Bowtie2 nao foi criado em '%s'. Execucao interrompida.",
      index_prefix
    ), call. = FALSE)
  }
  invisible(index_prefix)
}

## --- estatisticas de alinhamento a partir do log do bowtie2 ------------------------

#' Extrai total de reads e taxa de alinhamento (%) do log de texto que o
#' bowtie2 imprime ao final do alinhamento (capturado por bowtie2_samtools()).
parse_alignment_log <- function(bowtie2_log) {
  rate_line <- grep("overall alignment rate", bowtie2_log, value = TRUE)
  mapping_rate <- if (length(rate_line) > 0) {
    as.numeric(sub("^([0-9.]+)%.*", "\\1", rate_line[1]))
  } else {
    NA_real_
  }
  total_line <- grep("reads;", bowtie2_log, value = TRUE)
  total_reads <- if (length(total_line) > 0) {
    as.numeric(sub("^([0-9]+) reads;.*", "\\1", total_line[1]))
  } else {
    NA_real_
  }
  list(total_reads = total_reads, mapping_rate_pct = mapping_rate)
}

## --- orquestracao por amostra ---------------------------------------------------------

#' Alinha uma amostra do inicio ao fim: Rbowtie2 -> BAM -> ordenado+indexado
#' (Rsamtools). Devolve uma lista com o caminho do BAM final e as
#' estatisticas de alinhamento.
align_sample <- function(sample_id, fastq1, fastq2 = NULL, index_prefix, threads = 4) {
  validate_file_exists(fastq1, "FASTQ (mate 1 ou single-end)")
  if (!is.null(fastq2)) validate_file_exists(fastq2, "FASTQ (mate 2)")
  if (!file.exists(paste0(index_prefix, ".1.bt2")) && !file.exists(paste0(index_prefix, ".1.bt2l"))) {
    stop(sprintf(
      "Validacao falhou: indice Bowtie2 nao encontrado em '%s'. Rode build_bowtie2_index() primeiro.",
      index_prefix
    ), call. = FALSE)
  }
  log_message("04_alignment", sprintf("Alinhando amostra '%s' (Rbowtie2).", sample_id))

  raw_prefix <- file.path(tempdir(), sample_id)
  ## bowtie2_samtools() devolve (de forma invisivel) as linhas do log de
  ## alinhamento (stderr do bowtie2, gravado por ela em ".bowtie2.cerr.txt"
  ## no diretorio de trabalho) -- por isso o valor de retorno e' capturado
  ## diretamente, sem capture.output() (a saida vai para o console do SO,
  ## nao para a conexao de output do R).
  bowtie2_log <- bowtie2_samtools(
    bt2Index = index_prefix, output = raw_prefix, outputType = "bam",
    seq1 = fastq1, seq2 = fastq2, bamFile = NULL, overwrite = TRUE,
    sprintf("--threads %d", threads)
  )
  raw_bam <- paste0(raw_prefix, ".bam")
  validate_file_exists(raw_bam, "BAM gerado pelo Rbowtie2")

  bam_file <- file.path(PROJECT_DIRS$bam, sprintf("%s.sorted.bam", sample_id))
  ensure_dir(dirname(bam_file))
  sorted_prefix <- sub("\\.bam$", "", bam_file)
  sortBam(raw_bam, sorted_prefix, overwrite = TRUE)
  indexBam(bam_file)
  unlink(raw_bam)

  validate_file_exists(bam_file, "BAM ordenado")
  validate_file_exists(paste0(bam_file, ".bai"), "indice do BAM")

  stats <- parse_alignment_log(bowtie2_log)
  log_message("04_alignment", sprintf(
    "Amostra '%s': %s reads totais, taxa de alinhamento %s%%.",
    sample_id, stats$total_reads, stats$mapping_rate_pct
  ))
  c(list(sample_id = sample_id, bam = bam_file), stats)
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
               mapping_rate_pct = r$mapping_rate_pct, bam = r$bam,
               stringsAsFactors = FALSE)
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
