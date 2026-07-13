## ============================================================================
## Modulo 03 -- Controle de Qualidade (ShortRead)
## ============================================================================
## Descricao:
##   Roda controle de qualidade por amostra sobre os FASTQ brutos usando
##   ShortRead::qa() (substitui FastQC+MultiQC -- ver CLAUDE.md S5.2/S9: sem
##   dependencia de PATH, mesma finalidade). Gera um relatorio HTML por
##   amostra via ShortRead::report() e interpreta automaticamente os
##   resultados (duplicacao, conteudo de base/GC, contaminacao por
##   adaptador) a partir do objeto QA, sinalizando problemas sem interromper
##   a execucao (essa decisao fica para o Modulo 21).
##
## Entradas:
##   Dados/FASTQ/<SRR>_{1,2}.fastq.gz  (baixados pelo Modulo 01)
##
## Saidas:
##   Arquivos/qc/<amostra>/index.html         (relatorio ShortRead por amostra)
##   Arquivos/qc/qc_summary.csv               (interpretacao agregada por amostra)
##   Logs/03_qc.log
##
## Dependencias:
##   00_setup.R
##   ShortRead (Bioconductor)
##
## Funcoes definidas neste modulo:
##   run_qa(), interpret_qa(), run_module_03()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("ShortRead")
suppressMessages(library(ShortRead))

QC_DIR <- file.path(PROJECT_DIRS$arquivos, "qc")

## Limiares usados na interpretacao automatica (ver interpret_qa()). Nao sao
## "PASS/WARN/FAIL" do FastQC (ShortRead nao reproduz esse rotulo), mas
## limiares equivalentes checados diretamente sobre as metricas da QA.
QC_THRESHOLDS <- list(
  max_pct_n = 5,              # % maximo de bases 'N' aceitavel
  min_unique_seq_pct = 50,    # % minimo de sequencias unicas (duplicacao)
  max_adapter_pct = 5         # % maximo de reads com contaminacao de adaptador
)

## --- QA por amostra ------------------------------------------------------------

#' Roda ShortRead::qa() sobre um ou mais arquivos FASTQ (tratados como uma
#' unica amostra, ex.: mate 1 + mate 2) e grava o relatorio HTML em
#' Arquivos/qc/<sample_id>/.
run_qa <- function(fastq_files, sample_id, output_dir = QC_DIR) {
  missing_files <- fastq_files[!file.exists(fastq_files)]
  if (length(missing_files) > 0) {
    stop(sprintf(
      "Validacao falhou: %d arquivo(s) FASTQ ausente(s): %s. Execucao interrompida.",
      length(missing_files), paste(missing_files, collapse = ", ")
    ), call. = FALSE)
  }
  log_message("03_qc", sprintf("Rodando QA (ShortRead) para amostra '%s'.", sample_id))
  qa_result <- qa(fastq_files, type = "fastq")

  sample_dir <- file.path(output_dir, sample_id)
  ensure_dir(sample_dir)
  report(qa_result, dest = sample_dir)
  log_message("03_qc", sprintf("Relatorio QA salvo em '%s'.", sample_dir))

  qa_result
}

## --- interpretacao automatica -----------------------------------------------------

#' Interpreta um objeto QA do ShortRead (percentual de N, duplicacao,
#' contaminacao de adaptador) e sinaliza (log WARN) quando algum limiar de
#' QC_THRESHOLDS e' ultrapassado. Devolve uma linha de resumo (data.frame).
interpret_qa <- function(qa_result, sample_id) {
  base_calls <- qa_result[["baseCalls"]]
  total_bases <- sum(base_calls)
  pct_n <- if (total_bases > 0) 100 * sum(base_calls[, "N"]) / total_bases else NA_real_

  seq_dist <- qa_result[["sequenceDistribution"]]
  total_reads <- sum(seq_dist$nOccurrences * seq_dist$nReads)
  unique_reads <- sum(seq_dist$nReads[seq_dist$nOccurrences == 1])
  pct_unique <- if (total_reads > 0) 100 * unique_reads / total_reads else NA_real_

  adapter <- qa_result[["adapterContamination"]]
  pct_adapter <- if (!is.null(adapter) && nrow(adapter) > 0) {
    100 * mean(adapter$contamination, na.rm = TRUE)
  } else {
    NA_real_
  }

  if (!is.na(pct_n) && pct_n > QC_THRESHOLDS$max_pct_n) {
    log_message("03_qc", sprintf("Amostra '%s': %.1f%% de bases N (limite %d%%).",
                                  sample_id, pct_n, QC_THRESHOLDS$max_pct_n), level = "WARN")
  }
  if (!is.na(pct_unique) && pct_unique < QC_THRESHOLDS$min_unique_seq_pct) {
    log_message("03_qc", sprintf("Amostra '%s': apenas %.1f%% de sequencias unicas (duplicacao alta).",
                                  sample_id, pct_unique), level = "WARN")
  }
  if (!is.na(pct_adapter) && pct_adapter > QC_THRESHOLDS$max_adapter_pct) {
    log_message("03_qc", sprintf("Amostra '%s': %.1f%% de contaminacao por adaptador (limite %d%%).",
                                  sample_id, pct_adapter, QC_THRESHOLDS$max_adapter_pct), level = "WARN")
  }

  data.frame(
    sample_id = sample_id, total_reads = total_reads, pct_n = pct_n,
    pct_unique_sequences = pct_unique, pct_adapter_contamination = pct_adapter,
    stringsAsFactors = FALSE
  )
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 03 para uma tabela de amostras (colunas: sample_id,
#' fastq1, fastq2 [opcional]), salvando o resumo interpretado em CSV.
run_module_03 <- function(samples_df) {
  log_message("03_qc", "Iniciando Modulo 03 -- Controle de Qualidade.")
  ensure_dir(QC_DIR)
  summaries <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    fastq_files <- if ("fastq2" %in% names(row) && !is.na(row$fastq2)) {
      c(row$fastq1, row$fastq2)
    } else {
      row$fastq1
    }
    qa_result <- run_qa(fastq_files, row$sample_id)
    interpret_qa(qa_result, row$sample_id)
  })
  summary_df <- do.call(rbind, summaries)
  out_file <- file.path(QC_DIR, "qc_summary.csv")
  write.csv(summary_df, out_file, row.names = FALSE)
  log_message("03_qc", sprintf("Resumo de QC salvo em '%s'.", out_file))
  log_message("03_qc", "Modulo 03 concluido.")
  invisible(summary_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_03(samples_df) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R), com uma tabela de amostras derivada do metadata
## padronizado do Modulo 02.
