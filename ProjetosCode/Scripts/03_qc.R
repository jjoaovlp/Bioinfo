## ============================================================================
## Modulo 03 -- Controle de Qualidade (FastQC + MultiQC)
## ============================================================================
## Descricao:
##   Roda FastQC por amostra sobre os FASTQ brutos, agrega os relatorios com
##   MultiQC, e interpreta automaticamente os resultados (qualidade por base,
##   conteudo de GC, duplicacao, adaptadores) a partir do summary.txt que o
##   FastQC gera por amostra -- sinalizando falhas (FAIL) em modulos criticos
##   sem interromper a execucao (essa decisao fica para o Modulo 21).
##
## Entradas:
##   Dados/FASTQ/<SRR>_{1,2}.fastq.gz  (baixados pelo Modulo 01)
##
## Saidas:
##   Arquivos/qc/fastqc/<amostra>_fastqc.{html,zip}
##   Arquivos/qc/multiqc/multiqc_report.html
##   Arquivos/qc/qc_summary.csv          (interpretacao agregada por amostra/modulo)
##   Logs/03_qc.log
##
## Dependencias:
##   00_setup.R
##   FastQC e MultiQC no PATH (verificados via check_external_tool())
##
## Funcoes definidas neste modulo:
##   run_fastqc(), run_multiqc(), parse_fastqc_summary(),
##   interpret_qc_flags(), run_module_03()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))

QC_DIR <- file.path(PROJECT_DIRS$arquivos, "qc")
QC_FASTQC_DIR <- file.path(QC_DIR, "fastqc")
QC_MULTIQC_DIR <- file.path(QC_DIR, "multiqc")

## Modulos do FastQC cujo FAIL/WARN e' relevante o suficiente para logar em
## nivel WARN (o restante e' logado em INFO). Ver CLAUDE.md para a lista
## completa reportada por amostra em qc_summary.csv.
CRITICAL_FASTQC_MODULES <- c(
  "Per base sequence quality", "Per sequence GC content",
  "Sequence Duplication Levels", "Adapter Content"
)

## --- FastQC --------------------------------------------------------------------

#' Roda FastQC sobre um vetor de arquivos FASTQ, gravando html+zip em
#' Arquivos/qc/fastqc/. Nao interrompe o modulo inteiro se FastQC nao
#' estiver instalado -- apenas avisa e devolve FALSE.
run_fastqc <- function(fastq_files, output_dir = QC_FASTQC_DIR) {
  if (!check_external_tool("fastqc")) {
    log_message("03_qc", "FastQC ausente do PATH -- pulando Modulo 03.", level = "WARN")
    return(FALSE)
  }
  missing_files <- fastq_files[!file.exists(fastq_files)]
  if (length(missing_files) > 0) {
    stop(sprintf(
      "Validacao falhou: %d arquivo(s) FASTQ ausente(s): %s. Execucao interrompida.",
      length(missing_files), paste(missing_files, collapse = ", ")
    ), call. = FALSE)
  }
  ensure_dir(output_dir)
  log_message("03_qc", sprintf("Rodando FastQC em %d arquivo(s).", length(fastq_files)))
  status <- system2("fastqc", args = c(fastq_files, "-o", output_dir, "--quiet"))
  if (status != 0) {
    log_message("03_qc", sprintf("FastQC terminou com codigo de saida %d.", status), level = "ERROR")
    return(FALSE)
  }
  TRUE
}

## --- MultiQC --------------------------------------------------------------------

#' Agrega os relatorios FastQC de um diretorio com MultiQC.
run_multiqc <- function(input_dir = QC_FASTQC_DIR, output_dir = QC_MULTIQC_DIR) {
  if (!check_external_tool("multiqc")) {
    log_message("03_qc", "MultiQC ausente do PATH -- pulando agregacao.", level = "WARN")
    return(FALSE)
  }
  validate_dir_exists(input_dir, "diretorio de relatorios FastQC")
  ensure_dir(output_dir)
  log_message("03_qc", sprintf("Rodando MultiQC sobre '%s'.", input_dir))
  status <- system2("multiqc", args = c(input_dir, "-o", output_dir, "--force", "--quiet"))
  if (status != 0) {
    log_message("03_qc", sprintf("MultiQC terminou com codigo de saida %d.", status), level = "ERROR")
    return(FALSE)
  }
  TRUE
}

## --- interpretacao automatica -----------------------------------------------------

#' Le o summary.txt (PASS/WARN/FAIL por modulo) de dentro do .zip que o
#' FastQC gera para uma amostra, sem exigir descompactacao manual previa.
parse_fastqc_summary <- function(fastqc_zip) {
  validate_file_exists(fastqc_zip, "arquivo .zip do FastQC")
  inner_dir <- sub("\\.zip$", "", basename(fastqc_zip))
  summary_path <- file.path(inner_dir, "summary.txt")
  con <- unz(fastqc_zip, summary_path)
  on.exit(close(con), add = TRUE)
  lines <- tryCatch(
    readLines(con),
    error = function(e) {
      log_message("03_qc",
        sprintf("Nao foi possivel ler summary.txt de '%s': %s", fastqc_zip, conditionMessage(e)),
        level = "WARN")
      character(0)
    }
  )
  if (length(lines) == 0) {
    return(data.frame(Status = character(0), Module = character(0), Sample = character(0)))
  }
  parts <- strsplit(lines, "\t")
  data.frame(
    Status = vapply(parts, `[`, character(1), 1),
    Module = vapply(parts, `[`, character(1), 2),
    Sample = vapply(parts, `[`, character(1), 3),
    stringsAsFactors = FALSE
  )
}

#' Agrega o summary.txt de todos os .zip de Arquivos/qc/fastqc/ em uma unica
#' tabela e sinaliza (log WARN) amostras com FAIL em modulos criticos.
interpret_qc_flags <- function(fastqc_dir = QC_FASTQC_DIR) {
  zips <- list.files(fastqc_dir, pattern = "_fastqc\\.zip$", full.names = TRUE)
  if (length(zips) == 0) {
    log_message("03_qc", "Nenhum .zip de FastQC encontrado para interpretar.", level = "WARN")
    return(invisible(NULL))
  }
  all_summaries <- do.call(rbind, lapply(zips, parse_fastqc_summary))

  critical_fails <- all_summaries[
    all_summaries$Status == "FAIL" & all_summaries$Module %in% CRITICAL_FASTQC_MODULES,
  ]
  if (nrow(critical_fails) > 0) {
    for (i in seq_len(nrow(critical_fails))) {
      log_message("03_qc", sprintf(
        "FAIL em modulo critico '%s' na amostra '%s'.",
        critical_fails$Module[i], critical_fails$Sample[i]
      ), level = "WARN")
    }
  }

  ensure_dir(QC_DIR)
  out_file <- file.path(QC_DIR, "qc_summary.csv")
  write.csv(all_summaries, out_file, row.names = FALSE)
  log_message("03_qc", sprintf("Interpretacao de QC salva em '%s'.", out_file))
  all_summaries
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 03 completo para um vetor de arquivos FASTQ: FastQC,
#' MultiQC e interpretacao automatica dos summaries.
run_module_03 <- function(fastq_files) {
  log_message("03_qc", "Iniciando Modulo 03 -- Controle de Qualidade.")
  ok <- run_fastqc(fastq_files)
  if (isTRUE(ok)) {
    run_multiqc()
    interpret_qc_flags()
  }
  log_message("03_qc", "Modulo 03 concluido.")
  invisible(ok)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_03(fastq_files) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R), passando os caminhos dos FASTQ a processar (tipicamente
## lidos do metadata padronizado do Modulo 02).
