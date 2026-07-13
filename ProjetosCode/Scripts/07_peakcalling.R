## ============================================================================
## Modulo 07 -- Peak Calling (MACS3 via WSL)
## ============================================================================
## Descricao:
##   Chama picos com MACS3 (broad para XPC, narrow para ELK1/STAT1/STAT2 --
##   ver Scripts/config/peakcalling_config.R) para cada amostra de ChIP do
##   metadata padronizado (Modulo 02). MACS3 roda dentro do WSL (ver
##   CLAUDE.md S5.2/S9 -- bioconda nao publica build para Windows), atraves
##   das funcoes win_to_wsl_path()/run_macs3() definidas em 00_setup.R.
##
##   Decisao critica (CLAUDE.md S9.1): amostras sem input disponivel (as 6
##   XPC-KO de GSE214182) rodam MACS3 SEM controle experimental
##   (--nolambda), nunca com input de outro genotipo e nunca substituindo
##   por H3K4me3 nesse caso especifico. A disponibilidade de input e'
##   sempre lida da coluna Input do metadata -- nunca assumida.
##
## Entradas:
##   Dados/BAM/<amostra>.filtered.bam   (Modulo 05)
##   Dados/Metadata/chipseq_metadata.csv (Modulo 02, para saber Protein/Input)
##
## Saidas:
##   Dados/Peaks/<amostra>_peaks.{broadPeak,narrowPeak}
##   Arquivos/peakcalling/peakcalling_log.csv  (parametros efetivamente usados por amostra)
##   Logs/07_peakcalling.log
##
## Dependencias:
##   00_setup.R (win_to_wsl_path, run_macs3, check_macs3_wsl)
##   Scripts/config/peakcalling_config.R
##   MACS3 3.0.4 dentro do WSL (ambiente conda "chipseq")
##
## Funcoes definidas neste modulo:
##   determine_macs3_args(), call_peaks_macs3(), run_module_07()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
source(here::here("Scripts", "config", "peakcalling_config.R"))

PEAKCALLING_DIR <- file.path(PROJECT_DIRS$arquivos, "peakcalling")

## --- decisao de parametros por amostra ----------------------------------------------

#' Determina os argumentos MACS3 para uma amostra: tipo de pico (broad/narrow)
#' vem de PEAKCALLING_CONFIG por proteina; uso de --nolambda vem
#' exclusivamente de `has_input` (nunca assumido) -- ver CLAUDE.md S9.1.
determine_macs3_args <- function(protein, has_input) {
  if (!protein %in% names(PEAKCALLING_CONFIG)) {
    stop(sprintf(
      "Validacao falhou: proteina '%s' sem configuracao em PEAKCALLING_CONFIG. Execucao interrompida.",
      protein
    ), call. = FALSE)
  }
  cfg <- PEAKCALLING_CONFIG[[protein]]
  args <- c("-f", "BAM", "-g", MACS3_GSIZE, "-q", cfg$qvalue)
  if (identical(cfg$peak_type, "broad")) {
    args <- c(args, "--broad", "--broad-cutoff", cfg$broad_cutoff)
  }
  if (!has_input) {
    args <- c(args, "--nolambda")
    log_message("07_peakcalling",
      "Amostra sem input -- rodando MACS3 com --nolambda (ver CLAUDE.md S9.1).",
      level = "WARN")
  }
  list(peak_type = cfg$peak_type, args = args)
}

## --- chamada do MACS3 (via WSL) -----------------------------------------------------

#' Roda `macs3 callpeak` (dentro do WSL) para uma amostra. `input_bam` deve
#' ser NA/NULL quando a amostra genuinamente nao tem controle disponivel
#' (nunca substituir silenciosamente por outro genotipo).
call_peaks_macs3 <- function(sample_id, treatment_bam, input_bam, protein, output_dir = PROJECT_DIRS$peaks) {
  validate_file_exists(treatment_bam, "BAM de tratamento (ChIP)")
  has_input <- !is.null(input_bam) && !is.na(input_bam)
  if (has_input) validate_file_exists(input_bam, "BAM de input/controle")

  decision <- determine_macs3_args(protein, has_input)
  ensure_dir(output_dir)

  ## shQuote() e' essencial aqui -- win_to_wsl_path() preserva espacos do
  ## caminho original do Windows (ex. "/mnt/c/Users/Joao - PC/..."), que
  ## quebrariam em tokens separados no bash -lc do WSL sem aspas.
  args <- c(
    "callpeak",
    "-t", shQuote(win_to_wsl_path(treatment_bam)),
    if (has_input) c("-c", shQuote(win_to_wsl_path(input_bam))) else NULL,
    decision$args,
    "-n", sample_id,
    "--outdir", shQuote(win_to_wsl_path(output_dir))
  )

  log_message("07_peakcalling", sprintf(
    "MACS3 callpeak (%s, %s controle) para '%s'.",
    decision$peak_type, ifelse(has_input, "com", "sem"), sample_id
  ))
  run_macs3(args)

  peak_ext <- if (identical(decision$peak_type, "broad")) "broadPeak" else "narrowPeak"
  peak_file <- file.path(output_dir, sprintf("%s_peaks.%s", sample_id, peak_ext))
  validate_file_exists(peak_file, sprintf("arquivo de picos (%s)", peak_ext))

  data.frame(
    sample_id = sample_id, protein = protein, peak_type = decision$peak_type,
    has_input = has_input, peak_file = peak_file,
    n_peaks = length(readLines(peak_file)),
    stringsAsFactors = FALSE
  )
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 07 para uma tabela de amostras (colunas: sample_id,
#' protein, bam, input_bam [NA quando nao ha controle]). Salva um log com os
#' parametros efetivamente usados por amostra (broad/narrow, com/sem
#' --nolambda, numero de picos).
run_module_07 <- function(samples_df) {
  log_message("07_peakcalling", "Iniciando Modulo 07 -- Peak Calling.")
  if (!check_macs3_wsl()) {
    stop("Validacao falhou: MACS3 (via WSL) indisponivel. Execucao interrompida.", call. = FALSE)
  }
  ensure_dir(PEAKCALLING_DIR)

  results <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    input_bam <- if ("input_bam" %in% names(row)) row$input_bam else NA
    call_peaks_macs3(row$sample_id, row$bam, input_bam, row$protein)
  })
  results_df <- do.call(rbind, results)

  out_file <- file.path(PEAKCALLING_DIR, "peakcalling_log.csv")
  write.csv(results_df, out_file, row.names = FALSE)
  log_message("07_peakcalling", sprintf("Log de peak calling salvo em '%s'.", out_file))
  log_message("07_peakcalling", "Modulo 07 concluido.")
  invisible(results_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_07(samples_df) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R). `samples_df` deve vir do metadata padronizado do
## Modulo 02 apos o alinhamento/filtragem (Modulos 04-05).
