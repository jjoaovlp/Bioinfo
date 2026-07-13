## ============================================================================
## Modulo 11 -- Conversao para GRanges
## ============================================================================
## Descricao:
##   Converte os picos de cada amostra para GRanges padronizados: nomes
##   unicos no formato "<amostra>_peak_<indice>", verificacao de
##   consistencia (sem ranges NA/largura zero, sem nomes duplicados, sem
##   conjunto vazio) antes de salvar. Este e' o formato canonico usado pelos
##   modulos 12+ (padronizacao de genoma, universo regulatorio, overlaps).
##
## Entradas:
##   Dados/Peaks/<amostra>_peaks.*Peak   (Modulo 07)
##
## Saidas:
##   Arquivos/granges/<amostra>.rds   (GRanges padronizado)
##   Logs/11_granges.log
##
## Dependencias:
##   00_setup.R (import_peaks())
##   GenomicRanges
##
## Funcoes definidas neste modulo:
##   standardize_peak_names(), validate_granges_consistency(),
##   peaks_to_granges(), run_module_11()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GenomicRanges")
suppressMessages(library(GenomicRanges))

GRANGES_DIR <- file.path(PROJECT_DIRS$arquivos, "granges")

## --- padronizacao de nomes ------------------------------------------------------

#' Atribui nomes unicos e legiveis a cada pico: "<amostra>_peak_0001", etc.
standardize_peak_names <- function(gr, sample_id) {
  names(gr) <- sprintf("%s_peak_%04d", sample_id, seq_along(gr))
  gr
}

## --- validacao de consistencia --------------------------------------------------

#' Verifica que o GRanges de uma amostra e' internamente consistente:
#' nao-vazio, sem ranges NA/largura zero, sem nomes duplicados. Interrompe a
#' execucao (nao so avisa) quando a inconsistencia compromete a validade
#' cientifica da etapa seguinte (universo regulatorio, overlaps).
validate_granges_consistency <- function(gr, sample_id) {
  if (length(gr) == 0) {
    stop(sprintf(
      "Validacao falhou: GRanges de '%s' esta vazio (nenhum pico). Execucao interrompida.",
      sample_id
    ), call. = FALSE)
  }
  if (any(is.na(start(gr))) || any(is.na(end(gr))) || any(width(gr) <= 0)) {
    stop(sprintf(
      "Validacao falhou: GRanges de '%s' tem ranges invalidos (NA ou largura <= 0). Execucao interrompida.",
      sample_id
    ), call. = FALSE)
  }
  if (any(duplicated(names(gr)))) {
    stop(sprintf(
      "Validacao falhou: nomes de pico duplicados em '%s'. Execucao interrompida.", sample_id
    ), call. = FALSE)
  }
  invisible(TRUE)
}

## --- conversao por amostra ---------------------------------------------------------

#' Importa os picos de uma amostra, padroniza os nomes e valida a
#' consistencia antes de devolver o GRanges final.
peaks_to_granges <- function(peak_file, sample_id) {
  gr <- import_peaks(peak_file)
  gr <- standardize_peak_names(gr, sample_id)
  validate_granges_consistency(gr, sample_id)
  log_message("11_granges", sprintf("'%s': %d pico(s) convertido(s) e validado(s).", sample_id, length(gr)))
  gr
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 11 para uma tabela de amostras (colunas: sample_id,
#' peak_file), salvando o GRanges padronizado de cada uma em
#' Arquivos/granges/<amostra>.rds.
run_module_11 <- function(samples_df) {
  log_message("11_granges", "Iniciando Modulo 11 -- Conversao para GRanges.")
  ensure_dir(GRANGES_DIR)
  granges_list <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    gr <- peaks_to_granges(row$peak_file, row$sample_id)
    saveRDS(gr, file.path(GRANGES_DIR, sprintf("%s.rds", row$sample_id)))
    gr
  })
  names(granges_list) <- samples_df$sample_id
  log_message("11_granges", "Modulo 11 concluido.")
  invisible(granges_list)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_11(samples_df) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R).
