## ============================================================================
## Modulo 16 -- Sobreposicao Multipla
## ============================================================================
## Descricao:
##   Encontra as regioes compartilhadas simultaneamente por todas as
##   proteinas do projeto (XPC, ELK1, STAT1, STAT2), via
##   `Reduce(GenomicRanges::intersect, ...)` sobre um GRanges representativo
##   por proteina (tipicamente a condicao WT, que e' a unica presente para
##   as 4 -- ver CLAUDE.md S7).
##
## Entradas:
##   Arquivos/granges_hg38/<amostra>.rds   (Modulo 12)
##
## Saidas:
##   Arquivos/overlap/shared_regions_all_proteins.csv
##   Logs/16_multiple_overlap.log
##
## Dependencias:
##   00_setup.R
##   GenomicRanges (Bioconductor)
##
## Funcoes definidas neste modulo:
##   find_shared_regions(), run_module_16()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GenomicRanges")
suppressMessages(library(GenomicRanges))

OVERLAP_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "overlap")

## --- intersecao multipla --------------------------------------------------------

#' Reduz uma lista de GRanges (uma por proteina) a regiao compartilhada por
#' TODAS elas simultaneamente, via Reduce(GenomicRanges::intersect, ...).
find_shared_regions <- function(protein_granges) {
  if (length(protein_granges) < 2) {
    stop("Validacao falhou: e' preciso pelo menos 2 proteinas para intersecao multipla. Execucao interrompida.",
         call. = FALSE)
  }
  shared <- Reduce(GenomicRanges::intersect, protein_granges)
  log_message("16_multiple_overlap", sprintf(
    "%d regiao(oes) compartilhada(s) por todas as %d proteinas.",
    length(shared), length(protein_granges)
  ))
  shared
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 16: le o GRanges hg38 representativo de cada proteina
#' (`protein_sample_ids`, uma lista nomeada por proteina) e salva as regioes
#' compartilhadas por todas em CSV.
run_module_16 <- function(protein_sample_ids) {
  log_message("16_multiple_overlap", "Iniciando Modulo 16 -- Sobreposicao Multipla.")
  ensure_dir(OVERLAP_ARQ_DIR)
  granges_dir <- file.path(PROJECT_DIRS$arquivos, "granges_hg38")

  protein_granges <- lapply(protein_sample_ids, function(sids) {
    grs <- lapply(sids, function(sid) {
      gr_file <- file.path(granges_dir, sprintf("%s.rds", sid))
      validate_file_exists(gr_file, sprintf("GRanges hg38 de '%s' (Modulo 12)", sid))
      granges(readRDS(gr_file))
    })
    GenomicRanges::reduce(unlist(GRangesList(grs), use.names = FALSE))
  })

  shared <- find_shared_regions(protein_granges)
  out_file <- file.path(OVERLAP_ARQ_DIR, "shared_regions_all_proteins.csv")
  write.csv(as.data.frame(shared), out_file, row.names = FALSE)
  log_message("16_multiple_overlap", sprintf("Regioes compartilhadas salvas em '%s'.", out_file))
  log_message("16_multiple_overlap", "Modulo 16 concluido.")
  invisible(shared)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_16(protein_sample_ids) explicitamente (interativamente ou a
## partir de 22_master_pipeline.R). `protein_sample_ids` e' uma lista
## nomeada por proteina, ex. list(XPC=c("GSM..."), ELK1=c("GSM..."), ...).
