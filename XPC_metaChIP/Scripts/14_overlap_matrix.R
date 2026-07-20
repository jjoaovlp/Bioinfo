## ============================================================================
## Modulo 14 -- Matriz de Ocupacao Cromatinica
## ============================================================================
## Descricao:
##   Constroi a matriz binaria de ocupacao: uma linha por regiao do universo
##   regulatorio (Modulo 13), uma coluna por combinacao proteina+genotipo
##   efetivamente presente no projeto. O numero de colunas segue o desenho
##   real de cada proteina (CLAUDE.md S7/S9) -- nao ha coluna "ELK1 KO" ou
##   "STAT1 KO" porque esses genotipos nao existem nos dados (nunca inventar
##   colunas vazias para "completar" um layout hipotetico).
##
## Entradas:
##   Arquivos/regulatory_universe.rds        (Modulo 13)
##   Arquivos/granges_hg38/<amostra>.rds      (Modulo 12)
##
## Saidas:
##   Dados/Metadata/occupancy_matrix.csv   (regiao x proteina_genotipo, binario)
##   Logs/14_overlap_matrix.log
##
## Dependencias:
##   00_setup.R
##   GenomicRanges (Bioconductor)
##
## Funcoes definidas neste modulo:
##   build_occupancy_matrix(), run_module_14()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GenomicRanges")
suppressMessages(library(GenomicRanges))

## --- matriz de ocupacao --------------------------------------------------------

#' Constroi a matriz binaria regiao x condicao: `condition_granges` e' uma
#' lista nomeada (ex. "XPC_WT", "XPC_KO", "STAT2_WT", "STAT2_KO", "ELK1_WT",
#' "STAT1_WT") de GRanges hg38, uma por combinacao proteina+genotipo
#' efetivamente presente no projeto. Uma regiao e' "ocupada" (1) se qualquer
#' amostra daquela combinacao tem um pico sobreposto.
build_occupancy_matrix <- function(universe_gr, condition_granges) {
  occupancy <- vapply(condition_granges, function(gr) {
    as.integer(overlapsAny(universe_gr, gr))
  }, integer(length(universe_gr)))

  data.frame(
    Regiao = names(universe_gr),
    chr = as.character(seqnames(universe_gr)),
    start = start(universe_gr), end = end(universe_gr),
    occupancy,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 14: le o universo regulatorio (Modulo 13) e os GRanges
#' hg38 por combinacao proteina+genotipo (`condition_sample_ids`, uma lista
#' nomeada de vetores de sample_id -- os sample_ids de uma mesma condicao sao
#' combinados via union antes do overlap), salvando a matriz binaria em CSV.
run_module_14 <- function(condition_sample_ids) {
  log_message("14_overlap_matrix", "Iniciando Modulo 14 -- Matriz de Ocupacao Cromatinica.")
  universe_file <- file.path(PROJECT_DIRS$arquivos, "regulatory_universe.rds")
  validate_file_exists(universe_file, "universo regulatorio (Modulo 13)")
  universe_gr <- readRDS(universe_file)

  granges_dir <- file.path(PROJECT_DIRS$arquivos, "granges_hg38")
  condition_granges <- lapply(condition_sample_ids, function(sids) {
    grs <- lapply(sids, function(sid) {
      gr_file <- file.path(granges_dir, sprintf("%s.rds", sid))
      validate_file_exists(gr_file, sprintf("GRanges hg38 de '%s' (Modulo 12)", sid))
      granges(readRDS(gr_file))
    })
    ## unlist(GRangesList(...)), nao do.call(c, ...) -- ver nota em
    ## 13_regulatory_universe.R sobre o bug de dispatch S4.
    GenomicRanges::reduce(unlist(GRangesList(grs), use.names = FALSE))
  })

  matrix_df <- build_occupancy_matrix(universe_gr, condition_granges)
  out_file <- file.path(PROJECT_DIRS$metadata, "occupancy_matrix.csv")
  write.csv(matrix_df, out_file, row.names = FALSE)
  log_message("14_overlap_matrix", sprintf(
    "Matriz de ocupacao (%d regioes x %d condicoes) salva em '%s'.",
    nrow(matrix_df), length(condition_granges), out_file
  ))
  log_message("14_overlap_matrix", "Modulo 14 concluido.")
  invisible(matrix_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_14(condition_sample_ids) explicitamente (interativamente ou a
## partir de 22_master_pipeline.R). `condition_sample_ids` e' uma lista
## nomeada, ex. list(XPC_WT = c("GSM..."), XPC_KO = c("GSM..."), ...).
