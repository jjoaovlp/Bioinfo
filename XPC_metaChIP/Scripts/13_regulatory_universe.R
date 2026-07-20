## ============================================================================
## Modulo 13 -- Universo Regulatorio
## ============================================================================
## Descricao:
##   Combina os GRanges padronizados em hg38 (Modulo 12) de todas as
##   amostras/proteinas num unico conjunto de regioes nao sobrepostas
##   (`GenomicRanges::reduce()`), com um identificador unico por regiao.
##   Este e' o "universo regulatorio" usado pelos Modulos 14+ (matriz de
##   ocupacao, overlaps, hotspots, redes).
##
## Entradas:
##   Arquivos/granges_hg38/<amostra>.rds   (Modulo 12)
##
## Saidas:
##   Arquivos/regulatory_universe.rds
##   Arquivos/regulatory_universe.bed
##   Logs/13_regulatory_universe.log
##
## Dependencias:
##   00_setup.R
##   GenomicRanges, rtracklayer (Bioconductor)
##
## Funcoes definidas neste modulo:
##   build_regulatory_universe(), run_module_13()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GenomicRanges")
install_if_missing("rtracklayer")
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

#' Combina uma lista de GRanges (uma por amostra, todas ja em hg38) num
#' unico conjunto de regioes reduzidas e nao sobrepostas, com identificador
#' unico ("region_00001", etc.) por regiao.
build_regulatory_universe <- function(granges_list) {
  ## unlist(GRangesList(...)) -- nao do.call(c, ...) -- e' o idioma correto
  ## para combinar uma lista de GRanges em um so objeto GRanges: do.call(c,
  ## lapply(...)) nao dispara o dispatch S4 de c() para GRanges e devolve
  ## silenciosamente uma "list" comum, quebrando reduce() a seguir.
  all_gr <- unlist(GRangesList(lapply(granges_list, granges)), use.names = FALSE)
  universe <- GenomicRanges::reduce(all_gr)
  if (length(universe) == 0) {
    stop("Validacao falhou: universo regulatorio ficou vazio. Execucao interrompida.", call. = FALSE)
  }
  names(universe) <- sprintf("region_%05d", seq_along(universe))
  log_message("13_regulatory_universe", sprintf(
    "Universo regulatorio: %d regiao(oes) a partir de %d amostra(s).",
    length(universe), length(granges_list)
  ))
  universe
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 13: le todos os GRanges hg38 de Arquivos/granges_hg38/
#' (uma por sample_id em `sample_ids`), constroi o universo regulatorio e
#' salva como RDS e BED.
run_module_13 <- function(sample_ids) {
  log_message("13_regulatory_universe", "Iniciando Modulo 13 -- Universo Regulatorio.")
  granges_dir <- file.path(PROJECT_DIRS$arquivos, "granges_hg38")
  granges_list <- lapply(sample_ids, function(sid) {
    gr_file <- file.path(granges_dir, sprintf("%s.rds", sid))
    validate_file_exists(gr_file, sprintf("GRanges hg38 de '%s' (Modulo 12)", sid))
    readRDS(gr_file)
  })
  names(granges_list) <- sample_ids

  universe <- build_regulatory_universe(granges_list)

  out_rds <- file.path(PROJECT_DIRS$arquivos, "regulatory_universe.rds")
  out_bed <- file.path(PROJECT_DIRS$arquivos, "regulatory_universe.bed")
  saveRDS(universe, out_rds)
  rtracklayer::export(universe, out_bed, format = "BED")
  log_message("13_regulatory_universe", sprintf("Universo regulatorio salvo em '%s' e '%s'.", out_rds, out_bed))

  log_message("13_regulatory_universe", "Modulo 13 concluido.")
  invisible(universe)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_13(sample_ids) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R).
