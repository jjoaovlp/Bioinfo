## ============================================================================
## Modulo 17 -- Hotspots Regulatorios
## ============================================================================
## Descricao:
##   Calcula o "occupancy score" de cada regiao do universo regulatorio
##   (numero de condicoes proteina+genotipo que a ocupam, a partir da matriz
##   binaria do Modulo 14), identifica as regioes ocupadas por multiplas
##   proteinas simultaneamente, anota genomicamente (ChIPseeker) e salva o
##   ranking.
##
## Entradas:
##   Dados/Metadata/occupancy_matrix.csv   (Modulo 14)
##   Arquivos/regulatory_universe.rds       (Modulo 13)
##
## Saidas:
##   Arquivos/hotspots/hotspots_ranking.csv
##   Logs/17_hotspots.log
##
## Dependencias:
##   00_setup.R
##   GenomicRanges, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene (Bioconductor)
##
## Funcoes definidas neste modulo:
##   compute_occupancy_score(), annotate_hotspots(), run_module_17()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GenomicRanges")
install_if_missing("ChIPseeker")
install_if_missing("TxDb.Hsapiens.UCSC.hg38.knownGene")
suppressMessages({
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

HOTSPOTS_DIR <- file.path(PROJECT_DIRS$arquivos, "hotspots")

## --- occupancy score --------------------------------------------------------------

#' Calcula o occupancy score de cada regiao (soma das colunas binarias da
#' matriz de ocupacao do Modulo 14 -- numero de condicoes proteina+genotipo
#' que ocupam aquela regiao).
compute_occupancy_score <- function(occupancy_df, condition_cols) {
  occupancy_df$occupancy_score <- rowSums(occupancy_df[, condition_cols, drop = FALSE])
  occupancy_df
}

## --- anotacao genomica dos hotspots -------------------------------------------------

#' Anota genomicamente as regioes com occupancy_score >= min_score (ocupadas
#' por 2+ proteinas/condicoes) via ChIPseeker::annotatePeak(). O universo
#' regulatorio (Modulo 13) so tem o identificador da regiao em names(gr),
#' nao numa coluna BED "name" -- por isso e' copiado para uma metadata
#' column explicita (`region_id`) antes de anotar, ja que names() se perde
#' ao converter o resultado do annotatePeak() para data.frame.
annotate_hotspots <- function(universe_gr, hotspot_ids, tss_region = c(-3000, 3000)) {
  hotspot_gr <- universe_gr[hotspot_ids]
  mcols(hotspot_gr)$region_id <- names(hotspot_gr)
  annotatePeak(hotspot_gr, tssRegion = tss_region,
               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
               annoDb = "org.Hs.eg.db", verbose = FALSE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 17: le a matriz de ocupacao (Modulo 14) e o universo
#' regulatorio (Modulo 13), calcula o occupancy score, anota genomicamente
#' as regioes com score >= min_score e salva o ranking completo em CSV.
run_module_17 <- function(min_score = 2) {
  log_message("17_hotspots", "Iniciando Modulo 17 -- Hotspots Regulatorios.")
  ensure_dir(HOTSPOTS_DIR)

  matrix_file <- file.path(PROJECT_DIRS$metadata, "occupancy_matrix.csv")
  validate_file_exists(matrix_file, "matriz de ocupacao (Modulo 14)")
  occupancy_df <- read.csv(matrix_file, stringsAsFactors = FALSE)

  reserved_cols <- c("Regiao", "chr", "start", "end")
  condition_cols <- setdiff(names(occupancy_df), reserved_cols)
  occupancy_df <- compute_occupancy_score(occupancy_df, condition_cols)
  occupancy_df <- occupancy_df[order(-occupancy_df$occupancy_score), ]

  universe_file <- file.path(PROJECT_DIRS$arquivos, "regulatory_universe.rds")
  validate_file_exists(universe_file, "universo regulatorio (Modulo 13)")
  universe_gr <- readRDS(universe_file)

  hotspot_ids <- occupancy_df$Regiao[occupancy_df$occupancy_score >= min_score]
  if (length(hotspot_ids) > 0) {
    hotspot_anno <- annotate_hotspots(universe_gr, hotspot_ids)
    anno_df <- as.data.frame(hotspot_anno)[, c("region_id", "SYMBOL", "annotation", "distanceToTSS")]
    names(anno_df)[1] <- "Regiao"
    occupancy_df <- merge(occupancy_df, anno_df, by = "Regiao", all.x = TRUE)
    occupancy_df <- occupancy_df[order(-occupancy_df$occupancy_score), ]
  } else {
    log_message("17_hotspots",
      sprintf("Nenhuma regiao com occupancy_score >= %d.", min_score), level = "WARN")
  }

  out_file <- file.path(HOTSPOTS_DIR, "hotspots_ranking.csv")
  write.csv(occupancy_df, out_file, row.names = FALSE)
  log_message("17_hotspots", sprintf(
    "Ranking de hotspots salvo em '%s' (%d regiao(oes) com score >= %d).",
    out_file, length(hotspot_ids), min_score
  ))
  log_message("17_hotspots", "Modulo 17 concluido.")
  invisible(occupancy_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_17() explicitamente (interativamente ou a partir de
## 22_master_pipeline.R).
