## ============================================================================
## Modulo 15 -- Sobreposicao entre Pares
## ============================================================================
## Descricao:
##   Calcula todas as intersecoes par a par entre as condicoes
##   (proteina+genotipo) do projeto: numero de regioes sobrepostas e indice
##   de Jaccard (intersecao/uniao, em pb), via `GenomicRanges::findOverlaps()`.
##   Gera tabela completa e heatmap.
##
## Entradas:
##   Arquivos/granges_hg38/<amostra>.rds   (Modulo 12)
##
## Saidas:
##   Arquivos/overlap/pairwise_overlap_table.csv
##   Figuras/overlap/jaccard_heatmap.png
##   Logs/15_pairwise_overlap.log
##
## Dependencias:
##   00_setup.R
##   GenomicRanges, ggplot2 (Bioconductor/CRAN)
##
## Funcoes definidas neste modulo:
##   compute_jaccard(), build_pairwise_overlap_table(), plot_jaccard_heatmap(),
##   run_module_15()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GenomicRanges")
suppressMessages(library(GenomicRanges))

OVERLAP_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "overlap")
OVERLAP_FIG_DIR <- file.path(PROJECT_DIRS$figuras, "overlap")

## --- indice de Jaccard --------------------------------------------------------------

#' Indice de Jaccard (em pares de base) entre dois GRanges: bp da
#' intersecao dividido pelo bp da uniao.
compute_jaccard <- function(gr1, gr2) {
  inter_bp <- sum(width(GenomicRanges::intersect(gr1, gr2)))
  union_bp <- sum(width(GenomicRanges::reduce(unlist(GRangesList(gr1, gr2), use.names = FALSE))))
  if (union_bp == 0) return(NA_real_)
  inter_bp / union_bp
}

## --- tabela de sobreposicao par a par -------------------------------------------------

#' Constroi a tabela completa de sobreposicao par a par: para cada par de
#' condicoes, numero de regioes sobrepostas e indice de Jaccard.
build_pairwise_overlap_table <- function(condition_granges) {
  names_v <- names(condition_granges)
  pairs <- combn(names_v, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    gr1 <- condition_granges[[pair[1]]]
    gr2 <- condition_granges[[pair[2]]]
    n_overlap <- length(subsetByOverlaps(gr1, gr2))
    data.frame(
      condition_A = pair[1], condition_B = pair[2],
      n_regions_A = length(gr1), n_regions_B = length(gr2),
      n_overlapping = n_overlap, jaccard = compute_jaccard(gr1, gr2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

## --- heatmap ------------------------------------------------------------------------

#' Monta a matriz simetrica de Jaccard e salva um heatmap (PNG, 300 dpi).
plot_jaccard_heatmap <- function(overlap_table, condition_names, output_dir = OVERLAP_FIG_DIR) {
  ensure_dir(output_dir)
  install_if_missing("ggplot2")
  n <- length(condition_names)
  jac_mat <- matrix(1, n, n, dimnames = list(condition_names, condition_names))
  for (i in seq_len(nrow(overlap_table))) {
    a <- overlap_table$condition_A[i]; b <- overlap_table$condition_B[i]
    jac_mat[a, b] <- overlap_table$jaccard[i]
    jac_mat[b, a] <- overlap_table$jaccard[i]
  }
  df <- as.data.frame(as.table(jac_mat))
  names(df) <- c("A", "B", "Jaccard")
  p <- ggplot2::ggplot(df, ggplot2::aes(A, B, fill = Jaccard)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Indice de Jaccard entre condicoes", x = NULL, y = NULL)
  ggplot2::ggsave(file.path(output_dir, "jaccard_heatmap.png"), p, width = 7, height = 6, dpi = 300)
  invisible(TRUE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 15: le os GRanges hg38 por condicao (proteina+genotipo)
#' e calcula todas as sobreposicoes par a par.
run_module_15 <- function(condition_sample_ids) {
  log_message("15_pairwise_overlap", "Iniciando Modulo 15 -- Sobreposicao entre Pares.")
  ensure_dir(OVERLAP_ARQ_DIR)
  granges_dir <- file.path(PROJECT_DIRS$arquivos, "granges_hg38")

  condition_granges <- lapply(condition_sample_ids, function(sids) {
    grs <- lapply(sids, function(sid) {
      gr_file <- file.path(granges_dir, sprintf("%s.rds", sid))
      validate_file_exists(gr_file, sprintf("GRanges hg38 de '%s' (Modulo 12)", sid))
      granges(readRDS(gr_file))
    })
    GenomicRanges::reduce(unlist(GRangesList(grs), use.names = FALSE))
  })

  overlap_table <- build_pairwise_overlap_table(condition_granges)
  out_file <- file.path(OVERLAP_ARQ_DIR, "pairwise_overlap_table.csv")
  write.csv(overlap_table, out_file, row.names = FALSE)
  log_message("15_pairwise_overlap", sprintf("Tabela de sobreposicao par a par salva em '%s'.", out_file))

  plot_jaccard_heatmap(overlap_table, names(condition_granges))
  log_message("15_pairwise_overlap", "Modulo 15 concluido.")
  invisible(overlap_table)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_15(condition_sample_ids) explicitamente (interativamente ou a
## partir de 22_master_pipeline.R).
