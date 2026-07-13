## ============================================================================
## Modulo 10 -- Enriquecimento Funcional
## ============================================================================
## Descricao:
##   Roda enriquecimento funcional (GO, KEGG, Reactome, Hallmark) sobre os
##   genes anotados no Modulo 09 para cada proteina, usando
##   `clusterProfiler` (GO/KEGG/Hallmark via enricher) e `ReactomePA`
##   (Reactome). Todos os 4 datasets sao humanos (CLAUDE.md S9).
##
## Entradas:
##   Arquivos/annotation/<amostra>_annotation.csv   (Modulo 09; usa coluna geneId, o Entrez ID)
##
## Saidas:
##   Arquivos/enrichment/<proteina>_{go,kegg,reactome,hallmark}.csv
##   Figuras/enrichment/<proteina>_{go,kegg,reactome,hallmark}_dotplot.png
##   Logs/10_enrichment.log
##
## Dependencias:
##   00_setup.R
##   clusterProfiler, ReactomePA, msigdbr, org.Hs.eg.db (Bioconductor)
##
## Funcoes definidas neste modulo:
##   get_hallmark_term2gene(), run_enrichment_go(), run_enrichment_kegg(),
##   run_enrichment_reactome(), run_enrichment_hallmark(),
##   save_enrichment_result(), run_module_10()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("clusterProfiler")
install_if_missing("ReactomePA")
install_if_missing("msigdbr")
install_if_missing("org.Hs.eg.db")
suppressMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(org.Hs.eg.db)
})

ENRICHMENT_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "enrichment")
ENRICHMENT_FIG_DIR <- file.path(PROJECT_DIRS$figuras, "enrichment")

## --- conjuntos de genes Hallmark (MSigDB, via msigdbr) --------------------------------

#' Monta a tabela TERM2GENE (gene set -> ENTREZID) dos conjuntos Hallmark
#' (categoria "H" do MSigDB) para humano, usada por clusterProfiler::enricher().
get_hallmark_term2gene <- function() {
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
  unique(hallmark[, c("gs_name", "entrez_gene")])
}

## --- enriquecimentos ---------------------------------------------------------------

#' GO (Biological Process, por padrao) via clusterProfiler::enrichGO().
run_enrichment_go <- function(entrez_ids, ont = "BP") {
  enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = ont,
           pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
}

#' KEGG via clusterProfiler::enrichKEGG() (organismo humano, "hsa").
run_enrichment_kegg <- function(entrez_ids) {
  enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
}

#' Reactome via ReactomePA::enrichPathway().
run_enrichment_reactome <- function(entrez_ids) {
  ReactomePA::enrichPathway(gene = entrez_ids, organism = "human",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
}

#' Hallmark (MSigDB) via clusterProfiler::enricher() com TERM2GENE customizado.
run_enrichment_hallmark <- function(entrez_ids, term2gene = get_hallmark_term2gene()) {
  enricher(gene = entrez_ids, TERM2GENE = term2gene, pAdjustMethod = "BH", pvalueCutoff = 0.05)
}

## --- salvar resultado (tabela + dotplot) ----------------------------------------------

#' Salva a tabela completa de um resultado de enriquecimento e, se houver
#' termos significativos, um dotplot (PNG, 300 dpi).
save_enrichment_result <- function(enrich_result, protein, source_label,
                                    arq_dir = ENRICHMENT_ARQ_DIR, fig_dir = ENRICHMENT_FIG_DIR) {
  ensure_dir(arq_dir)
  out_file <- file.path(arq_dir, sprintf("%s_%s.csv", protein, source_label))
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    log_message("10_enrichment", sprintf("%s/%s: nenhum termo significativo.", protein, source_label),
                level = "WARN")
    write.csv(data.frame(), out_file, row.names = FALSE)
    return(invisible(FALSE))
  }
  write.csv(as.data.frame(enrich_result), out_file, row.names = FALSE)

  ensure_dir(fig_dir)
  install_if_missing("enrichplot")
  p <- tryCatch(enrichplot::dotplot(enrich_result, showCategory = 15), error = function(e) NULL)
  if (!is.null(p)) {
    install_if_missing("ggplot2")
    ggplot2::ggsave(file.path(fig_dir, sprintf("%s_%s_dotplot.png", protein, source_label)),
                     p, width = 8, height = 7, dpi = 300)
  }
  log_message("10_enrichment", sprintf("%s/%s salvo em '%s'.", protein, source_label, out_file))
  invisible(TRUE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 10 para uma proteina, a partir dos Entrez ID lidos da
#' tabela de anotacao do Modulo 09 (coluna `geneId` -- e' assim que
#' ChIPseeker::annotatePeak() nomeia o Entrez ID do gene mais proximo,
#' mesmo com annoDb="org.Hs.eg.db" preenchendo SYMBOL/ENSEMBL/GENENAME a
#' parte). Roda os 4 enriquecimentos (GO, KEGG, Reactome, Hallmark) e salva
#' tabela + dotplot de cada um.
run_module_10 <- function(protein, annotation_csv) {
  log_message("10_enrichment", sprintf("Iniciando Modulo 10 -- Enriquecimento (%s).", protein))
  validate_file_exists(annotation_csv, "tabela de anotacao (Modulo 09)")
  annotation_df <- read.csv(annotation_csv, stringsAsFactors = FALSE)
  if (!"geneId" %in% names(annotation_df)) {
    stop("Validacao falhou: coluna geneId (Entrez ID) ausente na tabela de anotacao. Execucao interrompida.",
         call. = FALSE)
  }
  entrez_ids <- unique(na.omit(as.character(annotation_df$geneId)))
  if (length(entrez_ids) == 0) {
    stop("Validacao falhou: nenhum ENTREZID valido para enriquecimento. Execucao interrompida.",
         call. = FALSE)
  }

  results <- list(
    go = run_enrichment_go(entrez_ids),
    kegg = run_enrichment_kegg(entrez_ids),
    reactome = run_enrichment_reactome(entrez_ids),
    hallmark = run_enrichment_hallmark(entrez_ids)
  )
  for (source_label in names(results)) {
    save_enrichment_result(results[[source_label]], protein, source_label)
  }
  log_message("10_enrichment", "Modulo 10 concluido.")
  invisible(results)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_10(protein, annotation_csv) explicitamente (interativamente ou
## a partir de 22_master_pipeline.R), uma vez por proteina.
