## ============================================================================
## Modulo 19 -- Redes Funcionais (PPI, GO, Reactome)
## ============================================================================
## Descricao:
##   Constroi tres redes funcionais a partir dos genes alvo/enriquecimentos
##   do projeto: rede de interacao proteina-proteina (PPI) via `STRINGdb`
##   para os genes anotados como alvo/proximos aos picos, e redes de
##   similaridade de termos GO/Reactome (`enrichplot::emapplot()`) a partir
##   dos resultados do Modulo 10.
##
## Entradas:
##   Arquivos/annotation/<amostra>_annotation.csv   (Modulo 09, coluna SYMBOL)
##   Resultados de enrichGO()/enrichPathway() do Modulo 10 (em memoria)
##
## Saidas:
##   Arquivos/network/ppi_network.graphml
##   Figuras/network/ppi_network.png
##   Figuras/network/go_network.png
##   Figuras/network/reactome_network.png
##   Logs/19_pathway_network.log
##
## Dependencias:
##   00_setup.R
##   STRINGdb, enrichplot, igraph (Bioconductor/CRAN)
##
## Funcoes definidas neste modulo:
##   build_ppi_network(), build_term_network(), run_module_19()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("STRINGdb")
install_if_missing("enrichplot")
install_if_missing("igraph")
suppressMessages({
  library(STRINGdb)
  library(igraph)
})

NETWORK_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "network")
NETWORK_FIG_DIR <- file.path(PROJECT_DIRS$figuras, "network")

## --- rede PPI (STRINGdb) --------------------------------------------------------

#' Constroi a rede de interacao proteina-proteina (STRING, especie humana)
#' para um vetor de SYMBOLs, devolvendo um objeto igraph. Genes sem
#' correspondencia no STRING sao descartados (nunca assumidos presentes).
build_ppi_network <- function(gene_symbols, score_threshold = 400) {
  sdb <- STRINGdb$new(version = "12", species = 9606, score_threshold = score_threshold)
  genes_df <- data.frame(gene = unique(gene_symbols), stringsAsFactors = FALSE)
  mapped <- sdb$map(genes_df, "gene", removeUnmappedRows = TRUE)
  if (nrow(mapped) == 0) {
    stop("Validacao falhou: nenhum gene mapeado no STRINGdb. Execucao interrompida.", call. = FALSE)
  }
  log_message("19_pathway_network", sprintf(
    "%d de %d genes mapeados no STRINGdb.", nrow(mapped), nrow(genes_df)
  ))
  g <- sdb$get_subnetwork(mapped$STRING_id)
  ## renomeia os nos de STRING_id para o SYMBOL original, para leitura humana
  id_to_symbol <- setNames(mapped$gene, mapped$STRING_id)
  igraph::V(g)$name <- id_to_symbol[igraph::V(g)$name]
  g
}

## --- redes de similaridade de termos (GO/Reactome) --------------------------------

#' Constroi e salva a rede de similaridade de termos (emapplot) de um
#' resultado de enriquecimento (enrichGO()/enrichPathway()), se houver
#' termos significativos suficientes (emapplot precisa de 2+ termos).
build_term_network <- function(enrich_result, label, output_dir = NETWORK_FIG_DIR) {
  install_if_missing("enrichplot")
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) < 2) {
    log_message("19_pathway_network",
      sprintf("%s: menos de 2 termos significativos -- pulando rede.", label), level = "WARN")
    return(invisible(NULL))
  }
  ensure_dir(output_dir)
  sim <- enrichplot::pairwise_termsim(enrich_result)
  p <- enrichplot::emapplot(sim, showCategory = 30)
  install_if_missing("ggplot2")
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_network.png", label)), p, width = 9, height = 9, dpi = 300)
  log_message("19_pathway_network", sprintf("Rede de %s salva em '%s'.", label, output_dir))
  invisible(p)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 19: constroi a rede PPI para os genes fornecidos e as
#' redes de similaridade de termos GO/Reactome a partir dos resultados do
#' Modulo 10 (ja calculados em memoria -- nao repete o enriquecimento aqui).
run_module_19 <- function(gene_symbols, go_result = NULL, reactome_result = NULL) {
  log_message("19_pathway_network", "Iniciando Modulo 19 -- Redes Funcionais.")
  ensure_dir(NETWORK_ARQ_DIR)
  ensure_dir(NETWORK_FIG_DIR)

  ppi_graph <- build_ppi_network(gene_symbols)
  igraph::write_graph(ppi_graph, file.path(NETWORK_ARQ_DIR, "ppi_network.graphml"), format = "graphml")

  install_if_missing("ggraph")
  install_if_missing("tidygraph")
  p <- ggraph::ggraph(tidygraph::as_tbl_graph(ppi_graph), layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.3) +
    ggraph::geom_node_point() +
    ggraph::geom_node_text(ggplot2::aes(label = name), size = 2, repel = TRUE) +
    ggplot2::theme_void() +
    ggplot2::ggtitle("Rede PPI (STRINGdb)")
  ggplot2::ggsave(file.path(NETWORK_FIG_DIR, "ppi_network.png"), p, width = 9, height = 9, dpi = 300)

  build_term_network(go_result, "go")
  build_term_network(reactome_result, "reactome")

  log_message("19_pathway_network", "Modulo 19 concluido.")
  invisible(ppi_graph)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_19(gene_symbols, go_result, reactome_result) explicitamente
## (interativamente ou a partir de 22_master_pipeline.R), com `go_result`/
## `reactome_result` vindos do Modulo 10.
