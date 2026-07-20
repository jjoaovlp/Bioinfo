## ============================================================================
## Modulo 18 -- Rede Bipartida Regulatoria (Proteina -> Regiao -> Gene)
## ============================================================================
## Descricao:
##   Constroi uma rede de tres niveis (Proteina -> Regiao Cromatinica ->
##   Gene) a partir da matriz de ocupacao (Modulo 14) e da anotacao
##   genomica dos hotspots (Modulo 17): uma aresta Proteina-Regiao para cada
##   condicao que ocupa a regiao, e uma aresta Regiao-Gene para o gene mais
##   proximo anotado. Calcula grau, betweenness, closeness e comunidades
##   (`igraph`/`tidygraph`), e exporta em formato Cytoscape (GraphML).
##
## Entradas:
##   Dados/Metadata/occupancy_matrix.csv     (Modulo 14)
##   Arquivos/hotspots/hotspots_ranking.csv  (Modulo 17, colunas Regiao/SYMBOL)
##
## Saidas:
##   Arquivos/network/bipartite_edges.csv
##   Arquivos/network/bipartite_metrics.csv
##   Arquivos/network/bipartite_network.graphml   (para importar no Cytoscape)
##   Figuras/network/bipartite_network.png
##   Logs/18_bipartite_network.log
##
## Dependencias:
##   00_setup.R
##   igraph, tidygraph, ggraph (CRAN)
##
## Funcoes definidas neste modulo:
##   build_bipartite_edges(), build_network_graph(), compute_network_metrics(),
##   export_cytoscape(), run_module_18()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("igraph")
install_if_missing("tidygraph")
install_if_missing("ggraph")
suppressMessages({
  library(igraph)
  library(tidygraph)
  library(ggraph)
})

NETWORK_ARQ_DIR <- file.path(PROJECT_DIRS$arquivos, "network")
NETWORK_FIG_DIR <- file.path(PROJECT_DIRS$figuras, "network")

## --- construcao das arestas -----------------------------------------------------

#' Monta a tabela de arestas de tres niveis: Proteina-Regiao (a partir da
#' matriz binaria de ocupacao) e Regiao-Gene (a partir da anotacao dos
#' hotspots, coluna SYMBOL). Regioes sem gene anotado (SYMBOL == NA) nao
#' geram aresta Regiao-Gene, mas continuam no grafo via a aresta
#' Proteina-Regiao.
build_bipartite_edges <- function(occupancy_df, gene_annotation_df) {
  reserved_cols <- c("Regiao", "chr", "start", "end", "occupancy_score",
                      "SYMBOL", "annotation", "distanceToTSS")
  condition_cols <- setdiff(names(occupancy_df), reserved_cols)

  protein_region_edges <- do.call(rbind, lapply(condition_cols, function(cond) {
    occupied <- occupancy_df[occupancy_df[[cond]] == 1, "Regiao"]
    if (length(occupied) == 0) return(NULL)
    data.frame(from = cond, to = occupied, edge_type = "Protein-Region", stringsAsFactors = FALSE)
  }))

  region_gene <- unique(gene_annotation_df[!is.na(gene_annotation_df$SYMBOL), c("Regiao", "SYMBOL")])
  region_gene_edges <- if (nrow(region_gene) > 0) {
    data.frame(from = region_gene$Regiao, to = region_gene$SYMBOL,
               edge_type = "Region-Gene", stringsAsFactors = FALSE)
  } else {
    NULL
  }

  rbind(protein_region_edges, region_gene_edges)
}

## --- grafo e metricas --------------------------------------------------------------

#' Constroi o grafo igraph (nao direcionado) a partir da tabela de arestas.
build_network_graph <- function(edges_df) {
  igraph::graph_from_data_frame(edges_df[, c("from", "to")], directed = FALSE)
}

#' Calcula grau, betweenness, closeness e comunidade (Louvain) de cada no.
compute_network_metrics <- function(g) {
  communities <- igraph::cluster_louvain(g)
  data.frame(
    node = igraph::V(g)$name,
    degree = igraph::degree(g),
    betweenness = igraph::betweenness(g),
    closeness = igraph::closeness(g),
    community = igraph::membership(communities),
    stringsAsFactors = FALSE
  )
}

#' Exporta o grafo em GraphML, formato lido nativamente pelo Cytoscape.
export_cytoscape <- function(g, path) {
  igraph::write_graph(g, path, format = "graphml")
  invisible(path)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 18 completo: monta as arestas, constroi o grafo, calcula
#' metricas, exporta para Cytoscape e salva uma figura da rede.
run_module_18 <- function(occupancy_csv, hotspots_csv) {
  log_message("18_bipartite_network", "Iniciando Modulo 18 -- Rede Bipartida Regulatoria.")
  ensure_dir(NETWORK_ARQ_DIR)
  validate_file_exists(occupancy_csv, "matriz de ocupacao (Modulo 14)")
  validate_file_exists(hotspots_csv, "ranking de hotspots (Modulo 17)")

  occupancy_df <- read.csv(occupancy_csv, stringsAsFactors = FALSE)
  hotspots_df <- read.csv(hotspots_csv, stringsAsFactors = FALSE)

  edges_df <- build_bipartite_edges(occupancy_df, hotspots_df)
  if (is.null(edges_df) || nrow(edges_df) == 0) {
    stop("Validacao falhou: nenhuma aresta gerada para a rede bipartida. Execucao interrompida.",
         call. = FALSE)
  }
  write.csv(edges_df, file.path(NETWORK_ARQ_DIR, "bipartite_edges.csv"), row.names = FALSE)

  g <- build_network_graph(edges_df)
  metrics_df <- compute_network_metrics(g)
  write.csv(metrics_df, file.path(NETWORK_ARQ_DIR, "bipartite_metrics.csv"), row.names = FALSE)

  export_cytoscape(g, file.path(NETWORK_ARQ_DIR, "bipartite_network.graphml"))

  ensure_dir(NETWORK_FIG_DIR)
  install_if_missing("ggplot2")
  ## theme_void() deixa o fundo do painel transparente (element_blank()) --
  ## em visualizadores sem fundo branco proprio isso aparece preto e torna
  ## os rotulos illegiveis (texto preto sobre "preto"). Fundo branco
  ## explicito no tema + bg="white" no ggsave corrige (bug real encontrado
  ## 2026-07-18). Com todos os hotspots (score>=2) a rede tende a virar um
  ## "hairball" ilegivel por volume -- para uma visao legivel de um
  ## subconjunto especifico, ver a rede focada da metanalise
  ## (Figuras/metanalise/rede_focada_xpc.png).
  p <- ggraph::ggraph(tidygraph::as_tbl_graph(g), layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.3) +
    ggraph::geom_node_point() +
    ggraph::geom_node_text(ggplot2::aes(label = name), size = 2, repel = TRUE, color = "black") +
    ggplot2::theme_void() +
    ggplot2::ggtitle("Rede bipartida: Proteina -> Regiao -> Gene") +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))
  ggplot2::ggsave(file.path(NETWORK_FIG_DIR, "bipartite_network.png"), p, width = 10, height = 10, dpi = 300, bg = "white")

  log_message("18_bipartite_network", sprintf(
    "Rede com %d no(s) e %d aresta(s) salva em '%s'.",
    igraph::vcount(g), igraph::ecount(g), NETWORK_ARQ_DIR
  ))
  log_message("18_bipartite_network", "Modulo 18 concluido.")
  invisible(list(graph = g, metrics = metrics_df, edges = edges_df))
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_18(occupancy_csv, hotspots_csv) explicitamente (interativamente
## ou a partir de 22_master_pipeline.R).
