## ============================================================================
## Modulo 20 -- Padronizacao de Visualizacoes
## ============================================================================
## Descricao:
##   Define a paleta de cores, o tema ggplot2 e a funcao de salvamento
##   padrao usados por todo o projeto: toda figura deve ser salva em
##   PNG+PDF+SVG, resolucao minima 300 dpi, com titulo, legenda de eixo e
##   fonte dos dados. Os modulos 03-19 ja salvam figuras individuais em PNG
##   a 300 dpi (via `ggsave()` direto); este modulo fornece `save_figure()`
##   como a forma padronizada de salvar novas figuras (usada pelo Modulo 22
##   no relatorio final) e pode ser adotado retroativamente pelos modulos
##   anteriores sem mudar a logica de cada um (ver CLAUDE.md S10 -- TODO).
##
## Entradas:
##   Um objeto ggplot (gerado por qualquer modulo)
##
## Saidas:
##   Figuras/<subpasta>/<nome>.{png,pdf,svg}
##   Logs/20_visualization.log
##
## Dependencias:
##   00_setup.R
##   ggplot2 (CRAN)
##
## Funcoes definidas neste modulo:
##   project_theme(), save_figure(), run_module_20()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("ggplot2")
suppressMessages(library(ggplot2))

## Paleta categorica (colorblind-safe, Okabe-Ito) usada para diferenciar
## proteinas/condicoes nas figuras do projeto.
PROJECT_PALETTE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

## --- tema e paleta --------------------------------------------------------------

#' Tema ggplot2 padrao do projeto: fundo limpo, texto legivel, adequado
#' para publicacao cientifica.
project_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(size = base_size * 0.7, color = "grey40"),
      legend.position = "right"
    )
}

## --- salvamento padronizado ------------------------------------------------------

#' Salva uma figura ggplot em PNG+PDF+SVG (300 dpi minimo), aplicando titulo,
#' subtitulo e fonte dos dados quando fornecidos. `output_dir` deve ser um
#' subdiretorio de Figuras/ (ex. "qc", "annotation", "network").
save_figure <- function(plot, filename, output_dir, title = NULL, subtitle = NULL,
                          source_label = NULL, width = 8, height = 6, dpi = 300,
                          formats = c("png", "pdf", "svg")) {
  if (dpi < 300) {
    stop("Validacao falhou: dpi deve ser >= 300 para figuras do projeto. Execucao interrompida.",
         call. = FALSE)
  }
  if (!is.null(title) || !is.null(subtitle) || !is.null(source_label)) {
    plot <- plot +
      ggplot2::labs(
        title = title, subtitle = subtitle,
        caption = if (!is.null(source_label)) sprintf("Fonte: %s", source_label) else NULL
      ) +
      project_theme()
  }
  full_dir <- file.path(PROJECT_DIRS$figuras, output_dir)
  ensure_dir(full_dir)
  for (fmt in formats) {
    out_file <- file.path(full_dir, sprintf("%s.%s", filename, fmt))
    ggplot2::ggsave(out_file, plot, width = width, height = height, dpi = dpi)
    log_message("20_visualization", sprintf("Figura salva em '%s'.", out_file))
  }
  invisible(TRUE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 20: verifica a arvore de Figuras/ e gera um manifesto
#' (contagem de figuras por subpasta) para o relatorio final (Modulo 22).
run_module_20 <- function() {
  log_message("20_visualization", "Iniciando Modulo 20 -- Padronizacao de Visualizacoes.")
  ensure_dir(PROJECT_DIRS$figuras)
  subdirs <- list.dirs(PROJECT_DIRS$figuras, recursive = FALSE)
  manifest <- data.frame(
    subpasta = basename(subdirs),
    n_figuras = vapply(subdirs, function(d) length(list.files(d)), integer(1)),
    stringsAsFactors = FALSE
  )
  out_file <- file.path(PROJECT_DIRS$arquivos, "figures_manifest.csv")
  write.csv(manifest, out_file, row.names = FALSE)
  log_message("20_visualization", sprintf("Manifesto de figuras salvo em '%s'.", out_file))
  log_message("20_visualization", "Modulo 20 concluido.")
  invisible(manifest)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## save_figure() a partir de qualquer modulo, ou run_module_20() para gerar
## o manifesto (interativamente ou a partir de 22_master_pipeline.R).
