## ============================================================================
## Modulo 22 -- Pipeline Mestre
## ============================================================================
## Descricao:
##   Orquestra os Modulos 01-21 na ordem correta, com mensagens de progresso
##   claras, tempo de execucao por modulo, e tolerancia a falhas
##   configuravel (`stop_on_error`): por padrao, um modulo que falhar e'
##   registrado como "FALHOU" e o pipeline continua para os modulos
##   seguintes cujos insumos nao dependam do que falhou (cada `run_module_*`
##   ja valida suas proprias entradas e para sozinho se algo essencial
##   faltar). Ao final, gera relatorio HTML (sempre) e PDF (se houver
##   LaTeX/pandoc disponivel), tabela resumo, sessionInfo(), tempo total de
##   execucao e a arvore de arquivos produzidos.
##
##   `config` e' uma lista nomeada com os argumentos especificos de cada
##   modulo (ver PIPELINE_STEPS) -- um modulo cujo argumento nao esteja em
##   `config` e' pulado (status "PULADO"), nao falha o pipeline inteiro.
##
## Entradas:
##   Depende do que estiver em `config` -- tipicamente `samples_df`,
##   `index_prefix`, `xpc_samples`, `stat2_sample_sheet`, `condition_sample_ids`,
##   `protein_sample_ids`, `occupancy_csv`, `hotspots_csv`, `gene_symbols`,
##   `metadata_df`, etc. (um por modulo que precisa deles).
##
## Saidas:
##   Arquivos/pipeline_summary.csv
##   Arquivos/file_tree.txt
##   Logs/sessioninfo_<timestamp>.txt   (via save_session_info(), Modulo 00)
##   Arquivos/relatorio_final.html
##   Arquivos/relatorio_final.pdf   (se disponivel)
##   Logs/22_master_pipeline.log
##
## Dependencias:
##   00_setup.R e Scripts/01_download.R .. Scripts/21_validation.R
##
## Funcoes definidas neste modulo:
##   run_master_pipeline(), generate_file_tree(), generate_html_report(),
##   generate_pdf_report()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
for (n in sprintf("%02d", 1:21)) {
  script_file <- Sys.glob(here::here("Scripts", sprintf("%s_*.R", n)))
  if (length(script_file) == 1) source(script_file)
}

## Sequencia fixa de etapas do pipeline. Cada `fn` recebe o `config` inteiro
## e devolve invisibly o resultado do run_module_* correspondente; se um
## argumento necessario nao estiver em `config`, a propria funcao de cada
## modulo ja para com validate_file_exists()/stop() -- aqui so decidimos se
## o modulo tem o minimo de config para ser tentado.
PIPELINE_STEPS <- list(
  list(id = "01", name = "Download",               requires = character(0), fn = function(cfg) run_module_01()),
  list(id = "02", name = "Metadata",                requires = character(0), fn = function(cfg) run_module_02()),
  list(id = "03", name = "QC",                      requires = "samples_df", fn = function(cfg) run_module_03(cfg$samples_df)),
  list(id = "04", name = "Alignment",                requires = c("samples_df", "index_prefix"), fn = function(cfg) run_module_04(cfg$samples_df, cfg$index_prefix)),
  list(id = "05", name = "Filtering",                requires = c("samples_df", "genome_build"), fn = function(cfg) run_module_05(cfg$samples_df, cfg$genome_build)),
  list(id = "06", name = "ChIP-QC",                  requires = "samples_df", fn = function(cfg) run_module_06(cfg$samples_df)),
  list(id = "07", name = "Peak Calling",             requires = "samples_df", fn = function(cfg) run_module_07(cfg$samples_df)),
  list(id = "08", name = "Differential Binding",      requires = c("xpc_samples", "stat2_sample_sheet"), fn = function(cfg) run_module_08(cfg$xpc_samples, cfg$stat2_sample_sheet)),
  list(id = "09", name = "Annotation",                requires = "samples_df", fn = function(cfg) run_module_09(cfg$samples_df)),
  list(id = "10", name = "Enrichment",                requires = c("proteins", "annotation_csvs"), fn = function(cfg) lapply(cfg$proteins, function(p) run_module_10(p, cfg$annotation_csvs[[p]]))),
  list(id = "11", name = "GRanges",                   requires = "samples_df", fn = function(cfg) run_module_11(cfg$samples_df)),
  list(id = "12", name = "Genome Standardization",     requires = "samples_df", fn = function(cfg) run_module_12(cfg$samples_df)),
  list(id = "13", name = "Regulatory Universe",        requires = "sample_ids", fn = function(cfg) run_module_13(cfg$sample_ids)),
  list(id = "14", name = "Overlap Matrix",              requires = "condition_sample_ids", fn = function(cfg) run_module_14(cfg$condition_sample_ids)),
  list(id = "15", name = "Pairwise Overlap",            requires = "condition_sample_ids", fn = function(cfg) run_module_15(cfg$condition_sample_ids)),
  list(id = "16", name = "Multiple Overlap",            requires = "protein_sample_ids", fn = function(cfg) run_module_16(cfg$protein_sample_ids)),
  list(id = "17", name = "Hotspots",                    requires = character(0), fn = function(cfg) run_module_17(if (is.null(cfg$min_score)) 2 else cfg$min_score)),
  list(id = "18", name = "Bipartite Network",           requires = c("occupancy_csv", "hotspots_csv"), fn = function(cfg) run_module_18(cfg$occupancy_csv, cfg$hotspots_csv)),
  list(id = "19", name = "Pathway Network",             requires = "gene_symbols", fn = function(cfg) run_module_19(cfg$gene_symbols, cfg$go_result, cfg$reactome_result)),
  list(id = "20", name = "Visualization",               requires = character(0), fn = function(cfg) run_module_20()),
  list(id = "21", name = "Validation",                  requires = character(0), fn = function(cfg) run_module_21(metadata_df = cfg$metadata_df))
)

## --- relatorio final ----------------------------------------------------------------

#' Gera a arvore de arquivos produzidos pelo projeto (Dados/Figuras/Arquivos),
#' salva em texto simples.
generate_file_tree <- function() {
  dirs <- c(PROJECT_DIRS$dados, PROJECT_DIRS$figuras, PROJECT_DIRS$arquivos)
  tree <- unlist(lapply(dirs, function(d) list.files(d, recursive = TRUE, full.names = TRUE)))
  out_file <- file.path(PROJECT_DIRS$arquivos, "file_tree.txt")
  writeLines(tree, out_file)
  out_file
}

#' Gera o relatorio final em HTML (sempre funciona, sem depender de
#' pandoc/LaTeX -- HTML puro via base R).
generate_html_report <- function(summary_df, total_time_min) {
  rows <- paste(sprintf(
    "<tr><td>%s</td><td>%s</td><td>%s</td><td>%.1f</td></tr>",
    summary_df$modulo, summary_df$nome, summary_df$status, summary_df$tempo_seg
  ), collapse = "\n")
  html <- sprintf('<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Relatorio Final -- Pipeline ChIP-seq</title>
<style>body{font-family:sans-serif;margin:2em;} table{border-collapse:collapse;width:100%%;}
td,th{border:1px solid #ccc;padding:6px 10px;text-align:left;}
th{background:#f0f0f0;} .OK{color:green;} .FALHOU{color:red;} .PULADO{color:#999;}</style>
</head><body>
<h1>Relatorio Final -- Pipeline Integrativo de ChIP-seq</h1>
<p>Gerado em %s. Tempo total: %.1f minuto(s).</p>
<table><tr><th>Modulo</th><th>Nome</th><th>Status</th><th>Tempo (s)</th></tr>
%s
</table>
<p>Ver CLAUDE.md para o registro completo de decisoes metodologicas, dependencias e
historico do projeto.</p>
</body></html>', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), total_time_min, rows)

  out_file <- file.path(PROJECT_DIRS$arquivos, "relatorio_final.html")
  writeLines(html, out_file)
  out_file
}

#' Tenta gerar o relatorio final em PDF via rmarkdown (precisa de
#' pandoc+LaTeX/tinytex). Se nao disponivel, avisa e pula -- nao interrompe
#' o pipeline por causa disso.
generate_pdf_report <- function(summary_df) {
  if (!rmarkdown::pandoc_available()) {
    log_message("22_master_pipeline", "Pandoc indisponivel -- pulando relatorio PDF.", level = "WARN")
    return(invisible(NULL))
  }
  rmd_content <- c(
    "---", "title: \"Relatorio Final -- Pipeline ChIP-seq\"", "output: pdf_document", "---",
    "", "```{r, echo=FALSE}", "knitr::kable(summary_df)", "```"
  )
  rmd_file <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, rmd_file)
  out_file <- file.path(PROJECT_DIRS$arquivos, "relatorio_final.pdf")
  tryCatch({
    rmarkdown::render(rmd_file, output_format = "pdf_document", output_file = out_file, quiet = TRUE,
                       envir = list2env(list(summary_df = summary_df)))
    out_file
  }, error = function(e) {
    log_message("22_master_pipeline",
      sprintf("Falha ao gerar PDF (provavelmente falta LaTeX): %s", conditionMessage(e)), level = "WARN")
    NULL
  })
}

## --- execucao do pipeline mestre -----------------------------------------------------

#' Executa o Modulo 22: roda os Modulos 01-21 na ordem, pulando os que nao
#' tem o config minimo, registrando status/tempo de cada um, e gera o
#' relatorio final (HTML sempre, PDF se possivel) + sessionInfo + arvore de
#' arquivos. `stop_on_error = TRUE` interrompe o pipeline inteiro no
#' primeiro modulo que falhar; o padrao (FALSE) registra a falha e segue
#' para os proximos.
run_master_pipeline <- function(config = list(), stop_on_error = FALSE) {
  log_message("22_master_pipeline", "Iniciando Pipeline Mestre.")
  start_time <- Sys.time()
  rows <- list()

  for (step in PIPELINE_STEPS) {
    missing_req <- setdiff(step$requires, names(config))
    if (length(missing_req) > 0) {
      log_message("22_master_pipeline", sprintf(
        "Modulo %s (%s): PULADO -- falta config: %s.",
        step$id, step$name, paste(missing_req, collapse = ", ")
      ), level = "WARN")
      rows[[step$id]] <- data.frame(modulo = step$id, nome = step$name, status = "PULADO",
                                     tempo_seg = 0, stringsAsFactors = FALSE)
      next
    }
    step_start <- Sys.time()
    status <- tryCatch({
      step$fn(config)
      "OK"
    }, error = function(e) {
      log_message("22_master_pipeline", sprintf(
        "Modulo %s (%s): FALHOU -- %s", step$id, step$name, conditionMessage(e)
      ), level = "ERROR")
      if (stop_on_error) stop(e)
      "FALHOU"
    })
    elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
    log_message("22_master_pipeline", sprintf("Modulo %s (%s): %s (%.1fs).", step$id, step$name, status, elapsed))
    rows[[step$id]] <- data.frame(modulo = step$id, nome = step$name, status = status,
                                   tempo_seg = round(elapsed, 1), stringsAsFactors = FALSE)
  }

  summary_df <- do.call(rbind, rows)
  total_time_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  write.csv(summary_df, file.path(PROJECT_DIRS$arquivos, "pipeline_summary.csv"), row.names = FALSE)
  save_session_info()
  generate_file_tree()
  html_path <- generate_html_report(summary_df, total_time_min)
  pdf_path <- generate_pdf_report(summary_df)

  log_message("22_master_pipeline", sprintf(
    "Pipeline concluido em %.1f minuto(s). Relatorio: '%s'.", total_time_min, html_path
  ))
  invisible(list(summary = summary_df, html_report = html_path, pdf_report = pdf_path))
}
