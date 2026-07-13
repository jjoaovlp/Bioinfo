## ============================================================================
## Modulo 00 -- Setup do Ambiente
## ============================================================================
## Descricao:
##   Prepara o ambiente de analise: instala/verifica pacotes R (via
##   BiocManager), verifica disponibilidade de softwares externos no PATH,
##   define funcoes auxiliares reutilizadas por todos os demais modulos
##   (logging, validacao defensiva de arquivos/diretorios, caminhos
##   relativos via `here`) e grava sessionInfo() para reprodutibilidade.
##
## Entradas:
##   Nenhuma. Parte do zero de uma instalacao limpa.
##
## Saidas:
##   Logs/00_setup.log
##   Logs/sessioninfo_<timestamp>.txt
##
## Dependencias:
##   here, BiocManager (instalados automaticamente se ausentes)
##
## Funcoes definidas neste modulo (fonte unica -- os demais modulos fazem
## source(here::here("Scripts", "00_setup.R")) para reaproveita-las):
##   ensure_dir(), project_path(), log_message(), validate_file_exists(),
##   validate_dir_exists(), install_if_missing(), install_project_packages(),
##   check_external_tool(), check_external_tools(), win_to_wsl_path(),
##   run_wsl_command(), check_macs3_wsl(), run_macs3(), save_session_info()
## ============================================================================

set.seed(1234)

if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

library(here)

## --- diretorios do projeto --------------------------------------------------

PROJECT_DIRS <- list(
  dados    = here("Dados"),
  geo      = here("Dados", "GEO"),
  fastq    = here("Dados", "FASTQ"),
  bam      = here("Dados", "BAM"),
  peaks    = here("Dados", "Peaks"),
  bigwig   = here("Dados", "BigWig"),
  metadata = here("Dados", "Metadata"),
  figuras  = here("Figuras"),
  arquivos = here("Arquivos"),
  scripts  = here("Scripts"),
  config   = here("Scripts", "config"),
  logs     = here("Logs")
)

#' Cria um diretorio recursivamente caso ainda nao exista. Idempotente.
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

invisible(lapply(PROJECT_DIRS, ensure_dir))

#' Atalho para here() restrito as pastas do projeto (nunca caminhos absolutos).
project_path <- function(...) here(...)

## --- logging ------------------------------------------------------------------

#' Grava uma mensagem de log com timestamp, modulo e nivel, no console e em
#' Logs/<modulo>.log. Nao interrompe a execucao -- use stop() explicitamente
#' (via validate_*()) quando uma falha for critica o suficiente para isso.
log_message <- function(module, msg, level = c("INFO", "WARN", "ERROR")) {
  level <- match.arg(level)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [%s] %s", timestamp, level, module, msg)
  message(line)
  log_file <- file.path(PROJECT_DIRS$logs, paste0(module, ".log"))
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

## --- validacao defensiva --------------------------------------------------------

#' Interrompe a execucao com mensagem clara caso o arquivo nao exista.
#' Usada nos pontos em que a ausencia do arquivo compromete a validade
#' cientifica da etapa seguinte (ex.: BAM final antes do peak calling).
validate_file_exists <- function(path, description = "arquivo") {
  if (!file.exists(path)) {
    stop(sprintf(
      "Validacao falhou: %s nao encontrado em '%s'. Execucao interrompida.",
      description, path
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Interrompe a execucao com mensagem clara caso o diretorio nao exista.
validate_dir_exists <- function(path, description = "diretorio") {
  if (!dir.exists(path)) {
    stop(sprintf(
      "Validacao falhou: %s nao encontrado em '%s'. Execucao interrompida.",
      description, path
    ), call. = FALSE)
  }
  invisible(TRUE)
}

## --- pacotes R ------------------------------------------------------------------

#' Instala um pacote via BiocManager caso ainda nao esteja disponivel.
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_message("00_setup", sprintf("Instalando pacote ausente: %s", pkg))
    tryCatch(
      BiocManager::install(pkg, update = FALSE, ask = FALSE),
      error = function(e) {
        log_message(
          "00_setup",
          sprintf("Falha ao instalar '%s': %s", pkg, conditionMessage(e)),
          level = "ERROR"
        )
        stop(e)
      }
    )
  }
  invisible(TRUE)
}

## Lista fixada de pacotes usados em algum ponto do pipeline (ver CLAUDE.md, S5.1).
## Desde 2026-07-13, ShortRead/Rsamtools/GenomicAlignments/Rbowtie2/ChIPQC/csaw
## substituem FastQC+MultiQC/samtools/bedtools/Bowtie2/deepTools como pacotes R
## nativos (ver CLAUDE.md S5.2 e S9) -- MACS3 continua externo (Python).
REQUIRED_PACKAGES <- c(
  "here", "BiocManager", "GEOquery", "Biobase",
  "GenomicRanges", "GenomicFeatures", "rtracklayer",
  "ShortRead", "Rsamtools", "GenomicAlignments", "Rbowtie2", "ChIPQC", "csaw",
  "ChIPseeker", "org.Hs.eg.db",
  "DiffBind", "clusterProfiler", "ReactomePA", "msigdbr",
  "STRINGdb", "igraph", "tidygraph", "ggraph",
  "data.table", "dplyr", "tidyr", "purrr", "stringr", "readr",
  "ggplot2", "yaml", "rmarkdown", "knitr"
)

#' Instala/verifica todas as dependencias do projeto.
install_project_packages <- function(packages = REQUIRED_PACKAGES) {
  invisible(lapply(packages, install_if_missing))
  log_message("00_setup", "Verificacao/instalacao de pacotes R concluida.")
}

## --- softwares externos (fora do R) ----------------------------------------------

## Ver CLAUDE.md S5.2 para o uso de cada um e o modulo que o requer. FastQC,
## MultiQC, Bowtie2, samtools e bedtools foram substituidos por pacotes R
## nativos (ShortRead, Rbowtie2, Rsamtools, GenomicRanges) -- so
## prefetch/fasterq-dump (download SRA) continuam precisando estar no PATH do
## Windows. MACS3 nao tem build bioconda para Windows e roda dentro do WSL
## (ver check_macs3_wsl()/run_macs3() abaixo), por isso nao entra em
## EXTERNAL_TOOLS (que assume checagem via Sys.which() no PATH do Windows).
EXTERNAL_TOOLS <- c("prefetch", "fasterq-dump")

#' Verifica se um software externo esta no PATH. Nao interrompe a execucao
#' do Modulo 00 (esses binarios so sao necessarios a partir do modulo
#' correspondente) -- apenas registra um aviso.
check_external_tool <- function(tool) {
  path <- Sys.which(tool)
  if (identical(unname(path), "")) {
    log_message("00_setup",
      sprintf("Software externo nao encontrado no PATH: %s", tool),
      level = "WARN")
    return(FALSE)
  }
  log_message("00_setup", sprintf("Software externo encontrado: %s -> %s", tool, path))
  TRUE
}

#' Roda check_external_tool() para todas as ferramentas de EXTERNAL_TOOLS.
check_external_tools <- function(tools = EXTERNAL_TOOLS) {
  status <- vapply(tools, check_external_tool, logical(1))
  names(status) <- tools
  status
}

## --- MACS3 via WSL (sem build bioconda para Windows) ------------------------------

## Distribuicao WSL e caminho do binario MACS3 dentro do ambiente conda
## "chipseq" criado especificamente para isso (ver CLAUDE.md S5.2/S9).
WSL_DISTRO <- "Ubuntu-24.04"
MACS3_WSL_BIN <- "~/miniconda3/envs/chipseq/bin/macs3"

#' Converte um caminho absoluto do Windows (ex. "C:/Users/x/y.bam") para o
#' caminho equivalente dentro do WSL (ex. "/mnt/c/Users/x/y.bam"), para poder
#' passar arquivos do projeto (sempre resolvidos via here()) para comandos
#' rodados dentro do WSL.
win_to_wsl_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!grepl("^[A-Za-z]:/", path)) {
    stop(sprintf(
      "Validacao falhou: '%s' nao parece um caminho absoluto do Windows (esperado 'C:/...').",
      path
    ), call. = FALSE)
  }
  drive <- tolower(substr(path, 1, 1))
  resto <- substring(path, 3)
  paste0("/mnt/", drive, resto)
}

#' Roda uma linha de comando dentro da distribuicao WSL configurada, via
#' `wsl -d <distro> -- bash -lc "<comando>"`. Devolve o codigo de saida.
run_wsl_command <- function(command_string, distro = WSL_DISTRO) {
  system2("wsl", args = c("-d", distro, "--", "bash", "-lc", shQuote(command_string)))
}

#' Verifica se o MACS3 esta disponivel dentro do WSL (nao interrompe a
#' execucao -- apenas avisa, seguindo o mesmo padrao de check_external_tool()).
check_macs3_wsl <- function() {
  status <- suppressWarnings(run_wsl_command(sprintf("%s --version", MACS3_WSL_BIN)))
  if (!identical(status, 0L)) {
    log_message("00_setup",
      "MACS3 (via WSL) nao encontrado ou WSL indisponivel -- Modulo 07 nao podera rodar.",
      level = "WARN")
    return(FALSE)
  }
  log_message("00_setup", "MACS3 (via WSL) disponivel.")
  TRUE
}

#' Roda o MACS3 (dentro do WSL) com os argumentos fornecidos. Os caminhos de
#' arquivo em `args` devem ja vir convertidos via win_to_wsl_path() -- esta
#' funcao nao faz a conversao sozinha, pois nem todo argumento e' um caminho
#' (ex. "-q", "0.01").
run_macs3 <- function(args) {
  command_string <- paste(MACS3_WSL_BIN, paste(args, collapse = " "))
  log_message("00_setup", sprintf("wsl macs3 %s", paste(args, collapse = " ")))
  status <- run_wsl_command(command_string)
  if (!identical(status, 0L)) {
    stop(sprintf("Falha ao rodar MACS3 via WSL (codigo %s). Execucao interrompida.", status),
         call. = FALSE)
  }
  invisible(status)
}

## --- sessionInfo -------------------------------------------------------------------

#' Salva sessionInfo() em Logs/ com timestamp, para reprodutibilidade.
save_session_info <- function() {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_file <- file.path(PROJECT_DIRS$logs, sprintf("sessioninfo_%s.txt", timestamp))
  writeLines(capture.output(sessionInfo()), out_file)
  log_message("00_setup", sprintf("sessionInfo() salvo em '%s'.", out_file))
  invisible(out_file)
}

## --- execucao principal --------------------------------------------------------------
## Guard: quando este script eh apenas source()ado por outro modulo para
## reaproveitar as funcoes acima, nao repetimos instalacao/checagem/sessionInfo
## a cada chamada dentro da mesma sessao de R.
if (!exists(".setup_00_done", envir = .GlobalEnv)) {
  log_message("00_setup", "Iniciando Modulo 00 -- Setup do ambiente.")
  install_project_packages()
  check_external_tools()
  check_macs3_wsl()
  save_session_info()
  assign(".setup_00_done", TRUE, envir = .GlobalEnv)
  log_message("00_setup", "Modulo 00 concluido.")
}
