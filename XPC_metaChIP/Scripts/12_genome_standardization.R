## ============================================================================
## Modulo 12 -- Padronizacao do Genoma (liftOver hg19 -> hg38)
## ============================================================================
## Descricao:
##   Padroniza todas as amostras para hg38. Le o `genome_build` de cada
##   amostra a partir do metadata/Modulo 11; quando ja for hg38, devolve o
##   GRanges inalterado; quando for hg19 (caso do ELK1, que a pagina GEO do
##   GSE91923 disponibiliza em hg19 e GRCh38 -- ver CLAUDE.md S10), aplica
##   `rtracklayer::liftOver()` com a chain oficial hg19->hg38 da UCSC. Nunca
##   assume a montagem -- interrompe a execucao se `genome_build` nao for
##   "hg19" nem "hg38" (nao ha chain cadastrada para outra montagem/especie,
##   coerente com a decisao de nao misturar especies no universo regulatorio
##   -- CLAUDE.md S9).
##
## Entradas:
##   Arquivos/granges/<amostra>.rds   (Modulo 11)
##
## Saidas:
##   Arquivos/granges_hg38/<amostra>.rds   (GRanges padronizado em hg38)
##   Arquivos/genome/hg19ToHg38.over.chain (baixada uma unica vez)
##   Logs/12_genome_standardization.log
##
## Dependencias:
##   00_setup.R
##   rtracklayer (Bioconductor)
##
## Funcoes definidas neste modulo:
##   download_liftover_chain(), liftover_to_hg38(), run_module_12()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("rtracklayer")
suppressMessages(library(rtracklayer))

GENOME_DIR <- file.path(PROJECT_DIRS$arquivos, "genome")
GRANGES_HG38_DIR <- file.path(PROJECT_DIRS$arquivos, "granges_hg38")

## URL oficial da UCSC para a chain de liftOver hg19 -> hg38.
LIFTOVER_HG19_TO_HG38_URL <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

## --- chain de liftOver --------------------------------------------------------------

#' Baixa (uma unica vez) a chain oficial hg19->hg38 da UCSC e devolve o
#' objeto Chain pronto para uso em rtracklayer::liftOver().
download_liftover_chain <- function(dest_dir = GENOME_DIR) {
  ensure_dir(dest_dir)
  chain_gz <- file.path(dest_dir, "hg19ToHg38.over.chain.gz")
  chain_path <- file.path(dest_dir, "hg19ToHg38.over.chain")
  if (!file.exists(chain_path)) {
    tryCatch({
      log_message("12_genome_standardization", "Baixando chain de liftOver hg19->hg38 (UCSC).")
      utils::download.file(LIFTOVER_HG19_TO_HG38_URL, chain_gz, quiet = TRUE, mode = "wb")
      con_in <- gzfile(chain_gz, "rt")
      writeLines(readLines(con_in), chain_path)
      close(con_in)
      unlink(chain_gz)
    }, error = function(e) {
      log_message("12_genome_standardization",
        sprintf("Falha ao baixar/descompactar chain de liftOver: %s", conditionMessage(e)),
        level = "ERROR")
      stop(e)
    })
  } else {
    log_message("12_genome_standardization", sprintf("Chain de liftOver ja existe em '%s'.", chain_path))
  }
  validate_file_exists(chain_path, "chain de liftOver hg19->hg38")
  rtracklayer::import.chain(chain_path)
}

## --- padronizacao por amostra ---------------------------------------------------------

#' Padroniza um GRanges para hg38: devolve inalterado se `genome_build` ja
#' for "hg38"; aplica liftOver se for "hg19"; interrompe a execucao para
#' qualquer outro valor (nunca assume a montagem -- CLAUDE.md S10).
liftover_to_hg38 <- function(gr, genome_build, sample_id) {
  if (identical(genome_build, "hg38")) {
    log_message("12_genome_standardization", sprintf("'%s' ja esta em hg38.", sample_id))
    return(gr)
  }
  if (!identical(genome_build, "hg19")) {
    stop(sprintf(
      "Validacao falhou: genome_build '%s' (amostra '%s') sem chain de liftOver cadastrada. Execucao interrompida.",
      genome_build, sample_id
    ), call. = FALSE)
  }
  chain <- download_liftover_chain()
  lifted_list <- rtracklayer::liftOver(gr, chain)
  n_input <- length(gr)
  n_mapped <- sum(lengths(lifted_list) > 0)
  if (n_mapped < n_input) {
    log_message("12_genome_standardization", sprintf(
      "'%s': %d de %d regioes nao mapearam de hg19 para hg38 e foram descartadas.",
      sample_id, n_input - n_mapped, n_input
    ), level = "WARN")
  }
  lifted <- unlist(lifted_list)
  if (length(lifted) == 0) {
    stop(sprintf(
      "Validacao falhou: nenhuma regiao de '%s' mapeou para hg38. Execucao interrompida.", sample_id
    ), call. = FALSE)
  }
  lifted
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 12 para uma tabela de amostras (colunas: sample_id,
#' genome_build), lendo o GRanges de Arquivos/granges/ (Modulo 11) e
#' salvando a versao padronizada em hg38 em Arquivos/granges_hg38/.
run_module_12 <- function(samples_df) {
  log_message("12_genome_standardization", "Iniciando Modulo 12 -- Padronizacao do genoma.")
  ensure_dir(GRANGES_HG38_DIR)
  granges_dir <- file.path(PROJECT_DIRS$arquivos, "granges")

  result <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    gr_file <- file.path(granges_dir, sprintf("%s.rds", row$sample_id))
    validate_file_exists(gr_file, sprintf("GRanges de '%s' (Modulo 11)", row$sample_id))
    gr <- readRDS(gr_file)
    gr_hg38 <- liftover_to_hg38(gr, row$genome_build, row$sample_id)
    saveRDS(gr_hg38, file.path(GRANGES_HG38_DIR, sprintf("%s.rds", row$sample_id)))
    gr_hg38
  })
  names(result) <- samples_df$sample_id
  log_message("12_genome_standardization", "Modulo 12 concluido.")
  invisible(result)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_12(samples_df) explicitamente (interativamente ou a partir de
## 22_master_pipeline.R).
