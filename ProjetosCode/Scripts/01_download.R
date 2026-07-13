## ============================================================================
## Modulo 01 -- Download dos dados do GEO
## ============================================================================
## Descricao:
##   Baixa metadata de serie (SOFT/series matrix) do GEO para os datasets do
##   projeto via GEOquery, lista arquivos suplementares disponiveis (peaks,
##   bigWig, contagens ja processadas) para permitir reutilizacao quando
##   compativeis, e baixa FASTQ bruto via SRA Toolkit quando necessario.
##   Todo download e registrado em Logs/01_download.log.
##
##   Ver CLAUDE.md secao 7 e 9 para o racional cientifico por tras da escolha
##   e composicao final de cada dataset (inclui a troca de GSE217805 por
##   GSE222667 + GSE247724 por incompatibilidade de especie).
##
## Entradas:
##   Nenhuma (baixa da internet). Requer conexao ativa.
##
## Saidas:
##   Dados/GEO/<GSE>_series_matrix.txt / .soft  (metadata bruta do GEO)
##   Arquivos/geo_supp_files_<GSE>.rds           (listagem de arquivos supl.)
##   Dados/FASTQ/<SRR>_{1,2}.fastq.gz            (apenas se download_sra_fastq() for chamada explicitamente)
##   Logs/01_download.log
##
## Dependencias:
##   00_setup.R (funcoes auxiliares, diretorios do projeto)
##   GEOquery (metadata + arquivos suplementares)
##   prefetch / fasterq-dump do SRA Toolkit (apenas se baixar FASTQ bruto)
##
## Funcoes definidas neste modulo:
##   get_dataset_registry(), download_series_metadata(), list_supplementary_files(),
##   download_supplementary_files(), download_sra_fastq(), run_module_01()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("GEOquery")
library(GEOquery)

## --- registro dos datasets do projeto -----------------------------------------
## Estado final pos relatorio de compatibilidade (CLAUDE.md S7/S9). Cada linha
## eh uma serie GEO a baixar; um mesmo `dataset_id` pode ter mais de uma serie
## (caso do STAT2, que combina GSE222667 [WT] + GSE247724 [STAT1-KO]).

#' Retorna o registro de datasets do projeto como data.frame.
#'
#' @return data.frame com colunas: dataset_id, protein, gse, organism,
#'   cell_line, role (papel da serie dentro do dataset combinado), diffbind
#'   (se este dataset participa do Modulo 08 de differential binding).
get_dataset_registry <- function() {
  data.frame(
    dataset_id = c("XPC", "ELK1", "STAT1", "STAT2", "STAT2"),
    protein     = c("XPC", "ELK1", "STAT1", "STAT2", "STAT2"),
    gse         = c("GSE214182", "GSE91923", "GSE222667", "GSE222667", "GSE247724"),
    organism    = c("Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens"),
    cell_line   = c("U2OS", "A549", "Huh7.5", "Huh7.5", "Huh7.5"),
    role        = c("primary", "primary", "primary_wt_only", "primary_wt", "primary_ko"),
    diffbind    = c(TRUE, FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

## --- metadata de serie (SOFT / series matrix) ----------------------------------

#' Baixa a metadata de uma serie GEO (series matrix + GPL) para Dados/GEO/ e
#' devolve o objeto ExpressionSet (ou lista de ExpressionSets) do GEOquery.
#' Nunca sobrescreve dados brutos ja baixados: GEOquery::getGEO usa cache por
#' padrao quando o arquivo ja existe em `destdir`.
download_series_metadata <- function(gse_id) {
  validate_dir_exists(PROJECT_DIRS$geo, "diretorio Dados/GEO")
  log_message("01_download", sprintf("Baixando metadata de serie: %s", gse_id))
  gset <- tryCatch(
    GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE,
                      destdir = PROJECT_DIRS$geo, returnType = "ExpressionSet"),
    error = function(e) {
      log_message("01_download",
        sprintf("Falha ao baixar metadata de %s: %s", gse_id, conditionMessage(e)),
        level = "ERROR")
      NULL
    }
  )
  if (is.null(gset)) {
    return(invisible(NULL))
  }
  log_message("01_download", sprintf("Metadata de %s baixada com sucesso.", gse_id))
  gset
}

## --- arquivos suplementares (peaks/bigWig/contagens ja processados) -----------

#' Lista (sem baixar) os arquivos suplementares disponiveis para uma serie,
#' para decidir se ha arquivos processados compativeis reaproveitaveis antes
#' de refazer alinhamento/peak calling do zero.
list_supplementary_files <- function(gse_id) {
  log_message("01_download", sprintf("Listando arquivos suplementares de %s", gse_id))
  files <- tryCatch(
    GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, fetch_files = FALSE),
    error = function(e) {
      log_message("01_download",
        sprintf("Nao foi possivel listar arquivos suplementares de %s: %s",
                gse_id, conditionMessage(e)),
        level = "WARN")
      NULL
    }
  )
  out_file <- file.path(PROJECT_DIRS$arquivos, sprintf("geo_supp_files_%s.rds", gse_id))
  if (!is.null(files)) {
    saveRDS(files, out_file)
    log_message("01_download", sprintf("Listagem salva em '%s' (%d arquivo(s)).",
                                        out_file, nrow(files)))
  }
  files
}

#' Baixa de fato os arquivos suplementares de uma serie para Dados/GEO/.
#' So deve ser chamada explicitamente (nao roda automaticamente em
#' run_module_01()) para evitar downloads grandes nao intencionais.
download_supplementary_files <- function(gse_id, pattern = NULL) {
  dest <- file.path(PROJECT_DIRS$geo, gse_id)
  ensure_dir(dest)
  log_message("01_download", sprintf("Baixando arquivos suplementares de %s para '%s'", gse_id, dest))
  tryCatch(
    GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = dest,
                                filter_regex = pattern, fetch_files = TRUE),
    error = function(e) {
      log_message("01_download",
        sprintf("Falha ao baixar arquivos suplementares de %s: %s",
                gse_id, conditionMessage(e)),
        level = "ERROR")
      NULL
    }
  )
}

## --- FASTQ bruto via SRA Toolkit -------------------------------------------------

#' Baixa e converte para FASTQ um vetor de SRR accessions via `prefetch` +
#' `fasterq-dump`. Nunca sobrescreve um FASTQ ja existente. Interrompe apenas
#' o download do SRR problematico (via tryCatch), nao o modulo inteiro.
download_sra_fastq <- function(srr_ids, dest_dir = PROJECT_DIRS$fastq) {
  if (!check_external_tool("prefetch") || !check_external_tool("fasterq-dump")) {
    log_message("01_download",
      "SRA Toolkit ausente do PATH -- pulando download de FASTQ bruto.",
      level = "WARN")
    return(invisible(NULL))
  }
  ensure_dir(dest_dir)
  for (srr in srr_ids) {
    existing <- Sys.glob(file.path(dest_dir, paste0(srr, "*.fastq.gz")))
    if (length(existing) > 0) {
      log_message("01_download", sprintf("FASTQ de %s ja existe, pulando download.", srr))
      next
    }
    tryCatch({
      log_message("01_download", sprintf("prefetch %s", srr))
      system2("prefetch", args = c(srr, "--output-directory", dest_dir))
      log_message("01_download", sprintf("fasterq-dump %s", srr))
      system2("fasterq-dump", args = c(file.path(dest_dir, srr), "--split-files",
                                        "--outdir", dest_dir))
    }, error = function(e) {
      log_message("01_download",
        sprintf("Falha ao baixar %s: %s", srr, conditionMessage(e)),
        level = "ERROR")
    })
  }
  invisible(TRUE)
}

## --- execucao do modulo ----------------------------------------------------------

#' Executa o Modulo 01 para todos os datasets do registro: baixa metadata de
#' serie e lista (sem baixar) arquivos suplementares, permitindo decidir
#' depois se ha dados processados reaproveitaveis. Download de FASTQ bruto
#' e de arquivos suplementares completos e feito sob demanda (chamada
#' explicita de download_sra_fastq()/download_supplementary_files()), nao
#' automaticamente aqui, dado o volume potencial de dados.
run_module_01 <- function(registry = get_dataset_registry()) {
  log_message("01_download", "Iniciando Modulo 01 -- Download.")
  series_ids <- unique(registry$gse)
  results <- lapply(series_ids, function(gse_id) {
    metadata <- download_series_metadata(gse_id)
    supp <- list_supplementary_files(gse_id)
    list(gse = gse_id, metadata = metadata, supp_files = supp)
  })
  names(results) <- series_ids
  log_message("01_download", "Modulo 01 concluido.")
  invisible(results)
}

## Este modulo nao baixa nada automaticamente ao ser source()ado -- apenas
## define as funcoes acima. A execucao (download de fato) e disparada
## explicitamente chamando run_module_01(), seja interativamente, seja a
## partir de 22_master_pipeline.R.
