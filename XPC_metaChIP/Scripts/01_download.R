## ============================================================================
## Modulo 01 -- Download dos dados do GEO
## ============================================================================
## Descricao:
##   Baixa metadata de serie (SOFT/series matrix) do GEO para os datasets do
##   projeto via GEOquery, lista arquivos suplementares disponiveis (peaks,
##   bigWig, contagens ja processadas) para permitir reutilizacao quando
##   compativeis, resolve SRX->SRR/URL de FASTQ via a API publica da ENA
##   (European Nucleotide Archive -- sem precisar de SRA Toolkit, ver
##   resolve_srr_table()) e baixa o FASTQ bruto diretamente do espelho HTTPS
##   da ENA. Todo download e registrado em Logs/01_download.log.
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
##   Dados/FASTQ/<SRR>.fastq.gz (ou _1/_2)       (apenas se download_fastq_ena() for chamada explicitamente)
##   Logs/01_download.log
##
## Dependencias:
##   00_setup.R (funcoes auxiliares, diretorios do projeto)
##   GEOquery (metadata + arquivos suplementares)
##   Nenhum software externo necessario para o FASTQ -- resolve_srr_table()/
##   download_fastq_ena() usam so a API HTTP da ENA e download.file() do R.
##
## Funcoes definidas neste modulo:
##   get_dataset_registry(), download_series_metadata(), list_supplementary_files(),
##   download_supplementary_files(), extract_srx_from_pdata(), get_ena_fastq_info(),
##   resolve_srr_table(), download_fastq_ena(), run_module_01()
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

## --- resolucao de SRR e download de FASTQ via ENA (sem SRA Toolkit) -------------------

## GEO so expoe o accession SRX (nao o SRR de fato) no campo `relation` da
## pData de cada amostra. A API publica da ENA (European Nucleotide Archive)
## resolve SRX -> SRR e ja devolve a URL direta do FASTQ (pronto, sem
## precisar de prefetch/fasterq-dump/SRA Toolkit) -- ela hospeda uma copia
## dos mesmos dados do SRA em formato FASTQ.gz direto.

#' Extrai o accession SRX do(s) campo(s) `relation*` da pData do GEO
#' (formato "SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRXxxxxxx").
#' Devolve NA para amostras sem link de SRA (nao deveria acontecer para
#' ChIP-seq, mas nunca assumido).
extract_srx_from_pdata <- function(pdata) {
  relation_cols <- grep("^relation", colnames(pdata), value = TRUE)
  if (length(relation_cols) == 0) return(rep(NA_character_, nrow(pdata)))
  apply(pdata[, relation_cols, drop = FALSE], 1, function(row) {
    m <- regmatches(row, regexpr("SRX[0-9]+", row))
    m <- m[nzchar(m)]
    if (length(m) == 0) NA_character_ else m[1]
  })
}

#' Consulta a API publica da ENA (Portal API, sem autenticacao) para
#' resolver um accession SRX (ou SRR) em corrida(s) SRR, URL(s) direta(s) de
#' FASTQ, tamanho e library layout. Devolve NULL (com log WARN) se a
#' consulta falhar -- nunca interrompe o modulo inteiro por uma amostra.
get_ena_fastq_info <- function(srx_id) {
  url <- sprintf(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,library_layout",
    srx_id
  )
  tryCatch(
    utils::read.delim(url, stringsAsFactors = FALSE),
    error = function(e) {
      log_message("01_download",
        sprintf("Falha ao consultar ENA para %s: %s", srx_id, conditionMessage(e)),
        level = "WARN")
      NULL
    }
  )
}

#' Resolve SRX->SRR/fastq_ftp/tamanho para cada GSM de uma pData do GEO,
#' devolvendo uma tabela (GSM, SRX, SRR, fastq_ftp, fastq_bytes,
#' library_layout) pronta para juntar ao metadata padronizado (Modulo 02).
resolve_srr_table <- function(pdata) {
  srx_ids <- extract_srx_from_pdata(pdata)
  gsms <- pdata$geo_accession
  rows <- lapply(seq_along(srx_ids), function(i) {
    if (is.na(srx_ids[i])) {
      return(data.frame(GSM = gsms[i], SRX = NA_character_, SRR = NA_character_,
                         fastq_ftp = NA_character_, fastq_bytes = NA_character_,
                         LibraryLayout = NA_character_, stringsAsFactors = FALSE))
    }
    ena <- get_ena_fastq_info(srx_ids[i])
    if (is.null(ena) || nrow(ena) == 0) {
      return(data.frame(GSM = gsms[i], SRX = srx_ids[i], SRR = NA_character_,
                         fastq_ftp = NA_character_, fastq_bytes = NA_character_,
                         LibraryLayout = NA_character_, stringsAsFactors = FALSE))
    }
    ## uma amostra pode ter mais de uma corrida (SRR) -- concatena com ";"
    data.frame(
      GSM = gsms[i], SRX = srx_ids[i],
      SRR = paste(ena$run_accession, collapse = ";"),
      fastq_ftp = paste(ena$fastq_ftp, collapse = ";"),
      fastq_bytes = paste(ena$fastq_bytes, collapse = ";"),
      LibraryLayout = ena$library_layout[1],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Verifica se um .fastq.gz esta integro (stream gzip completo, nao
#' truncado) tentando ler ate o fim via uma conexao gzfile(). Um download
#' truncado por timeout nao gera erro em download.file() por padrao (so um
#' arquivo incompleto silencioso) -- por isso essa checagem e' essencial
#' apos cada download, nunca assumir que "download.file() nao lancou erro"
#' significa "arquivo integro".
verify_gzip_integrity <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  result <- tryCatch({
    repeat {
      chunk <- readBin(con, "raw", n = 1e7)
      if (length(chunk) == 0) break
    }
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  result
}

#' Baixa os FASTQ de uma corrida diretamente do espelho HTTPS da ENA (sem
#' SRA Toolkit). `fastq_ftp` pode ter 1 URL (single-end) ou 2 separadas por
#' ";" (paired-end, mate 1/2); `expected_bytes` (opcional, mesmo formato,
#' vindo de fastq_bytes do Modulo 02/resolve_srr_table()) e' comparado ao
#' tamanho final para detectar truncamento. Nunca sobrescreve um FASTQ ja
#' existente E integro; um arquivo existente mas corrompido/incompleto e'
#' removido e baixado de novo. O timeout de download e' elevado para esta
#' chamada (arquivos de ChIP-seq facilmente passam de 1GB e nao cabem no
#' timeout padrao de 60s do R) e restaurado ao sair.
download_fastq_ena <- function(fastq_ftp, expected_bytes = NA_character_,
                                 dest_dir = PROJECT_DIRS$fastq, timeout_sec = 3600) {
  if (is.na(fastq_ftp) || !nzchar(fastq_ftp)) {
    return(invisible(character(0)))
  }
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout_sec, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  ensure_dir(dest_dir)
  urls <- strsplit(fastq_ftp, ";")[[1]]
  expected <- suppressWarnings(as.numeric(strsplit(as.character(expected_bytes), ";")[[1]]))
  vapply(seq_along(urls), function(i) {
    u <- urls[i]
    full_url <- if (startsWith(u, "http")) u else paste0("https://", u)
    dest <- file.path(dest_dir, basename(u))
    exp_size <- if (length(expected) >= i) expected[i] else NA_real_

    already_ok <- file.exists(dest) &&
      (is.na(exp_size) || abs(file.size(dest) - exp_size) < 1024) &&
      verify_gzip_integrity(dest)
    if (already_ok) {
      log_message("01_download", sprintf("FASTQ ja existe e esta integro: '%s'.", dest))
      return(dest)
    }
    if (file.exists(dest)) {
      log_message("01_download", sprintf("FASTQ existente incompleto/corrompido, refazendo: '%s'.", dest),
                  level = "WARN")
      unlink(dest)
    }
    tryCatch({
      log_message("01_download", sprintf("Baixando '%s'.", full_url))
      utils::download.file(full_url, dest, quiet = TRUE, mode = "wb")
      if (!is.na(exp_size) && abs(file.size(dest) - exp_size) > 1024) {
        stop(sprintf("tamanho baixado (%.0f) difere do esperado (%.0f) -- download truncado.",
                      file.size(dest), exp_size))
      }
      if (!verify_gzip_integrity(dest)) {
        stop("arquivo .gz corrompido/truncado (falhou verificacao de integridade).")
      }
      log_message("01_download", sprintf("Download de '%s' concluido e verificado.", dest))
    }, error = function(e) {
      log_message("01_download", sprintf("Falha ao baixar '%s': %s", full_url, conditionMessage(e)),
                  level = "ERROR")
      unlink(dest)
    })
    dest
  }, character(1))
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
