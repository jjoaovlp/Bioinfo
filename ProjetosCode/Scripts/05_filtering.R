## ============================================================================
## Modulo 05 -- Filtragem do BAM
## ============================================================================
## Descricao:
##   Filtra o BAM alinhado em uma unica passada via `Rsamtools::filterBam()`:
##   remove duplicatas de PCR (mesma posicao+fita, sem depender de
##   samtools/Picard), aplica o filtro de MAPQ (nativo via `ScanBamParam`) e
##   remove reads sobrepostos a regioes da blacklist ENCODE (via
##   `GenomicRanges::overlapsAny()` -- substitui `bedtools intersect -v`, ver
##   CLAUDE.md S5.2/S9). Produz o BAM final usado a partir do peak calling
##   (Modulo 07).
##
## Entradas:
##   Dados/BAM/<amostra>.sorted.bam(.bai)   (Modulo 04)
##
## Saidas:
##   Dados/BAM/<amostra>.filtered.bam(.bai)   -- BAM final
##   Arquivos/filtering/filtering_stats.csv
##   Arquivos/genome/<genoma>-blacklist.bed   -- blacklist ENCODE baixada
##   Logs/05_filtering.log
##
## Dependencias:
##   00_setup.R
##   Rsamtools, GenomicRanges, rtracklayer (Bioconductor)
##
## Funcoes definidas neste modulo:
##   download_encode_blacklist(), build_bam_filter(), filter_sample(),
##   run_module_05()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("Rsamtools")
install_if_missing("GenomicRanges")
install_if_missing("rtracklayer")
suppressMessages({
  library(Rsamtools)
  library(GenomicRanges)
  library(rtracklayer)
})

FILTERING_DIR <- file.path(PROJECT_DIRS$arquivos, "filtering")
GENOME_DIR <- file.path(PROJECT_DIRS$arquivos, "genome")

## URL padrao do repositorio oficial Boyle-Lab/Blacklist (ENCODE), a mesma
## fonte normalmente usada em pipelines de ChIP-seq para excluir regioes de
## sinal artefactual conhecidas.
ENCODE_BLACKLIST_URLS <- c(
  hg38 = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
  hg19 = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz"
)

## --- blacklist ENCODE --------------------------------------------------------------

#' Baixa a blacklist ENCODE oficial (Boyle-Lab/Blacklist) para um genoma, se
#' ainda nao existir localmente em Arquivos/genome/, e devolve como GRanges.
download_encode_blacklist <- function(genome_build = "hg38", dest_dir = GENOME_DIR) {
  if (!genome_build %in% names(ENCODE_BLACKLIST_URLS)) {
    stop(sprintf(
      "Validacao falhou: sem URL de blacklist ENCODE cadastrada para '%s'. Execucao interrompida.",
      genome_build
    ), call. = FALSE)
  }
  ensure_dir(dest_dir)
  bed_gz <- file.path(dest_dir, sprintf("%s-blacklist.bed.gz", genome_build))
  bed_path <- file.path(dest_dir, sprintf("%s-blacklist.bed", genome_build))
  if (!file.exists(bed_path)) {
    tryCatch({
      log_message("05_filtering", sprintf("Baixando blacklist ENCODE de %s.", genome_build))
      utils::download.file(ENCODE_BLACKLIST_URLS[[genome_build]], bed_gz, quiet = TRUE, mode = "wb")
      con_in <- gzfile(bed_gz, "rt")
      writeLines(readLines(con_in), bed_path)
      close(con_in)
      unlink(bed_gz)
    }, error = function(e) {
      log_message("05_filtering",
        sprintf("Falha ao baixar/descompactar blacklist de %s: %s", genome_build, conditionMessage(e)),
        level = "ERROR")
      stop(e)
    })
  } else {
    log_message("05_filtering", sprintf("Blacklist de %s ja existe em '%s'.", genome_build, bed_path))
  }
  validate_file_exists(bed_path, "blacklist ENCODE descompactada")
  rtracklayer::import(bed_path, format = "BED")
}

## --- filtro combinado (dedup + blacklist) ------------------------------------------

#' Constroi a FilterRules usada por filterBam(): remove reads duplicados
#' (mesmo cromossomo+posicao+fita) e reads sobrepostos a blacklist_gr. O
#' filtro de MAPQ e' aplicado separadamente via ScanBamParam(mapqFilter=),
#' de forma nativa e mais eficiente.
build_bam_filter <- function(blacklist_gr) {
  FilterRules(list(
    dedup_and_blacklist = function(x) {
      keep_dup <- !duplicated(paste(x$rname, x$pos, x$strand))
      read_gr <- GRanges(seqnames = x$rname, ranges = IRanges(start = x$pos, width = pmax(x$qwidth, 1)))
      keep_blacklist <- !overlapsAny(read_gr, blacklist_gr)
      keep_dup & keep_blacklist
    }
  ))
}

## --- orquestracao por amostra ---------------------------------------------------------

#' Aplica a filtragem completa a um BAM (MAPQ + dedup + blacklist em uma
#' unica passada de filterBam()), indexa o BAM final e registra o numero de
#' reads antes/depois.
filter_sample <- function(sample_id, bam_in, blacklist_gr, min_mapq = 30) {
  validate_file_exists(bam_in, "BAM de entrada")
  log_message("05_filtering", sprintf("Filtrando amostra '%s'.", sample_id))

  final_bam <- file.path(PROJECT_DIRS$bam, sprintf("%s.filtered.bam", sample_id))
  ensure_dir(dirname(final_bam))

  n_input <- countBam(bam_in)$records
  what <- c("rname", "pos", "strand", "qwidth")
  filterBam(
    bam_in, final_bam,
    filter = build_bam_filter(blacklist_gr),
    param = ScanBamParam(what = what, mapqFilter = min_mapq)
  )
  indexBam(final_bam)
  n_final <- countBam(final_bam)$records

  validate_file_exists(final_bam, "BAM final filtrado")
  validate_file_exists(paste0(final_bam, ".bai"), "indice do BAM final")

  log_message("05_filtering", sprintf(
    "Amostra '%s': %d -> %d reads apos MAPQ>=%d + dedup + blacklist.",
    sample_id, n_input, n_final, min_mapq
  ))
  data.frame(sample_id = sample_id, reads_input = n_input, reads_final = n_final,
             bam = final_bam, stringsAsFactors = FALSE)
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 05 para uma tabela de amostras (colunas: sample_id, bam),
#' baixando a blacklist ENCODE do genoma indicado uma unica vez e salvando o
#' resumo de reads perdidos por amostra em CSV.
run_module_05 <- function(samples_df, genome_build = "hg38", min_mapq = 30) {
  log_message("05_filtering", "Iniciando Modulo 05 -- Filtragem.")
  ensure_dir(FILTERING_DIR)
  blacklist_gr <- download_encode_blacklist(genome_build)
  results <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    filter_sample(row$sample_id, row$bam, blacklist_gr, min_mapq = min_mapq)
  })
  stats_df <- do.call(rbind, results)
  out_file <- file.path(FILTERING_DIR, "filtering_stats.csv")
  write.csv(stats_df, out_file, row.names = FALSE)
  log_message("05_filtering", sprintf("Estatisticas de filtragem salvas em '%s'.", out_file))
  log_message("05_filtering", "Modulo 05 concluido.")
  invisible(stats_df)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_05(samples_df, genome_build) explicitamente (interativamente ou
## a partir de 22_master_pipeline.R).
