## ============================================================================
## Modulo 05 -- Filtragem do BAM
## ============================================================================
## Descricao:
##   Marca e remove duplicatas de PCR (via samtools fixmate/markdup, sem
##   depender de Picard/Java), filtra por qualidade de mapeamento (MAPQ) e
##   remove reads sobrepostos a regioes da blacklist ENCODE (via bedtools),
##   produzindo o BAM final usado a partir do peak calling (Modulo 07).
##
## Entradas:
##   Dados/BAM/<amostra>.sorted.bam(.bai)   (Modulo 04)
##
## Saidas:
##   Dados/BAM/<amostra>.filtered.bam(.bai)   -- BAM final
##   Arquivos/filtering/<amostra>_filtering_stats.csv
##   Arquivos/genome/<genoma>-blacklist.bed   -- blacklist ENCODE baixada
##   Logs/05_filtering.log
##
## Dependencias:
##   00_setup.R
##   samtools e bedtools no PATH (verificados via check_external_tool())
##
## Funcoes definidas neste modulo:
##   mark_duplicates(), remove_marked_duplicates(), filter_mapq(),
##   download_encode_blacklist(), remove_blacklist_regions(),
##   filter_sample(), run_module_05()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))

FILTERING_DIR <- file.path(PROJECT_DIRS$arquivos, "filtering")
GENOME_DIR <- file.path(PROJECT_DIRS$arquivos, "genome")

## URL padrao do repositorio oficial Boyle-Lab/Blacklist (ENCODE), a mesma
## fonte normalmente usada em pipelines de ChIP-seq para excluir regioes de
## sinal artefactual conhecidas.
ENCODE_BLACKLIST_URLS <- c(
  hg38 = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
  hg19 = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz"
)

## --- duplicatas (samtools, sem Picard/Java) --------------------------------------------

#' Marca duplicatas de PCR em um BAM via samtools (name-sort -> fixmate ->
#' coord-sort -> markdup), sem remove-las. Usa arquivos temporarios
#' intermediarios (nunca pipes de shell) para portabilidade entre SOs.
mark_duplicates <- function(bam_in, bam_out, threads = 4) {
  if (!check_external_tool("samtools")) {
    stop("Validacao falhou: 'samtools' ausente do PATH. Execucao interrompida.", call. = FALSE)
  }
  validate_file_exists(bam_in, "BAM de entrada")
  ensure_dir(dirname(bam_out))

  namesorted <- tempfile(fileext = ".bam")
  fixmated <- tempfile(fileext = ".bam")
  coordsorted <- tempfile(fileext = ".bam")
  on.exit(unlink(c(namesorted, fixmated, coordsorted)), add = TRUE)

  log_message("05_filtering", sprintf("samtools sort -n '%s'", bam_in))
  st1 <- system2("samtools", args = c("sort", "-@", threads, "-n", "-o", namesorted, bam_in))
  if (st1 != 0) stop(sprintf("Falha em 'samtools sort -n' (codigo %d).", st1), call. = FALSE)

  log_message("05_filtering", "samtools fixmate -m")
  st2 <- system2("samtools", args = c("fixmate", "-m", namesorted, fixmated))
  if (st2 != 0) stop(sprintf("Falha em 'samtools fixmate' (codigo %d).", st2), call. = FALSE)

  log_message("05_filtering", "samtools sort (coordenada)")
  st3 <- system2("samtools", args = c("sort", "-@", threads, "-o", coordsorted, fixmated))
  if (st3 != 0) stop(sprintf("Falha em 'samtools sort' (codigo %d).", st3), call. = FALSE)

  log_message("05_filtering", sprintf("samtools markdup -> '%s'", bam_out))
  st4 <- system2("samtools", args = c("markdup", coordsorted, bam_out))
  if (st4 != 0) stop(sprintf("Falha em 'samtools markdup' (codigo %d).", st4), call. = FALSE)

  validate_file_exists(bam_out, "BAM com duplicatas marcadas")
  invisible(bam_out)
}

#' Remove (nao so marca) os reads flagados como duplicata (flag SAM 1024).
remove_marked_duplicates <- function(bam_in, bam_out, threads = 4) {
  if (!check_external_tool("samtools")) {
    stop("Validacao falhou: 'samtools' ausente do PATH. Execucao interrompida.", call. = FALSE)
  }
  validate_file_exists(bam_in, "BAM com duplicatas marcadas")
  ensure_dir(dirname(bam_out))
  log_message("05_filtering", sprintf("Removendo duplicatas marcadas -> '%s'", bam_out))
  status <- system2("samtools", args = c("view", "-@", threads, "-b", "-F", "1024",
                                          "-o", bam_out, bam_in))
  if (status != 0) stop(sprintf("Falha ao remover duplicatas (codigo %d).", status), call. = FALSE)
  invisible(bam_out)
}

## --- filtro de qualidade de mapeamento -------------------------------------------------

#' Mantem apenas reads com MAPQ >= min_mapq (padrao 30, recomendacao comum
#' para ChIP-seq de fator de transcricao/histona).
filter_mapq <- function(bam_in, bam_out, min_mapq = 30, threads = 4) {
  if (!check_external_tool("samtools")) {
    stop("Validacao falhou: 'samtools' ausente do PATH. Execucao interrompida.", call. = FALSE)
  }
  validate_file_exists(bam_in, "BAM de entrada")
  ensure_dir(dirname(bam_out))
  log_message("05_filtering", sprintf("Filtrando MAPQ >= %d -> '%s'", min_mapq, bam_out))
  status <- system2("samtools", args = c("view", "-@", threads, "-b", "-q", min_mapq,
                                          "-o", bam_out, bam_in))
  if (status != 0) stop(sprintf("Falha no filtro de MAPQ (codigo %d).", status), call. = FALSE)
  invisible(bam_out)
}

## --- blacklist ENCODE --------------------------------------------------------------

#' Baixa a blacklist ENCODE oficial (Boyle-Lab/Blacklist) para um genoma, se
#' ainda nao existir localmente em Arquivos/genome/.
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
  if (file.exists(bed_path)) {
    log_message("05_filtering", sprintf("Blacklist de %s ja existe em '%s'.", genome_build, bed_path))
    return(invisible(bed_path))
  }
  tryCatch({
    log_message("05_filtering", sprintf("Baixando blacklist ENCODE de %s.", genome_build))
    utils::download.file(ENCODE_BLACKLIST_URLS[[genome_build]], bed_gz, quiet = TRUE, mode = "wb")
    R.utils_available <- requireNamespace("R.utils", quietly = TRUE)
    if (R.utils_available) {
      R.utils::gunzip(bed_gz, bed_path, remove = TRUE)
    } else {
      con_in <- gzfile(bed_gz, "rt")
      writeLines(readLines(con_in), bed_path)
      close(con_in)
      unlink(bed_gz)
    }
  }, error = function(e) {
    log_message("05_filtering",
      sprintf("Falha ao baixar/descompactar blacklist de %s: %s", genome_build, conditionMessage(e)),
      level = "ERROR")
    stop(e)
  })
  validate_file_exists(bed_path, "blacklist ENCODE descompactada")
  invisible(bed_path)
}

#' Remove reads que se sobrepoem a regioes da blacklist via
#' `bedtools intersect -v`.
remove_blacklist_regions <- function(bam_in, blacklist_bed, bam_out) {
  if (!check_external_tool("bedtools")) {
    stop("Validacao falhou: 'bedtools' ausente do PATH. Execucao interrompida.", call. = FALSE)
  }
  validate_file_exists(bam_in, "BAM de entrada")
  validate_file_exists(blacklist_bed, "BED da blacklist ENCODE")
  ensure_dir(dirname(bam_out))
  log_message("05_filtering", sprintf("Removendo regioes da blacklist -> '%s'", bam_out))
  status <- system2("bedtools",
    args = c("intersect", "-v", "-abam", bam_in, "-b", blacklist_bed),
    stdout = bam_out)
  if (status != 0) stop(sprintf("Falha no bedtools intersect (codigo %d).", status), call. = FALSE)
  invisible(bam_out)
}

## --- orquestracao por amostra ---------------------------------------------------------

#' Aplica a sequencia completa de filtragem a um BAM: marca+remove
#' duplicatas, filtra MAPQ, remove blacklist, indexa o BAM final e registra
#' o numero de reads perdido em cada etapa.
filter_sample <- function(sample_id, bam_in, blacklist_bed, min_mapq = 30, threads = 4) {
  log_message("05_filtering", sprintf("Filtrando amostra '%s'.", sample_id))
  ensure_dir(FILTERING_DIR)

  marked <- tempfile(fileext = ".bam")
  deduped <- tempfile(fileext = ".bam")
  mapq_filtered <- tempfile(fileext = ".bam")
  final_bam <- file.path(PROJECT_DIRS$bam, sprintf("%s.filtered.bam", sample_id))
  on.exit(unlink(c(marked, deduped, mapq_filtered)), add = TRUE)

  count_reads <- function(bam) {
    as.numeric(system2("samtools", args = c("view", "-c", bam), stdout = TRUE))
  }

  n_input <- count_reads(bam_in)
  mark_duplicates(bam_in, marked, threads = threads)
  remove_marked_duplicates(marked, deduped, threads = threads)
  n_dedup <- count_reads(deduped)
  filter_mapq(deduped, mapq_filtered, min_mapq = min_mapq, threads = threads)
  n_mapq <- count_reads(mapq_filtered)
  remove_blacklist_regions(mapq_filtered, blacklist_bed, final_bam)
  n_final <- count_reads(final_bam)

  status <- system2("samtools", args = c("index", final_bam))
  if (status != 0) stop(sprintf("Falha ao indexar BAM final (codigo %d).", status), call. = FALSE)
  validate_file_exists(paste0(final_bam, ".bai"), "indice do BAM final")

  stats <- data.frame(
    sample_id = sample_id, reads_input = n_input, reads_pos_dedup = n_dedup,
    reads_pos_mapq = n_mapq, reads_final = n_final, bam = final_bam,
    stringsAsFactors = FALSE
  )
  log_message("05_filtering", sprintf(
    "Amostra '%s': %d -> %d (dedup) -> %d (MAPQ) -> %d (blacklist).",
    sample_id, n_input, n_dedup, n_mapq, n_final
  ))
  stats
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 05 para uma tabela de amostras (colunas: sample_id, bam),
#' baixando a blacklist ENCODE do genoma indicado uma unica vez e salvando o
#' resumo de reads perdidos por etapa em CSV.
run_module_05 <- function(samples_df, genome_build = "hg38", min_mapq = 30, threads = 4) {
  log_message("05_filtering", "Iniciando Modulo 05 -- Filtragem.")
  ensure_dir(FILTERING_DIR)
  blacklist_bed <- download_encode_blacklist(genome_build)
  results <- lapply(seq_len(nrow(samples_df)), function(i) {
    row <- samples_df[i, ]
    filter_sample(row$sample_id, row$bam, blacklist_bed, min_mapq = min_mapq, threads = threads)
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
