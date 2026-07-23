## ============================================================================
## Modulo 08 -- Differential Binding (DiffBind / csaw)
## ============================================================================
## Descricao:
##   Compara WT vs deficiente para as duas proteinas que tem esse desenho
##   genuino (ver CLAUDE.md S7/S9): XPC (WT vs XPC-KO) e STAT2 (WT vs
##   STAT1-KO). ELK1 e STAT1 NAO entram neste modulo (sem braco deficiente
##   disponivel -- decisao registrada em CLAUDE.md S9).
##
##   Dois caminhos distintos, conforme a disponibilidade de input:
##   - STAT2 (input pareado nos dois genotipos): fluxo padrao do `DiffBind`
##     (dba() -> dba.count() -> dba.normalize() -> dba.contrast() -> dba.analyze()).
##   - XPC (XPC-KO sem input -- CLAUDE.md S9.1): NAO se compara so por
##     presenca/ausencia de pico. Constroi-se uma regiao consenso
##     (WT peaks UNIAO XPC-KO peaks) via GenomicRanges::reduce(), conta-se
##     reads de ChIP nessa regiao consenso com `csaw::regionCounts()`, e a
##     diferenca e' testada por contagem (edgeR, via `csaw::asDGEList()`).
##
## Entradas:
##   Dados/BAM/<amostra>.filtered.bam       (Modulo 05)
##   Dados/Peaks/<amostra>_peaks.*Peak       (Modulo 07)
##
## Saidas:
##   Arquivos/differential/<proteina>_diffbind_full.csv   (tabela completa)
##   Arquivos/differential/<proteina>_peaks_gained.csv
##   Arquivos/differential/<proteina>_peaks_lost.csv
##   Arquivos/differential/<proteina>_peaks_stable.csv
##   Logs/08_diffbind.log
##
## Dependencias:
##   00_setup.R
##   DiffBind, csaw, edgeR, GenomicRanges, rtracklayer (Bioconductor)
##
## Funcoes definidas neste modulo (import_peaks() vem de 00_setup.R):
##   build_consensus_regions(), count_reads_csaw(), classify_regions(),
##   run_diffbind_consensus_count(), run_diffbind_standard(), run_module_08()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
install_if_missing("DiffBind")
install_if_missing("csaw")
install_if_missing("edgeR")
install_if_missing("rtracklayer")
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(csaw)
  library(edgeR)
})

DIFFBIND_DIR <- file.path(PROJECT_DIRS$arquivos, "differential")

## Limiar de FDR usado para classificar uma regiao como ganha/perdida (nao
## estavel). Ver CLAUDE.md para justificativa (0.05 e' o padrao aceito para
## ChIP-seq differential binding).
FDR_THRESHOLD <- 0.05

## --- regiao consenso (import_peaks() vem de 00_setup.R) ---------------------------

#' Constroi a regiao consenso (uniao reduzida) de dois ou mais conjuntos de
#' picos -- usada quando um dos grupos nao tem input proprio e a comparacao
#' precisa ser feita por contagem de reads na mesma regiao, nao por
#' presenca/ausencia de pico (ver CLAUDE.md S9.1).
build_consensus_regions <- function(peak_files) {
  peaks_list <- lapply(peak_files, import_peaks)
  ## unlist(GRangesList(...)), nao do.call(c, ...) -- este ultimo pode nao
  ## disparar o dispatch S4 de c() para GRanges e devolver uma "list" comum
  ## em vez de um GRanges (bug encontrado e corrigido em 2026-07-13, ver
  ## CLAUDE.md).
  all_peaks <- unlist(GRangesList(lapply(peaks_list, granges)), use.names = FALSE)
  consensus <- GenomicRanges::reduce(all_peaks)
  log_message("08_diffbind", sprintf("Regiao consenso: %d regioes.", length(consensus)))
  consensus
}

## --- contagem por regiao (csaw) -----------------------------------------------------

#' Conta reads de ChIP em cada regiao consenso para uma lista de BAMs via
#' csaw::regionCounts().
count_reads_csaw <- function(bam_files, regions) {
  for (b in bam_files) validate_file_exists(b, "BAM filtrado")
  log_message("08_diffbind", sprintf("Contando reads em %d regiao(oes) para %d amostra(s).",
                                      length(regions), length(bam_files)))
  csaw::regionCounts(bam.files = bam_files, regions = regions)
}

## --- classificacao das regioes --------------------------------------------------

#' Classifica cada regiao como "Gained" (mais sinal no grupo deficiente),
#' "Lost" (menos sinal no deficiente) ou "Stable", com base em FDR e no
#' sinal do log fold-change. `logfc_col`/`fdr_col` identificam as colunas
#' relevantes na tabela de resultados (nomes variam entre DiffBind e edgeR).
classify_regions <- function(results_df, logfc_col, fdr_col, fdr_threshold = FDR_THRESHOLD) {
  results_df$status <- ifelse(
    results_df[[fdr_col]] >= fdr_threshold, "Stable",
    ifelse(results_df[[logfc_col]] > 0, "Gained", "Lost")
  )
  results_df
}

## --- caminho A: contagem em regiao consenso (XPC, XPC-KO sem input) -------------------

#' Fluxo de differential binding por contagem (nao por presenca/ausencia de
#' pico) -- usado para XPC porque o XPC-KO nao tem input (CLAUDE.md S9.1).
#' `wt_bams`/`ko_bams` sao os BAMs filtrados; `wt_peaks`/`ko_peaks` os
#' arquivos de picos correspondentes (usados so para construir o consenso).
run_diffbind_consensus_count <- function(wt_bams, ko_bams, wt_peaks, ko_peaks) {
  consensus <- build_consensus_regions(c(wt_peaks, ko_peaks))
  se <- count_reads_csaw(c(wt_bams, ko_bams), consensus)

  group <- factor(c(rep("WT", length(wt_bams)), rep("KO", length(ko_bams))), levels = c("WT", "KO"))
  y <- csaw::asDGEList(se, group = group)
  ## TMM (trimmed mean of M-values) -- corrige diferencas de composicao entre
  ## bibliotecas, alem do escalonamento por tamanho de biblioteca que
  ## asDGEList() ja faz sozinho (ver CLAUDE.md S9.1/plano tarefa 1b).
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  test <- edgeR::glmQLFTest(fit, coef = 2)

  tbl <- as.data.frame(SummarizedExperiment::rowRanges(se))
  tbl <- cbind(tbl, edgeR::topTags(test, n = Inf, sort.by = "none")$table)
  classify_regions(tbl, logfc_col = "logFC", fdr_col = "FDR")
}

## --- caminho B: DiffBind padrao (STAT2, input pareado nos dois genotipos) -------------

#' Fluxo padrao de DiffBind -- usado para STAT2, que tem input em WT e em
#' STAT1-KO (desenho completo, sem a limitacao do XPC-KO).
run_diffbind_standard <- function(sample_sheet) {
  install_if_missing("DiffBind")
  suppressMessages(library(DiffBind))
  dbObj <- dba(sampleSheet = sample_sheet)
  dbObj <- dba.count(dbObj)
  ## TMM em vez do default (so library-size) -- corrige diferencas de
  ## composicao entre bibliotecas (ver CLAUDE.md S9.1/plano tarefa 1b).
  dbObj <- dba.normalize(dbObj, normalize = DBA_NORM_TMM)
  dbObj <- dba.contrast(dbObj, categories = DBA_CONDITION, minMembers = 2)
  dbObj <- dba.analyze(dbObj)
  report <- dba.report(dbObj)
  tbl <- as.data.frame(report)
  classify_regions(tbl, logfc_col = "Fold", fdr_col = "FDR")
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 08 completo para XPC (via contagem em consenso) e STAT2
#' (via DiffBind padrao), salvando as tabelas completas e por categoria
#' (ganhos/perdidos/estaveis) em Arquivos/differential/.
run_module_08 <- function(xpc_samples, stat2_sample_sheet) {
  log_message("08_diffbind", "Iniciando Modulo 08 -- Differential Binding.")
  ensure_dir(DIFFBIND_DIR)

  ## XPC: WT vs XPC-KO por contagem em regiao consenso (CLAUDE.md S9.1)
  xpc_wt <- xpc_samples[xpc_samples$genotype == "WT", ]
  xpc_ko <- xpc_samples[xpc_samples$genotype == "XPC-KO", ]
  xpc_results <- run_diffbind_consensus_count(
    xpc_wt$bam, xpc_ko$bam, xpc_wt$peak_file, xpc_ko$peak_file
  )
  save_diffbind_results("XPC", xpc_results)

  ## STAT2: WT vs STAT1-KO via DiffBind padrao (input pareado nos dois lados)
  stat2_results <- run_diffbind_standard(stat2_sample_sheet)
  save_diffbind_results("STAT2", stat2_results)

  log_message("08_diffbind", "Modulo 08 concluido.")
  invisible(list(XPC = xpc_results, STAT2 = stat2_results))
}

#' Salva a tabela completa e as 3 tabelas por categoria (ganhos/perdidos/
#' estaveis) de uma proteina em Arquivos/differential/.
#'
#' Cria DIFFBIND_DIR se ainda nao existir -- necessario quando esta funcao e
#' chamada diretamente (fora de run_module_08(), que ja garante o
#' ensure_dir() antes), como nos scripts de retomada do XPC.
save_diffbind_results <- function(protein, results_df) {
  ensure_dir(DIFFBIND_DIR)
  write.csv(results_df, file.path(DIFFBIND_DIR, sprintf("%s_diffbind_full.csv", protein)),
            row.names = FALSE)
  for (status in c("Gained", "Lost", "Stable")) {
    subset_df <- results_df[results_df$status == status, ]
    label <- switch(status, Gained = "gained", Lost = "lost", Stable = "stable")
    write.csv(subset_df, file.path(DIFFBIND_DIR, sprintf("%s_peaks_%s.csv", protein, label)),
              row.names = FALSE)
  }
  log_message("08_diffbind", sprintf(
    "%s: %d regioes (%d ganhas, %d perdidas, %d estaveis).",
    protein, nrow(results_df),
    sum(results_df$status == "Gained"), sum(results_df$status == "Lost"),
    sum(results_df$status == "Stable")
  ))
  invisible(TRUE)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_08(xpc_samples, stat2_sample_sheet) explicitamente
## (interativamente ou a partir de 22_master_pipeline.R). `xpc_samples`
## precisa das colunas sample_id/genotype/bam/peak_file; `stat2_sample_sheet`
## segue o formato de sample sheet do DiffBind (SampleID, Condition,
## bamReads, Peaks, PeakCaller, ...).
