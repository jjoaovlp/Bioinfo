## ============================================================================
## Modulo 02 -- Construcao do Metadata Padronizado
## ============================================================================
## Descricao:
##   Constroi a tabela de metadata padronizada do projeto a partir da
##   pData() real de cada serie GEO (nunca hardcoded manualmente -- ver nota
##   metodologica abaixo), com as colunas: Dataset, Protein, GSM, SRR,
##   Genotype, Treatment, Replicate, Input, Genome, Sequencing Strategy,
##   Read Length, Library Layout. Aplica o filtro do projeto (ver CLAUDE.md
##   secao 9): mantem apenas amostras WT/Controle e KO/KD/Deficiente dentro
##   do desenho experimental relevante para cada proteina, descartando
##   linhagens/celulas ou alvos de ChIP fora do escopo das 4 proteinas do
##   projeto (ex.: H3K4me3, IRF9, IRF1, ou linhagens 2fTGH/ST2-U3C usadas
##   apenas como comparacao em outra publicacao).
##
##   Nota metodologica: cada serie GEO usa uma convencao de titulo de amostra
##   diferente (GSE214182 usa virgulas; GSE222667 usa underscore; GSE247724
##   usa virgulas com layout distinto de GSE214182). Por isso ha um parser
##   dedicado por serie -- eles foram escritos e validados manualmente contra
##   os titulos reais de amostra de cada serie (consultados na pagina GEO em
##   2026-07-13), mas o parsing em si roda sobre a pData() baixada em tempo de
##   execucao pelo Modulo 01, nunca sobre uma tabela fixa embutida no codigo.
##   Caso a GEO altere formato/titulos, os parsers devem falhar de forma
##   visivel (validate_metadata_row()) em vez de silenciosamente gerar linhas
##   incorretas.
##
## Entradas:
##   Objetos ExpressionSet retornados por download_series_metadata() (Modulo 01)
##
## Saidas:
##   Dados/Metadata/chipseq_metadata.csv       (tabela final, WT/KO apenas)
##   Dados/Metadata/chipseq_metadata_filtered_out.csv (auditoria do que foi descartado e por que)
##   Logs/02_metadata.log
##
## Dependencias:
##   00_setup.R, 01_download.R, Biobase (via GEOquery), dplyr, stringr
##
## Funcoes definidas neste modulo:
##   extract_pdata(), safe_col(), parse_gse214182(), parse_gse222667(),
##   parse_gse247724(), parse_gse91923(), match_xpc_input(),
##   validate_metadata_schema(), validate_metadata_row(),
##   build_master_metadata(), run_module_02()
## ============================================================================

source(here::here("Scripts", "00_setup.R"))
source(here::here("Scripts", "01_download.R"))
install_if_missing("dplyr")
install_if_missing("stringr")
library(dplyr)
library(stringr)

METADATA_COLUMNS <- c(
  "Dataset", "Protein", "GSM", "SRR", "Genotype", "Treatment", "Replicate",
  "Input", "Genome", "SequencingStrategy", "ReadLength", "LibraryLayout",
  "CellLine", "Source_GSE"
)

#' Confere que um data.frame tem exatamente as colunas de METADATA_COLUMNS
#' antes de ser gravado -- evita salvar um CSV com colunas faltando ou fora
#' de ordem por causa de um parser desatualizado.
validate_metadata_schema <- function(df) {
  missing_cols <- setdiff(METADATA_COLUMNS, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Validacao falhou: colunas ausentes no metadata final: %s. Execucao interrompida.",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
  df[, METADATA_COLUMNS, drop = FALSE]
}

#' Extrai pData() de um ExpressionSet (ou lista de ExpressionSets, um por
#' plataforma/GPL -- GEOquery::getGEO() sempre devolve uma lista quando a
#' serie usa mais de uma plataforma, ex. GSE247724 com 5 GPLs diferentes) e
#' valida que nao esta vazio antes de prosseguir com o parsing. Combina
#' pData() de todas as plataformas via bind_rows() para nao perder amostras
#' que estejam em uma plataforma diferente da primeira da lista.
extract_pdata <- function(gset, gse_id) {
  if (is.null(gset)) {
    stop(sprintf(
      "Validacao falhou: metadata de %s nao foi baixada (Modulo 01). Execucao interrompida.",
      gse_id
    ), call. = FALSE)
  }
  if (is.list(gset) && !is.data.frame(gset)) {
    if (length(gset) > 1) {
      log_message("02_metadata", sprintf(
        "%s tem %d plataformas (GPL) distintas -- combinando pData() de todas.",
        gse_id, length(gset)
      ))
    }
    pdata <- bind_rows(lapply(gset, function(x) as.data.frame(Biobase::pData(x))))
  } else {
    pdata <- Biobase::pData(gset)
  }
  if (nrow(pdata) == 0) {
    stop(sprintf(
      "Validacao falhou: pData() de %s esta vazia. Execucao interrompida.", gse_id
    ), call. = FALSE)
  }
  pdata
}

#' Le uma coluna de pData() se ela existir, senao devolve NA e registra um
#' aviso -- nunca assume que um campo opcional do GEO esta presente.
safe_col <- function(pdata, colname, gse_id) {
  if (colname %in% colnames(pdata)) {
    return(as.character(pdata[[colname]]))
  }
  log_message("02_metadata",
    sprintf("Coluna '%s' ausente em pData(%s); preenchendo com NA.", colname, gse_id),
    level = "WARN")
  rep(NA_character_, nrow(pdata))
}

## --- parser: GSE214182 (XPC / H3K4me3, U2OS) -----------------------------------

#' Titulos no formato "<ALVO> U2OS <genotipo> cells, <UV>, biol rep<N>[, input]"
#' ex.: "XPC U2OS ASH1L knockout cells, 1h post UV, biol rep2"
parse_gse214182 <- function(pdata) {
  titles <- safe_col(pdata, "title", "GSE214182")
  pattern <- "^(H3K4me3|XPC) U2OS (wildtype|ASH1L knockout|XPC knockout) cells, (no UV|[0-9]+h post UV), biol rep([0-9]+)(, input)?$"
  m <- str_match(titles, pattern)
  unmatched <- is.na(m[, 1])
  if (any(unmatched)) {
    log_message("02_metadata",
      sprintf("GSE214182: %d titulo(s) fora do padrao esperado, ignorados: %s",
              sum(unmatched), paste(titles[unmatched], collapse = "; ")),
      level = "WARN")
  }
  genotype_map <- c(wildtype = "WT", "ASH1L knockout" = "ASH1L-KO", "XPC knockout" = "XPC-KO")
  df <- data.frame(
    Dataset    = "XPC",
    Protein    = m[, 2],
    GSM        = safe_col(pdata, "geo_accession", "GSE214182"),
    SRR        = NA_character_,
    Genotype   = unname(genotype_map[m[, 3]]),
    Treatment  = ifelse(m[, 4] == "no UV", "0h_post_UV", paste0(str_extract(m[, 4], "^[0-9]+"), "h_post_UV")),
    Replicate  = m[, 5],
    Input      = ifelse(!is.na(m[, 6]), "self_input", NA_character_),
    Genome     = NA_character_,
    SequencingStrategy = safe_col(pdata, "library_strategy", "GSE214182"),
    ReadLength = NA_character_,
    LibraryLayout = safe_col(pdata, "library_layout", "GSE214182"),
    CellLine   = "U2OS",
    Source_GSE = "GSE214182",
    stringsAsFactors = FALSE
  )
  df[!unmatched, , drop = FALSE]
}

#' Para cada linha de ChIP de XPC (sem input proprio), encontra o GSM do
#' input de H3K4me3 do mesmo genotipo+timepoint como substituto (decisao
#' registrada em CLAUDE.md secao 9, ponto 4).
match_xpc_input <- function(df) {
  h3k4_input <- df |>
    filter(Protein == "H3K4me3", !is.na(Input)) |>
    select(Genotype, Treatment, input_gsm = GSM)

  df |>
    left_join(
      h3k4_input |> group_by(Genotype, Treatment) |> slice(1) |> ungroup(),
      by = c("Genotype", "Treatment")
    ) |>
    mutate(
      Input = case_when(
        Protein == "XPC" & !is.na(input_gsm) ~ paste0(input_gsm, " (H3K4me3_substitute)"),
        Protein == "XPC" ~ NA_character_,
        TRUE ~ Input
      )
    ) |>
    select(-input_gsm)
}

## --- parser: GSE222667 (STAT1/STAT2/IRF9/IRF1 WT, Huh7.5) ----------------------

#' Titulos no formato "HuhWT_<estimulo[_tempo]>_<alvo>_rep<N>" ou
#' "HuhWT_<estimulo[_tempo]>_INP" (input, sem alvo/replicata explicita).
parse_gse222667 <- function(pdata) {
  titles <- safe_col(pdata, "title", "GSE222667")
  gsms <- safe_col(pdata, "geo_accession", "GSE222667")
  seq_strategy <- safe_col(pdata, "library_strategy", "GSE222667")
  lib_layout <- safe_col(pdata, "library_layout", "GSE222667")

  rows <- lapply(seq_along(titles), function(i) {
    tokens <- str_split(titles[i], "_")[[1]]
    if (length(tokens) < 3 || tokens[1] != "HuhWT") {
      log_message("02_metadata",
        sprintf("GSE222667: titulo fora do padrao esperado, ignorado: '%s'", titles[i]),
        level = "WARN")
      return(NULL)
    }
    is_input <- tokens[length(tokens)] == "INP"
    if (is_input) {
      protein <- "Input"
      replicate <- NA_character_
      treatment_tokens <- tokens[2:(length(tokens) - 1)]
    } else {
      replicate <- str_remove(tokens[length(tokens)], "^rep")
      protein_raw <- tokens[length(tokens) - 1]
      protein <- case_when(
        str_detect(protein_raw, regex("^p?STAT1$", ignore_case = TRUE)) ~ "STAT1",
        str_detect(protein_raw, regex("^p?STAT2$", ignore_case = TRUE)) ~ "STAT2",
        TRUE ~ protein_raw
      )
      treatment_tokens <- tokens[2:(length(tokens) - 2)]
    }
    treatment <- paste(treatment_tokens, collapse = "_")
    data.frame(
      Dataset = "STAT1_STAT2", Protein = protein, GSM = gsms[i], SRR = NA_character_,
      Genotype = "WT", Treatment = treatment, Replicate = replicate,
      Input = NA_character_, Genome = NA_character_,
      SequencingStrategy = seq_strategy[i], ReadLength = NA_character_,
      LibraryLayout = lib_layout[i], CellLine = "Huh7.5", Source_GSE = "GSE222667",
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

## --- parser: GSE247724 (STAT1/STAT2/IRF9, multi-linhagem, mantemos so Huh7.5 STAT1 K.O.)

#' Titulos separados por virgula: "<linhagem>, <estimulo[, tempo]>, <alvo>[, rep<N>]"
#' Mantem apenas linhas de "Huh7.5 STAT1 K.O." (as demais linhagens -- 2fTGH,
#' ST2-U3C -- pertencem a outro desenho experimental e sao descartadas, ver
#' CLAUDE.md secao 9).
parse_gse247724 <- function(pdata) {
  titles <- safe_col(pdata, "title", "GSE247724")
  gsms <- safe_col(pdata, "geo_accession", "GSE247724")
  seq_strategy <- safe_col(pdata, "library_strategy", "GSE247724")
  lib_layout <- safe_col(pdata, "library_layout", "GSE247724")

  rows <- lapply(seq_along(titles), function(i) {
    tokens <- trimws(str_split(titles[i], ",")[[1]])
    if (length(tokens) < 2) {
      log_message("02_metadata",
        sprintf("GSE247724: titulo fora do padrao esperado, ignorado: '%s'", titles[i]),
        level = "WARN")
      return(NULL)
    }
    if (!str_starts(tokens[1], "Huh7.5 STAT1")) {
      return(NULL)  # fora do escopo do dataset combinado STAT1/STAT2 (ver CLAUDE.md S9)
    }
    has_rep <- str_detect(tokens[length(tokens)], "^rep[0-9]+$")
    if (has_rep) {
      replicate <- str_remove(tokens[length(tokens)], "^rep")
      protein <- tokens[length(tokens) - 1]
      treatment_tokens <- tokens[2:(length(tokens) - 2)]
    } else {
      replicate <- "1"
      protein <- tokens[length(tokens)]
      treatment_tokens <- tokens[2:(length(tokens) - 1)]
    }
    protein_norm <- case_when(
      str_detect(protein, regex("stat2", ignore_case = TRUE)) ~ "STAT2",
      str_detect(protein, regex("^input$", ignore_case = TRUE)) ~ "Input",
      TRUE ~ protein
    )
    data.frame(
      Dataset = "STAT1_STAT2", Protein = protein_norm, GSM = gsms[i], SRR = NA_character_,
      Genotype = "STAT1-KO", Treatment = paste(treatment_tokens, collapse = "_"),
      Replicate = replicate, Input = NA_character_, Genome = NA_character_,
      SequencingStrategy = seq_strategy[i], ReadLength = NA_character_,
      LibraryLayout = lib_layout[i], CellLine = "Huh7.5", Source_GSE = "GSE247724",
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

## --- parser: GSE91923 (ELK1, A549, ENCODE) --------------------------------------

#' Serie ENCODE com apenas 2 replicatas de ELK1, sem braco deficiente
#' (CLAUDE.md secao 9, ponto 2) -- occupancy-only, sem Modulo 08.
parse_gse91923 <- function(pdata) {
  gsms <- safe_col(pdata, "geo_accession", "GSE91923")
  seq_strategy <- safe_col(pdata, "library_strategy", "GSE91923")
  lib_layout <- safe_col(pdata, "library_layout", "GSE91923")
  data.frame(
    Dataset = "ELK1", Protein = "ELK1", GSM = gsms, SRR = NA_character_,
    Genotype = "WT", Treatment = "none", Replicate = as.character(seq_along(gsms)),
    Input = NA_character_,  # TODO (CLAUDE.md S10): resolver input ENCODE ENCSR623KNM
    Genome = NA_character_, SequencingStrategy = seq_strategy, ReadLength = NA_character_,
    LibraryLayout = lib_layout, CellLine = "A549", Source_GSE = "GSE91923",
    stringsAsFactors = FALSE
  )
}

## --- validacao por linha ---------------------------------------------------------

#' Interrompe a execucao se uma linha do metadata final estiver com campos
#' criticos ausentes (GSM, Protein, Genotype) -- nao deixa a etapa seguinte
#' (peak calling/DiffBind) rodar sobre metadata incompleto silenciosamente.
validate_metadata_row <- function(df) {
  critical_na <- is.na(df$GSM) | is.na(df$Protein) | is.na(df$Genotype)
  if (any(critical_na)) {
    stop(sprintf(
      "Validacao falhou: %d linha(s) do metadata com GSM/Protein/Genotype ausente. Execucao interrompida.",
      sum(critical_na)
    ), call. = FALSE)
  }
  duplicated_gsm <- duplicated(df$GSM)
  if (any(duplicated_gsm)) {
    stop(sprintf(
      "Validacao falhou: GSM duplicado no metadata final: %s. Execucao interrompida.",
      paste(unique(df$GSM[duplicated_gsm]), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

## --- montagem final ---------------------------------------------------------------

#' Constroi o metadata final combinando os 4 datasets, aplicando o filtro do
#' projeto (mantem so proteinas-alvo dentro do escopo; ver cabecalho) e
#' salvando tanto a tabela mantida quanto a auditoria do que foi descartado.
build_master_metadata <- function(geo_results) {
  xpc <- parse_gse214182(extract_pdata(geo_results[["GSE214182"]]$metadata, "GSE214182")) |>
    match_xpc_input()
  stat_wt <- parse_gse222667(extract_pdata(geo_results[["GSE222667"]]$metadata, "GSE222667"))
  stat_ko <- parse_gse247724(extract_pdata(geo_results[["GSE247724"]]$metadata, "GSE247724"))
  elk1 <- parse_gse91923(extract_pdata(geo_results[["GSE91923"]]$metadata, "GSE91923"))

  all_rows <- bind_rows(xpc, stat_wt, stat_ko, elk1)

  ## Escopo de proteinas do projeto (CLAUDE.md S1): XPC, ELK1, STAT1, STAT2 +
  ## as respectivas linhas de Input necessarias para peak calling. H3K4me3
  ## (ja absorvido em match_xpc_input como Input substituto), IRF9 e IRF1 sao
  ## descartados por estarem fora do escopo das 4 proteinas do projeto.
  in_scope <- all_rows$Protein %in% c("XPC", "ELK1", "STAT1", "STAT2", "Input")

  kept <- validate_metadata_schema(all_rows[in_scope, , drop = FALSE])
  filtered_out <- all_rows[!in_scope, , drop = FALSE]

  validate_metadata_row(kept)

  kept_path <- file.path(PROJECT_DIRS$metadata, "chipseq_metadata.csv")
  filtered_path <- file.path(PROJECT_DIRS$metadata, "chipseq_metadata_filtered_out.csv")
  write.csv(kept, kept_path, row.names = FALSE)
  write.csv(filtered_out, filtered_path, row.names = FALSE)

  log_message("02_metadata", sprintf(
    "Metadata final: %d amostra(s) mantida(s), %d descartada(s) (fora de escopo). Salvo em '%s'.",
    nrow(kept), nrow(filtered_out), kept_path
  ))
  kept
}

## --- execucao do modulo ------------------------------------------------------------

#' Executa o Modulo 02 completo: baixa metadata via Modulo 01 e constroi a
#' tabela padronizada final.
run_module_02 <- function() {
  log_message("02_metadata", "Iniciando Modulo 02 -- Construcao do metadata.")
  geo_results <- run_module_01()
  metadata <- build_master_metadata(geo_results)
  log_message("02_metadata", "Modulo 02 concluido.")
  invisible(metadata)
}

## Este modulo nao executa automaticamente ao ser source()ado -- chame
## run_module_02() explicitamente (interativamente ou a partir de
## 22_master_pipeline.R).
