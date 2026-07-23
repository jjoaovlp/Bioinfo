# XPC/diffbind — differential binding WT vs XPC-KO

Duas subpastas, mesmo desenho (WT vs XPC-KO, consenso via `GenomicRanges::
reduce()` + contagem de reads via `csaw`/`edgeR`), **normalização diferente**:

## `atual_input_pareado_TMM/` — versão vigente (2026-07-22)

- **Peak calling do XPC-KO**: input H3K4me3 do WT, pareado por timepoint
  (0h→GSM6600680, 1h→GSM6600686, 3h→GSM6600692) — não `--nolambda`.
- **Normalização**: TMM (`edgeR::calcNormFactors()`).
- **Consenso**: 1576 regiões (WT∪KO, já sem o background artefatual do KO).
- **Resultado**: 0 Gained / 0 Lost / 1576 Stable — nenhuma diferença
  sobrevive a FDR<0.05. Resultado nulo honesto (reflexo da heterogeneidade do
  ChIP-WT, não um bug), substituindo os 3283 "Gained" artefatuais da versão
  anterior.

## `ANTES_nolambda_semTMM/` — versão anterior (legado, até 2026-07-13)

- **Peak calling do XPC-KO**: `--nolambda` (sem controle) → 250 mil+ picos de
  background por amostra.
- **Sem TMM** (só escalonamento por tamanho de biblioteca via
  `csaw::asDGEList()`).
- **Consenso**: 353.522 regiões (inflado pelo background do KO).
- **Resultado**: 3283 "Gained" (mais sinal no KO que no WT — biologicamente
  implausível, artefato do peak-calling assimétrico, não biologia real).

Preservado para comparação/auditoria — ver `CLAUDE.md` §9.1-REVISÃO e
`Analises/qc_comparativo/peak_count_KO_prepos.png` /
`diffbind_categorias_prepos.png` para as figuras antes/depois.
