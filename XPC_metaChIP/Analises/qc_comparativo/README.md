# qc_comparativo — figuras ANTES/DEPOIS de cada normalização + métricas por amostra

Não é uma análise em si — são as figuras/tabelas que **comparam** os efeitos
das normalizações aplicadas em 2026-07-22, usadas para validar que cada
correção fez o esperado:

- **`peak_count_KO_prepos.{png,csv}`** — nº de picos do XPC-KO, `--nolambda`
  (antes) vs input H3K4me3 pareado (depois). Ver `../XPC/diffbind/README.md`.
- **`diffbind_{categorias,MA,logCPM_boxplot}_prepos.png`** — categorias
  Gained/Lost/Stable, MA-plot e distribuição de logCPM do diffbind do XPC,
  antes (`--nolambda`, sem TMM) vs depois (input pareado + TMM).
- **`gene_sets_prepos_topN.{png,csv}`**, **`jaccard_prepos_topN_
  {nearest,promotor}.png`** — tamanho dos gene-sets e Jaccard, `Metanalise/principal_sem_normalizacao/`
  (full) vs `Metanalise/principal_normalizado_topN/` (top-1000).
- **`fragmentsize_comparativo.png`**, **`qc_metricas_painel.png`**,
  **`qc_metricas_todas_amostras.csv`** — comprimento de fragmento e painel
  multi-métrica (SSD, fragment length, nº de picos, profundidade pós-filtro,
  fração removida na filtragem) por amostra, de **todas** as ~87 amostras
  processadas no projeto (não só as usadas na metanálise principal).
- **`PENDENCIA_FRiP.txt`** — FRiP continua `NA` para a maioria das amostras
  (ChIPQC rodou sem picos supridos); regenerar exigiria rodar
  `ChIPQCsample()` com picos para 80+ amostras — estimado em 50-90h de
  execução serial, não feito nesta revisão.
