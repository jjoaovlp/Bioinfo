# Metanalise/design_anterior_pre_revisao — design ORIGINAL (legado, anterior a 2026-07-22)

**Não é a análise principal atual** — preservado como registro histórico/
auditoria, não apagado. A análise principal vigente é `../Metanalise/principal_sem_normalizacao/` (ou
`../Metanalise/principal_normalizado_topN/` para a versão normalizada).

**Normalização:** nenhuma (sem TMM no diffbind associado, sem top-N).

**Amostras/estado por proteína (design antigo):**
- **XPC** = pooled WT, 10 amostras (0h+1h+3h) — igual ao atual
- **STAT1/STAT2/IRF9** = **untreated + IFNα2h pooled** (4 réplicas cada) —
  DIFERENTE do atual, que usa só IFNα2h
- **ELK1** = constitutivo — igual ao atual

**Por que foi substituído (2026-07-22):** o "untreated" quase não contribuía
sinal (STAT2 UN = 5 regiões no genoma inteiro, ver `../Metanalise/baseline_controle/`) e só
diluía o desenho sem ganho biológico — a metanálise principal passou a usar
IFNα2h puro (estado ativo do ISGF3). O núcleo XPC∩STAT1∩STAT2∩IRF9 (8 genes)
é **idêntico** entre este design antigo e o atual — a mudança não alterou o
achado central, só a interpretação do que "STAT1/STAT2/IRF9" representa no
pool.

**Também aqui:** `jaccard_regiao_pairwise.csv` / `jaccard_heatmap_regiao.png`
/ `regioes_compartilhadas_todas_proteinas.csv` — Jaccard calculado no **nível
de região genômica** (não gene), do antigo Módulo 15/16
(`pairwise_overlap_table`/`multiple_overlap`), usando o mesmo pool "WT"
UN+IFNα2h. Explica por que XPC parece "não ter nada a ver" com as demais
proteínas nesse arquivo (Jaccard 0.001-0.003) — os picos de XPC raramente
caem exatamente sobre a mesma coordenada de STAT/IRF9/ELK1; a relação real só
aparece no nível de gene mais próximo (arquivos `gene_sets_*`/`genes_comuns_*`
deste mesmo diretório).

**Conteúdo:** estrutura igual a `Metanalise/principal_sem_normalizacao/` (gene_sets, jaccard de gene,
upset, venn, combinações + `interseccoes/<combo>/`), mais os 3 arquivos de
Jaccard por região citados acima.
