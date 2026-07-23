# Pipeline Integrativo de ChIP-seq (XPC, ELK1, STAT1, STAT2, IRF9)

Pipeline modular em R para análise integrativa de ChIP-seq, do download bruto no GEO
até redes regulatórias multi-fator de transcrição e metanálise cross-proteína. Ver
[`CLAUDE.md`](CLAUDE.md) para o registro completo do projeto (histórico, decisões
metodológicas, dependências, checklist de validação e pendências) e
[`Analises/RESUMO_METANALISE.md`](Analises/RESUMO_METANALISE.md) para os resultados
da metanálise e [`Analises/README.md`](Analises/README.md) para o índice de pastas
(análise do XPC vs metanálise cross-proteína).

## O experimento

O projeto integra ChIP-seq de 5 proteínas humanas (U2OS/Huh7.5/A549) de 4 estudos
públicos do GEO, buscando a interseção entre o reparo de DNA acoplado à transcrição
(XPC) e a resposta a interferon (eixo STAT1/STAT2/IRF9, complexo ISGF3), com ELK1
(fator de transcrição ubíquo, ENCODE) como referência de "ruído de fundo" esperado.

| Proteína(s) | GSE | Linhagem | Desenho | Peak type | Differential binding (Módulo 08)? |
|---|---|---|---|---|---|
| **XPC** | GSE214182 | U2OS | WT vs **XPC-KO**, 0h/1h/3h pós-UV (ASH1L-KO como comparação secundária) | Broad | Sim (WT vs XPC-KO) |
| **ELK1** | GSE91923 | A549 (ENCODE ENCSR623KNM) | 2 réplicas, sem KO/KD, constitutivo | Narrow | Não (sem braço deficiente) |
| **STAT1** | GSE222667 | Huh7.5 | WT, untreated + IFNα/IFNγ em múltiplos timepoints | Narrow | Não (sem STAT1-ChIP em fundo deficiente) |
| **STAT2** | GSE222667 (WT) + GSE247724 (STAT1-KO) | Huh7.5 | WT vs **STAT1-KO**, timepoints casados (UN/2h/24h/72h IFNα) | Narrow | Sim (WT vs STAT1-KO) |
| **IRF9** | GSE222667 | Huh7.5 | WT, untreated + IFNα 2h (mesmo desenho de STAT1/STAT2) | Narrow | Não |

Nenhuma das 5 proteínas tem input próprio nativo. Adaptação metodológica: XPC/ASH1L-KO
usam H3K4me3 (mesma linhagem/timepoint) como substituto de input; XPC-KO usa esse
**mesmo substituto**, pareado por timepoint (ver `CLAUDE.md` §9.1-REVISÃO — decisão
corrigida em 2026-07-22 após o `--nolambda` original gerar background massivo); STAT1/
STAT2/IRF9/ELK1 usam input próprio do estudo original.

## Desenho experimental detalhado — amostras por proteína e por etapa

Cada proteína passou por um subconjunto diferente de etapas (nem toda amostra
baixada chegou a ser peak-called; nem todo peak-calling entrou na metanálise
principal). A tabela abaixo mostra exatamente onde cada amostra parou.

### XPC — GSE214182 (U2OS, WT vs XPC-KO vs ASH1L-KO, 0h/1h/3h pós-UV)

| Genótipo | Timepoint | GSMs | Alinhado/filtrado? | Peak-called? | N picos | Usado em |
|---|---|---|---|---|---|---|
| WT | 0h | 715,716,717,718 | Sim | Sim (broad, com input) | 1462, 0, 0, 7 | Consenso pooled (metanálise) + `XPC_0h` (timepoints) |
| WT | 1h | 724,725,726 | Sim | Sim | 0, 0, 0 | Consenso pooled — **contribui 0 regiões** |
| WT | 3h | 732,733,734 | Sim | Sim | 0, 109, 0 | Consenso pooled + `XPC_3h` (timepoints) |
| **XPC-KO** | 0h | 722,723 | Sim | Sim (revisado 2026-07-22) | 0, 4 | Diffbind (WT vs KO) |
| **XPC-KO** | 1h | 730,731 | Sim | Sim (revisado) | 0, 5 | Diffbind |
| **XPC-KO** | 3h | 738,739 | Sim | Sim (revisado) | 0, 8 | Diffbind |
| ASH1L-KO | 0h/1h/3h | 719-721,727-729,735-737 | FASTQ baixado (0h chegou a ser alinhado, nunca filtrado; 1h/3h nem alinhados) | Não | — | Não usado (citado no metadata como comparação secundária, nunca executado); FASTQ/BAM removidos em 2026-07-23 (liberação de espaço) |
| Input H3K4me3 (substituto) | 0h/1h/3h | 680, 686, 692 | Sim | — (usado só como `-c` do MACS3) | — | Input pareado por timepoint para WT e XPC-KO |

**Ponto de atenção nº1 — heterogeneidade do ChIP-WT**: das 10 réplicas WT, só
**3 produzem picos** (715, 718, 733); as outras 7 dão 0, inclusive **as 3 de 1h
inteiras**. Confirmado com teste de sanidade (q=0.5, muito mais permissivo) que
não é limiar — é ausência real de estrutura de enriquecimento. Consequência:
o consenso pooled usado na metanálise principal é **dominado por GSM6600715**
(0h). Ver `Analises/RESUMO_METANALISE.md` §8 para o diagnóstico completo.

**Ponto de atenção nº2 — peak calling do XPC-KO revisado**: até 2026-07-13 o
XPC-KO rodava com `--nolambda` (sem controle), gerando 250 mil+ picos de
background por amostra. Revisado em 2026-07-22 para usar o mesmo input
H3K4me3 do WT, pareado por timepoint — ver `CLAUDE.md` §9.1-REVISÃO.

### STAT1 — GSE222667 (Huh7.5, WT, sem KO)

| Condição | GSMs | Alinhado/filtrado? | Peak-called? | Usado em |
|---|---|---|---|---|
| Untreated (UN) | 563, 564 | Sim | Sim | `Analises/Metanalise/baseline_controle/` (controle) |
| IFNα 2h | 567, 568 | Sim | Sim | `Analises/Metanalise/principal_sem_normalizacao/` + `Metanalise/principal_normalizado_topN/` (**estado ativo da metanálise principal**) |
| IFNα 0.5h/8h/24h/72h | 565,566,569-574 | Sim | **Não** | Não usado — disponível para uma futura análise de cinética |
| IFNγ 0.5h/4h/24h/72h | 575-582 | Sim | **Não** | Não usado — estímulo biologicamente distinto (GAS vs ISGF3), fora do escopo atual |

### STAT2 — GSE222667 (WT) + GSE247724 (STAT1-KO)

| Condição | GSMs | Alinhado/filtrado? | Peak-called? | Usado em |
|---|---|---|---|---|
| WT untreated | 583, 584 | Sim | Sim | `Metanalise/baseline_controle/` |
| WT IFNα 2h | 587, 588 | Sim | Sim | `Metanalise/principal_sem_normalizacao/` + `Metanalise/principal_normalizado_topN/` (**estado ativo**) |
| WT IFNα 0.5h/8h/24h/72h | 585,586,589-594 | Sim | **Não** | Não usado |
| **STAT1-KO** (untreated/2h/24h/72h) | GSM7899570-591 (8 ChIP + 4 input) | Sim | **Não** | **Não usado — diffbind WT-vs-STAT1-KO nunca executado neste projeto**, apesar de `run_diffbind_standard()` (Módulo 08) já suportar esse fluxo |

### IRF9 — GSE222667 (Huh7.5, WT)

| Condição | GSMs | Alinhado/filtrado? | Peak-called? | Usado em |
|---|---|---|---|---|
| Untreated (UN) | 595, 596 | Sim | Sim | `Metanalise/baseline_controle/` |
| IFNα 2h | 599, 600 | Sim | Sim | `Metanalise/principal_sem_normalizacao/` + `Metanalise/principal_normalizado_topN/` (**estado ativo**) |
| IFNα 0.5h/8h/24h/72h, IFNγ (todos) | 597,598,601-614 | **Não** (nunca baixado) | — | Não usado |

Nota: IRF9 está listado em `chipseq_metadata_filtered_out.csv`, não no metadata
principal (filtrado por critério de QC do Módulo 02) — usado manualmente na
metanálise mesmo assim, pois os picos UN/IFNα2h existem e são de boa qualidade.

### ELK1 — GSE91923 (A549, ENCODE ENCSR623KNM)

| Condição | GSMs | Usado em |
|---|---|---|
| Constitutivo (única condição) | 2423754, 2423755 | Todas as análises (`Metanalise/principal_sem_normalizacao`, `Metanalise/principal_normalizado_topN`, `Metanalise/rede_regulatoria`) — nunca varia entre elas |

Genoma nativo hg19 → liftOver para hg38 (Módulo 12), único dataset que precisou
dessa conversão.

## Metodologia por etapa — decisões e pontos de atenção

### Peak calling (Módulo 07 — MACS3 via WSL)

- **Broad** (XPC) vs **narrow** (STAT1/STAT2/IRF9/ELK1) — `q=0.01`, broad-cutoff
  `0.1` quando aplicável (`Scripts/config/peakcalling_config.R`).
- **Input**: todas as amostras rodam **com** controle (`-c`), nunca `--nolambda`
  por omissão — a coluna `Input` do metadata decide, nunca é assumida
  (`determine_macs3_args()`). Única exceção histórica: XPC-KO, corrigida em
  2026-07-22 (ver acima).
- **Fallback `--nomodel`**: usado quando o MACS3 não consegue construir o
  modelo de picos pareados (sinal fraco demais) — usa o comprimento de
  fragmento estimado no Módulo 06 (cross-correlation) como `--extsize`.
  Aplicado em várias réplicas fracas de XPC-WT e nos 6 XPC-KO re-chamados.

### Differential binding (Módulo 08 — `csaw`/`edgeR` para XPC, `DiffBind` para STAT2)

- **XPC (WT vs XPC-KO)**: região consenso (`GenomicRanges::reduce()` da união
  WT∪KO) + contagem de reads via `csaw::regionCounts()` + teste por edgeR
  (`glmQLFTest`), **não** por presença/ausência de pico — necessário porque o
  KO tem ChIP mais fraco por design (sem proteína pra imunoprecipitar).
  **TMM adicionado em 2026-07-22** (`edgeR::calcNormFactors()`).
  Resultado atual: 1576 regiões consenso, **0 Gained / 0 Lost / 1576 Stable**
  (nenhuma diferença sobrevive a FDR<0.05) — resultado nulo honesto, não bug;
  reflexo direto da heterogeneidade do ChIP-WT (ponto de atenção acima).
- **STAT2 (WT vs STAT1-KO)**: fluxo padrão `DiffBind` (`dba.count()` →
  `dba.normalize(normalize=DBA_NORM_TMM)` → `dba.analyze()`) — **implementado
  mas nunca executado** neste projeto (os BAMs do STAT1-KO existem, os picos
  não foram gerados).
- **ELK1 e STAT1**: fora do Módulo 08 por design — nenhum dos dois tem braço
  deficiente/KO disponível.

### Anotação (Módulo 09 — ChIPseeker)

Duas visões por proteína, sempre `TxDb.Hsapiens.UCSC.hg38.knownGene` +
`org.Hs.eg.db`, `tssRegion=c(-3000,3000)`:
- **Nearest gene**: gene mais próximo de cada pico, qualquer distância.
- **Promotor**: subconjunto restrito a `annotation` começando com
  `"Promoter"` (±3kb do TSS) — mais específico, usado como visão principal
  nas figuras/tabelas de interseção.

### Metanálise cross-proteína (`Analises/Metanalise/principal_sem_normalizacao/` e `Metanalise/principal_normalizado_topN/`)

- **Seleção de estado por proteína** (revisada 2026-07-22): XPC = pooled WT
  (10 amostras, 3 timepoints); STAT1/STAT2/IRF9 = **só IFNα 2h** (2 réplicas
  cada); ELK1 = constitutivo (2 réplicas). O "untreated" foi removido do pool
  principal porque contribui sinal quase nulo (STAT2 UN = 5 regiões no genoma
  inteiro) — fica isolado em `Metanalise/baseline_controle/` como controle.
- **Duas normalizações de tamanho de gene-set**: `Metanalise/principal_sem_normalizacao/` usa todos os
  picos do estado ativo; `Metanalise/principal_normalizado_topN/` restringe a **top-1000 picos por
  proteína** (rankeados por `signalValue` do MACS3) antes de anotar — nivela
  STAT1/STAT2 (que saturam com ~19-20 mil genes-alvo no nearest gene) ao
  mesmo patamar de XPC (817)/IRF9 (921)/ELK1 (212, inalterado por já ter <1000
  picos).
- **Ponto de atenção — desbalanço STAT vs demais**: mesmo restrito a IFNα2h,
  STAT1/STAT2 continuam gerando dezenas de milhares de picos (ativação
  massiva do complexo ISGF3), muito mais que XPC/IRF9/ELK1. Interseções
  amplas envolvendo STAT (sem restringir a promotor ou usar o topN) tendem a
  ser dominadas pelo volume de picos do STAT, não necessariamente por
  co-ligação biológica real — por isso a versão `Metanalise/principal_normalizado_topN/` existe, como
  checagem de robustez.
- **Núcleo XPC∩STAT1∩STAT2∩IRF9 (8 genes)**: estável entre o design antigo
  (UN+IFNα2h pooled) e o novo (só IFNα2h) — não é um artefato da seleção de
  estado.

### XPC por timepoint (`Analises/XPC/timepoints/`)

Metanálise separada por timepoint (0h e 3h como âncoras distintas contra o
eixo interferon IFNα2h) em vez do consenso pooled único. 1h fica documentado
como vazio (0 picos nas 3 réplicas — mesma causa do peak calling por amostra,
não é um artefato da metanálise). 0h tem interseção muito maior com o eixo
interferon que 3h simplesmente porque o consenso pooled do XPC é dominado por
GSM6600715 (0h) — não é evidência de que XPC∩interferon seja específico do
timepoint 0h biologicamente, é um reflexo de qual réplica teve ChIP funcional.

### Rede regulatória (`Analises/Metanalise/rede_regulatoria/`)

Rede bipartida Proteína→Região→Gene sobre a matriz de ocupação binária
(Módulo 14, `XPC_WT`/`STAT1_WT`/`STAT2_WT`/`IRF9_WT`/`ELK1_WT` — pool "WT"
histórico do Módulo 13/14, **anterior** à revisão de estado da metanálise
principal, ver nota abaixo). Duas versões: completa (score≥2, ilegível como
imagem estática por volume) e focada em `occupancy_score≥3` (1462 regiões,
rótulos só nos 9 genes de score 4 + top-5 de score 3 por `signalValue`).

**Nota de consistência**: as colunas `STAT1_WT`/`STAT2_WT`/`IRF9_WT` da matriz
de ocupação (Módulos 13/14/17) usam o pool **UN+IFNα2h** (4 réplicas cada),
não a seleção "só IFNα2h" adotada na metanálise principal em 2026-07-22 — a
rede não foi re-gerada com o novo pool porque o pedido explícito da revisão
era sobre a metanálise de gene-sets, não sobre a rede. Isso é uma
inconsistência conhecida entre as duas análises, documentada aqui para quem
for comparar números diretamente entre `Analises/Metanalise/rede_regulatoria/` e `Analises/Metanalise/principal_sem_normalizacao/`.

## Resumo do pipeline (22 módulos)

`Scripts/00_setup.R` → `22_master_pipeline.R`, cada módulo com uma responsabilidade
única (download, metadata, alinhamento via WSL/bowtie2, filtragem, ChIP-QC, peak
calling via WSL/MACS3, differential binding, anotação, enriquecimento funcional,
GRanges/liftOver, universo regulatório, matriz de ocupação, hotspots, rede bipartida,
rede de pathways, validação). Substituições R-nativas usadas onde os pacotes
Bioconductor originalmente previstos (`Rsubread`, `Rbowtie2`) não funcionavam no
Windows sob Smart App Control — ver `CLAUDE.md` §5/§9 para o histórico completo dessas
adaptações.

## Como rodar

Todos os scripts assumem que o diretório de trabalho do R é a raiz deste projeto,
resolvido via `here::here()` — nunca caminhos absolutos.

```r
# 1. Instalar/checar dependências
source("Scripts/00_setup.R")

# 2. Baixar dados do GEO (metadata, FASTQ ou arquivos processados)
source("Scripts/01_download.R")

# 3. Construir metadata padronizado (WT vs deficiente, sem amostras tratadas)
source("Scripts/02_metadata.R")

# ... demais módulos em Scripts/, executados em ordem numérica ...

# Pipeline completo:
source("Scripts/22_master_pipeline.R")
```

Os scripts que geram a metanálise revisada (2026-07-22) — re-call de picos do
XPC-KO, diffbind com TMM, metanálise principal full/top-N, metanálise por
timepoint, rede bipartida score≥3, figuras comparativas — foram escritos como
scripts de orquestração ad-hoc (não fazem parte do `22_master_pipeline.R` numerado,
por reutilizarem as funções dos módulos 07/08/09/10/18 com parâmetros específicos
dessa revisão); ver `CLAUDE.md` §9.1-REVISÃO e `Analises/RESUMO_METANALISE.md` §0
para o que cada um faz.

## Estrutura de pastas

```
XPC_metaChIP/
├── Dados/
│   ├── GEO/         # metadata bruto do GEO (SOFT/series matrix), nunca editado
│   ├── FASTQ/        # reads brutos (não versionado)
│   ├── BAM/          # alinhamentos (não versionado)
│   ├── Peaks/        # picos MACS3
│   ├── BigWig/       # tracks de cobertura (não versionado)
│   └── Metadata/     # metadata padronizado (CSV, versionado)
├── Analises/          # metanálise/enriquecimento/rede consolidados (ver abaixo)
├── Figuras/           # figuras dos módulos individuais (QC, annotation, alignment...)
├── Arquivos/          # arquivos intermediários dos módulos individuais (RDS, tabelas)
├── Scripts/           # 00_setup.R ... 22_master_pipeline.R
├── Logs/              # logs de execução por módulo
├── CLAUDE.md
└── README.md
```

### Estrutura de `Analises/` (reorganizada em 2026-07-22, refinada em seguida)

Separada em dois grandes ramos — **`XPC/`** (análises de uma proteína só,
não cruzam com as demais) e **`Metanalise/`** (tudo que cruza 2+ proteínas)
— mais `qc_comparativo/`, que é transversal aos dois. Cada subpasta tem seu
próprio `README.md` dizendo se tem normalização (TMM/top-N) e quais
amostras/estado usa:

```
Analises/
├── RESUMO_METANALISE.md
│
├── XPC/                              # SÓ XPC (nenhuma interseção com outras proteínas)
│   ├── individual/                    # anotação + enriquecimento por timepoint e geral
│   ├── timepoints/                    # XPC 0h vs 3h, cada um contra o eixo interferon
│   └── diffbind/                      # WT vs XPC-KO
│       ├── atual_input_pareado_TMM/   # vigente: input pareado + TMM
│       └── ANTES_nolambda_semTMM/     # legado: --nolambda, sem TMM
│
├── Metanalise/                        # TUDO que cruza 2+ proteínas
│   ├── principal_sem_normalizacao/    # estado ativo (IFNα2h), todos os picos
│   │   ├── gene_sets_*, jaccard_*, upset_*.png, venn_*.png
│   │   └── interseccoes/<combo>/      # 1 pasta por combinação (2-5 proteínas):
│   │       genes_{nearest,promotor}.txt + {nearest,promotor}_{go,kegg,reactome,hallmark}.csv + dotplots
│   ├── principal_normalizado_topN/    # igual acima, top-1000 picos/proteína por signalValue
│   ├── baseline_controle/             # controle: só amostras untreated
│   ├── design_anterior_pre_revisao/   # design ORIGINAL (pré 2026-07-22, UN+IFNα2h pooled) — legado
│   │   └── interseccoes/<combo>/
│   ├── nucleo_XPC_interferon/         # destaque: os 8 genes XPC∩STAT1∩STAT2∩IRF9 + enriquecimento
│   └── rede_regulatoria/              # rede bipartida Proteína→Região→Gene (score≥2 e score≥3)
│
└── qc_comparativo/                    # figuras pré/pós normalização + painel de métricas por amostra
                                        # (cruza XPC e Metanalise, por isso fica fora dos dois ramos)
```

Dados brutos grandes (FASTQ/BAM/BigWig) nunca são versionados — ver `.gitignore` na
raiz do repositório.

**Limpeza de espaço (2026-07-23):** removidos ~154GB de dado bruto que não
entra em nenhuma análise atual — FASTQ/BAM das 54 amostras STAT1/STAT2/
STAT1-KO nunca peak-called (as mesmas listadas como "Não usado" nas tabelas
acima), o FASTQ/BAM parcial do ASH1L-KO, os `.sorted.bam` (pré-filtro,
Módulo 05) de todas as amostras que **são** usadas — só o `.filtered.bam`
final é lido por qualquer módulo a partir do 06 — e as 2 réplicas cruas de
input do ELK1 já mescladas em `ELK1input_pooled.filtered.bam`. Nada usado em
análise atual foi tocado; recuperar os dados removidos exigiria baixar/
realinhar do zero (várias horas por amostra).
