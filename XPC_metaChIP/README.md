# Pipeline Integrativo de ChIP-seq (XPC, ELK1, STAT1, STAT2, IRF9)

Pipeline modular em R para análise integrativa de ChIP-seq, do download bruto no GEO
até redes regulatórias multi-fator de transcrição e metanálise cross-proteína. Ver
[`CLAUDE.md`](CLAUDE.md) para o registro completo do projeto (histórico, decisões
metodológicas, dependências, checklist de validação e pendências) e
[`Analises/RESUMO_METANALISE.md`](Analises/RESUMO_METANALISE.md) para os resultados
da metanálise.

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

## Resumo do pipeline (22 módulos)

`Scripts/00_setup.R` → `22_master_pipeline.R`, cada módulo com uma responsabilidade
única (download, metadata, alinhamento via WSL/bowtie2, filtragem, ChIP-QC, peak
calling via WSL/MACS3, differential binding, anotação, enriquecimento funcional,
GRanges/liftOver, universo regulatório, matriz de ocupação, hotspots, rede bipartida,
rede de pathways, validação). Substituições R-nativas usadas onde os pacotes
Bioconductor originalmente previstos (`Rsubread`, `Rbowtie2`) não funcionavam no
Windows sob Smart App Control — ver `CLAUDE.md` §5/§9 para o histórico completo dessas
adaptações.

## Metanálise cross-proteína

A metanálise (`Analises/`) cruza os 5 conjuntos de picos por gene-alvo (gene mais
próximo e promotor ±3kb), Jaccard, interseções completas (todas as 26 combinações de
2 a 5 proteínas) e enriquecimento funcional (GO/KEGG/Reactome/Hallmark) por
combinação. Revisada em 2026-07-22 (ver seção "Estrutura de `Analises/`" abaixo e
`Analises/RESUMO_METANALISE.md` §0 para o que mudou e por quê):

- **Estado ativo por proteína**: XPC = pooled WT (10 amostras, 0h+1h+3h pós-UV);
  STAT1/STAT2/IRF9 = **só IFNα 2h** (estado ativo puro — o "untreated" foi removido
  do pool principal por ter sinal quase nulo, ver `meta_baseline/`); ELK1 =
  constitutivo (2 ENCODE).
- **Duas normalizações**: `meta_geral/` (todos os picos do estado ativo) e
  `meta_topN/` (top-1000 picos por `signalValue` MACS3 — nivela STAT1/STAT2, que
  saturam com 40-62 mil picos no IFNα2h, ao mesmo patamar de XPC/IRF9/ELK1).
- **XPC por timepoint** (`Analises/XPC/timepoints/`): 0h e 3h como âncoras
  separadas contra o eixo interferon (1h documentado como vazio — 0 picos nas 3
  réplicas, causa no peak calling por amostra, não na metanálise).
- **Baseline/controle** (`Analises/meta_baseline/`): só amostras não-estimuladas —
  confirma que STAT1/STAT2/IRF9 dependem de ativação por interferon para ocupar
  cromatina (STAT2 untreated = 5 regiões no genoma inteiro).

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

### Estrutura de `Analises/` (reorganizada em 2026-07-22)

Uma pasta por análise, com todos os arquivos e figuras daquela análise juntos
(sem separar por Arquivos/Figuras):

```
Analises/
├── RESUMO_METANALISE.md         # relatório completo da metanálise
├── meta_geral/                   # metanálise principal (IFNα2h, full)
│   ├── gene_sets_*, jaccard_*, upset_*.png, venn_*.png
│   └── interseccoes/<combo>/     # 1 pasta por combinação (2-5 proteínas):
│       genes_{nearest,promotor}.txt + {nearest,promotor}_{go,kegg,reactome,hallmark}.csv + dotplots
├── meta_topN/                    # igual acima, top-1000 picos/proteína por signalValue
├── meta_baseline/                 # controle: só amostras untreated
├── meta_geral_ANTES_revisao/      # análise ORIGINAL (pré 2026-07-22, UN+IFNα2h pooled) — legado, não apagado
│   └── interseccoes/<combo>/
├── XPC/                           # tudo que é análise individual do XPC
│   ├── individual/                # enriquecimento/anotação por timepoint e geral
│   ├── timepoints/                 # metanálise XPC 0h/3h vs eixo interferon
│   ├── nucleo/                     # núcleo XPC∩STAT1∩STAT2∩IRF9 (8 genes) + enriquecimento
│   └── diffbind/                   # WT vs XPC-KO (re-call + TMM) + backups pré-revisão
├── rede/                           # rede bipartida (score≥2 completa + score≥3 focada) + GraphML
└── qc_comparativo/                 # figuras pré/pós normalização + painel de métricas por amostra
```

Dados brutos grandes (FASTQ/BAM/BigWig) nunca são versionados — ver `.gitignore` na
raiz do repositório.
