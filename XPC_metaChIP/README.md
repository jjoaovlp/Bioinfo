# Pipeline Integrativo de ChIP-seq (XPC, ELK1, STAT1, STAT2)

Pipeline modular em R para análise integrativa de ChIP-seq, do download bruto no GEO
até redes regulatórias multi-fator de transcrição. Ver [`CLAUDE.md`](CLAUDE.md) para o
registro completo do projeto (histórico, decisões metodológicas, dependências, checklist
de validação e pendências).

## Como rodar

Todos os scripts assumem que o diretório de trabalho do R é a raiz deste projeto
(`ProjetosCode/`), resolvido via `here::here()` — nunca caminhos absolutos.

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

## Estrutura de pastas

```
ProjetosCode/
├── Dados/
│   ├── GEO/        # metadata bruto do GEO (SOFT/series matrix), nunca editado
│   ├── FASTQ/       # reads brutos (não versionado)
│   ├── BAM/         # alinhamentos (não versionado)
│   ├── Peaks/       # picos MACS3
│   ├── BigWig/       # tracks de cobertura (não versionado)
│   └── Metadata/    # metadata padronizado (CSV, versionado)
├── Figuras/          # todas as figuras (PNG/PDF/SVG, 300 dpi+)
├── Arquivos/         # arquivos intermediários (RDS, GRanges, tabelas)
├── Scripts/          # 00_setup.R ... 22_master_pipeline.R
├── Logs/             # logs de execução por módulo
├── CLAUDE.md
└── README.md
```

Dados brutos grandes (FASTQ/BAM/BigWig) nunca são versionados — ver `.gitignore` na
raiz do repositório.
