# CLAUDE.md — Pipeline Integrativo de ChIP-seq

Registro permanente do projeto. Atualizado a cada alteração estrutural, novo módulo,
correção de bug ou mudança metodológica.

---

## 1. Objetivo do projeto

Pipeline modular, reprodutível e documentado em R para análise integrativa de ChIP-seq
de quatro proteínas — **XPC**, **ELK1**, **STAT1**, **STAT2** — cobrindo desde o download
de dados públicos do GEO até a construção de um universo regulatório multi-proteína,
matriz de ocupação cromatínica, redes bipartidas (proteína → região → gene) e redes
funcionais (PPI/GO/Reactome). Todo o código é R, dividido em 22 módulos (`Scripts/00_*.R`
a `Scripts/22_*.R`), sem scripts monolíticos.

## 2. Organização das pastas

```
ProjetosCode/
├── Dados/
│   ├── GEO/        # metadata bruto do GEO (SOFT/series matrix) — nunca editado
│   ├── FASTQ/       # reads brutos — NÃO versionado
│   ├── BAM/         # alinhamentos — NÃO versionado
│   ├── Peaks/       # picos MACS3
│   ├── BigWig/       # tracks de cobertura — NÃO versionado
│   └── Metadata/    # metadata padronizado (CSV) — versionado
├── Figuras/          # todas as figuras (PNG/PDF/SVG, ≥300 dpi)
├── Arquivos/         # intermediários (RDS, GRanges, tabelas, sessionInfo)
├── Scripts/
│   └── config/       # arquivos de configuração (parâmetros de peak calling etc.)
├── Logs/             # logs de execução por módulo
├── CLAUDE.md
└── README.md
```

Caminhos sempre relativos via `here::here()`. Dados brutos do GEO nunca são alterados
in-place; apenas lidos.

## 3. Lista dos módulos implementados

| # | Script | Status | Descrição |
|---|--------|--------|-----------|
| 00 | `00_setup.R` | ✅ implementado e testado | Instalação de dependências (BiocManager), verificação de versões, `sessionInfo()`, funções auxiliares (log, validação, diretórios) |
| 01 | `01_download.R` | ✅ implementado e testado | Download de metadata/SOFT/FASTQ/processados do GEO para os 4 datasets |
| 02 | `02_metadata.R` | ✅ implementado e testado | Metadata padronizado; filtro para manter apenas WT/Controle e KO/KD/Deficiente (remove amostras tratadas com agente que não seja o próprio desenho WT-vs-KO). Rodado de ponta a ponta em 2026-07-13: 89 amostras mantidas, 84 descartadas — ver histórico. |
| 03 | `03_qc.R` | ✅ implementado e testado de ponta a ponta | `ShortRead::qa()`/`report()` + interpretação automática (% N, duplicação, adaptador) |
| 04 | `04_alignment.R` | ✅ implementado e testado de ponta a ponta | `Rbowtie2` (build+align+bam) + `Rsamtools` (sort/index), parsing do log de alinhamento |
| 05 | `05_filtering.R` | ✅ implementado e testado de ponta a ponta | `Rsamtools::filterBam()` (MAPQ nativo + dedup + blacklist ENCODE via `GenomicRanges`) em uma única passada |
| 06 | `06_chip_qc.R` | ✅ implementado (ChIPQCsample testado; batch não testado) | `ChIPQC`: fingerprint (SSD), fragment size, coverage; correlação/PCA em lote com 2+ amostras |
| 07 | `07_peakcalling.R` | ✅ implementado e testado (mecanismo WSL/MACS3 validado) | MACS3 via WSL (broad para XPC; narrow para ELK1/STAT1/STAT2); `--nolambda` automático quando `Input` do metadata está ausente (CLAUDE.md S9.1) |
| 08 | `08_diffbind.R` | ✅ implementado (mecânica testada; teste estatístico requer réplicas reais) | XPC: consenso (`GenomicRanges::reduce`) + contagem via `csaw`/`edgeR` (sem presença/ausência, CLAUDE.md S9.1). STAT2: `DiffBind` padrão (input pareado). ELK1/STAT1 fora (sem braço deficiente) |
| 09 | `09_annotation.R` | ✅ implementado e testado de ponta a ponta | `ChIPseeker::annotatePeak()` (TxDb hg38 + org.Hs.eg.db) |
| 10 | `10_enrichment.R` | ✅ implementado e testado de ponta a ponta | GO/KEGG (`clusterProfiler`), Reactome (`ReactomePA`), Hallmark (`msigdbr`+`enricher`) |
| 11 | `11_granges.R` | ✅ implementado e testado de ponta a ponta | Conversão/padronização para GRanges + validação de consistência |
| 12 | `12_genome_standardization.R` | ✅ implementado e testado com liftOver real | liftOver hg19→hg38 (ELK1) via chain oficial UCSC |
| 13 | `13_regulatory_universe.R` | ✅ implementado e testado de ponta a ponta | `GenomicRanges::reduce()` sobre todas as amostras hg38 |
| 14 | `14_overlap_matrix.R` | ✅ implementado e testado de ponta a ponta | Matriz binária de ocupação (região × proteína_genótipo) |
| 15 | `15_pairwise_overlap.R` | ✅ implementado e testado de ponta a ponta | `findOverlaps()`, índice de Jaccard, heatmap |
| 16 | `16_multiple_overlap.R` | ✅ implementado e testado de ponta a ponta | Interseção múltipla (`Reduce(intersect, ...)`) entre as 4 proteínas |
| 17 | `17_hotspots.R` | ✅ implementado e testado de ponta a ponta | Occupancy score, anotação genômica (ChIPseeker) e ranking |
| 18 | `18_bipartite_network.R` | ✅ implementado e testado de ponta a ponta | Rede Proteína→Região→Gene (igraph/tidygraph/ggraph), grau/betweenness/closeness/comunidades, export Cytoscape (GraphML) |
| 19 | `19_pathway_network.R` | ✅ implementado e testado de ponta a ponta (PPI real) | Rede PPI via `STRINGdb` (real, testada com GAPDH/TP53/MYC/EGFR/STAT1) + redes de similaridade GO/Reactome (`enrichplot::emapplot`) |
| 20 | `20_visualization.R` | ⬜ pendente | Padronização de figuras |
| 21 | `21_validation.R` | ⬜ pendente | Validação científica e técnica, com interrupção em falha crítica |
| 22 | `22_master_pipeline.R` | ⬜ pendente | Orquestração + relatório final |

## 4. Histórico de alterações

- **2026-07-13** — Estrutura de pastas (`Dados/`, `Figuras/`, `Arquivos/`, `Scripts/`,
  `Logs/`) e `.gitignore` criados (sessão anterior).
- **2026-07-13** — Relatório de compatibilidade dos 4 datasets executado (ver §9).
  Detectadas 3 incompatibilidades científicas antes de qualquer código de análise ser
  escrito:
  1. **STAT1/STAT2 eram de camundongo** (GSE217805) enquanto XPC/ELK1 são humanos →
     dataset substituído por combinação de duas subséries humanas do mesmo grupo
     (GSE222667 + GSE247724, células Huh7.5). Motivo: `liftOver` só reconcilia versões
     de montagem dentro da mesma espécie; misturar camundongo/humano no universo
     regulatório (Módulos 13-19) invalidaria as sobreposições. Referência: ver §9.
  2. **ELK1 (GSE91923) não tem braço KO/KD** (apenas 2 réplicas ENCODE) → Módulo 08
     (DiffBind) não será executado para ELK1; ELK1 permanece nos módulos de
     ocupação/anotação/enriquecimento/rede.
  3. **XPC (GSE214182) não tem input pareado** para o ChIP de XPC → Módulo 07 usará o
     input de H3K4me3 do mesmo genótipo/timepoint como substituto, com ressalva
     metodológica registrada.
  Decisões tomadas e confirmadas pelo usuário; ver §9 para o relatório completo.
- **2026-07-13** — README.md e CLAUDE.md criados.
- **2026-07-13** — `00_setup.R`, `01_download.R`, `02_metadata.R` implementados.
- **2026-07-13** — Descoberta uma instalação real de R (4.5.2 e 4.6.0, em
  `C:\Program Files\R\`) e da extensão R do VS Code (`reditorsupport.r`) na
  máquina do usuário — ao contrário do que constava em §5.2, o ambiente
  **tem** R e ~515 pacotes já instalados (incluindo `here`, `GEOquery`,
  `Biobase`, `dplyr`, `stringr`, `GenomicRanges`, `ChIPseeker`, `DiffBind`).
  Os módulos 00-02 foram executados de verdade (não só revisados) contra os
  4 datasets reais do GEO, e 4 bugs genuínos foram encontrados e corrigidos:
  1. **`here::here()` resolvia para a raiz do repositório Git** (`Bioinfo/`),
     não para `ProjetosCode/`, porque `.git` fica em `Bioinfo/` e
     `ProjetosCode/` não tinha marcador de projeto próprio — isso criava uma
     árvore `Dados/Figuras/Arquivos/Logs` duplicada e errada na raiz do repo.
     Corrigido criando `ProjetosCode/.here` (arquivo vazio, convenção do
     pacote `here`) para ancorar o root corretamente.
  2. **`GEOquery::getGEO()` na versão instalada retorna `SummarizedExperiment`
     por padrão**, não mais `ExpressionSet`, quebrando `Biobase::pData()`.
     Corrigido passando `returnType = "ExpressionSet"` explicitamente em
     `download_series_metadata()` (01_download.R).
  3. **`match_xpc_input()` colidia a coluna `Replicate`** no `left_join()`
     (virava `Replicate.x`/`Replicate.y`, quebrando o schema do metadata).
     Corrigido removendo `Replicate` do `select()` do lado do input (não é
     necessário para o match por Genotype+Treatment).
  4. **`extract_pdata()` só lia a primeira plataforma (GPL)** quando uma
     série GEO usa múltiplas plataformas — é o caso de GSE247724, que tem 5
     GPLs diferentes e cujas amostras "Huh7.5 STAT1 K.O." ficavam justamente
     fora da primeira. Corrigido combinando `pData()` de todas as
     plataformas via `bind_rows()`.
  Após as correções, os 4 parsers (`parse_gse214182`, `parse_gse222667`,
  `parse_gse247724`, `parse_gse91923`) foram validados linha a linha contra
  os dados reais baixados (61, 82, 89 e 2 amostras respectivamente) e
  batem exatamente com os títulos verificados manualmente na página do GEO.
  **Achado adicional:** as amostras **XPC-KO não têm nenhum input
  disponível**, nem mesmo o substituto de H3K4me3 (H3K4me3 só foi feito em
  WT e ASH1L-KO) — ver §10 (pendência: peak calling do XPC-KO no Módulo 07
  terá que rodar sem controle pareado).
- **2026-07-13** — Implementados `03_qc.R` (FastQC/MultiQC + interpretação
  automática), `04_alignment.R` (Bowtie2/samtools) e `05_filtering.R`
  (dedup via samtools, MAPQ, blacklist ENCODE via bedtools). Como as
  ferramentas externas (FastQC, Bowtie2, samtools, bedtools) ainda não estão
  no PATH, essas funções não podem rodar de ponta a ponta, mas toda a lógica
  de parsing pura (`parse_flagstat()`, `parse_fastqc_summary()`,
  `interpret_qc_flags()`) foi testada com dados sintéticos no R real e um
  bug real foi encontrado e corrigido: `interpret_qc_flags()` tentava gravar
  `qc_summary.csv` sem garantir que `Arquivos/qc/` existisse.
- **2026-07-13** — Decisão de arquitetura: substituir FastQC/MultiQC, samtools e
  bedtools por pacotes R nativos já instalados (`ShortRead`, `Rsamtools`+
  `GenomicAlignments`, `GenomicRanges`) nos Módulos 03-05, eliminando a dependência de
  PATH para essas três ferramentas sem mudar a metodologia. Bowtie2 (Módulo 04) passa a
  usar `Rbowtie2` (Bioconductor) — mesmo algoritmo, só empacotado como pacote R, evitando
  configuração manual de PATH (decisão do usuário; alternativa `Rsubread::align()`
  rejeitada por mudar o algoritmo de alinhamento). deepTools (Módulo 06) substituído por
  `ChIPQC`. MACS3 (Módulo 07) não tem equivalente R real — decisão do usuário: instalar
  Python + MACS3 de verdade e adicionar ao PATH, em vez de trocar por um peak caller R
  nativo (ex. `mosaics`), para preservar a metodologia pedida originalmente. Decisão
  detalhada sobre o peak calling do XPC-KO sem input registrada em §9.1.
- **2026-07-13** — Execução real da migração de ferramentas e reescrita completa dos
  Módulos 03-05, com testes de ponta a ponta (não só sintaxe) no R real:
  - Instalado `Rbowtie2` via BiocManager.
  - Instalado **Python 3.12** (winget) — necessário porque `Rbowtie2::bowtie2_build()`
    chama `python3` internamente; o instalador oficial do python.org só cria
    `python.exe`, então foi criada uma cópia `python3.exe` ao lado (bug conhecido do
    Windows). Sem isso, `bowtie2_build()`/`bowtie2_samtools()` falhavam silenciosamente
    (o stub "App Execution Alias" da Microsoft Store intercepta `python3`, imprime uma
    mensagem e sai sem gerar erro no R, então nem o fallback interno do Rbowtie2 para
    `python` disparava).
  - Instalado **Miniconda** — primeiro tentado no Windows nativo, mas o `bioconda` **não
    publica build de MACS3 para `win-64`**. Miniconda foi então instalado **dentro do
    WSL** (Ubuntu-24.04, já presente na máquina), onde o ambiente `chipseq` com MACS3
    3.0.4 foi criado com sucesso (sem precisar de `sudo` — Miniconda fica no `$HOME` do
    usuário do WSL). Funções `win_to_wsl_path()`, `run_wsl_command()`,
    `check_macs3_wsl()` e `run_macs3()` adicionadas a `00_setup.R` para o Módulo 07
    poder chamar o MACS3 de dentro do R.
  - Python, Miniconda e o ambiente `chipseq` foram adicionados ao **PATH do usuário**
    (persistente, via `[Environment]::SetEnvironmentVariable`), não só à sessão atual.
  - `03_qc.R`, `04_alignment.R` e `05_filtering.R` foram **reescritos do zero** para
    usar `ShortRead::qa()`, `Rbowtie2::bowtie2_build()`/`bowtie2_samtools()` e
    `Rsamtools::filterBam()` (MAPQ nativo + dedup + blacklist ENCODE numa única
    passada), respectivamente, no lugar dos wrappers de linha de comando anteriores.
  - Testados de ponta a ponta com dados reais (FASTQ de teste do próprio pacote
    `Rbowtie2`, genoma lambda): QA real gerado, alinhamento com 95% de taxa,
    ordenação/indexação, filtro combinado MAPQ+dedup+blacklist (2000→1838 reads),
    download real da blacklist ENCODE hg38 (636 regiões). Isso revelou e corrigiu mais
    4 bugs reais:
    1. `ShortRead::qa()[["adapterContamination"]]` é um `data.frame` (coluna
       `contamination`), não um vetor — `interpret_qa()` tentava `mean(as.numeric(.))`
       direto em um `data.frame` e falhava.
    2. `Rbowtie2::bowtie2_samtools()` tem um parâmetro `bamFile` posicional entre
       `seq2` e `...` — passar um argumento extra do bowtie2 (`--threads N`) sem nome
       no final da chamada era absorvido por `bamFile` em vez de cair em `...`.
       Corrigido nomeando `bamFile = NULL` explicitamente.
    3. A saída/log do alinhamento (taxa de alinhamento, total de reads) não é
       capturável via `capture.output()` — `bowtie2_samtools()` escreve direto no
       console do SO e devolve o log (via `.bowtie2.cerr.txt`) como **valor de
       retorno invisível**; `parse_alignment_log()` só funcionou depois de capturar
       o retorno da função diretamente.
    4. `interpret_qc_flags()`/diretório `Arquivos/qc/` (já corrigido antes) — mantido
       aqui por completude do padrão "sempre `ensure_dir()` antes de `write.csv()`".
- **2026-07-13** — Implementados `06_chip_qc.R` (ChIPQC: fingerprint/SSD, fragment
  size, coverage por amostra; correlação/PCA em lote) e `07_peakcalling.R` (MACS3 via
  WSL, broad/narrow por proteína, `--nolambda` automático quando a coluna `Input` do
  metadata está ausente) + `Scripts/config/peakcalling_config.R` com os parâmetros
  científicos por proteína. `ChIPQCsample()` e o mecanismo completo de chamada do
  MACS3 via WSL (conversão de caminho, quoting, invocação real) foram testados de
  ponta a ponta. Bug real encontrado e corrigido: caminhos do Windows com espaço
  (comuns em nomes de usuário) quebravam em múltiplos tokens no `bash -lc` do WSL
  porque `win_to_wsl_path()` não era envolvida em `shQuote()` antes de compor os
  argumentos do MACS3 em `call_peaks_macs3()`.
- **2026-07-13** — Implementado `08_diffbind.R` (differential binding — ver §3/§9.1) e,
  em seguida, `09_annotation.R` (ChIPseeker), `10_enrichment.R` (GO/KEGG/Reactome/
  Hallmark) e `11_granges.R` (padronização/validação de GRanges). `import_peaks()` foi
  promovida de `08_diffbind.R` para `00_setup.R` (função utilitária reaproveitada por
  08/09/11 e futuros módulos que leem picos). Testados de ponta a ponta com um
  narrowPeak sintético usando coordenadas **reais** de hg38 (promotor de GAPDH, entre
  outras) — Módulo 09 anotou corretamente o GAPDH (Entrez 2597), Módulo 10 rodou os 4
  enriquecimentos sem erro (sem termos significativos, esperado para 4 genes
  aleatórios), Módulo 11 converteu e validou o GRanges. Um bug real foi corrigido:
  `annotatePeak()` não cria uma coluna chamada `ENTREZID` — o Entrez ID do gene mais
  próximo fica na coluna `geneId` mesmo com `annoDb="org.Hs.eg.db"` (que só adiciona
  `SYMBOL`/`ENSEMBL`/`GENENAME`); `10_enrichment.R` esperava `ENTREZID` e falhava.
  Também trocado `msigdbr(category=)` (depreciado) por `msigdbr(collection=)`.
- **2026-07-13** — Implementados `12_genome_standardization.R` (liftOver hg19→hg38 via
  chain oficial da UCSC), `13_regulatory_universe.R` (`GenomicRanges::reduce()` sobre
  todas as amostras) e `14_overlap_matrix.R` (matriz binária região×proteína_genótipo).
  Testados de ponta a ponta com GRanges sintéticos hg38 e, no caso do Módulo 12, com
  **liftOver real** de uma coordenada hg19 conhecida (promotor do GAPDH,
  chr12:6.643.093-6.647.537 em hg19) — o resultado (chr12:6.533.927-6.538.371 em hg38)
  bate com a coordenada hg38 real do GAPDH, confirmando a chain e a conversão.
  Dois bugs reais encontrados e corrigidos:
  1. **`do.call(c, lapply(granges_list, granges))` não dispara o dispatch S4** de
     `c()` para `GRanges` de forma confiável e pode devolver uma `list` comum em vez
     de um `GRanges`, quebrando `reduce()` a seguir (`"unable to find an inherited
     method for function 'reduce' for signature 'x = list'"`). Corrigido em **três
     lugares** (`08_diffbind.R`, `13_regulatory_universe.R`, `14_overlap_matrix.R`)
     trocando para o idioma correto do Bioconductor: `unlist(GRangesList(lista), use.names=FALSE)`.
  2. **`S4Vectors::lengths` não existe** (não é exportado) — `lengths()` é o genérico
     base do R que `GenomicRanges`/`IRanges` estendem via método S4; corrigido em
     `12_genome_standardization.R` chamando `lengths()` sem qualificar o pacote.
- **2026-07-13** — Implementados `15_pairwise_overlap.R` (Jaccard + heatmap),
  `16_multiple_overlap.R` (interseção múltipla entre as 4 proteínas) e
  `17_hotspots.R` (occupancy score + anotação + ranking). Testados de ponta a ponta
  com um cenário sintético de 4 "proteínas" (3 sobrepostas no promotor do GAPDH,
  1 em outro cromossomo): Módulo 16 corretamente devolveu 0 regiões compartilhadas
  por todas as 4 (já que a 4ª não se sobrepõe às outras — comportamento esperado,
  não bug); Módulo 17 corretamente identificou a região do GAPDH como hotspot
  (occupancy_score=3) e a anotou. Bug real corrigido: o universo regulatório
  (Módulo 13) só tem o identificador da região em `names(gr)`, não numa coluna BED
  "name" — como `names()` se perde ao converter o resultado de `annotatePeak()`
  para `data.frame`, `annotate_hotspots()` foi corrigida para copiar
  `names(hotspot_gr)` para uma metadata column explícita (`region_id`) antes de
  anotar.
- **2026-07-13** — Implementados `18_bipartite_network.R` (rede Proteína→Região→Gene,
  métricas de grau/betweenness/closeness/comunidades via `igraph`, export Cytoscape em
  GraphML) e `19_pathway_network.R` (rede PPI via `STRINGdb` + redes de similaridade
  GO/Reactome via `enrichplot::emapplot`). Testados de ponta a ponta: Módulo 18 com o
  cenário sintético dos módulos 13-17 (rede de 7 nós/5 arestas, 2 comunidades
  corretamente detectadas); Módulo 19 com uma **rede PPI real** via STRINGdb
  (GAPDH/TP53/MYC/EGFR/STAT1 — 5/5 genes mapeados, 20 arestas). Nenhum bug novo
  encontrado nesta dupla de módulos.

## 5. Dependências

### 5.1 Pacotes do R (via BiocManager/CRAN)

| Pacote | Uso |
|---|---|
| `here` | caminhos relativos |
| `BiocManager` | instalação de pacotes Bioconductor |
| `GEOquery`, `Biobase` | download de metadata/SOFT/FASTQ processados do GEO; acesso a `pData()` |
| `GenomicRanges`, `GenomicFeatures`, `rtracklayer` | GRanges, liftOver, e (desde 2026-07-13) remoção de regiões da blacklist ENCODE via `findOverlaps()`/`subsetByOverlaps()` — substitui `bedtools` |
| `ShortRead` | QC nativo em R (`qa()`/`report()`) — substitui FastQC/MultiQC no Módulo 03 |
| `Rsamtools`, `GenomicAlignments` | manipulação de BAM (sort/index/flagstat-equivalente) — substitui `samtools` nos Módulos 04-05 |
| `Rbowtie2` | alinhamento (mesmo Bowtie2, empacotado como pacote R) — Módulo 04 |
| `ChIPQC` | fingerprint/PCA/coverage/correlação — substitui `deepTools` no Módulo 06 |
| `csaw` | differential binding por contagem de reads em regiões consenso — apoio ao Módulo 08 (ver §9.1) |
| `ChIPseeker`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db` | anotação |
| `DiffBind` | binding diferencial |
| `clusterProfiler`, `ReactomePA`, `msigdbr` | enriquecimento funcional |
| `STRINGdb`, `igraph`, `tidygraph`, `ggraph` | redes |
| `data.table`, `dplyr`, `tidyr`, `purrr`, `stringr`, `readr` | manipulação de dados |
| `ggplot2` | visualização |
| `yaml` | arquivos de configuração |
| `rmarkdown`, `knitr` | relatório final |

Lista completa fixada e versões efetivas serão gravadas em `Logs/` por `00_setup.R`
(`sessionInfo()`), a cada execução.

### 5.2 Softwares externos necessários (fora do R)

Revisado em 2026-07-13: sempre que existir um pacote R nativo equivalente já
instalado, ele substitui a ferramenta de linha de comando (elimina dependência de
PATH, sem mudar a metodologia). Ver histórico de decisões abaixo.

| Software | Uso | Observação |
|---|---|---|
| SRA Toolkit (`prefetch`, `fasterq-dump`) | download de FASTQ do SRA | Módulo 01 — sem equivalente R nativo, continua externo (já presente no PATH desta máquina em `Bioinfo/SRAtoolkit/bin`) |
| MACS3 (3.0.4) | peak calling | Módulo 07 — instalado dentro de um ambiente conda (`chipseq`) **no WSL** (Ubuntu-24.04), não no PATH nativo do Windows — ver detalhes abaixo |

**Substituídos por pacotes R nativos (não precisam mais estar no PATH):**

| Ferramenta externa (spec original) | Substituto R nativo | Módulo |
|---|---|---|
| FastQC + MultiQC | `ShortRead::qa()`/`report()` | 03 — testado de ponta a ponta |
| Bowtie2 | `Rbowtie2` (mesmo algoritmo, empacotado) | 04 — testado de ponta a ponta |
| samtools (view/sort/index/flagstat/fixmate/markdup) | `Rsamtools` + `GenomicAlignments` | 04-05 — testado de ponta a ponta |
| bedtools (`intersect -v` para blacklist) | `GenomicRanges::overlapsAny()` dentro de `Rsamtools::filterBam()` | 05 — testado de ponta a ponta |
| deepTools (fingerprint/PCA/coverage/correlação) | `ChIPQC` (já instalado) | 06 — a testar quando o módulo for escrito |

**R está instalado** (versões 4.5.2 e 4.6.0 em `C:\Program Files\R\`, ~515 pacotes já
presentes incluindo todo `REQUIRED_PACKAGES` de §5.1) e os módulos 00-05 já foram
executados de ponta a ponta com dados reais/sintéticos reais (ver histórico abaixo).

**MACS3 no Windows — por que via WSL:** o `bioconda` não publica build de MACS3 para
Windows (`win-64`); compilar via `pip` no Windows falha por faltar o Visual C++ Build
Tools (dependência transitiva `cykhash`). A solução foi instalar o Miniconda **dentro
do WSL** (Ubuntu-24.04, já presente na máquina) — lá o bioconda funciona normalmente,
sem precisar de `sudo` (Miniconda instala no `$HOME` do usuário). O ambiente `chipseq`
com MACS3 3.0.4 fica em `~/miniconda3/envs/chipseq/bin/macs3` dentro do WSL. O Módulo
07 (a implementar) deve chamar o MACS3 via as funções auxiliares já criadas em
`00_setup.R`: `win_to_wsl_path()` (converte `C:/...` para `/mnt/c/...`) e `run_macs3()`
(roda `wsl -d Ubuntu-24.04 -- bash -lc "~/miniconda3/envs/chipseq/bin/macs3 ..."`).
`check_macs3_wsl()` verifica a disponibilidade sem interromper o Módulo 00.

### 5.3 Versões

`sessionInfo()` já foi salvo em `Logs/sessioninfo_<timestamp>.txt` na execução real de
2026-07-13 (R 4.6.0). Novas execuções regravam um snapshot atualizado a cada rodada.

## 6. Fluxograma resumido do pipeline

```
00 Setup → 01 Download → 02 Metadata → 03 QC → 04 Alignment → 05 Filtering
   → 06 ChIP-QC → 07 Peak Calling → 08 DiffBind (XPC, STAT2 apenas)
   → 09 Annotation → 10 Enrichment → 11 GRanges → 12 Genome Standardization (hg38)
   → 13 Regulatory Universe → 14 Overlap Matrix → 15 Pairwise Overlap
   → 16 Multiple Overlap → 17 Hotspots → 18 Bipartite Network → 19 Pathway Network
   → 20 Visualization → 21 Validation (gate a cada etapa) → 22 Master Pipeline
```

## 7. Datasets utilizados (estado final, pós relatório de compatibilidade)

| Proteína(s) | GSE | Organismo | Linhagem | Desenho | Peak type | DiffBind (M08)? |
|---|---|---|---|---|---|---|
| XPC | GSE214182 | *H. sapiens* | U2OS | WT vs **XPC-KO** (ASH1L-KO disponível como comparação secundária/exploratória) | Broad | Sim (WT vs XPC-KO) |
| ELK1 | GSE91923 | *H. sapiens* | A549 (ENCODE, ENCSR623KNM) | apenas 2 réplicas, sem KO/KD | Narrow | Não (sem braço deficiente) |
| STAT1 | GSE222667 | *H. sapiens* | Huh7.5 | apenas WT (IFNα/IFNγ, múltiplos timepoints) | Narrow | Não (sem STAT1-ChIP em fundo deficiente — não faz sentido biológico) |
| STAT2 | GSE222667 (WT) + GSE247724 (Huh7.5 STAT1-KO) | *H. sapiens* | Huh7.5 | WT vs **STAT1-KO**, timepoints casados (0h/2h/24h/72h IFNα) | Narrow | Sim (WT vs STAT1-KO) |

**Dataset originalmente especificado e substituído:** STAT1/STAT2 = GSE217805
(*Mus musculus*, MEFs, WT vs Arid1a-KO) — descartado por incompatibilidade de espécie
com XPC/ELK1 (decisão do usuário em 2026-07-13, ver §4).

## 8. Checklist de validação científica

- [ ] Organismo/montagem/tipo de sequenciamento/réplicas confirmados por dataset (Módulo 02)
- [ ] Relatório de compatibilidade entre datasets (feito manualmente nesta sessão — ver §9; Módulo 02 deve automatizar)
- [ ] FRiP, NSC, RSC calculados e dentro de limites aceitáveis (Módulo 06/21)
- [ ] Taxa de alinhamento e nº de reads validados (Módulo 04/21)
- [ ] Ausência de input tratada e documentada (XPC — ver §9)
- [ ] Picos vazios / GRanges inconsistentes verificados (Módulo 21)
- [ ] Nenhuma amostra tratada (fora do desenho WT-vs-KO) incluída (Módulo 02)
- [ ] Coordenadas padronizadas em hg38 para todas as proteínas antes do universo regulatório (Módulo 12)
- [ ] Comparações de differential binding restritas a WT vs deficiente genuínos (XPC, STAT2) — ELK1/STAT1 excluídos do Módulo 08 por design

## 9. Relatório de compatibilidade (2026-07-13)

Resumo das inconsistências encontradas e decisões tomadas antes de iniciar o
desenvolvimento dos módulos analíticos (evita construir 22 módulos sobre premissa
inválida):

1. **Espécie.** GSE217805 (STAT1/STAT2 conforme especificação original) é
   *Mus musculus*; XPC (GSE214182) e ELK1 (GSE91923) são *Homo sapiens*. `liftOver`
   reconcilia apenas montagens da mesma espécie — não converte camundongo→humano de
   forma biologicamente válida para sobreposição de picos. **Decisão do usuário:**
   trocar o dataset. Substituído por GSE222667 (ChIP-seq WT, Huh7.5,
   pSTAT1/pSTAT2/IRF9/IRF1 + input, IFNα/IFNγ) combinado com GSE247724 (ChIP-seq
   Huh7.5 STAT1-KO, IRF9/STAT2/pSTAT2 + input, IFNα), ambos humanos, mesma linhagem
   parental, timepoints IFNα casados (0h/2h/24h/72h).
2. **ELK1 sem braço deficiente.** GSE91923 tem apenas 2 réplicas ENCODE de A549, sem
   KO/KD. **Decisão do usuário:** pular Módulo 08 (DiffBind) para ELK1; manter nos
   demais módulos (ocupação, anotação, enriquecimento, rede).
3. **XPC-KO com ChIP de XPC.** GSE214182 tem "XPC ChIP-seq em células XPC-knockout"
   (GSM6600722/723/730/731/738/739) — no desenho original do estudo-fonte, isso é um
   controle negativo de especificidade de anticorpo (não deveria haver proteína XPC
   para imunoprecipitar). **Decisão do usuário:** usar mesmo assim como o braço
   "deficiente" do Módulo 08 (WT vs XPC-KO), documentando que o resultado esperado é
   perda quase total de picos específicos — inclusive, isso serve como validação
   positiva de sensibilidade do próprio pipeline (Módulo 21). ASH1L-KO permanece
   disponível como comparação secundária/exploratória (testa o papel de ASH1L na
   ocupação de XPC), mas não é a comparação primária do Módulo 08.
4. **XPC sem input pareado.** Amostras de ChIP de XPC em GSE214182 não têm input
   próprio (só as amostras de H3K4me3 têm). **Decisão do usuário:** reaproveitar o
   input de H3K4me3 do mesmo genótipo/timepoint como controle substituto no MACS3,
   com ressalva metodológica registrada nos metadados e no relatório do Módulo 07.
5. **STAT1 sem braço deficiente com ChIP do próprio STAT1.** Por construção biológica
   (não há amostra "STAT1 ChIP-seq em célula STAT1-KO" em nenhuma das séries — faria
   tão pouco sentido quanto XPC-ChIP-em-XPC-KO, mas aqui nem foi tentado pelos autores),
   STAT1 é tratado como occupancy-only, mesma lógica aplicada ao ELK1. Decisão tomada
   por consistência metodológica com a decisão nº2 (não é necessário perguntar de novo
   ao usuário — mesmo princípio já aprovado).

### 9.1 Decisão de Peak Calling para XPC-KO sem Input (GSE214182)

**Problema:** as 6 amostras XPC-KO (GSM6600722/723/730/731/738/739) não têm nenhum
input/controle pareado — nem próprio, nem substituto de H3K4me3 (que só existe para
WT/ASH1L-KO, ver decisão nº4 acima).

**Decisão (2026-07-13):** peak calling do XPC-KO com **MACS3 sem controle
experimental** (`--nolambda`), usando apenas o modelo de Poisson baseado na
profundidade de sequenciamento da própria amostra.

**Alternativas rejeitadas e por quê:**
- *Input de WT como controle* — rejeitado: introduziria viés dependente de genótipo
  (o background de cromatina de WT não é o mesmo do XPC-KO).
- *H3K4me3 como substituto* — rejeitado: já é usado para WT/ASH1L-KO (decisão nº4);
  usá-lo também para XPC-KO misturaria o mesmo controle em genótipos com braços de
  comparação diferentes.
- *Excluir XPC-KO da análise* — rejeitado: preserva a comparação biológica
  WT-vs-deficiente que é o objetivo central do Módulo 08 para esta proteína.

**Estratégia por condição (a implementar no Módulo 07):**
```
if (Input_available == TRUE) {
  macs3 callpeak -t CHIP.bam -c INPUT.bam -f BAM --broad -q 0.01   # XPC-WT, ASH1L-KO
} else {
  macs3 callpeak -t CHIP.bam -f BAM --broad --nolambda -q 0.01     # XPC-KO
}
```
`Input_available` deve ser derivado automaticamente da coluna `Input` do metadata
(`chipseq_metadata.csv`, Módulo 02) — nunca assumir que todas as amostras têm input.

**Regras proibidas nesta decisão:** não usar input de WT para XPC-KO; não usar
H3K4me3 como substituto para XPC-KO; não excluir XPC-KO automaticamente da análise;
não comparar XPC-WT vs XPC-KO apenas por presença/ausência de pico.

**Mitigação a jusante (Módulo 08):** como o peak calling do XPC-KO não tem controle,
a análise diferencial **não pode depender só dos arquivos `.broadPeak`**. O Módulo 08
deve: (1) construir uma região consenso XPC-WT ∪ XPC-KO via `GenomicRanges::reduce()`;
(2) quantificar reads de ChIP em todas as regiões consenso; (3) comparar WT vs XPC-KO
por **contagem de reads na mesma região** (`DiffBind` e/ou `csaw`), não por
presença/ausência binária de pico.

**Validação adicional obrigatória para XPC-KO (Módulo 21):** FRiP score, número total
de picos, largura média dos picos, distribuição de fold enrichment, correlação entre
réplicas, PCA baseado na contagem de regiões consenso. Se uma réplica XPC-KO for
outlier extremo em relação às demais, sinalizar no relatório final.

**Impacto esperado:** FRiP/FDR menos confiáveis para XPC-KO que para XPC-WT (sem
modelagem de background local); espera-se também um número de picos "específicos"
muito menor no XPC-KO, o que serve como validação positiva de sensibilidade do
antocorpo/pipeline (não há proteína XPC para imunoprecipitar nesse genótipo).

## 10. Pendências futuras (TODO)

- [x] **XPC-KO sem input/controle disponível** — decisão registrada em §9.1 e
      **implementada** no Módulo 07 (`determine_macs3_args()` lê `has_input`
      diretamente, nunca assume). Falta a parte de differential binding por
      contagem de reads em regiões consenso no Módulo 08.
- [x] ~~Implementar Módulo 06 (ChIP-QC)~~ — `06_chip_qc.R` implementado com `ChIPQC`;
      `ChIPQCsample()` testado de ponta a ponta (fingerprint/SSD, fragment size,
      coverage, com `peaks=NULL` antes do peak calling). `run_chipqc_batch()`
      (correlação/PCA com 2+ amostras) não testado por falta de 2 amostras reais
      nesta sessão — validar quando houver BAMs reais de mais de uma amostra.
- [x] ~~Implementar Módulo 07 (Peak Calling)~~ — `07_peakcalling.R` implementado;
      mecanismo completo testado de ponta a ponta via WSL (conversão de caminho,
      quoting, chamada real do MACS3, leitura do arquivo de picos). Um bug real de
      quoting foi encontrado e corrigido: caminhos do Windows com espaço (ex. nomes
      de usuário) quebravam em tokens separados no `bash -lc` do WSL sem `shQuote()`
      ao redor de cada `win_to_wsl_path()`. O teste com dados sintéticos não gerou
      picos reais (MACS3 precisa de ≥100 pares de picos +/- para construir o modelo
      de fragmento, e o FASTQ de teste do `Rbowtie2` tem poucochíssimas reads) — isso
      é uma limitação do dado de teste, não do pipeline; validar com FASTQ real
      quando disponível.
- [x] ~~Implementar Módulo 08 (Differential Binding)~~ — `08_diffbind.R` implementado
      com os dois fluxos (consenso+contagem via `csaw`/`edgeR` para XPC; `DiffBind`
      padrão para STAT2). A mecânica (`build_consensus_regions()`,
      `count_reads_csaw()`) foi testada de ponta a ponta com dados reais/sintéticos;
      o teste estatístico final (`edgeR::glmQLFit`) falha com 1 amostra por grupo
      (sem graus de liberdade para estimar dispersão) — isso é esperado e só
      funciona com réplicas reais (o metadata real tem XPC WT n=10, XPC-KO n=6),
      não é um bug do código.
- [x] ~~Instalar `Rbowtie2`~~ — instalado via BiocManager em 2026-07-13, testado de
      ponta a ponta (build de índice + alinhamento real, 95% de taxa de alinhamento).
- [x] ~~Instalar Python + MACS3~~ — Python 3.12 instalado (winget); MACS3 3.0.4
      instalado dentro de um ambiente conda no **WSL** (bioconda não publica build
      para Windows) — ver §5.2 para detalhes e `win_to_wsl_path()`/`run_macs3()` em
      `00_setup.R`. Falta só usar essas funções ao escrever o Módulo 07.
- [ ] Confirmar montagem de genoma real (hg38 vs hg19) para XPC e STAT1/STAT2 a partir
      dos metadados SOFT completos no Módulo 01 (não estava explícita nas páginas GEO
      consultadas) — validar antes do Módulo 12.
- [ ] Buscar/baixar o input/controle ENCODE de ELK1 (ENCSR623KNM) via ENCODE portal, já
      que o input não aparece diretamente entre as 2 GSMs do GSE91923.
- [ ] Implementar Módulos 06–22 (ChIP-QC via ChIPQC, peak calling, DiffBind, anotação,
      enriquecimento, GRanges, padronização de genoma, universo regulatório, matriz de
      ocupação, overlaps, hotspots, redes, visualização, validação, master pipeline).
- [ ] Considerar IRF9 (ChIP'd em ambas as séries de STAT1/STAT2) como 5ª proteína
      opcional em versão futura do projeto (fora do escopo atual de 4 proteínas).
- [ ] Criar `Scripts/config/peakcalling_config.R` junto com o Módulo 07.
- [ ] Criar primeira tag (`v0.1.0` — estrutura inicial) quando o usuário autorizar.
