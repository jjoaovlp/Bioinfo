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
| 01 | `01_download.R` | ✅ implementado e testado (SRR real resolvido para as 89 amostras) | Download de metadata/SOFT/processados do GEO; resolução SRX→SRR e download de FASTQ via API pública da ENA (sem SRA Toolkit) |
| 02 | `02_metadata.R` | ✅ implementado e testado | Metadata padronizado; filtro para manter apenas WT/Controle e KO/KD/Deficiente (remove amostras tratadas com agente que não seja o próprio desenho WT-vs-KO). Rodado de ponta a ponta em 2026-07-13: 89 amostras mantidas, 84 descartadas — ver histórico. |
| 03 | `03_qc.R` | ✅ implementado e testado de ponta a ponta | `ShortRead::qa()`/`report()` + interpretação automática (% N, duplicação, adaptador) |
| 04 | `04_alignment.R` | ✅ implementado e corrigido para dados reais (ver §4, bug do .fastq.gz) | `Rbowtie2` (build+align+bam) + `Rsamtools` (sort/index), parsing do log de alinhamento; descompacta FASTQ .gz antes de alinhar |
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
| 20 | `20_visualization.R` | ✅ implementado e testado | Paleta/tema padrão, `save_figure()` multi-formato (PNG+PDF+SVG, ≥300dpi), manifesto de figuras |
| 21 | `21_validation.R` | ✅ implementado e testado com metadata real | Bateria de checks PASS/WARN/FAIL (duplicatas, espécie, genoma, réplicas, alinhamento, picos, GRanges, strand) — interrompe em qualquer FAIL |
| 22 | `22_master_pipeline.R` | ✅ implementado e testado de ponta a ponta (incl. download real) | Orquestra 01-21 com status OK/FALHOU/PULADO por módulo, tempo de execução, relatório HTML (+ PDF se houver pandoc), sessionInfo, árvore de arquivos |

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
- **2026-07-13** — **Pipeline completo (22/22 módulos) implementado.** Últimos três:
  `20_visualization.R` (paleta Okabe-Ito colorblind-safe, tema ggplot2 padrão,
  `save_figure()` multi-formato PNG+PDF+SVG ≥300dpi, manifesto de figuras),
  `21_validation.R` (bateria de checks PASS/WARN/FAIL — duplicata de amostra, espécie
  mista, genoma inválido, réplicas insuficientes, BAM sem índice, taxa de alinhamento,
  picos vazios, FRiP, GRanges inválido, orientação de strand — interrompe a execução em
  qualquer FAIL) e `22_master_pipeline.R` (orquestra 01-21, tolerante a falha por
  padrão: um módulo sem o `config` necessário é marcado "PULADO", um que falhar é
  marcado "FALHOU" e o pipeline segue; gera relatório HTML sempre, PDF se houver
  pandoc, `sessionInfo()`, árvore de arquivos e tabela resumo com tempo por módulo).
  **Testado de ponta a ponta de verdade**, incluindo download real do GEO: rodando
  `run_master_pipeline()` só com `metadata_df` no config, os Módulos 01 e 02 executaram
  de verdade (download real, ~138s e ~133s), os Módulos 03-16/18-19 foram corretamente
  pulados (faltava config), o Módulo 17 falhou com uma mensagem clara e específica
  (dependia do Módulo 14, que tinha sido pulado) sem derrubar o pipeline, e os Módulos
  20/21 rodaram OK — relatório HTML gerado corretamente, PDF pulado (sem pandoc nesta
  máquina). Nenhum bug encontrado nesta etapa final.
- **2026-07-13** — **Resolução real de SRR implementada** (`01_download.R`): GEO só
  expõe o accession SRX no campo `relation` da pData, não o SRR de fato. Adicionadas
  `extract_srx_from_pdata()` (regex sobre `relation*`), `get_ena_fastq_info()` (consulta
  a API pública da ENA — European Nucleotide Archive — por `run_accession`,
  `fastq_ftp`, `fastq_bytes`, `library_layout`) e `resolve_srr_table()` (junta tudo por
  GSM). `download_fastq_ena()` baixa o FASTQ diretamente do espelho HTTPS da ENA via
  `download.file()` — **sem precisar de SRA Toolkit/prefetch/fasterq-dump em lugar
  nenhum** (nem Windows nem WSL). `build_master_metadata()` (Módulo 02) agora popula a
  coluna `SRR` de verdade (antes sempre `NA`) e adiciona `fastq_ftp`/`fastq_bytes` ao
  metadata. **Testado de ponta a ponta com as 89 amostras reais**: 89/89 SRR resolvidos
  com sucesso; tamanho total do FASTQ calculado a partir de `fastq_bytes`: **~95,3 GB**
  (XPC 64,8GB, STAT2 13,5GB, STAT1 8,8GB, Input 6,6GB, ELK1 1,5GB). Decisão do usuário:
  baixar as 89 amostras completas (95GB) em vez de um piloto menor.
  (Nota: durante a investigação desta solução, um ambiente `sra-tools` chegou a ser
  instalado no ambiente conda `chipseq` do WSL como alternativa — acabou não sendo
  necessário, mas fica disponível caso o download via ENA falhe para alguma amostra.)
- **2026-07-13** — **Bug crítico encontrado e corrigido em `download_fastq_ena()`**
  durante a primeira tentativa real de baixar as 89 amostras: o `timeout` padrão do R
  (`getOption("timeout")` = 60s) cortava o download no meio de arquivos grandes
  (ChIP-seq facilmente passa de 1-4GB por amostra), e `download.file()` **não lança um
  erro capturável** nesse caso — o arquivo fica silenciosamente truncado e o
  `tryCatch` não detectava nada de errado. As primeiras 21 amostras baixadas nessa
  tentativa (2,9GB) estavam **todas corrompidas** (`gzip -t` falhou com "unexpected
  end of file"); o processo foi interrompido manualmente assim que isso foi
  descoberto (via `Stop-Process` pelo PID específico, não por nome — matar por nome
  poderia afetar outras sessões de R do usuário) e os arquivos truncados apagados.
  Correção em `01_download.R`: `download_fastq_ena()` agora (1) eleva
  `options(timeout=)` para 3600s durante o download (restaurado com `on.exit()`), (2)
  compara o tamanho baixado com `fastq_bytes` esperado (vindo da ENA via
  `resolve_srr_table()`), e (3) verifica a integridade do stream gzip por completo
  (`verify_gzip_integrity()`) — um arquivo existente mas incompleto/corrompido é
  apagado e baixado de novo, nunca aceito silenciosamente. Testado de ponta a ponta
  com um arquivo real (SRR23051534, 346.653.269 bytes) — tamanho bateu exatamente e a
  integridade foi confirmada. **Velocidade real medida: ~2,7 MB/s**, o que projeta
  **~10 horas** para as 89 amostras (95GB) — bem mais que a estimativa inicial de
  "algumas horas". Decisão do usuário: prosseguir mesmo assim com as 89 amostras
  completas, rodando em processo desacoplado (`Start-Process` no PowerShell,
  independente da sessão de ferramentas) monitorado periodicamente.
- **2026-07-14** — Download das 89 amostras concluído (89GB, todas íntegras). Genoma
  hg38 (UCSC) baixado e índice Bowtie2 construído com sucesso. Módulo 03 (QC) rodado
  de ponta a ponta nas 89 amostras reais. Ao iniciar o Módulo 04 (alinhamento) real,
  dois bugs adicionais apareceram:
  1. O processo desacoplado do índice hg38 falhou na primeira tentativa com o mesmo
     bug do `python3`/PATH já conhecido (§ acima) — o processo `Start-Process` não
     herdou a atualização de PATH feita numa sessão PowerShell anterior. Corrigido
     definindo `$env:PATH` explicitamente na mesma chamada que lança o processo
     desacoplado, e relançado com sucesso.
  2. **Bug crítico real, encontrado só com dados reais**: o binário do Bowtie2 usado
     por esta instalação do `Rbowtie2` no Windows **não lê `.fastq.gz` corretamente**
     — processa só uma fração minúscula do arquivo e para silenciosamente, sem
     lançar erro. Só apareceu porque os testes anteriores sempre usaram FASTQ de
     exemplo do próprio pacote `Rbowtie2` (não comprimidos); esta foi a primeira vez
     alinhando um `.fastq.gz` real. Sintoma: a amostra GSM6600715 (111.135.914 reads
     reais, confirmado via `zcat | wc -l`) "alinhou" em 15 segundos reportando **22
     reads processados** — impossível para um arquivo real desse tamanho. Teste
     isolado confirmou: uma fatia de 10.000 reads comprimida produzia só 10 reads
     alinhados; a mesma fatia descomprimida produzia os 10.000 corretamente
     (94,49% de taxa de alinhamento, idêntico nos dois casos após a correção).
     **69 das 89 amostras já tinham sido "alinhadas" incorretamente** antes de o
     problema ser notado (o processo rodava ~15s por amostra em vez de minutos) —
     todos os 69 BAMs errados foram apagados. Corrigido em `04_alignment.R`:
     `decompress_fastq_if_needed()` descompacta qualquer FASTQ `.gz` para um arquivo
     temporário antes de `bowtie2_samtools()`, removido logo depois (via `on.exit()`).
     Validado com a mesma fatia de teste: 10.000/10.000 reads, 94,49% de alinhamento,
     idêntico ao arquivo descomprimido de referência. Módulo 04 relançado do zero
     para as 89 amostras com a correção.
- **2026-07-15** — **Decisão de priorização do usuário**: XPC é a proteína mais
  importante do projeto. Dentro do XPC, o **ASH1L-KO fica de fora por enquanto**
  (é a comparação secundária/exploratória — CLAUDE.md §9 ponto 3 —, não a principal
  do Módulo 08, que é WT vs XPC-KO). Isso reduz o XPC de 64,8GB (25 amostras) para
  **39,5GB (16 amostras, WT+XPC-KO)** — uma redução de ~39% dentro do XPC, e de ~27%
  no volume total do projeto (95,3GB → ~70GB nesta rodada). O ASH1L-KO (9 amostras,
  25,4GB) será processado depois, numa rodada separada, quando "sobrar tempo".
  Na hora da decisão, o Módulo 04 já tinha processado 9 amostras (GSM6600715-723):
  4 WT + 3 ASH1L-KO + 2 XPC-KO — os 3 ASH1L-KO já feitos não foram descartados
  (ficam prontos para quando essa parte for retomada), só não se processou mais
  nenhum ASH1L-KO daqui em diante. A fila foi reconstruída a partir do metadata,
  removendo as 9 amostras já concluídas e todas as `Genotype == "ASH1L-KO"`
  restantes, e reordenada por prioridade: **XPC (WT+XPC-KO) → Input → STAT2 → STAT1
  → ELK1** (74 amostras nesta fila). O processo anterior foi parado de forma limpa
  (`Stop-Process` pelo PID específico, nunca por nome) depois que a amostra em
  andamento terminou, para não desperdiçar trabalho já feito; um processo órfão do
  `bowtie2-align-s` (cujo processo R pai já tinha sido encerrado) e ~46GB de
  arquivos temporários da amostra que estava em andamento no momento da parada
  (que não foram limpos por `Stop-Process -Force` pular o `on.exit()` do R) também
  foram removidos.
- **2026-07-15** — **Interrupção não planejada**: o computador reiniciou sozinho
  (`LastBootUpTime` às 06:57, sem ação do usuário — provavelmente atualização
  automática do Windows) enquanto o Módulo 04 alinhava `GSM6600725`, matando o
  processo silenciosamente (sem erro capturável, já que o SO inteiro reiniciou).
  Retomado sem perda: os 10 BAMs já concluídos permaneceram intactos, nenhum
  arquivo parcial/corrompido foi deixado para trás, e o script (idempotente por
  design — sempre filtra `Dados/BAM/*.sorted.bam` já existentes) foi relançado e
  retomou exatamente de `GSM6600725` em diante.
- **2026-07-15** — **Auto-orquestração da análise do XPC**: decisão do usuário de
  que, assim que XPC (WT+XPC-KO) e seus inputs terminarem de alinhar, a análise
  (filtragem→ChIP-QC→peak calling→differential binding) deve começar
  automaticamente, sem esperar STAT1/STAT2/ELK1. **Achado durante o planejamento**:
  as 3 amostras de H3K4me3 que seriam o input substituto do XPC-WT
  (GSM6600680/686/692, uma por timepoint — decisão registrada em §9 ponto 4) tinham
  sido excluídas do metadata final no Módulo 02 (ficaram só em
  `chipseq_metadata_filtered_out.csv`, fora do escopo das 4 proteínas do projeto) e
  por isso nunca tinham sido baixadas. O SRR/`fastq_ftp` já estava resolvido para
  elas (6,1GB no total). Decisão do usuário: baixar e alinhar essas 3 amostras
  também antes de considerar "XPC + inputs" completo (preserva a metodologia
  original de usar H3K4me3 como substituto de input para XPC-WT).

  Implementado um script orquestrador único (`run_xpc_priority_pipeline.R`) que:
  (1) termina o alinhamento do XPC (WT+XPC-KO) restante + as 3 H3K4me3; (2) roda
  automaticamente os Módulos 05 (filtragem), 06 (ChIP-QC), 07 (peak calling — WT
  usa o H3K4me3 correspondente como input via a coluna `Input` do metadata,
  XPC-KO usa `--nolambda`, sem intervenção manual) e 08 (differential binding,
  chamando diretamente `build_consensus_regions()`/`run_diffbind_consensus_count()`/
  `save_diffbind_results()` de `08_diffbind.R` — não o `run_module_08()` completo,
  que exige também dados de STAT2 ainda não alinhados); (3) só então continua para
  o restante da fila (Input→STAT2→STAT1→ELK1).

  Durante o download dos 3 H3K4me3, a rede caiu repetidamente no meio de
  transferências grandes (`download.file()` do R e `curl` falharam várias vezes
  com "Failure when receiving data from the peer"/erro de TLS do schannel,
  mesmo com `--retry`). Resolvido usando o **BITS** (Background Intelligent
  Transfer Service, serviço nativo do Windows) via `Start-BitsTransfer`/
  `Resume-BitsTransfer` — mais resiliente a quedas de conexão que
  `download.file()`/`curl` para arquivos grandes nesta rede. Os 3 arquivos foram
  baixados e verificados (tamanho exato + integridade gzip) com sucesso.
- **2026-07-16** — **Bug real encontrado por inspeção visual das figuras de QC**:
  `run_module_06()` (`06_chip_qc.R`) era chamado sem `blacklist_gr` em
  `run_xpc_priority_pipeline.R`, então `ChIPQCsample(..., blacklist = NULL)` nunca
  calculava o ponto "Post_Blacklist" — o gráfico `plotSSD()` (fingerprint) mostrava
  a legenda com duas categorias (Pre/Post_Blacklist) mas só o ponto "Pre_Blacklist"
  aparecia, em todas as amostras conferidas (GSM6600715, GSM6600722). O gráfico de
  fragment size/cross-coverage (`plotCC()`) estava correto (pico em shift ~150-160bp,
  esperado para ChIP-seq). Corrigido em `06_chip_qc.R`: `run_module_06()` agora
  carrega a blacklist ENCODE automaticamente via `download_encode_blacklist()`
  (reaproveitada de `05_filtering.R`, agora `source()`ada como dependência) quando
  `blacklist_gr` não é informado. **Decisão do usuário**: aplicar a correção só daqui
  para frente (STAT2/STAT1/ELK1) — as 9 figuras de fingerprint de XPC já geradas
  antes da correção não serão reprocessadas (o dado em si não muda, já que o
  Módulo 05 já filtra os BAMs contra a blacklist antes do ChIP-QC rodar; o efeito é
  só cosmético/informativo no gráfico).
- **2026-07-17** — **Crash real do pipeline**: `run_xpc_priority_pipeline.R` (PID
  2688) completou o ChIP-QC individual das 19 amostras de XPC+H3K4me3 com sucesso
  (`Arquivos/chip_qc/chipqc_metrics.csv` salvo), mas a etapa em lote
  (`run_chipqc_batch()` → `ChIPQC()`, correlação+PCA entre as 19 amostras) falhou e
  **derrubou o processo R inteiro**: `ChIPQC()` tentou paralelizar via BiocParallel
  (10 workers SNOW, default no Windows), e os processos-filho não carregam
  automaticamente o namespace de `GenomeInfoDb` — `seqlevels<-` (usada
  internamente por `ChIPQC()`) não foi encontrada em nenhum worker, gerando
  "BiocParallel errors... Execução interrompida". Nenhum dado científico foi
  perdido (BAMs filtrados e métricas individuais das 19 amostras intactos), só a
  etapa de comparação em lote não chegou a rodar, e por isso o script nunca
  chegou no peak calling (Módulo 07).

  **Corrigido** em `run_chipqc_batch()` (`06_chip_qc.R`): força
  `BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)` antes de
  chamar `ChIPQC()`, evitando os workers SNOW (mais lento, mas roda no processo
  principal que já tem o namespace carregado).

  **Decisão do usuário**: como re-rodar `ChIPQC()` em lote recomputaria
  `ChIPQCsample()` para as 19 amostras do zero (~24h, mesmo tempo do ChIP-QC
  individual já feito), pular a correlação/PCA em lote por agora e seguir direto
  para peak calling (Módulo 07) + differential binding (Módulo 08) — nenhum dos
  dois depende do resultado do `ChIPQC()` em lote (usam só os BAMs filtrados e as
  métricas individuais). O heatmap de correlação e o PCA entre amostras podem ser
  gerados depois, separadamente, agora que o bug está corrigido.

  Retomado com um novo script (`run_xpc_resume_after_chipqc.R`, novo PID 25560)
  que pula direto para peak calling+diffbind e depois continua para o restante da
  fila (Input→STAT2→STAT1→ELK1), igual ao script original.

  **Bug secundário encontrado no processo de retomada**: a máquina tem duas
  instalações de R (`R-4.5.2` e `R-4.6.0` em `C:\Program Files\R\`) — só a
  `R-4.6.0` tem a biblioteca de ~515 pacotes do usuário
  (`AppData\Local\R\win-library\4.6`); `R-4.5.2` só tem os ~30 pacotes base.
  Resolver a versão do R automaticamente (ex.: `Get-ChildItem ... | Select-Object
  -First 1`) pode pegar a errada por ordem alfabética. Usar sempre o caminho
  explícito `C:\Program Files\R\R-4.6.0\bin\Rscript.exe` ao lançar scripts deste
  projeto.
- **2026-07-17** — **Segundo crash real (Módulo 07, monitoramento autônomo)**: o
  MACS3 não conseguiu construir o modelo de picos pareados para `GSM6600718`
  ("MACS3 needs at least 100 paired peaks... but can only find 82! Process for
  pairing-model is terminated!") — amostra de sinal mais fraco. Derrubou o
  processo inteiro de novo (o `run_macs3()` propaga o código de saída não-zero
  do WSL como `stop()`). **Corrigido** em `call_peaks_macs3()`
  (`07_peakcalling.R`): agora tenta de novo automaticamente com `--nomodel
  --extsize <fragment_length> --shift 0` quando isso acontece, usando o
  comprimento de fragmento já estimado pelo Módulo 06 (cross-correlation, lido
  de `Arquivos/chip_qc/chipqc_metrics.csv` via nova função
  `lookup_chipqc_fragment_length()`) — é a alternativa que o próprio MACS3
  sugere no aviso de erro. Retomado com `run_xpc_resume2_peakcalling.R` (PID
  11440), reaproveitando os peaks já prontos (GSM6600715/716/717) e testado com
  sucesso: as 16 amostras de XPC terminaram peak calling (fallback funcionou
  para GSM6600718).

  **Terceiro crash real (Módulo 08, monitoramento autônomo)**: com peak calling
  100% completo, `save_diffbind_results()` falhou ao escrever o CSV
  (`não é possível abrir a conexão`) porque `Arquivos/differential/` nunca tinha
  sido criada — `ensure_dir(DIFFBIND_DIR)` só existia dentro de
  `run_module_08()` (o pipeline completo, que exige também STAT2), não na
  própria `save_diffbind_results()`, que é chamada diretamente pelos scripts de
  retomada do XPC (decisão de arquitetura registrada acima). O consenso+contagem
  de reads (`run_diffbind_consensus_count()`, 353522 regiões, 16 amostras) já
  tinha rodado antes do crash mas não foi salvo, precisando ser recalculado.
  **Corrigido**: `save_diffbind_results()` agora chama `ensure_dir(DIFFBIND_DIR)`
  no próprio corpo, tornando-a independente de quem a chama. Retomado com
  `run_xpc_resume3_diffbind.R` (PID 19960), reaproveitando os 16 peaks prontos e
  refazendo só o diffbind.

  Padrão notado nos 3 crashes: cada um revelou uma lacuna genuína de robustez
  (bug de paralelismo, limite do MACS3 em amostras de sinal fraco, diretório de
  saída assumido mas nunca criado) — nenhum foi um falso alarme de monitoramento.
  Todos foram corrigidos no código-fonte (não só contornados no script de
  retomada), então protegem também as próximas proteínas (STAT2/STAT1/ELK1).

  **Marco: análise de XPC completa** (2026-07-17, 06:06:09) — 16 amostras de
  XPC (10 WT + 6 XPC-KO) com peaks MACS3 (broad), differential binding via
  DiffBind (consenso de 353.522 regiões) salvo em
  `Arquivos/differential/XPC_*.csv`: 3.283 regiões ganhas, 2.499 perdidas,
  347.740 estáveis (WT vs XPC-KO). O pipeline seguiu sozinho, sem intervenção,
  para o restante da fila (64 amostras — Input/STAT2/STAT1/ELK1, Módulo 04).

  Aproveitando que o alinhamento do restante estava rodando (não compete por
  recursos com uma etapa nativa em R que só depende de BAMs já prontos),
  relançada separadamente a etapa de ChIPQC em lote (correlação+PCA das 19
  amostras de XPC+H3K4me3) que tinha crashado antes — agora com o fix de
  `BiocParallel::SerialParam()` já aplicado (`run_xpc_chipqc_batch.R`, PID
  12832).
- **2026-07-17** — **Incidente de processo órfão do BiocParallel**: ao parar o
  primeiro script de reprocessamento do fingerprint do XPC (PID 24120, para
  adicionar a regeneração do `chipqc_metrics.csv`), o worker-filho que o
  `ChIPQCsample()` havia spawnado via BiocParallel (PID 24748) **ficou órfão** e
  continuou rodando, queimando um core inteiro e a caminho de colidir na escrita
  dos mesmos arquivos de saída do run relançado. Detectado porque a CPU do
  launcher legítimo estava ~0 enquanto deveria haver trabalho pesado; mapeada a
  árvore de processos (`Get-CimInstance Win32_Process`) e o órfão foi morto por
  PID específico antes de sobrescrever qualquer arquivo. **Lição**: ao parar um
  processo que usa ChIPQC/BiocParallel, matar a árvore inteira (launcher +
  worker-filho), nunca só o launcher. O launcher fica com CPU ~0; é o worker que
  faz o trabalho pesado — para monitorar saúde, olhar a CPU do worker.
- **2026-07-17** — **Bug real: `save_batch_qc_plots()` gerava imagens em branco**.
  O ChIPQC em lote terminou com sucesso (SerialParam funcionou), mas
  `correlation_heatmap.png` e `pca.png` sairam **byte-a-byte idênticos e
  totalmente brancos**. Causa: `plotCorHeatmap()`/`plotPrincomp()` do ChIPQC
  desenham em BASE GRAPHICS (não devolvem objeto ggplot), mas o código usava
  `ggplot2::ggsave()`, que captura o dispositivo ggplot vazio → imagem branca,
  a mesma nas duas chamadas. **Corrigido** em `save_batch_qc_plots()`
  (`06_chip_qc.R`): abre dispositivo `grDevices::png()`, chama a função de plot,
  fecha com `dev.off()` (via helper `save_base_plot()` com `on.exit()` e
  `tryCatch()`). Também: `run_chipqc_batch()` agora salva o objeto
  `ChIPQCexperiment` como `Arquivos/chip_qc/chipqc_experiment.rds`, para permitir
  regenerar os gráficos sem recomputar (~10h) tudo de novo. **Pendência**: os
  plots de correlação/PCA do XPC especificamente ficaram em branco (o run que os
  gerou já tinha terminado com a versão bugada e não salvou o RDS) — serão
  regenerados depois, quando não competir com o reprocessamento do fingerprint
  em andamento. A correção já protege os plots em lote das próximas proteínas
  (STAT2/STAT1/ELK1), que ainda vão rodar o Módulo 06.
- **2026-07-18** — **RE-ESCOPO do downstream** (aprovado pelo usuário; plano em
  `.claude/plans/cozy-percolating-wind.md`). XPC recebe análise individual
  completa (peak calling ✓, diffbind ✓, **anotação ✓** — ocupação WT 1577
  regiões + regiões diferenciais ganhas/perdidas/estáveis separadas em
  `Arquivos/annotation/XPC_*`; achado: perdidas 45% em promotor). As demais
  (STAT1, STAT2, ELK1, **IRF9** como 5ª proteína) entram **só na metanálise
  (Módulos 13–19), apenas WT**, nos timepoints **untreated + IFNα 2h** (2
  réplicas), para achar o que há em comum — âncora no XPC (genes em comum
  XPC∩combinações + Venn/UpSet, camada adicional). **Sem** diffbind das outras,
  **sem** anotação individual das outras (a associação a genes das regiões
  comuns vem dos Módulos 17/18), **sem** processar KO downstream. Filtragem
  (Módulo 05) mantida completa (MAPQ+dedup+blacklist; a sobreposição Pre/Post nas
  figuras é só o sub-passo da blacklist, marginal). ELK1 usa o input ENCODE
  (controle ENCSR949BZP do experimento ELK1 ENCSR623KNM: ENCFF002ECM/ECL, A549,
  input library, single-end 36nt). IRF9-WT (GSM6928595/596 UN, GSM6928599/600
  IFNα2h) baixado da ENA (estava em `chipseq_metadata_filtered_out.csv`).
- **2026-07-18** — **BLOQUEIO: Smart App Control (SAC) passou a bloquear o
  bowtie2 nativo.** Ao alinhar as 6 amostras novas (4 IRF9 + 2 input ELK1) o
  `bowtie2-align-s.exe` (não-assinado) falhou com "Permission denied" (status
  126). Causa: o **Smart App Control do Windows 11 virou ENFORCEMENT**
  (`HKLM:\SYSTEM\CurrentControlSet\Control\CI\Policy\VerifiedAndReputablePolicyState=1`)
  e passou a bloquear executáveis não-assinados — funcionou nas 86 amostras
  anteriores porque o SAC estava em avaliação/desligado, e flipou para
  enforcement após um reboot/atualização. **Não afeta** ChIPQC/DiffBind/
  ChIPseeker (pacotes R, rodam dentro do Rscript.exe assinado) nem o MACS3 (roda
  no WSL/Linux). **Solução escolhida pelo usuário: alinhar via WSL** (não desligar
  o SAC, que é irreversível no Win11). Instalados `bowtie2 2.5.5` + `samtools
  1.24` no env conda `chipseq` do WSL (mesmo usado p/ MACS3); o índice hg38
  (`.bt2`, independente de plataforma) e os FASTQ são lidos via `/mnt/c`; o
  bowtie2 do Linux lê `.fastq.gz` direto (sem o bug do bowtie2 nativo). Script
  `scratchpad/align_new_wsl.sh` (bowtie2 -p4 | samtools sort). **Qualquer
  alinhamento futuro deve usar o WSL enquanto o SAC estiver em enforcement.**
- **2026-07-18** — **Títulos das figuras de ChIP-QC agora incluem o nome
  legível da amostra.** `06_chip_qc.R`: novas `sample_label()`/`sample_title()`
  resolvem "Proteína Genótipo" (ex. "XPC WT") a partir de `Protein`/`Genotype`
  em `chipseq_metadata.csv` (com fallback para `chipseq_metadata_filtered_out.csv`
  — cobre IRF9/H3K4me3) por `sample_id` (GSM); amostras sintéticas sem GSM
  (inputs ENCODE combinados) usam `MANUAL_SAMPLE_LABELS`. `save_sample_qc_plots()`
  agora titula os gráficos de fragment size/cross-coverage e fingerprint/SSD
  como `"GSM... (XPC WT)"` em vez de só `"GSM..."`. Vale só para figuras
  regeneradas a partir de agora (não retroativo às já geradas).
- **2026-07-18** — **Cache do objeto ChIPQCsample (RDS) + CSV de métricas por
  amostra.** `ChIPQCsample()` é a etapa cara do Módulo 06 (minutos a horas por
  amostra) e o objeto `qc` não era salvo — trocar o título/estilo de uma figura
  já processada exigia recomputar tudo. Agora `run_chipqc_sample()`
  (`06_chip_qc.R`) salva automaticamente: (1) o objeto em
  `Arquivos/chip_qc/qc_objects/<sample_id>_chipqc.rds` (permite re-plotar/
  re-extrair sem recomputar), e (2) as métricas (frip/ssd/fragment_length) em
  CSV próprio por amostra em `Arquivos/chip_qc/metrics/<sample_id>_metrics.csv`
  — elimina a dependência do `chipqc_metrics.csv` em lote, que é sobrescrito a
  cada grupo processado e já precisou de backup manual entre rodadas. Novo
  helper `extract_sample_metrics()` centraliza a extração (reusado por
  `run_module_06()`).
- **2026-07-18** — **Metanálise WT concluída (XPC como referência).** Peak
  calling WT (narrow, MACS3/WSL) das 14 amostras: STAT1/STAT2/IRF9 (UN+IFNα2h) +
  ELK1 (com input ENCODE). **Achado biológico de validação**: STAT1/STAT2/IRF9
  têm pouquíssimos picos em untreated (2-67) e dezenas de milhares em IFNα2h
  (40-62k) — salto que confirma a ativação do complexo ISGF3 por interferon
  (timepoints UN+IFNα2h acertados). Rodados os Módulos 13-18 (universo
  regulatório, matriz de ocupação, Jaccard, interseção das 5 [0 regiões ao nível
  de coordenada exata — esperado dada a diversidade], hotspots anotados, rede
  bipartida Proteína→Região→Gene). **Camada XPC-âncora (nível de genes)**:
  XPC∩combinações, em 2 versões — (a) *nearest-gene*: STAT1/STAT2 saturam
  (~19-20k genes, quase o genoma, por terem 40-62k picos); (b) *promotor ±3kb*
  (recomendada): reduz a saturação mas STAT1/STAT2 seguem amplos (~14k).
  **Resultado interpretável** (limitado pelos conjuntos específicos): **XPC ∩
  eixo interferon (STAT1∩STAT2∩IRF9) = 8 genes** (promotor), incluindo **RIGI
  (DDX58/RIG-I)** e **IFI44L** — ISGs canônicos que validam o achado. XPC∩ELK1=0.
  Saídas em `Arquivos/metanalise/` (CSVs ambas versões), `Figuras/metanalise/`
  (Venn/UpSet ambas versões), `Arquivos/{overlap,hotspots,network}/`,
  `granges_hg38/`. **Ressalva**: amplitude de STAT1/STAT2 é parte biologia real
  (STATs induzidos ligam massivamente às 2h), parte possível sobre-sensibilidade
  do peak calling — interseções confiáveis são as com IRF9/ELK1/XPC.
- **2026-07-18/19** — **Extras da metanálise + enriquecimento + correções de
  visualização** (pedidos do usuário). **Enriquecimento funcional** dos 8
  genes-núcleo (`Scripts/10_enrichment.R`, GO/KEGG/Reactome/Hallmark): só
  Reactome achou termo significativo — **"Modulation of host responses by
  IFN-stimulated genes"** (R-HSA-9909505, p.adjust=0.022, genes RIGI/IFI44L,
  fold enrichment 332×) — validação estatística independente de que o núcleo
  XPC∩interferon está na via biológica esperada. **Bug real corrigido**: o
  `msigdbr` renomeou a coluna do Entrez ID (`entrez_gene`→`ncbi_gene` a partir
  da 7.5.1), quebrando `get_hallmark_term2gene()`; agora resolve o nome da
  coluna dinamicamente. **Tabela de interseções exclusivas** (estrutura do
  UpSet, cada gene em uma só combinação) salva em
  `metanalise_interseccoes_completas_promoter.csv`. **Rede focada legível**
  (`rede_focada_xpc.png`, só os 8 genes-núcleo + XPC∩IRF9/ELK1) substituindo a
  necessidade de usar a rede completa do Módulo 18 (que com todos os hotspots
  vira um "hairball" ilegível por volume — não é bug, é o preço de plotar
  milhares de nós). **Bug real corrigido**: `theme_void()` deixa o fundo do
  painel transparente (`element_blank()`), que aparece preto em visualizadores
  sem fundo próprio, tornando rótulos de texto preto invisíveis — corrigido com
  fundo branco explícito no tema + `ggsave(..., bg="white")`, aplicado também
  em `18_bipartite_network.R` (`build_network_graph`/rede completa). **Novos**:
  correlation heatmap dos perfis de genes-alvo entre proteínas
  (`correlation_heatmap_genes.png`) e Jaccard com valores numéricos anotados
  (`jaccard_heatmap_valores.png`). **Figuras de ChIP-QC melhoradas**
  (`save_sample_qc_plots()`): fragment size agora marca o comprimento de
  fragmento estimado com linha tracejada + anotação; fingerprint ganhou
  subtítulo com o valor de SSD e interpretação; ambas com `theme_minimal()`.
  Regeneradas via cache RDS (sem recomputar) para as 14 amostras WT.
  **`RESUMO_METANALISE.md`** (`Arquivos/metanalise/`, versionado
  explicitamente mesmo com a pasta no `.gitignore` — é conteúdo escrito à mão,
  não regenerável por script) consolida todos os achados em um único
  documento: genes-alvo por proteína, Jaccard, interseções (cumulativas e
  exclusivas), os 8 genes-núcleo, enriquecimento, e o que cada figura mostra.
  **Pendência**: `Figuras/chip_qc/correlation_heatmap.png` e `pca.png` do
  batch original de XPC+H3K4me3 (19 amostras, rodado em 2026-07-17) estavam
  **em branco** (a execução antecedeu a correção do bug de `theme_void()`
  transparente); relançado em background (`run_xpc_chipqc_batch.R`, ~10h,
  decisão do usuário) — desta vez `run_chipqc_batch()` salva o
  `ChIPQCexperiment` como `chipqc_experiment.rds`, então nunca mais precisa
  recomputar para replotar.
- **2026-07-19** — **Rede bipartida legível dos hotspots de ocupação máxima.**
  A rede completa do Módulo 18 (`Figuras/network/bipartite_network.png`) usa
  todos os hotspots com `occupancy_score>=2` (46.477 regiões — o threshold do
  Módulo 17 é baixo demais para uma rede legível, quase toda região tem 2
  proteínas por acaso) e continua um "hairball" mesmo após o fix do fundo
  transparente. Investigado: só **9 regiões** têm `occupancy_score` máximo (=4,
  ocupadas por 4 proteínas simultaneamente) — geradas como
  `Figuras/network/bipartite_network_top_hotspots.png`
  (`scratchpad/run_bipartite_top_hotspots.R`), totalmente legível. **Achado de
  consistência**: 2 desses 9 genes hotspot (PHACTR4, IFI44L) são exatamente 2
  dos 8 genes-núcleo já identificados na camada XPC-âncora (seção 5 do
  `RESUMO_METANALISE.md`) — os dois métodos independentes (nível de região via
  Módulo 17, nível de gene via ChIPseeker) convergem. A figura também mostra
  visualmente que ELK1 nunca compartilha hotspot com XPC (confirma Jaccard=0).
  **Usuário reportou que `bipartite_network.png` (saída padrão do Módulo 18)
  continuava ilegível** ("tudo preto, sem fundo") mesmo após o fix do
  `theme_void()`. Investigado a fundo: o fundo já era branco de fato — o
  efeito "tudo preto" é visual, causado pela densidade de ~46.477 nós/arestas
  sobrepostos (mesmo filtrando para score≥3 ainda sobram 1.342 genes únicos).
  **Confirmado: é puramente um problema de escala, não de estilo/renderização**
  — nenhum PNG estático consegue representar essa quantidade de nós de forma
  legível. **Ação**: `Figuras/network/bipartite_network.png` foi
  **sobrescrito com o conteúdo de `bipartite_network_top_hotspots.png`**
  (as 9 regiões de ocupação máxima) — agora ambos os nomes de arquivo mostram
  a mesma versão legível. A rede completa (46 mil+ nós) permanece disponível
  em `Arquivos/network/bipartite_network.graphml` para exploração interativa
  em Cytoscape/Gephi, que é a ferramenta correta para essa escala — não um
  PNG. `RESUMO_METANALISE.md` atualizado para refletir isso.
- **2026-07-19** — **Causa raiz real do correlation heatmap/PCA vazios do
  ChIP-QC.** Não era bug de renderização (o fix de `theme_void()`/fundo
  transparente do dia anterior era real mas não era a causa deste problema).
  `plotCorHeatmap()`/`plotPrincomp()` do `ChIPQC` **retornam `NULL` com a
  mensagem "No peaks to plot."** quando a `ChIPQCexperiment` não tem picos —
  e `run_xpc_chipqc_batch.R` chamava `ChIPQC()` sem a coluna `Peaks` (desenho
  original: essa etapa de QC roda antes do peak calling). O batch rodou de
  verdade (~4h40, mais rápido que a estimativa de 10h) mas produziu PNGs
  vazios porque não havia dado nenhum para plotar — confirmado inspecionando
  a classe do retorno (`NULL`) diretamente no R. **Decisão do usuário**: com
  o peak calling do XPC (16 amostras) já pronto, rodar peak calling rápido
  (broad, `--nolambda`) também para as 3 amostras H3K4me3 (nunca usadas como
  alvo científico, só input substituto do XPC-WT — nova entrada
  `H3K4me3` em `PEAKCALLING_CONFIG`) e relançar o ChIPQC em lote das 19 agora
  COM picos para todas (`run_h3k4me3_peaks_and_chipqc_batch.R`).

  **Resultado**: com picos para as 19 amostras, `correlation_heatmap.png` e
  `pca.png` renderizaram conteúdo real pela primeira vez (~20min, bem mais
  rápido que os batches sem picos). Achados: réplicas XPC-WT têm correlação
  quase nula entre si (blocos isolados na diagonal — condizente com o
  consenso WT pequeno e específico, 1.577 regiões); XPC-KO e H3K4me3 mostram
  mais correlação cruzada entre si. No PCA, GSM6600715 (maior amostra do
  lote, 111M reads) é um outlier isolado no PC1 (62% da variância).

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
| MACS3 (3.0.4) | peak calling | Módulo 07 — instalado dentro de um ambiente conda (`chipseq`) **no WSL** (Ubuntu-24.04), não no PATH nativo do Windows — ver detalhes abaixo |

**Substituídos por pacotes R nativos (não precisam mais estar no PATH):**

| Ferramenta externa (spec original) | Substituto R nativo | Módulo |
|---|---|---|
| SRA Toolkit (`prefetch`/`fasterq-dump`) | API pública da ENA (`resolve_srr_table()`/`download_fastq_ena()`, `download.file()`) — download direto do FASTQ já pronto, sem precisar converter `.sra` | 01 — testado com as 89 amostras reais |
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
