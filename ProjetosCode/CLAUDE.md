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
| 03 | `03_qc.R` | ✅ implementado (funções puras testadas) | FastQC + MultiQC + interpretação automática (PASS/WARN/FAIL) |
| 04 | `04_alignment.R` | ✅ implementado (parsing testado) | Bowtie2 + samtools, sort/index, `flagstat` |
| 05 | `05_filtering.R` | ✅ implementado (não executável sem ferramentas) | MarkDuplicates (via samtools, sem Picard), MAPQ, blacklist ENCODE (bedtools) |
| 06 | `06_chip_qc.R` | ⬜ pendente | DeepTools (fingerprint, PCA, coverage, correlation, fragment size) |
| 07 | `07_peakcalling.R` | ⬜ pendente | MACS3 (broad para XPC; narrow para ELK1/STAT1/STAT2) |
| 08 | `08_diffbind.R` | ⬜ pendente | DiffBind — apenas XPC (WT vs XPC-KO) e STAT2 (WT vs STAT1-KO); ELK1 e STAT1 ficam fora (sem braço deficiente disponível) |
| 09 | `09_annotation.R` | ⬜ pendente | ChIPseeker::annotatePeak() |
| 10 | `10_enrichment.R` | ⬜ pendente | clusterProfiler, ReactomePA, msigdbr |
| 11 | `11_granges.R` | ⬜ pendente | Conversão/padronização para GRanges |
| 12 | `12_genome_standardization.R` | ⬜ pendente | liftOver hg19→hg38 (ELK1) |
| 13 | `13_regulatory_universe.R` | ⬜ pendente | `GenomicRanges::reduce()` |
| 14 | `14_overlap_matrix.R` | ⬜ pendente | Matriz binária de ocupação |
| 15 | `15_pairwise_overlap.R` | ⬜ pendente | `findOverlaps()`, Jaccard, heatmaps |
| 16 | `16_multiple_overlap.R` | ⬜ pendente | Interseção múltipla (Reduce/intersect) |
| 17 | `17_hotspots.R` | ⬜ pendente | Occupancy score, ranking |
| 18 | `18_bipartite_network.R` | ⬜ pendente | igraph/tidygraph/ggraph, export Cytoscape |
| 19 | `19_pathway_network.R` | ⬜ pendente | STRINGdb, GO/Reactome network |
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

## 5. Dependências

### 5.1 Pacotes do R (via BiocManager/CRAN)

| Pacote | Uso |
|---|---|
| `here` | caminhos relativos |
| `BiocManager` | instalação de pacotes Bioconductor |
| `GEOquery`, `Biobase` | download de metadata/SOFT/FASTQ processados do GEO; acesso a `pData()` |
| `GenomicRanges`, `GenomicFeatures`, `rtracklayer` | GRanges, liftOver |
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

| Software | Uso | Observação |
|---|---|---|
| SRA Toolkit (`prefetch`, `fasterq-dump`) | download de FASTQ do SRA | Módulo 01 |
| FastQC | QC por amostra | Módulo 03 |
| MultiQC | agregação de QC | Módulo 03 |
| Bowtie2 | alinhamento | Módulo 04 |
| samtools | sort/index/stats | Módulos 04-05 |
| samtools (`fixmate`/`markdup`, sem Picard) | remoção de duplicatas | Módulo 05 |
| bedtools | remoção de regiões da blacklist ENCODE (`intersect -v`) | Módulo 05 |
| deepTools (`bamCoverage`, `plotFingerprint`, `multiBamSummary`, `plotCorrelation`, `plotPCA`) | fingerprint, PCA, coverage, correlação | Módulo 06 |
| MACS3 | peak calling | Módulo 07 |

**R está instalado** (versões 4.5.2 e 4.6.0 em `C:\Program Files\R\`, ~515 pacotes já
presentes incluindo todo `REQUIRED_PACKAGES` de §5.1) e os módulos 00-02 já foram
executados de verdade contra os 4 datasets reais (ver histórico acima). As ferramentas
de linha de comando (SRA Toolkit, FastQC, MultiQC, Bowtie2, samtools, MACS3, deepTools)
**continuam ausentes do PATH** (verificado em 2026-07-13) — necessárias a partir do
Módulo 03. `00_setup.R` verifica a presença de cada uma via `Sys.which()` e avisa (sem
interromper) quando ausente. Os módulos 03+ seguem sendo desenvolvidos como código
revisado (sintaxe validada via `Rscript`) mas não podem ser executados de ponta a ponta
até essas ferramentas serem instaladas.

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

## 10. Pendências futuras (TODO)

- [ ] **XPC-KO não tem nenhum input/controle disponível** (nem mesmo o substituto de
      H3K4me3, que só existe para WT/ASH1L-KO) — confirmado na execução real do
      Módulo 02 em 2026-07-13. O Módulo 07 terá que rodar MACS3 sem controle pareado
      para as 6 amostras XPC-KO (opção "peak calling sem controle", ex. `--nolambda`),
      com FRiP/FDR esperados menos confiáveis — documentar isso explicitamente no
      relatório do Módulo 07 e tratar como limitação conhecida na validação (Módulo 21).
- [ ] Confirmar montagem de genoma real (hg38 vs hg19) para XPC e STAT1/STAT2 a partir
      dos metadados SOFT completos no Módulo 01 (não estava explícita nas páginas GEO
      consultadas) — validar antes do Módulo 12.
- [ ] Buscar/baixar o input/controle ENCODE de ELK1 (ENCSR623KNM) via ENCODE portal, já
      que o input não aparece diretamente entre as 2 GSMs do GSE91923.
- [ ] Implementar Módulos 03–22 (QC, alinhamento, filtragem, ChIP-QC, peak calling,
      DiffBind, anotação, enriquecimento, GRanges, padronização de genoma, universo
      regulatório, matriz de ocupação, overlaps, hotspots, redes, visualização,
      validação, master pipeline).
- [ ] Nenhuma ferramenta externa (Bowtie2/MACS3/samtools/FastQC/MultiQC/deepTools) está
      instalada neste ambiente de desenvolvimento — execução real dos módulos 03+
      requer ambiente com essas dependências.
- [ ] Considerar IRF9 (ChIP'd em ambas as séries de STAT1/STAT2) como 5ª proteína
      opcional em versão futura do projeto (fora do escopo atual de 4 proteínas).
- [ ] Criar `Scripts/config/peakcalling_config.R` junto com o Módulo 07.
- [ ] Primeiro commit e primeira tag (`v0.1.0` — estrutura inicial) pendente de
      autorização do usuário para push ao GitHub.
