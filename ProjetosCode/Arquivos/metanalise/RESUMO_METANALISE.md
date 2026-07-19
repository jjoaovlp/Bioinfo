# Resumo da Metanálise WT — XPC como Referência

Gerado em 2026-07-19. Pipeline ChIP-seq (5 proteínas WT: XPC, STAT1, STAT2, IRF9, ELK1).
Todos os dados brutos (CSV/TXT) estão em `Arquivos/metanalise/`; figuras em `Figuras/metanalise/`.

## 1. Amostras e genes-alvo por proteína

| Proteína | Timepoints (WT) | Genes-alvo (promotor ±3kb) | Jaccard máx. com outra proteína |
|---|---|---|---|
| XPC | 0h/1h/3h pós-UV (10 amostras, broad) | 264 | 0.003 (STAT2) — praticamente ortogonal |
| STAT1 | untreated + IFNα 2h (4 amostras, narrow) | 14.840 | 0.616 (STAT2) |
| STAT2 | untreated + IFNα 2h (4 amostras, narrow) | 13.896 | 0.616 (STAT1) |
| IRF9 | untreated + IFNα 2h (4 amostras, narrow) | 501 | 0.016 (STAT2) |
| ELK1 | 2 réplicas ENCODE, com input (2 amostras, narrow) | 189 | 0.002 (STAT1) |

**"Genes-alvo (promotor)"** = gene com pelo menos um pico da proteína dentro de ±3kb do TSS
(ChIPseeker, `TxDb.Hsapiens.UCSC.hg38.knownGene`). STAT1/STAT2 são amplos porque, no timepoint
IFNα2h, têm 40–62 mil picos (ativação massiva do complexo ISGF3).

## 2. Índice de Jaccard entre proteínas (nível de REGIÃO, não gene)

Ver `Figuras/metanalise/jaccard_heatmap_valores.png`. Jaccard = |A∩B| / |A∪B| das coordenadas
de pico. Só mede sobreposição *física* de picos — por isso XPC parece "não ter nada a ver" com
as demais aqui (0.000–0.003): os picos de XPC (reparo de DNA, pós-UV) raramente caem exatamente
sobre um pico de STAT/IRF/ELK1 (resposta a interferon / fator ubíquo). A relação real aparece só
quando se olha o **gene mais próximo** de cada pico (seção 3), não a coordenada exata.

| | XPC | STAT1 | STAT2 | IRF9 | ELK1 |
|---|---|---|---|---|---|
| **XPC** | 1.000 | 0.002 | 0.003 | 0.001 | 0.000 |
| **STAT1** | | 1.000 | 0.616 | 0.013 | 0.001 |
| **STAT2** | | | 1.000 | 0.016 | 0.002 |
| **IRF9** | | | | 1.000 | 0.001 |
| **ELK1** | | | | | 1.000 |

## 3. Interseções de genes ancoradas no XPC (todas as combinações)

Fonte: `xpc_genes_promotor_todas_combinacoes.csv` (versão promotor) e
`xpc_genes_comuns_todas_combinacoes.csv` (versão nearest-gene, mais permissiva).
Aqui a contagem é **cumulativa**: "XPC∩STAT1∩STAT2" inclui os genes que também estão em IRF9/ELK1.

| Combinação | Nº proteínas | Genes em comum com XPC (promotor) |
|---|---|---|
| XPC ∩ STAT1 ∩ STAT2 ∩ IRF9 ∩ ELK1 | 5 | 0 |
| **XPC ∩ STAT1 ∩ STAT2 ∩ IRF9** | **4** | **8** ← núcleo do eixo interferon |
| XPC ∩ STAT1 ∩ STAT2 ∩ ELK1 | 4 | 0 |
| XPC ∩ STAT1 ∩ IRF9 ∩ ELK1 | 4 | 0 |
| XPC ∩ STAT2 ∩ IRF9 ∩ ELK1 | 4 | 0 |
| XPC ∩ STAT1 ∩ STAT2 | 3 | 157 |
| XPC ∩ STAT1 ∩ IRF9 | 3 | 8 |
| XPC ∩ STAT2 ∩ IRF9 | 3 | 8 |
| XPC ∩ STAT1 ∩ ELK1 | 3 | 0 |
| XPC ∩ STAT2 ∩ ELK1 | 3 | 0 |
| XPC ∩ IRF9 ∩ ELK1 | 3 | 0 |
| XPC ∩ STAT1 | 2 | 170 |
| XPC ∩ STAT2 | 2 | 164 |
| XPC ∩ IRF9 | 2 | 8 |
| XPC ∩ ELK1 | 2 | 0 |

## 4. Interseções EXCLUSIVAS (estrutura do UpSet plot)

Fonte: `metanalise_interseccoes_completas_promoter.csv`. Diferente da seção 3, aqui cada gene
conta em **apenas uma linha** — a combinação exata de proteínas que o têm como alvo (nem mais,
nem menos). É a mesma lógica de barras do `Figuras/metanalise/upset_promotor_5_proteinas.png`.

| Combinação exclusiva | Nº proteínas | Nº genes |
|---|---|---|
| XPC & STAT1 & STAT2 & IRF9 | 4 | 8 |
| STAT1 & STAT2 & IRF9 & ELK1 | 4 | 7 |
| STAT1 & STAT2 & IRF9 | 3 | 462 |
| STAT1 & STAT2 & ELK1 | 3 | 155 |
| XPC & STAT1 & STAT2 | 3 | 149 |
| *(demais combinações — ver CSV completo)* | | |

O grupo **STAT1&STAT2&IRF9 exclusivo (462 genes, sem XPC)** é rico em ISGs canônicos
(CXCL10, IFIT1/2/3/5, ISG15, IFI6, DDX60, GBP1, IRF1, CD274/PD-L1) — confirma que a assinatura de
interferon foi capturada corretamente pelo peak calling, mesmo fora do núcleo comum com XPC.

## 5. Núcleo XPC ∩ eixo interferon completo — os 8 genes centrais

`genes_nucleo_XPC_interferon.txt`:

**ARHGAP29-AS1, IFI44L, LOC105371874, MEF2A, MIR21, PHACTR4, RIGI, SLC1A3**

Dois são ISGs (Interferon-Stimulated Genes) canônicos:
- **RIGI** (DDX58/RIG-I) — sensor citoplasmático de RNA viral, ativa a resposta antiviral inata.
- **IFI44L** — Interferon-Induced Protein 44-Like, marcador clássico de assinatura de interferon.

### Enriquecimento funcional (Reactome, GO, KEGG, Hallmark)

Só o **Reactome** encontrou termo significativo (esperado — 8 genes é uma lista pequena;
GO/KEGG/Hallmark ficaram vazios por falta de poder estatístico):

| Via | p.adjust | Genes | Fold enrichment |
|---|---|---|---|
| **Modulation of host responses by IFN-stimulated genes** (R-HSA-9909505) | **0.022** | RIGI, IFI44L | 332× |

Isso **valida biologicamente** a análise: os genes que o XPC compartilha com STAT1/STAT2/IRF9
caem exatamente na via de resposta a interferon — não é um artefato do método.

## 6. Figuras geradas (e o que cada uma mostra)

| Arquivo | Conteúdo |
|---|---|
| `venn_xpc_stat1_stat2_irf9.png` | Venn (nearest-gene) — 4 conjuntos, números inflados por STAT1/2 |
| `venn_promotor_xpc_stat1_stat2_irf9.png` | Venn (promotor, recomendado) — mesma comparação, mais específica |
| `upset_5_proteinas.png` / `upset_promotor_5_proteinas.png` | UpSet das 5 proteínas — todas as combinações, incl. ELK1 |
| `rede_focada_xpc.png` | Rede proteína→gene **só** dos 8 genes-núcleo + XPC∩IRF9/ELK1 (legível) |
| `correlation_heatmap_genes.png` | Correlação de Pearson dos perfis binários gene×proteína |
| `jaccard_heatmap_valores.png` | Jaccard par-a-par com valores numéricos (nível de região) |

**Nota sobre `Figuras/network/bipartite_network.png`** (Módulo 18, gerada antes destas): usa
**todos** os hotspots (score ≥2 proteínas) e todos os genes — por volume vira um "hairball"
ilegível. Não é um bug de dados, é o preço de plotar milhares de nós; para uma visão legível,
usar `rede_focada_xpc.png` (subconjunto específico e pequeno).

## 7. Ressalva metodológica

A amplitude de STAT1/STAT2 (~14 mil genes-alvo de promotor) é parte biologia real (STATs
induzidos por IFN ligam massivamente ao DNA no pico da resposta, 2h) e parte possível
sobre-sensibilidade do peak calling (40–62 mil picos por amostra, bem mais que XPC/IRF9/ELK1).
Por isso as interseções mais confiáveis nesta análise são as que envolvem os conjuntos
específicos (XPC, IRF9, ELK1) — a seção 5 (núcleo de 8 genes) é o resultado central e mais
robusto da metanálise.
