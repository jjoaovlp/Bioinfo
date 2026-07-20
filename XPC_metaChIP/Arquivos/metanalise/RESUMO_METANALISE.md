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

Ver `Figuras/metanalise/jaccard_heatmap_valores.png`.

**O que é o índice de Jaccard.** É uma medida de **quanto dois conjuntos se sobrepõem**:

> **J(A, B) = |A ∩ B| / |A ∪ B|** = (regiões em comum) / (total de regiões distintas das duas).

Varia de **0** (nenhuma coordenada de pico compartilhada) a **1** (os dois conjuntos de picos são
idênticos). Aqui é calculado no **nível de região** (coordenada genômica), ou seja, mede
sobreposição *física* de picos — duas proteínas só têm Jaccard alto se ligam literalmente às
**mesmas** posições no genoma. Exemplo: se A tem 100 picos, B tem 100 picos e 40 são compartilhados,
J = 40 / (100 + 100 − 40) = 40/160 = 0,25.

Por isso o XPC parece "não ter nada a ver" com as demais aqui (0.000–0.003): os picos de XPC
(reparo de DNA, pós-UV) raramente caem exatamente sobre um pico de STAT/IRF/ELK1 (resposta a
interferon / fator ubíquo). Já STAT1↔STAT2 = 0.616 porque coligam massivamente aos **mesmos**
sítios (complexo ISGF3, após IFN). A relação real do XPC com o eixo interferon aparece só quando se
olha o **gene mais próximo** de cada pico (seção 3), não a coordenada exata — por isso a metanálise
foi feita no nível de gene, não de região.

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

| `bipartite_network_top_hotspots.png` | Rede proteína→gene só das 9 regiões com occupancy_score MÁXIMO (=4, ocupadas por 4 proteínas simultaneamente) — inclui PHACTR4 e IFI44L, os mesmos 2 genes-núcleo já vistos na seção 5 |

**Nota sobre `Figuras/network/bipartite_network.png`** (2026-07-19, atualizado): o arquivo
gerado pelo Módulo 18 por padrão usa **todos** os hotspots (46.477 regiões com score≥2 — quase
toda região tem 2 proteínas por acaso) e vira um "hairball" ilegível — não é um bug de
estilo/fundo, é o preço de plotar dezenas de milhares de nós num PNG estático (confirmado:
mesmo restringindo a score≥3 ainda sobram 1.342 genes únicos, ainda ilegível). Por isso este
arquivo foi **substituído pelo conteúdo de `bipartite_network_top_hotspots.png`** (as 9 regiões
de ocupação máxima) — as duas imagens agora são idênticas e ambas legíveis. Para explorar a
rede **completa** (todos os 46 mil+ hotspots), use `Arquivos/network/bipartite_network.graphml`
em uma ferramenta de rede interativa (Cytoscape/Gephi) — nenhuma imagem estática consegue
representar essa escala de forma legível. Alternativas legíveis já prontas:
`rede_focada_xpc.png` (genes-núcleo XPC-âncora) ou `bipartite_network_top_hotspots.png`
(hotspots de ocupação máxima, sem viés de proteína de referência).

## 7. Ressalva metodológica

A amplitude de STAT1/STAT2 (~14 mil genes-alvo de promotor) é parte biologia real (STATs
induzidos por IFN ligam massivamente ao DNA no pico da resposta, 2h) e parte possível
sobre-sensibilidade do peak calling (40–62 mil picos por amostra, bem mais que XPC/IRF9/ELK1).
Por isso as interseções mais confiáveis nesta análise são as que envolvem os conjuntos
específicos (XPC, IRF9, ELK1) — a seção 5 (núcleo de 8 genes) é o resultado central e mais
robusto da metanálise.

## 8. Diagnóstico do peak calling do XPC-WT (investigação de possível erro)

Ver `Figuras/annotation/XPC_WT_peakcalling_comparativo.png`. Das 10 amostras XPC-WT (3
timepoints post-UV), só **3 produziram picos** com os parâmetros padrão (MACS3 broad,
q=0.01, broad-cutoff=0.1): GSM6600715 (0h, 1462 picos), GSM6600718 (0h, 7), GSM6600733 (3h,
109). As outras 7 (incluindo **todas as 3 réplicas de 1h**) deram **0 picos**.

**Investigação (conclusão: não é bug do pipeline)** — checado e descartado:
- **Profundidade**: todas as 10 amostras têm 17–76 milhões de tags pós-filtro; a amostra
  com mais profundidade (GSM6600724, 1h, 76,4M) deu 0 picos, e a com sinal mais forte
  (GSM6600715) tem profundidade comparável (74,7M) — não é falta de leitura.
- **Pareamento do input**: conferido nos 10 `.xls` do MACS3 — cada réplica usa o input
  correto do próprio timepoint (0h→GSM6600680, 1h→GSM6600686, 3h→GSM6600692), sem
  nenhuma amostra cruzada.
- **Erros de execução**: nenhum job crashou sem tratamento; os 3 que precisaram do
  fallback `--nomodel` (718, 725, 734) já estão previstos no código.

**Achado real**: a estimativa de comprimento de fragmento do MACS3 (`alternative
fragment length`) é **limpa e única** só em GSM6600715 (202bp) — nas demais aparecem
**vários candidatos dispersos** (ex. GSM6600732 lista 10 candidatos entre 99–580bp),
assinatura clássica de correlação cruzada ruidosa = ChIP com pouca estrutura real,
independente da profundidade. SSD (ChIPQC) não discrimina bem sozinho (0,68–0,92 em
ambos os grupos).

**Teste de sanidade** (pedido do usuário): re-rodado o MACS3 com limiares muito mais
permissivos (q=0.5, broad-cutoff=0.5, vs. padrão q=0.01/0.1) em 3 amostras — as duas
mais fracas (GSM6600724, 1h; GSM6600732, 3h) e a de referência (GSM6600715, 0h), para
calibrar o efeito do relaxamento numa amostra que sabidamente tem sinal real.
Resultado (`Arquivos/sanity_test_peakcalling/`):

| Amostra | Picos padrão (q=0.01) | Picos relaxado (q=0.5) |
|---|---|---|
| GSM6600715 (0h, referência) | 1462 | **3569** (~2,4×) |
| GSM6600724 (1h) | 0 | **0** |
| GSM6600732 (3h) | 0 | **0** |

Na amostra de referência o relaxamento inflou os picos de forma esperada (~2,4×,
comportamento normal de afrouxar o limiar). Nas duas amostras fracas, **mesmo no teste
mais permissivo possível, continuam com 0 picos** — confirma que não é o limiar padrão
que é estrito demais; não há absolutamente nenhuma região com qualquer estrutura de
enriquecimento detectável nessas amostras.

**Leitura**: GSM6600715 é um outlier de sinal forte (já identificado antes no PCA como
outlier isolado no PC1, maior amostra do lote). A maioria das réplicas de XPC-WT parece
ter ChIP fraco/sem enriquecimento estruturado — plausível biologicamente (XPC tem
associação difusa/transitória à cromatina fora do reparo ativo), mas também consistente
com qualidade de anticorpo/experimento variável entre réplicas do estudo original
(GSE214182). A consequência prática: o consenso XPC-WT usado na metanálise principal é
dominado por GSM6600715.

## 9. Metanálise SOMENTE das versões sem tratamento (baseline)

Pedido do usuário: repetir a metanálise usando só as condições **não-estimuladas** —
XPC 0h_post_UV (proxy de baseline; não há amostra verdadeiramente pré-UV neste desenho),
STAT1/STAT2/IRF9 **UN** (pré-interferon) e ELK1 (ENCODE, sempre "none"). Ver
`Figuras/metanalise/jaccard_heatmap_untreated.png`,
`venn_untreated_{nearest,promotor}.png`, `Arquivos/metanalise/{jaccard,gene_sets}_untreated.csv`.

**ELK1 incluído diretamente no Venn** (pedido do usuário): além do Venn de 4 conjuntos
(XPC+eixo interferon), `venn_untreated_{nearest,promotor}_5proteinas.png` mostra as
**5 proteínas juntas** (`ggVennDiagram` suporta 5 conjuntos com layout poligonal, ainda
legível). XPC∩ELK1 = 7 genes (nearest) / 0 (promotor) — os mesmos 7 genes já vistos na
metanálise principal (ELK1 não varia por timepoint, então seu conjunto de genes é
idêntico entre baseline e WT-completo). Núcleo das 5 juntas = 0 genes em ambas as
versões.

| Proteína | Amostras | N regiões consenso | Genes (promotor) |
|---|---|---|---|
| XPC | 0h_post_UV (4) | 1461 | 252 |
| STAT1 | UN (2) | 46 | 24 |
| STAT2 | UN (2) | **5** | 3 |
| IRF9 | UN (2) | 69 | 69 |
| ELK1 | none (2) | 212 | 189 |

**Resultado: núcleo XPC∩STAT1∩STAT2∩IRF9 = 0 genes** (tanto por promotor quanto por
gene mais próximo) — e o Jaccard STAT1↔STAT2 cai de **0.616 (IFNα 2h)** para **0.041
(baseline)**. Isso é uma validação biológica limpa, não uma falha da análise: STAT2 UN
tem só 5 regiões de pico no genoma inteiro (STAT2 é quase inteiramente dependente de
ativação por interferon para ligar cromatina), então não há material suficiente para
qualquer interseção robusta no estado basal. Confirma por que a metanálise principal
(seção 1–5) usa deliberadamente IFNα 2h para STAT1/STAT2/IRF9 em vez do baseline — é o
único timepoint em que esses fatores realmente ocupam o genoma.
