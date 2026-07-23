# Resumo da Metanálise WT — XPC como Referência

Gerado em 2026-07-19. Pipeline ChIP-seq (5 proteínas WT: XPC, STAT1, STAT2, IRF9, ELK1).

## 0. REVISÃO (2026-07-22) — leia isto primeiro

As seções 1–10 abaixo descrevem a metanálise **original** (STAT1/STAT2/IRF9 =
untreated+IFNα2h pooled). Em 2026-07-22 essa análise foi **revisada e substituída**
como principal, por decisão do usuário, após diagnóstico de desbalanço severo entre
proteínas (STAT1/STAT2 saturando com 40-62 mil picos no IFNα2h, XPC/IRF9/ELK1 na
casa das centenas). O conteúdo original **não foi apagado** — foi movido para
`Analises/meta_geral_ANTES_revisao/` (arquivos + figuras + enriquecimento por
combinação), preservado como registro histórico/auditoria. Todas as referências de
caminho nas seções 1–10 abaixo (`Arquivos/metanalise/...`, `Figuras/metanalise/...`,
`Arquivos/enrichment/...`) apontam para essa pasta legada, não para a estrutura atual.

**O que mudou na análise principal (agora em `Analises/meta_geral/`):**
- STAT1/STAT2/IRF9 passam a usar **só IFNα 2h** (estado ativo puro), descartando o
  "untreated" que estava no pool. XPC continua pooled (10 WT, 0h+1h+3h); ELK1
  continua constitutivo (2 ENCODE). Motivo: o "untreated" quase não contribuía
  sinal (STAT2 UN = 5 regiões no genoma inteiro, ver seção 9) e só diluía o
  desenho sem ganho biológico.
- Nova variante **top-N** (`Analises/meta_topN/`, N=1000 picos por proteína
  rankeados por `signalValue` do MACS3): reduz STAT1 de ~19.9k→927 genes
  (nearest) e STAT2 de ~18.5k→919, deixando-os na mesma ordem de grandeza de
  XPC (817)/IRF9 (921)/ELK1 (212 — inalterado, já tinha <1000 picos). Serve
  para checar se as interseções XPC∩STAT/IRF9 são robustas à normalização de
  tamanho do gene-set, não só um artefato da amplitude do STAT.
- **Diffbind do XPC corrigido**: o XPC-KO era peak-called com `--nolambda` (sem
  input), gerando 250 mil+ picos de background que confundiam o consenso
  WT∪KO (3283 regiões "Gained" implausíveis). Corrigido usando o mesmo input
  H3K4me3 do WT, pareado por timepoint — ver CLAUDE.md §9.1-REVISÃO. Os picos
  do KO caíram para 0-8 por amostra (comportamento correto de controle
  negativo); refeito com **TMM**, o resultado é 0 Gained/0 Lost/1576 Stable
  (nenhuma diferença sobrevive a FDR<0.05 — ver `Analises/XPC/diffbind/`).
- **Metanálise por timepoint do XPC** (`Analises/XPC/timepoints/`): 0h e 3h
  como âncoras separadas contra o eixo interferon (1h documentado como vazio —
  0 picos nas 3 réplicas, mesma causa em qualquer análise, não é bug da
  metanálise). 0h∩interferon é muito maior que 3h∩interferon porque o
  consenso XPC-WT é dominado por GSM6600715 (0h) — ver seção 8 abaixo.
- **Rede bipartida em `occupancy_score≥3`** (1462 regiões, `Analises/rede/`):
  inclui todas as regiões, rotula só os 9 genes de score 4 + top-5 de score 3
  por `signalValue` (interferon-stimulated genes clássicos: MX1, IRF9, RIGI,
  HERC6, DDX60).
- **Nova organização de pastas**: tudo em `Analises/`, uma pasta por
  interseção de genes (`Analises/meta_geral/interseccoes/<combo>/`,
  `Analises/meta_topN/interseccoes/<combo>/`), tudo do XPC junto
  (`Analises/XPC/{individual,timepoints,nucleo,diffbind}/`), baseline em
  `Analises/meta_baseline/`, figuras comparativas pré/pós em
  `Analises/qc_comparativo/`.
- **Venns de 5 conjuntos com ELK1 removidos** (`venn_untreated_*_5proteinas.png`)
  — o usuário achou o resultado mais legível só com 4 conjuntos (ELK1
  permanece nos UpSet, que suportam melhor 5 conjuntos).
- **FRiP permanece pendente** — ChIPQC rodou sem picos supridos para a maioria
  das amostras; regenerar exigiria rodar `ChIPQCsample()` com picos para 30+
  amostras (custoso). Ver `Analises/qc_comparativo/PENDENCIA_FRiP.txt`.

**Organização adicional (mesmo dia, reorganização fina):** cada subpasta de
`Analises/` ganhou um `README.md` próprio indicando explicitamente se o
conteúdo tem normalização (TMM/top-N) ou não, e quais amostras/estado por
proteína — ver `meta_geral/README.md`, `meta_topN/README.md`,
`meta_baseline/README.md`, `meta_geral_ANTES_revisao/README.md`,
`XPC/README.md`, `XPC/diffbind/README.md`, `rede/README.md`,
`qc_comparativo/README.md`. `XPC/diffbind/` foi dividido em
`atual_input_pareado_TMM/` (vigente) vs `ANTES_nolambda_semTMM/` (legado) —
antes os dois conjuntos ficavam misturados na mesma pasta, só diferenciados
pelo sufixo do nome do arquivo. Também recuperados para
`meta_geral_ANTES_revisao/` os arquivos de Jaccard **por região** (não por
gene) do antigo Módulo 15/16 (`jaccard_regiao_pairwise.csv`,
`jaccard_heatmap_regiao.png`, `regioes_compartilhadas_todas_proteinas.csv`),
que tinham ficado para trás em `Arquivos/overlap/`/`Figuras/overlap/` na
reorganização anterior.

---

Todos os dados brutos (CSV/TXT) da versão ORIGINAL (pré-revisão) estão em
`Analises/meta_geral_ANTES_revisao/`.

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

## 10. TODAS as combinações de 2+ proteínas (não só ancoradas em XPC) + enriquecimento

Pedido do usuário: as seções 3–4 só cobriam combinações que incluíam XPC. Os dois UpSet
(`upset_5_proteinas.png` = gene mais próximo, `upset_promotor_5_proteinas.png` = só
promotor ±3kb) mostram **todas** as interseções possíveis entre as 5 proteínas, inclusive
as que não envolvem XPC (ex. STAT1∩STAT2 sozinhas). Gerado um CSV único e completo por
versão — `Arquivos/metanalise/genes_comuns_todas_combinacoes_{nearest,promotor}.csv` —
com todas as 26 combinações de 2 a 5 proteínas (contagem + lista de genes), e rodado
enriquecimento funcional (GO/KEGG/Reactome/Hallmark) para cada combinação com genes
suficientes (`Arquivos/enrichment/<combinação>_{nearest,promotor}_{go,kegg,reactome,hallmark}.csv`).

**Título das figuras** (pedido do usuário): os dois UpSet agora trazem, no próprio
gráfico, a que condições eles se referem — "WT apenas; XPC=todos os timepoints post-UV;
STAT1/STAT2/IRF9=untreated+IFNα 2h; ELK1=ENCODE" — e deixam explícito se é a versão
"gene mais próximo" ou "só promotor".

**Log de status por combinação** (`enriquecimento_log_{nearest,promotor}.csv`) — versão
promotor:

| Combinação | N genes | Status |
|---|---|---|
| STAT1_STAT2 | 13017 | enriquecido |
| STAT1_STAT2_IRF9 | 477 | enriquecido |
| STAT1_IRF9 / STAT2_IRF9 | 484 | enriquecido |
| STAT1_ELK1 / STAT2_ELK1 | 168-169 | enriquecido |
| XPC_STAT1 / XPC_STAT2 | 164-170 | enriquecido |
| XPC_STAT1_STAT2 | 157 | enriquecido |
| XPC_STAT1_IRF9 / XPC_STAT2_IRF9 / XPC_IRF9 | 8 | enriquecido (mesmo núcleo da seção 5) |
| STAT1_IRF9_ELK1 / STAT2_IRF9_ELK1 / IRF9_ELK1 | 7 | 0 termos significativos |
| **qualquer combinação com XPC + ELK1** | **0** | insuficiente — XPC nunca compartilha gene de promotor com ELK1 |

Confirma, de forma exaustiva (todas as 26 combinações, não só as com XPC), o padrão já
visto: XPC e ELK1 nunca convergem em nenhuma combinação; STAT1/STAT2/IRF9 convergem
fortemente entre si (esperado, complexo ISGF3); as combinações com XPC ficam
consistentemente em torno de 8 genes quando incluem todo o eixo interferon — o mesmo
núcleo de 8 genes da seção 5, agora confirmado como estável através de múltiplas
combinações parciais também.
