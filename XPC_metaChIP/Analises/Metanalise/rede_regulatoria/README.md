# rede — rede bipartida Proteína→Região→Gene

**Normalização:** nenhuma (rede binária de ocupação, não gene-set).

**Amostras/estado por proteína:** `XPC_WT`/`STAT1_WT`/`STAT2_WT`/`IRF9_WT`/
`ELK1_WT` — pool do Módulo 13/14 (**UN+IFNα2h**, o mesmo do design
`../Metanalise/design_anterior_pre_revisao/`, **não** a seleção "só IFNα2h" da metanálise
principal atual). Isso não foi atualizado na revisão de 2026-07-22 porque o
pedido era sobre os gene-sets da metanálise, não sobre a rede — ver nota de
consistência no `README.md` do projeto antes de comparar números entre
`rede/` e `Metanalise/principal_sem_normalizacao/`.

**Conteúdo:**
- `bipartite_network.graphml`/`.png` + `bipartite_edges.csv`/
  `bipartite_metrics.csv` — rede completa, `occupancy_score≥2` (46 mil+
  regiões — ilegível como imagem estática, usar o GraphML no Cytoscape/Gephi).
- `bipartite_network_score3plus.graphml`/`.png` +
  `hotspots_score3plus_com_signal.csv` — versão focada,
  `occupancy_score≥3` (1462 regiões: score3=1453, score4=9, score5=0),
  rótulos só nos 9 genes de score 4 + top-5 de score 3 por `signalValue`
  MACS3.
- `bipartite_network_top_hotspots.png` — as 9 regiões de ocupação máxima
  (score 4) isoladas, legível sem GraphML.
