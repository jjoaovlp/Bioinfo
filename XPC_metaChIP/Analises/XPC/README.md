# XPC — todas as análises individuais do XPC (não cross-proteína)

Sem normalização de gene-set (não se aplica — são análises de uma única
proteína, não interseções entre proteínas):

- **`individual/`** — anotação (`XPC_WT_{0h,3h}_annotation.csv`, ocupação,
  comparativo de peak-calling por réplica) e enriquecimento funcional
  (GO/KEGG/Reactome/Hallmark) do XPC-WT geral e por timepoint (0h/3h; 1h
  vazio — 0 picos).
- **`timepoints/`** — metanálise XPC 0h vs 3h (âncoras separadas) contra o
  eixo interferon no estado ativo (IFNα2h, mesma seleção de `../meta_geral/`).
  1h documentado como vazio.
- **`nucleo/`** — os 8 genes XPC∩STAT1∩STAT2∩IRF9 (núcleo, estável entre o
  design antigo e o atual — ver `../meta_geral_ANTES_revisao/README.md`) +
  enriquecimento funcional desse núcleo.
- **`diffbind/`** — WT vs XPC-KO; ver `diffbind/README.md` para a distinção
  entre a versão com TMM+input pareado (atual) e a versão antiga (`--nolambda`,
  sem TMM).

Todas as amostras usadas aqui são as 10 WT + 6 KO de GSE214182 (0h/1h/3h
pós-UV) — ver a tabela completa por GSM no `README.md` do projeto.
