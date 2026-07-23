# Metanalise/principal_sem_normalizacao — metanálise principal, SEM normalização de top-N

**Normalização:** nenhuma restrição de tamanho de gene-set — usa **todos** os
picos do estado ativo de cada proteína (o oposto de `Metanalise/principal_normalizado_topN/`, que
restringe a top-1000 por proteína).

**Amostras/estado por proteína** (ver `README.md` do projeto para a tabela
completa por GSM):
- **XPC** = pooled WT, 10 amostras (0h+1h+3h pós-UV)
- **STAT1** = só IFNα 2h (GSM6928567/568)
- **STAT2** = só IFNα 2h (GSM6928587/588)
- **IRF9** = só IFNα 2h (GSM6928599/600)
- **ELK1** = constitutivo, 2 réplicas ENCODE (GSM2423754/755)

**Ponto de atenção:** STAT1/STAT2 saturam (~19-20 mil genes-alvo nearest,
~14-15 mil promotor) por ativação massiva do ISGF3 no IFNα2h — interseções
amplas envolvendo STAT aqui tendem a ser dominadas pelo volume de picos do
STAT. Para checar robustez sem esse viés de volume, ver `Metanalise/principal_normalizado_topN/`.

**Conteúdo:** `gene_sets_por_proteina_{nearest,promotor}.csv`, `jaccard_
{nearest,promotor}.csv`, `upset_5_proteinas_{nearest,promotor}.png` (5
proteínas, inclui ELK1), `venn_{nearest,promotor}.png` (4 proteínas, sem
ELK1), `genes_comuns_todas_combinacoes_{nearest,promotor}.csv` (todas as 26
combinações de 2-5 proteínas) e `interseccoes/<combo>/` (1 pasta por
combinação com genes + enriquecimento GO/KEGG/Reactome/Hallmark + dotplots).

Ver também: `Metanalise/principal_normalizado_topN/` (top-N), `Metanalise/baseline_controle/` (controle untreated),
`Metanalise/design_anterior_pre_revisao/` (design anterior a 2026-07-22, UN+IFNα2h pooled).
