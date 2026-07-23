# Metanalise/baseline_controle — controle, amostras SEM tratamento (untreated)

**Normalização:** nenhuma (mesma lógica de `Metanalise/principal_sem_normalizacao/`, sem restrição
top-N) — a diferença aqui é a **seleção de amostras**, não a normalização.

**Amostras/estado por proteína:**
- **XPC** = só 0h pós-UV (4 réplicas, proxy de baseline — não há amostra
  verdadeiramente pré-UV neste desenho)
- **STAT1** = untreated (GSM6928563/564)
- **STAT2** = untreated (GSM6928583/584)
- **IRF9** = untreated (GSM6928595/596)
- **ELK1** = constitutivo (mesmas 2 réplicas — não varia por tratamento)

**Por que existe:** serve de controle negativo para a metanálise principal —
confirma que STAT1/STAT2/IRF9 dependem de ativação por interferon para
ocupar cromatina de forma robusta. Achado-chave: **STAT2 untreated tem só 5
regiões de pico no genoma inteiro** (vs. dezenas de milhares em IFNα2h).
Núcleo XPC∩STAT1∩STAT2∩IRF9 aqui = **0 genes** (tanto nearest quanto
promotor) — não é falha da análise, é o esperado dado que STAT2 quase não
liga cromatina sem estímulo.

**Conteúdo:** `gene_sets_untreated.csv`, `jaccard_untreated.csv`,
`venn_untreated_{nearest,promotor}.png`, `xpc_genes_untreated_*` (interseções
com XPC nesse estado).
