# meta_topN — metanálise principal, COM normalização top-N (N=1000)

**Normalização:** restringe cada proteína aos **top-1000 picos por
`signalValue`** (MACS3) antes de anotar — nivela o tamanho do gene-set entre
proteínas. Mesma seleção de estado/amostras que `meta_geral/` (XPC pooled WT,
STAT1/STAT2/IRF9 só IFNα2h, ELK1 constitutivo), só muda a normalização.

**Efeito da normalização** (nearest gene): STAT1 ~19.865→**927**, STAT2
~18.537→**919**, XPC 1.210→817 (caiu porque também é truncado a 1000, mas
tinha só 1577 picos), IRF9 1.012→921 (quase sem efeito, já tinha <1000
picos), ELK1 212→212 (inalterado, só 228 picos). Ver
`Analises/qc_comparativo/gene_sets_prepos_topN.png` para o comparativo visual.

**Uso recomendado:** checar se uma interseção que aparece em `meta_geral/`
sobrevive quando o volume de picos do STAT1/STAT2 é controlado — se sim, é
evidência mais forte de co-ligação real, não só efeito de escala.

**Conteúdo:** mesma estrutura de `meta_geral/` (gene_sets, jaccard, upset,
venn, combinações + `interseccoes/<combo>/`), mas com os gene-sets truncados.
