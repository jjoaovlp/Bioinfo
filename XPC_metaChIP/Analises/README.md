# Analises/ — índice

Tudo aqui se divide em dois ramos claros, mais uma pasta transversal:

## `XPC/` — análises de uma proteína só

Nenhum arquivo aqui cruza XPC com outra proteína. Ver `XPC/README.md`.
- `individual/` — anotação e enriquecimento do XPC-WT (geral e por timepoint)
- `timepoints/` — XPC 0h vs 3h, cada um comparado contra o eixo interferon
- `diffbind/` — WT vs XPC-KO (`atual_input_pareado_TMM/` vs
  `ANTES_nolambda_semTMM/` — ver `diffbind/README.md`)

## `Metanalise/` — tudo que cruza 2 ou mais proteínas

Ver `Metanalise/*/README.md` de cada subpasta para normalização e amostras.
- `principal_sem_normalizacao/` — análise principal vigente, IFNα2h, todos os picos
- `principal_normalizado_topN/` — igual, mas top-1000 picos/proteína (nivela STAT1/STAT2)
- `baseline_controle/` — controle: só amostras untreated
- `design_anterior_pre_revisao/` — design ORIGINAL (até 2026-07-22, UN+IFNα2h pooled) — legado, preservado
- `nucleo_XPC_interferon/` — destaque: os 8 genes XPC∩STAT1∩STAT2∩IRF9 (idêntico entre design antigo e atual)
- `rede_regulatoria/` — rede bipartida Proteína→Região→Gene (usa pool UN+IFNα2h — ver nota de consistência no `README.md` do projeto)

**Qual usar por padrão:** `Metanalise/principal_sem_normalizacao/`
(ou `principal_normalizado_topN/` para checar robustez ao volume de picos do
STAT). `design_anterior_pre_revisao/` só para auditoria histórica.

## `qc_comparativo/` — transversal (não é resultado de análise)

Figuras/tabelas que comparam **antes vs depois** de cada normalização
aplicada em 2026-07-22 (re-call do XPC-KO, TMM no diffbind, top-N nos
gene-sets) + painel de métricas de QC por amostra. Ver `qc_comparativo/README.md`.

---

Para o desenho experimental completo (amostras por GSM, o que foi baixado/
alinhado/peak-called/usado em cada etapa), ver o `README.md` na raiz do
projeto. Para o relato completo da metanálise (achados, ressalvas), ver
`RESUMO_METANALISE.md`.
