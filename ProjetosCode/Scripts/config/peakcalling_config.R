## ============================================================================
## Configuracao do Modulo 07 -- Peak Calling (MACS3)
## ============================================================================
## Parametros de peak calling por proteina, conforme especificado no projeto
## (CLAUDE.md S1/S7) e a decisao registrada em CLAUDE.md S9.1 (XPC-KO sem
## input -> --nolambda). Manter todos os parametros cientificos aqui, nunca
## hardcoded dentro de 07_peakcalling.R, para que qualquer mudanca de
## parametro fique registrada e facil de auditar (ver boas praticas do
## projeto: toda mudanca de parametro cientifico deve ser justificada e
## registrada no historico do CLAUDE.md).

PEAKCALLING_CONFIG <- list(
  XPC = list(peak_type = "broad", qvalue = 0.01, broad_cutoff = 0.1),
  ELK1 = list(peak_type = "narrow", qvalue = 0.01),
  STAT1 = list(peak_type = "narrow", qvalue = 0.01),
  STAT2 = list(peak_type = "narrow", qvalue = 0.01),
  ## IRF9: adicionada como 5a proteina da metanalise WT (re-escopo 2026-07-18,
  ## CLAUDE.md S4). Fator de transcricao (componente do complexo ISGF3 com
  ## STAT1/STAT2) -> picos narrow, mesmos parametros dos demais TFs.
  IRF9 = list(peak_type = "narrow", qvalue = 0.01),
  ## H3K4me3: marca histona de promotor ativo, usada nesta analise so como
  ## input/controle substituto do XPC-WT (CLAUDE.md S9 ponto 4), nunca como
  ## alvo cientifico proprio. Picos chamados aqui apenas para viabilizar
  ## correlation/PCA do ChIPQC em lote das 19 amostras XPC+H3K4me3 (essas 3
  ## amostras nao tem input proprio -> --nolambda, igual XPC-KO). Broad e' o
  ## tipo biologicamente correto p/ H3K4me3 (domanios largos em promotores
  ## ativos), mesmos parametros do XPC para consistencia.
  H3K4me3 = list(peak_type = "broad", qvalue = 0.01, broad_cutoff = 0.1)
)

## Genoma efetivo para MACS3 --gsize. "hs" = tamanho efetivo padrao do
## genoma humano usado pelo MACS3 (ver documentacao do MACS3); todos os 4
## datasets sao humanos apos a substituicao registrada em CLAUDE.md S9.
MACS3_GSIZE <- "hs"
