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
  STAT2 = list(peak_type = "narrow", qvalue = 0.01)
)

## Genoma efetivo para MACS3 --gsize. "hs" = tamanho efetivo padrao do
## genoma humano usado pelo MACS3 (ver documentacao do MACS3); todos os 4
## datasets sao humanos apos a substituicao registrada em CLAUDE.md S9.
MACS3_GSIZE <- "hs"
