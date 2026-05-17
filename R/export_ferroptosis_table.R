#' Export table of significant ferroptosis-related genes (ferroviz style)
#'
#' Cette fonction filtre les gènes significatifs selon un seuil de logFC et de
#' p-value ajustée, annote avec FerrDb, classe les gènes en Up/Down, puis
#' retourne un data.frame contenant uniquement les gènes ferroptosis
#' (driver / suppressor / both).
#'
#' @param resulths Data frame limma contenant logFC, adj.P.Val.
#' @param ferrdb Data frame FerrDb contenant au minimum "Symbol" et "cat".
#' @param fc Seuil absolu de logFC (défaut = 0.10).
#'
#' @return Un data.frame contenant les gènes ferroptosis significatifs.
#' @export
#'
#' @examples
#' tbl <- export_ferroptosis_table(resulths, ferrdb, fc = 0.1)
#' head(tbl)
#'
export_ferroptosis_table <- function(resulths, ferrdb, fc = 0.10) {

  library(dplyr)
  library(tidyr)

  # 1) Ajouter la colonne gene
  sig <- resulths %>%
    mutate(gene = row.names(.)) %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= fc) %>%
    mutate(regulation = ifelse(logFC >= 0, "up", "down"))

  # 2) Annotation FerrDb
  annotated <- sig %>%
    left_join(ferrdb, by = c("gene" = "Symbol")) %>%
    mutate(
      cat = replace_na(cat, "unclass"),
      ferroptosis = cat
    ) %>%
    filter(ferroptosis != "unclass")

  # 3) Table finale propre
  out <- annotated %>%
    select(
      gene,
      logFC,AveExpr,
      adj.P.Val,
      regulation,
      ferroptosis
    ) %>%
    arrange(ferroptosis, desc(logFC))
  out%>%distinct()->out
  return(out)
}
