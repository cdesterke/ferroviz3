#' Barplot ferroptosis style ferroviz (Up/Down proportions)
#'
#' Cette fonction filtre les gènes significatifs selon un seuil de logFC et de
#' p-value ajustée, annote avec FerrDb, classe les gènes en Up/Down, puis génère
#' un barplot proportionnel (position = "fill") montrant la proportion de gènes
#' ferroptosis (driver/suppressor/both) dans chaque direction de régulation.
#'
#' @param resulths Data frame limma contenant logFC, adj.P.Val.
#' @param ferrdbhs Data frame FerrDb contenant au minimum "hs.gene" et "class".
#' @param fc Seuil absolu de logFC (défaut = 0.10).
#' @param size Taille du texte du graphique (défaut = 16).
#'
#' @return Un objet ggplot2 représentant le barplot.
#' @export
#'
#' @examples
#' barploths(resulths, ferrdb, fc = 0.25, size = 16)
#'
barploths <- function(resulths, ferrdb, fc = 0.10, size = 16) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)

  # 1) Sélection des gènes significatifs
  sig <- resulths %>%
    mutate(gene = row.names(.)) %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= fc) %>%
    mutate(regulation = ifelse(logFC >= 0, "up", "down"))

  # 2) Annotation FerrDb (Symbol = gene)
  all <- sig %>%
    left_join(ferrdb, by = c("gene" = "Symbol")) %>%
    mutate(
      cat = replace_na(cat, "unclass"),
      ferroptosis = cat
    ) %>%
    filter(ferroptosis != "unclass")

  # 3) Barplot style ferroviz
  p <- ggplot(all, aes(x = factor(regulation), fill = factor(ferroptosis))) +
    geom_bar(position = "fill") +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_fill(vjust = 0.5),
      size = 8,
      colour = "white"
    ) +
    scale_fill_manual(values = c(
      "driver" = "tomato",
      "suppressor" = "royalblue",
      "both" = "seagreen3"
    )) +
    xlab(paste0("LogFC ≥ ", fc)) +
    ylab("Relative proportions") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    ) +
    theme(text = element_text(size = size))

  return(p)
}
