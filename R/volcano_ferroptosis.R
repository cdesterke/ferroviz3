#' Volcano plot annoté avec FerrDb (driver / suppressor / both)
#'
#' Cette fonction génère un volcano plot à partir d'un objet limma (resulths)
#' et d'une annotation FerrDb (driver, suppressor, both).  
#'
#' @param resulths Data frame issu de limma contenant logFC, P.Value, adj.P.Val.
#' @param ferrdb Data frame FerrDb contenant au minimum la colonne "Symbol" et "cat".
#' @param fc_thresh Seuil absolu de logFC (défaut = 0.5).
#' @param p_thresh Seuil de p-value ajustée (défaut = 0.05).
#' @param plot_title Titre du graphique (défaut = "DEGs and ferroptosis").
#' @param base_size Taille de base du thème ggplot2 (défaut = 14).
#' @param label_size Taille des labels ggrepel (défaut = 4).
#'
#' @return Un objet ggplot2 représentant le volcano plot.
#' @export
#'
#' @examples
#' load("ferrdb.rda")
#' volcano_ferroptosis(resulths, ferrdb, plot_title = "Mon Volcano",fc_thresh=0.25,base_size=16)
#'
volcano_ferroptosis <- function(resulths, ferrdb,
                                fc_thresh = 0.5,
                                p_thresh = 0.05,
                                plot_title = "DEGs and ferroptosis",
                                base_size = 14,
                                label_size = 4) {

  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  # Préparation des données
  res <- resulths %>%
    mutate(gene = row.names(resulths)) %>%
    relocate(gene) %>%
    mutate(significance = ifelse(abs(logFC) >= fc_thresh & adj.P.Val <= p_thresh,
                                 "YES", "no")) %>%
    left_join(ferrdb, by = c("gene" = "Symbol"))

  # Catégories ferroptose
  res2 <- res %>%
    mutate(
      negLogP = -log10(P.Value),
      cat2 = case_when(
        cat == "both" ~ "both",
        cat == "driver" ~ "driver",
        cat == "suppressor" ~ "suppressor",
        TRUE ~ "other"
      )
    )

  # Volcano plot
  p <- ggplot() +
    geom_point(
      data = res2 %>% filter(cat2 == "other"),
      aes(x = logFC, y = negLogP),
      color = "grey70",
      size = 1.5
    ) +
    geom_point(
      data = res2 %>% filter(cat2 != "other"),
      aes(x = logFC, y = negLogP, color = cat2),
      size = 4
    ) +
    scale_color_manual(
      values = c(
        "driver" = "tomato",
        "suppressor" = "royalblue",
        "both" = "seagreen3"
      )
    ) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_thresh),
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = res2 %>%
        filter(cat2 %in% c("driver", "suppressor", "both"),
               adj.P.Val <= p_thresh) %>%
        distinct(gene, .keep_all = TRUE),
      aes(x = logFC, y = negLogP, label = gene, color = cat2),
      size = label_size,
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.3,
      show.legend = FALSE
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    ) +
    labs(
      title = plot_title,
      x = "log2 Fold Change",
      y = "-log10(P-value)",
      color = "Category"
    )

  return(p)
}
