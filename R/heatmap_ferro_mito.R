#' Heatmap of Up/Down Ferroptosis Categories vs Mitochondrial Pathways
#'
#' This function integrates differential expression results (from limma) with a
#' ferroptosis–mitochondria annotation table. It filters significant genes based
#' on adjusted p-value and log-fold change thresholds, assigns Up/Down
#' regulation, and generates a heatmap showing the number of significant genes
#' for each combination of:
#' 
#' - ferroptosis category (driver / suppressor / both)
#' - mitochondrial pathway (e.g., Complex V, Lipid metabolism, etc.)
#' - regulation direction (Up / Down)
#'
#' The heatmap is faceted by regulation direction (Up vs Down).
#'
#' @param resulths A limma results table containing at least: logFC, adj.P.Val.
#' @param mito A tibble containing ferroptosis and mitochondrial annotations,
#'   with columns:
#'   \describe{
#'     \item{Symbol}{Gene symbol}
#'     \item{cat}{Ferroptosis category: driver / suppressor / both}
#'     \item{mitochondria}{Mitochondrial pathway annotation (factor or character)}
#'   }
#' @param fc Numeric. Minimum absolute log-fold change required to consider a
#'   gene significant. Default: 0.1.
#' @param size Numeric. Base font size for the heatmap theme. Default: 14.
#'
#' @return A ggplot2 heatmap object.
#'
#' @examples
#' \dontrun{
#' p <- heatmap_ferro_mito(resulths, mito, fc = 0.1)
#' print(p)
#' }
#'
#' @export
heatmap_ferro_mito <- function(resulths, mito, fc = 0.1, size = 14) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)

  # 1) Annotate resulths with mitochondrial + ferroptosis categories
  df <- resulths %>%
    mutate(gene = row.names(.)) %>%
    filter(adj.P.Val <= 0.05 & abs(logFC) >= fc) %>%
    mutate(regulation = ifelse(logFC >= 0, "Up", "Down")) %>%
    left_join(mito, by = c("gene" = "Symbol")) %>%
    mutate(
      cat = replace_na(cat, "unclass"),
      mitochondria = as.character(mitochondria),
      mitochondria = replace_na(mitochondria, "none")
    ) %>%
    filter(cat != "unclass", mitochondria != "none")

  # 2) Count combinations
  tab <- df %>%
    count(regulation, cat, mitochondria)

  # 3) Heatmap
  p <- ggplot(tab, aes(x = mitochondria, y = cat, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), color = "white", size = 5) +
    scale_fill_gradient(low = "steelblue", high = "tomato") +
    facet_wrap(~ regulation, ncol = 1) +
    theme_minimal(base_size = size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = size + 2, face = "bold")
    ) +
    labs(
      title = "Ferroptosis × Mitochondrial pathways (Up vs Down)",
      x = "Mitochondrial pathway",
      y = "Ferroptosis category",
      fill = "Count"
    )

  return(p)
}
