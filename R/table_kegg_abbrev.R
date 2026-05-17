#' KEGG Abbreviation Table
#'
#' This function extracts all KEGG metabolic pathways and their corresponding
#' abbreviations from a KEGG annotation table (kegg2). It removes duplicates,
#' sorts entries alphabetically, and returns a clean tibble suitable for
#' documentation, export, or downstream visualization.
#'
#' @param kegg2 A tibble containing at least the columns:
#'   \describe{
#'     \item{kegg}{Full KEGG pathway name}
#'     \item{kegg_abbrev}{Short KEGG pathway abbreviation}
#'   }
#'
#' @return A tibble with unique KEGG pathways and their abbreviations.
#'
#' @examples
#' \dontrun{
#' kegg_table <- as.data.frame(table_kegg_abbrev(kegg2))
#' print(kegg_table)
#' }
#'
#' @export
table_kegg_abbrev <- function(kegg2) {

  library(dplyr)

  kegg2 %>%
    select(kegg, kegg_abbrev) %>%
    mutate(
      kegg = as.character(kegg),
      kegg_abbrev = as.character(kegg_abbrev)
    ) %>%
    distinct() %>%
    arrange(kegg)
}
