# Chemistry colour values: ggmsa scheme_NT / scheme_AA, Artistic-2.0.
# Only palette values are adapted; no ggmsa rendering code is used. See inst/NOTICE.

#' Colour residue cells by chemical group
#'
#' A standard discrete fill scale for [geom_residue()]. DNA/RNA use separate
#' colours for A, C, G and T/U; AA groups residues by side-chain chemistry.
#' Gaps, stops and ambiguous symbols are neutral. Lowercase symbols receive
#' the same colours as uppercase symbols. The chemistry colours follow ggmsa;
#' this scale does not compute conservation or an alignment.
#'
#' @param alphabet Explicit DNA, RNA or AA palette. Defaults to DNA for direct
#'   calls. geom_residue() uses the data's alphabet column, otherwise AA.
#' @param ... Additional arguments to [ggplot2::scale_fill_manual()], including
#'   name, breaks and guide. Use scale_fill_manual() for custom colours.
#' @return A ggplot2 discrete fill scale. Explicit user scales take priority
#'   over a residue layer's inferred default.
#' @export
#' @examples
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), "DNA")
#' ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue() +
#'   scale_fill_residue("DNA", guide = "none")
scale_fill_residue <- function(alphabet = c("DNA", "RNA", "AA"), ...) {
  alphabet <- match.arg(alphabet)
  values <- stats::setNames(rep("#FFFFFF", 29), c(LETTERS, "-", ".", "*"))
  if (alphabet == "AA") {
    groups <- list("#F2BE3C" = c("A", "G", "I", "L", "M", "P", "V"),
      "#FF6D6D" = c("D", "E"), "#769DCC" = c("H", "K", "R"),
      "#74CE98" = c("C", "N", "Q", "S", "T"), "#FFFF66" = c("F", "W", "Y"))
    for (colour in names(groups)) values[groups[[colour]]] <- colour
  } else {
    values[c("A", "C", "G", "T", "U")] <- c("#FF6D6D", "#769DCC", "#F2BE3C", "#74CE98", "#74CE98")
  }
  values <- c(values, stats::setNames(values[LETTERS], tolower(LETTERS)))
  ggplot2::scale_fill_manual(values = values, ...)
}
