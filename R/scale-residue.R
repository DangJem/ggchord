#' Colour residue cells by alphabet
#'
#' @param alphabet DNA, RNA or AA.
#' @param name Fill legend title. The default is the alphabet name.
#' @param values Optional named colours overriding the defaults for those symbols.
#' @param na.value Colour for missing or unmapped values.
#' @param guide A standard discrete guide.
#' @param ... Additional options for a discrete fill scale.
#' @return An ordinary ggplot2 fill scale.
#' @details The DNA defaults are red A, blue C, yellow G and green T.
#'   Supply named `values` to replace individual symbol colours.
#' @export
#' @examples
#' library(ggplot2)
#'
#' cells <- as_ggchord_alignment(c(A = "AC-GT", B = "ACTGT"), alphabet = "DNA")
#' ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol)) +
#'   scale_fill_residue("DNA", values = c(G = "#D9A400"))
scale_fill_residue <- function(alphabet = "DNA", name = waiver(), values = NULL,
    na.value = "#B8C1C5", guide = "legend", ...) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function(), scale_fill_manual)
  alphabet <- ggchord_match_choice(alphabet, c("DNA", "RNA", "AA"), "alphabet")
  dots <- list(...)
  ggchord_check_param_names(dots, c("breaks", "labels", "limits", "drop", "na.translate",
    "expand", "call", "super"), "scale_fill_residue()")
  palette <- pal_residue(alphabet, values)
  if (!is.character(na.value) || length(na.value) != 1L || is.na(na.value))
    ggchord_stop("na.value must be one colour")
  tryCatch(col2rgb(c(palette, na.value)), error = function(e)
    ggchord_stop("scale_fill_residue(): ", conditionMessage(e)))
  if (is_waiver(name)) name <- alphabet
  scale_fill_manual(name = name, values = palette, na.value = na.value, guide = guide, ...)
}

#' Colours for residue symbols
#' @inheritParams scale_fill_residue
#' @return A named vector of colours for the selected alphabet, including gaps.
#' @export
pal_residue <- function(alphabet = "DNA", values = NULL) {
  ggchord_local_console_errors()
  alphabet <- ggchord_match_choice(alphabet, c("DNA", "RNA", "AA"), "alphabet")
  allowed <- if (alphabet == "AA") c(LETTERS, "*", "-", ".") else
    strsplit(if (alphabet == "RNA") "ACGURYSWKMBDHVN-." else "ACGTRYSWKMBDHVN-.", "")[[1L]]
  palette <- setNames(rep("#B8C1C5", length(allowed)), allowed)
  if (alphabet == "AA") {
    groups <- list("#E41A1C" = c("D", "E"), "#377EB8" = c("K", "R", "H"),
      "#4DAF4A" = c("S", "T", "N", "Q"), "#FFB000" = c("A", "V", "I", "L", "M", "F", "W", "Y"),
      "#984EA3" = "C", "#999999" = "G", "#A65628" = "P", "#222222" = "*")
    for (colour in names(groups)) palette[groups[[colour]]] <- colour
  } else if (alphabet == "RNA") {
    palette[c("A", "C", "G", "U")] <- c("#4DAF4A", "#377EB8", "#FFB000", "#E41A1C")
  } else {
    palette[c("A", "C", "G", "T")] <- c("#FF6D6D", "#769DCC", "#F2BE3C", "#74CE98")
  }
  palette[c("-", ".")] <- "transparent"
  if (!is.null(values)) {
    if (!is.character(values) || is.null(names(values)) || anyNA(values) ||
        anyNA(names(values)) || anyDuplicated(toupper(names(values))) ||
        any(!toupper(names(values)) %in% allowed))
      ggchord_stop("values must provide colours for distinct valid residue symbols")
    palette[toupper(names(values))] <- values
  }
  tryCatch(col2rgb(palette), error = function(e)
    ggchord_stop("pal_residue(): ", conditionMessage(e)))
  palette
}
