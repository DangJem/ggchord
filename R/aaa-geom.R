# Shared geom helper loaded before the public layer files.

#' Clone a geom and expose selected standard aesthetics under role names
#' @noRd
rename_geom_aes <- function(geom = ggplot2::GeomPolygon, renames) {
  new_geom <- ggplot2::ggproto(
    paste0("GeomChord", sub("^Geom", "", class(geom)[1])), geom
  )

  aes_names <- names(new_geom$default_aes)
  for (old in names(renames)) {
    aes_names[aes_names == old] <- renames[[old]]
  }
  names(new_geom$default_aes) <- aes_names

  old_handle_na <- geom$handle_na
  new_geom$handle_na <- function(self, data, params) {
    for (old in names(renames)) {
      colnames(data)[colnames(data) == renames[[old]]] <- old
    }
    old_handle_na(data, params)
  }
  old_draw_key <- geom$draw_key
  new_geom$draw_key <- function(data, params, size) {
    for (old in names(renames)) {
      colnames(data)[colnames(data) == renames[[old]]] <- old
    }
    old_draw_key(data, params, size)
  }
  new_geom
}
