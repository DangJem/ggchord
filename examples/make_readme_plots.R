# From the repository root, use the same standalone plotting examples as the website.
# The README circular flagship uses the dense pETDuet-1 record (139 enzyme rows).
library(ggchord)
library(ggplot2)

specs <- jsonlite::read_json("website/handbook-examples.json")
file_names <- c(chord = "showcasechord", collinear = "showcasecollinear",
  radial = "flagradial", sankey = "flagsankey", dotplot = "flagdotplot", fusion = "flagfusion")
for (view in names(file_names)) {
  id <- paste0("case-", view)
  plot <- source(file.path("website/examples", paste0(id, ".R")), local = new.env())$value
  size <- specs[[id]]
  ggsave(file.path("man/figures", paste0(file_names[[view]], ".png")), plot,
    width = size$width, height = size$height, dpi = size$dpi, units = "in", bg = size$background)
}

plasmid_example <- new.env()
source("website/examples/map-duet.R", local = plasmid_example)
plasmid <- plasmid_example$draw_reference_map(flagship_circular_example$dense)
ggsave("man/figures/showcaseplasmid.png", plasmid, width = 16.5, height = 8,
  dpi = 144, units = "in", bg = "#FFFCF5")
