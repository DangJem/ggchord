# 独立的原生 ggplot2 图例测量实验；不改变 ggchord 默认 Guide。
library(ggplot2)
out <- "tools/long-term-audit"
labels <- paste("Synthetic category", seq_len(20))
input <- data.frame(x = seq_along(labels), y = 1, category = factor(labels, levels = labels))
path <- file.path(tempdir(), "ggchord-guide-probe.png")
grDevices::png(path, width = 8, height = 6, units = "in", res = 144)
measure <- function(rows) {
  p <- ggplot(input, aes(x, y, colour = category)) + geom_point() +
    guides(colour = guide_legend(nrow = rows, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
  table <- ggplotGrob(p)
  box <- table$grobs[[which(table$layout$name == "guide-box-bottom")]]
  width <- grid::convertWidth(sum(box$widths), "mm", valueOnly = TRUE)
  height <- grid::convertHeight(sum(box$heights), "mm", valueOnly = TRUE)
  data.frame(rows = rows, width_mm = width, height_mm = height,
    available_width_mm = 8 * 25.4 - 12, available_height_mm = 6 * 25.4 / 3,
    fits_budget = width <= 8 * 25.4 - 12 && height <= 6 * 25.4 / 3)
}
candidates <- do.call(rbind, lapply(c(1, 2, 4, 5, 10), measure))
grDevices::dev.off()
eligible <- candidates[candidates$fits_budget, ]
chosen <- if (nrow(eligible)) eligible[which.min(eligible$height_mm), ] else data.frame()
utils::write.csv(candidates, file.path(out, "guide-candidates.csv"), row.names = FALSE)
jsonlite::write_json(list(canvas_inches = c(8, 6), input = "20 synthetic categories",
  candidate_selected = chosen, scope = "measured candidates only; no integration, multi-guide solver or locked user settings",
  impossible = !nrow(eligible)), file.path(out, "guide-results.json"), pretty = TRUE, digits = NA)

# 在原生标度里验证候选切割取值；这是目标规格实验，不是默认已修复。
cuts <- c(0, 50, 80, 90, 95, 100)
sample <- c(0, 20, 49.9, 50, 60, 79.9, 80, 85, 89.9, 90, 92, 94.9, 95, 97, 99.9, 100)
scale <- scale_fill_viridis_b(limits = c(0, 100), breaks = cuts)
scale$train(sample)
mapped <- data.frame(pident = sample, fill = scale$map(sample))
utils::write.csv(mapped, file.path(out, "pident-candidate-mapping.csv"), row.names = FALSE)
utils::write.csv(data.frame(lower = head(cuts, -1), upper = tail(cuts, -1),
  fraction = diff(cuts) / 100), file.path(out, "pident-guide-required-segments.csv"), row.names = FALSE)
cat("Guide candidates measured on one 8x6 PNG device; public integration pending\n")
