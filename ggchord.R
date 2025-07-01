## ggchord: Draw a chord diagram for two sequences and their BLAST alignment intervals
#' ggchord: Draw a chord diagram for two sequences and their BLAST alignment intervals
#'
#' This function reads BLASTN outfmt6 format results:
#' - The upper semicircle (0~π) represents the query sequence, mapped to radius = r_query;
#' - The lower semicircle (π~2π) represents the subject sequence, mapped to radius = r_subject;
#' - Alignment intervals are represented by ribbon polygons connecting corresponding regions on the query and subject arcs.
#'
#' @param blast_df Data frame containing columns: qstart, qend, sstart, send, qlen, slen, length
#' @param min_len Minimum alignment length filter, default 100
#' @param title Plot title, default "BLAST Chord Diagram"
#' @param ribbon_col Color for the connecting ribbons, default "steelblue"
#' @param gap_frac Total angular gap fraction between the two sequence arcs, default 0.02
#' @param r_query Radius for the query arc, default 1.0
#' @param r_subject Radius for the subject arc, default 0.8
#' @return A ggplot2 object
#' @import ggplot2
#' @export
ggchord <- function(
  blast_df,
  min_len    = 100,
  title      = "BLAST Chord Diagram",
  ribbon_col = "steelblue",
  gap_frac   = 0.02,
  r_query    = 1.0,
  r_subject  = 0.8
) {
  library(ggplot2)
  # 1. Filter alignment data
  df <- subset(blast_df, length >= min_len)
  if (nrow(df) == 0) stop("No alignments meet the minimum length requirement")
  # 2. Get sequence lengths
  qlen <- unique(df$qlen); slen <- unique(df$slen)
  if (length(qlen)!=1 || length(slen)!=1) stop("qlen/slen should be constant values")
  total_len <- qlen + slen
  # 3. Circle circumference and gaps
  total_circle <- 2*pi
  total_gap    <- total_circle * gap_frac
  gap_each     <- total_gap / 2
  usable_angle <- total_circle - total_gap
  # 4. Assign angular ranges
  theta_q_len <- usable_angle * (qlen  / total_len)
  theta_s_len <- usable_angle * (slen  / total_len)
  theta_q0    <- gap_each
  theta_q1    <- theta_q0 + theta_q_len
  theta_s0    <- theta_q1 + gap_each
  theta_s1    <- theta_s0 + theta_s_len
  # 5. Generate sequence arc data
  nseg <- 200
  arc_q_theta <- seq(theta_q0, theta_q1, length.out = nseg)
  arc_s_theta <- seq(theta_s0, theta_s1, length.out = nseg)
  arc_q <- data.frame(x = r_query   * cos(arc_q_theta), y = r_query   * sin(arc_q_theta))
  arc_s <- data.frame(x = r_subject * cos(arc_s_theta), y = r_subject * sin(arc_s_theta))
  # 6. Calculate angular boundaries for alignment intervals
  df$q0_ang <- theta_q0 + (df$qstart - 1)/qlen * theta_q_len
  df$q1_ang <- theta_q0 + (df$qend   - 1)/qlen * theta_q_len
  df$s0_ang <- theta_s0 + (df$sstart - 1)/slen * theta_s_len
  df$s1_ang <- theta_s0 + (df$send   - 1)/slen * theta_s_len
  # 7. Construct ribbon polygons
  ribbon_list <- lapply(seq_len(nrow(df)), function(i) {
    # Query-side arc, inner boundary
    tq <- seq(df$q0_ang[i], df$q1_ang[i], length.out = 50)
    xq <- r_query   * cos(tq); yq <- r_query   * sin(tq)
    # Subject-side arc, outer boundary, reversed order
    ts <- seq(df$s0_ang[i], df$s1_ang[i], length.out = 50)
    xs <- r_subject * cos(ts); ys <- r_subject * sin(ts)
    # Combine polygon points
    data.frame(
      x = c(xq, rev(xs)),
      y = c(yq, rev(ys)),
      group = i
    )
  })
  ribbons <- do.call(rbind, ribbon_list)
  # 8. Plotting
  p <- ggplot() +
    geom_path(data = arc_q,   aes(x, y), size = 1) +
    geom_path(data = arc_s,   aes(x, y), size = 1) +
    geom_polygon(data = ribbons, aes(x, y, group = group), fill = ribbon_col, color = NA, alpha = 0.6) +
    coord_equal() + theme_void() + ggtitle(title)
  return(p)
}

# Example usage:
blast_df <- read.table("vB_AbaM_CP14__PQ859668.1.o7", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(blast_df) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                        "qstart","qend","sstart","send","evalue","bitscore",
                        "%qcov_per_subject","qlen","slen","strand","title")
p <- ggchord(
  blast_df, min_len=500, title="vB_AbaM_CP14 vs PQ859668 Chord Diagram",
  gap_frac=0.02, r_query=1.2, r_subject=0.8
)
print(p)
