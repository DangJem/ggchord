## ggchord: A chord diagram function for pairwise sequences based on BLAST alignments (arc/line mode switch with precise ribbon alignment)
#' ggchord: Draw a chord diagram for two sequences and their BLAST alignment intervals
#'
#' This function reads BLASTN outfmt6 format results:
#' - Supports switching between "arc mode" and "line mode":
#'   - Arc mode: Trajectories are semi-circles or interpolated arcs, with curvature adjustable by `curvature`;
#'     Ribbon endpoints are strictly extracted from main arc point sets to ensure alignment;
#'   - Line mode: Trajectories are horizontal lines, with vertical distance between sequences adjustable by `line_gap_frac`;
#'
#' @param blast_df Data frame containing columns: qstart, qend, sstart, send, qlen, slen, length
#' @param min_len Minimum alignment length filter, default 100
#' @param title Plot title, default "BLAST Chord Diagram"
#' @param ribbon_col Color for the connecting ribbons, default "steelblue"
#' @param arc_mode Logical, TRUE for arc mode, FALSE for line mode, default TRUE
#' @param curvature Curvature of arcs in arc mode [0,1], 0=polyline (semi-circle polyline), 1=full arc, default 1
#' @param gap_frac_arc Horizontal gap fraction between sequence arcs in arc mode, default 0.02
#' @param line_gap_frac Vertical gap fraction between sequence lines in line mode, default 0.5
#' @param r_query Radius for query arc (arc mode) or Y-coordinate (ignored in line mode), default 1.0
#' @param r_subject Radius for subject arc (arc mode) or Y-coordinate (ignored in line mode), default 0.8
#' @return A ggplot2 object
#' @import ggplot2
#' @export
ggchord <- function(
    blast_df,
    min_len        = 100,
    title          = "BLAST Chord Diagram",
    ribbon_col     = "steelblue",
    arc_mode       = TRUE,
    curvature      = 1.0,
    gap_frac_arc   = 0.02,
    line_gap_frac  = 0.5,
    r_query        = 1.0,
    r_subject      = 0.8
) {
  library(ggplot2)
  df <- subset(blast_df, length >= min_len)
  if (nrow(df)==0) stop("No alignments meet the minimum length requirement")
  qlen <- unique(df$qlen); slen <- unique(df$slen)
  if (length(qlen)!=1||length(slen)!=1) stop("qlen/slen should be constant values")
  
  if (arc_mode) {
    # Arc mode...
    total_circle <- 2*pi
    total_gap    <- total_circle * gap_frac_arc
    gap_each     <- total_gap/2
    usable       <- total_circle - total_gap
    total_len    <- qlen + slen
    theta_q_len  <- usable * qlen/total_len
    theta_s_len  <- usable * slen/total_len
    th_q0 <- gap_each; th_q1 <- th_q0 + theta_q_len
    th_s0 <- th_q1 + gap_each; th_s1 <- th_s0 + theta_s_len
    nseg <- 500; u <- seq(0,1,length.out=nseg)
    # Generate main arcs
    gen_main <- function(theta0, theta1, radius) {
      ang <- seq(theta0, theta1, length.out=nseg)
      arc_xy  <- cbind(x=radius*cos(ang), y=radius*sin(ang))
      p0 <- c(radius*cos(theta0), radius*sin(theta0)); p1 <- c(radius*cos(theta1), radius*sin(theta1))
      lin_xy <- cbind((1-u)*p0[1]+u*p1[1], (1-u)*p0[2]+u*p1[2])
      main_xy <- (1-curvature)*lin_xy + curvature*arc_xy
      return(main_xy)
    }
    arc_q <- gen_main(th_q0, th_q1, r_query)
    arc_s <- gen_main(th_s0, th_s1, r_subject)
    # Index mapping for precise alignment
    idx_map <- function(angle, theta0, theta_len) {
      idx <- round((angle-theta0)/theta_len*(nseg-1))+1
      idx <- pmin(pmax(idx,1), nseg)
      return(idx)
    }
    ribbon_list <- lapply(seq_len(nrow(df)), function(i) {
      q0a <- th_q0 + (df$qstart[i]-1)/qlen*theta_q_len
      q1a <- th_q0 + (df$qend[i]  -1)/qlen*theta_q_len
      s0a <- th_s0 + (df$sstart[i]-1)/slen*theta_s_len
      s1a <- th_s0 + (df$send[i]  -1)/slen*theta_s_len
      i0q <- idx_map(q0a, th_q0, theta_q_len); i1q <- idx_map(q1a, th_q0, theta_q_len)
      i0s <- idx_map(s0a, th_s0, theta_s_len); i1s <- idx_map(s1a, th_s0, theta_s_len)
      seg_q <- arc_q[i0q:i1q,]; seg_s <- arc_s[i0s:i1s,]
      data.frame(x=c(seg_q[,1], rev(seg_s[,1])), y=c(seg_q[,2], rev(seg_s[,2])), group=i)
    })
    ribbons <- do.call(rbind, ribbon_list)
    p <- ggplot() +
      geom_path(data=data.frame(x=arc_q[,1],y=arc_q[,2]), aes(x,y), size=1) +
      geom_path(data=data.frame(x=arc_s[,1],y=arc_s[,2]), aes(x,y), size=1) +
      geom_polygon(data=ribbons, aes(x,y,group=group), fill=ribbon_col, alpha=0.6) +
      theme_void() + ggtitle(title)
  } else {
    # Line mode: Set two horizontal lines at y=+d/2 and y=-d/2
    sep <- line_gap_frac * (r_query + r_subject)
    yq <- sep/2; ys <- -sep/2
    arc_q <- data.frame(x=c(0, qlen), y=rep(yq,2))
    arc_s <- data.frame(x=c(0, slen), y=rep(ys,2))
    ribbon_list <- lapply(seq_len(nrow(df)), function(i) {
      xq0 <- df$qstart[i]-1; xq1 <- df$qend[i]-1
      xs0 <- df$sstart[i]-1; xs1 <- df$send[i]-1
      data.frame(x=c(xq0,xq1,xs1,xs0), y=c(rep(yq,2), rep(ys,2)), group=i)
    })
    ribbons <- do.call(rbind, ribbon_list)
    p <- ggplot() +
      geom_path(data=arc_q, aes(x,y), size=1) +
      geom_path(data=arc_s, aes(x,y), size=1) +
      geom_polygon(data=ribbons, aes(x,y,group=group), fill=ribbon_col, alpha=0.6) +
      theme_void() + ggtitle(title)
  }
  return(p)
}

# Complete example script
blast_df <- read.table("vB_AbaM_CP14__PQ859668.1.o7", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
colnames(blast_df) <- c("qaccver","saccver","pident","length","mismatch","gapopen",
                        "qstart","qend","sstart","send","evalue","bitscore",
                        "qcovs","qlen","slen","sstrand","stitle")
# Arc mode example
p_arc <- ggchord(blast_df, min_len=5000, title="Arc Mode", arc_mode=TRUE, curvature=1, gap_frac_arc=0, r_query=1, r_subject=1)
print(p_arc)
# Line mode example
p_line <- ggchord(blast_df, min_len=5000, title="Line Mode", arc_mode=FALSE, line_gap_frac=0, r_query=1, r_subject=1)
print(p_line)
