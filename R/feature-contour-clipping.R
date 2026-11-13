# Clip a straight mark to the union of actual feature polygons. Arrow heads,
# compressed bodies and curved contours may be narrower than the nominal band.
ggchord_feature_line_intervals <- function(a, b, polygons) {
  delta <- b-a
  cross <- function(x,y) x[1L]*y[2L]-x[2L]*y[1L]
  times <- c(0,1)
  for (p in polygons) {
    xy <- cbind(p$x,p$y); xy <- rbind(xy,xy[1L,])
    for (j in seq_len(nrow(xy)-1L)) {
      edge <- xy[j+1L,]-xy[j,]
      den <- cross(delta,edge)
      if (abs(den)<1e-12) next
      t <- cross(xy[j,]-a,edge)/den
      u <- cross(xy[j,]-a,delta)/den
      if (t>0 && t<1 && u>=0 && u<=1) times <- c(times,t)
    }
  }
  times <- sort(unique(times))
  if (length(times)<2L) return(matrix(numeric(),ncol=2L))
  intervals <- cbind(head(times,-1L),tail(times,-1L))
  tolerance <- 64*.Machine$double.eps*max(1,abs(c(a,b)))
  intervals <- intervals[(intervals[,2L]-intervals[,1L])*sqrt(sum(delta^2))>
    tolerance,,drop=FALSE]
  keep <- vapply(rowMeans(intervals),function(t) {
    point <- a+t*delta
    any(vapply(polygons,function(p) ggchord_point_in_polygon(
      point[1L],point[2L],p),logical(1)))
  },logical(1))
  intervals[keep,,drop=FALSE]
}

ggchord_clip_feature_boundary <- function(xy, polygons) {
  intervals <- ggchord_feature_line_intervals(xy[1L,],xy[2L,],polygons)
  lapply(seq_len(nrow(intervals)),function(i) {
    rbind(xy[1L,]+intervals[i,1L]*(xy[2L,]-xy[1L,]),
      xy[1L,]+intervals[i,2L]*(xy[2L,]-xy[1L,]))
  })
}
