# Least-squares packing with an explicit physical order. Adjacent labels have
# their required before/after clearance; only a run whose constraints bind is
# moved. Pool-adjacent-violators solves the equivalent isotonic regression.
ggchord_pack_ordered_axis <- function(preferred, order_value,
                                      before, after, gap) {
  n <- length(preferred)
  if (n < 2L) return(preferred)
  ord <- order(order_value, seq_len(n))
  distance <- after[ord][-n] + before[ord][-1L] + gap
  offset <- c(0, cumsum(distance))
  target <- preferred[ord] - offset
  means <- numeric(n)
  counts <- integer(n)
  first <- integer(n)
  last <- integer(n)
  blocks <- 0L
  for (i in seq_len(n)) {
    blocks <- blocks + 1L
    means[blocks] <- target[i]
    counts[blocks] <- 1L
    first[blocks] <- i
    last[blocks] <- i
    while (blocks > 1L && means[blocks - 1L] > means[blocks]) {
      count <- counts[blocks - 1L] + counts[blocks]
      means[blocks - 1L] <- (
        means[blocks - 1L] * counts[blocks - 1L] +
          means[blocks] * counts[blocks]
      ) / count
      counts[blocks - 1L] <- count
      last[blocks - 1L] <- last[blocks]
      blocks <- blocks - 1L
    }
  }
  fitted <- numeric(n)
  for (i in seq_len(blocks)) fitted[first[i]:last[i]] <- means[i]
  out <- numeric(n)
  out[ord] <- fitted + offset
  out
}
