Start with one complete circular backbone using bundled pUC19.

## Draw the backbone and centre label {#draw}

Run this lesson in a fresh R session after installing the version described on the [installation page](installation.html).

```r
library(ggchord)
```

{{code:firstcircle}}

{{figure:firstcircle}}

Ticks extend inward from the inner backbone edge, with labels farther inward. On a narrow canvas, the default removes colliding text while retaining the ticks.

## Read the ticks around the circle {#features}

One turn represents the full length of pUC19. The origin is at 12 o’clock by default, with coordinates increasing clockwise. Ticks show bp and the centre label identifies the sequence. The backbone and `coord_circular()` establish positions around the circle. The next section places a replication origin on it.

## Use your own sequence {#adapt}

Provide one sequence row with `accver` and full `length`, and use the same ID in the feature table. [Continue to features and labels](plasmid-build.html).
