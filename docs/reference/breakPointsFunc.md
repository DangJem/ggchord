# Generate major axis tick breakpoints

Generates uniform and visually appealing major tick positions based on
sequence length and target tick count (avoids excessively short end
ticks)

## Usage

``` r
breakPointsFunc(max_value, n = 5, tol = 0.5)
```

## Arguments

- max_value:

  Numeric, sequence length (maximum value)

- n:

  Integer, target number of ticks, default 5

- tol:

  Numeric (0-1), tolerance threshold for end tick length (proportion of
  the median length of other ticks), default 0.5

## Value

Numeric vector, major tick positions (including 0 and max_value)
