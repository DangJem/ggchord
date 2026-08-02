# Generate a default categorical palette

Returns the first `n` Set1 colors, interpolating with
[`colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html) when
`n` exceeds 9.

## Usage

``` r
chord_default_palette(n)
```

## Arguments

- n:

  Number of colors requested

## Value

A character vector of `n` colors
