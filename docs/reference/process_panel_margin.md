# Process panel margin parameters

Standardizes input margin parameters into a list containing t (top), r
(right), b (bottom), l (left). Supports single-value or list input.

## Usage

``` r
process_panel_margin(arg_list)
```

## Arguments

- arg_list:

  Numeric (single value) or list (named/unnamed), margin parameters

## Value

List containing four elements: t, r, b, l (numeric, margin sizes)
