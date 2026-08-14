# Normalise a named vector of group colours

Normalise a named vector of group colours

## Usage

``` r
resolve_ggchord_group_colors(colors, groups)
```

## Arguments

- colors:

  Named vector, or NULL. When NULL a default grey is used later.

- groups:

  Character vector of group names in display order.

## Value

Named character vector with names equal to \`groups\`; unnamed colour
vectors are recycled positionally.
