# Coerce a validation result to a flat data.frame

Combines the `errors` and `warnings` tables into a single data.frame and
adds a `severity` column, which is convenient for filtering, exporting
or printing the full report programmatically.

## Usage

``` r
# S3 method for class 'ggchord_validation'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `"ggchord_validation"` object.

- row.names:

  Ignored.

- optional:

  Ignored.

- ...:

  Ignored.

## Value

A data.frame with columns `table`, `category`, `row`, `column`,
`message` and `severity`.
