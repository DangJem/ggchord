# Missing value handling operator

Used to safely handle NULL values: returns y if x is NULL, otherwise
returns x

## Usage

``` r
if_null_else(x, y)
```

## Arguments

- x:

  Any R object (may be NULL)

- y:

  Default value to return when x is NULL

## Value

x if x is not NULL, otherwise y
