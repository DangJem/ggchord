# Process sequence-related parameters

Standardizes sequence parameters (e.g., radius, gap) from various input
formats (single value/vector/named vector) into a vector named by
sequence IDs

## Usage

``` r
process_sequence_param(
  param,
  seqs,
  param_name,
  default_value = NULL,
  allow_null = FALSE
)
```

## Arguments

- param:

  Input parameter (can be NULL, single value, vector, named vector)

- seqs:

  Character vector, list of sequence IDs

- param_name:

  Character, name of the parameter (used in error messages)

- default_value:

  Default value when param is NULL, default NULL

- allow_null:

  Logical, whether to allow param to be NULL, default FALSE

## Value

Named vector (names are seq_ids), standardized parameter values
