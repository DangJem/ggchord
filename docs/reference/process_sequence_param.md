# Process sequence-related parameters

Standardizes sequence parameters (e.g., radius, gap) from flexible input
formats into a vector named by sequence IDs. Supported formats:

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

  Input parameter (can be NULL, single value, vector, named vector, or
  list)

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

## Details

1\. A single value, recycled to every sequence. 2. A named vector by
sequence ID, e.g. \`c("MT108731.1" = 3, "OR222515.1" = 1)\`. 3. An
unnamed vector with length equal to the number of sequences (matched by
sequence order). 4. A list named by sequence ID. 5. A list named by
sequence order ("1", "2", ...). 6. An unnamed list matched by sequence
order (length 1 or equal to the number of sequences); a length-one list
recycles.
