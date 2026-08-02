# Process gene-related parameters

Standardizes gene parameters (e.g., offset, width) from various input
formats (single value/vector/list) into a list seperated by sequence and
strand

## Usage

``` r
process_gene_param(param, seqs, param_name, default_value, is_logical = FALSE)
```

## Arguments

- param:

  Input parameter (can be NULL, single value, vector, list)

- seqs:

  Character vector, list of sequence IDs

- param_name:

  Character, name of the parameter (used in error messages)

- default_value:

  Default value when param is NULL

- is_logical:

  Logical, whether the parameter is logical (TRUE/FALSE), default FALSE

## Value

List (named by seq_id), where each element is a vector with "+"/"-"
(parameter values for positive/negative strands)
