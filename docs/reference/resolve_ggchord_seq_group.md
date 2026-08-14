# Resolve a seq_group specification into a named character vector

Resolve a seq_group specification into a named character vector

## Usage

``` r
resolve_ggchord_seq_group(seq_data, seqs, seq_group)
```

## Arguments

- seq_data:

  data.frame containing at least \`seq_id\` and, optionally, a
  \`seq_group\` column.

- seqs:

  Character vector of sequence IDs in drawing order.

- seq_group:

  NULL, a column name in \`seq_data\`, or a parameter accepted by
  \[process_sequence_param()\] (single value, named vector, unnamed
  vector matching the sequences, or a list).

## Value

A named character vector with one element per sequence (names are
sequence IDs), or NULL when grouping is disabled.
