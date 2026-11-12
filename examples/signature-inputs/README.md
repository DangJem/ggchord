# Existing protein signature annotations

Canonical tiny cached inputs and explicit derived TSV files are packaged once
under `inst/extdata/`, with `signature-sources.json` URLs, date, CC0 attribution,
checksums and transformation rules. The local generator is
`data-raw/prepare-signature-inputs.py`; it verifies the raw snapshots before
re-encoding their actual observed location fragments. No annotation engine runs.
Missing MD5/status/run date/integrated description remain `-`. The derived TSVs
are coordinate round-trip controls, not newly run InterProScan outputs.

API-returned accessions are lowercase and remain verbatim. All actual source
locations in these two snapshots are single-fragment; general JSON reading
reuses native Feature segments for discontinuities. Model IDs may repeat and
are not unique feature-instance IDs. No common score metric or source characters
are inferred. Full cached objects remain the provenance for any richer grouping
that the simpler TSV cannot express.
