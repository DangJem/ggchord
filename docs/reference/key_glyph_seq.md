# Key glyph for sequence legends

Draws the path symbol only when the key data contains colour; otherwise
returns a blank (prevents ggplot2 4.x from mixing unrelated layers into
other legends with default grey/black symbols).

## Usage

``` r
key_glyph_seq(data, params, size)
```
