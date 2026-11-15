## Resubmission: ggchord 0.28.2

This revision follows the two rejected 0.28.1 submissions. Existing GitHub
release tags and assets are unchanged.

## Version and DESCRIPTION

The jump from CRAN 0.16.6 to 0.28.2 is intentional and reflects cumulative
project development; intermediate development versions were not all submitted
to CRAN. NEWS.md records the changes and migration information. DESCRIPTION
now uses "dot plots" and spells out "high-scoring segment pair" without the
redundant HSP abbreviation. "Synteny" is an intentional scientific term for
conserved relationships/order of genomic regions, not a spelling mistake.

## Optional reader examples

The VCF example previously read an entire compressed demonstration file before
selecting three records. It now stops after the original headers and the first
three observed records, verified identical to the previous selection. The
BigWig example retains one regional read of the existing packaged GC signal.
No API, scientific input values or complete website tutorial was changed.

Both optional Bioconductor decoder examples are in donttest blocks with an
explanation of their initialization cost. Complete runnable code and dependency/
platform guards are retained; no dontrun block was added. The blocks include
the dependency checks so that an ordinary example run does not initialize the
optional decoders. We understand that --as-cran also executes donttest code:
both final checks below include successful --run-donttest passes.

Independent fresh-process measurements found Windows VariantAnnotation loading
at 8.96 s and three-record parsing at 0.44 s; preparation and preview were below
the timer's resolution. Debian rtracklayer loading took 2.450 s, a regional read
0.055 s and preview output 0.002 s. The decoder initialization cost is retained,
not claimed to have been optimized away.

## R-devel validation

- Windows Server 2022 x64, R-devel 2026-11-14 r90890 ucrt:
  0 ERRORs, 0 WARNINGs, 1 NOTE (the intentional minor-version jump).
- Debian GNU/Linux forky/sid x86_64, R-devel 2026-11-08 r90859:
  0 ERRORs, 0 WARNINGs, 1 NOTE (the same version jump).

Both checks used R CMD build followed by R CMD check --as-cran, without
turning off examples, donttest examples or the PDF manual. Ordinary examples,
the --run-donttest pass, 34 testthat assertions, vignette rebuilding and PDF/HTML
manual checks passed. No examples timing or DESCRIPTION spelling NOTE remained
in these runs. Whole-topic times in the final recorded --run-donttest pass were
10.36 s for Windows VCF and 2.977 s for Debian BigWig. These sequential check
timings are distinct from the independent cold-load measurements above.

The two checks used identical package source; only CI environment setup changed
between commits. Initial runner failures (missing LaTeX/fonts on Windows and
libuv-dependent R packages on Debian) were repaired before the final results.
These GitHub runner results are not new CRAN incoming results or guarantees for
identical timings on CRAN hardware. No CRAN submission has been performed.

Windows run: https://github.com/DangJem/ggchord/actions/runs/43727617245
Debian run: https://github.com/DangJem/ggchord/actions/runs/43726759756
